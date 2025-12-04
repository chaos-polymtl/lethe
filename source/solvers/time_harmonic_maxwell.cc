// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/time_harmonic_maxwell.h>

using VectorType = GlobalVectorType;

template <int dim>
TimeHarmonicMaxwell<dim>::TimeHarmonicMaxwell(
  MultiphysicsInterface<dim>      *multiphysics_interface,
  const SimulationParameters<dim> &p_simulation_parameters,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  std::shared_ptr<SimulationControl> p_simulation_control)
  : AuxiliaryPhysics<dim, GlobalVectorType>()
  , multiphysics(multiphysics_interface)
  , computing_timer(p_triangulation->get_mpi_communicator(),
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , simulation_parameters(p_simulation_parameters)
  , triangulation(p_triangulation)
  , simulation_control(std::move(p_simulation_control))
  , dof_handler_trial_interior(
      std::make_shared<DoFHandler<dim>>(*triangulation))
  , dof_handler_trial_skeleton(
      std::make_shared<DoFHandler<dim>>(*triangulation))
  , dof_handler_test(std::make_shared<DoFHandler<dim>>(*triangulation))
  , extractor_E_real(0)
  , extractor_E_imag(dim)
  , extractor_H_real(2 * dim)
  , extractor_H_imag(3 * dim)
{
  if (simulation_parameters.mesh.simplex)
    {
      // for simplex meshes
      AssertThrow(
        false,
        ExcMessage(
          "TimeHarmonicMaxwell solver not yet implemented for simplex meshes."));
    }
  else
    {
      AssertThrow(dim == 3,
                  ExcMessage(
                    "TimeHarmonicMaxwell only implemented for 3D problems."));

      AssertThrow(
        simulation_parameters.fem_parameters.electromagnetics_trial_order <
          simulation_parameters.fem_parameters.electromagnetics_test_order,
        ExcMessage(
          "The DPG method requires the test space to be of higher order than the trial space."));

      // Usual case, for quad/hex meshes
      fe_trial_interior = std::make_shared<FESystem<dim>>(
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim,
        FE_DGQ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order) ^
          dim);
      fe_trial_skeleton = std::make_shared<FESystem<dim>>(
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_trial_order));
      fe_test = std::make_shared<FESystem<dim>>(
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order),
        FE_NedelecSZ<dim>(
          simulation_parameters.fem_parameters.electromagnetics_test_order));
      mapping = std::make_shared<MappingQ<dim>>(fe_trial_interior->degree);
      cell_quadrature = std::make_shared<QGauss<dim>>(fe_test->degree + 1);
      face_quadrature = std::make_shared<QGauss<dim - 1>>(fe_test->degree + 1);
    }

  // Initialize solution shared_ptr
  present_solution = std::make_shared<GlobalVectorType>();

  // Allocate solution transfer
  solution_transfer = std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(
    *dof_handler_trial_interior);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
TimeHarmonicMaxwell<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;

  // Interior output setup
  std::vector<std::string> solution_interior_names(dim, "E_real");
  for (unsigned int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("E_imag");
    }
  for (unsigned int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("H_real");
    }
  for (unsigned int i = 0; i < dim; ++i)
    {
      solution_interior_names.emplace_back("H_imag");
    }

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    solution_interior_data_component_interpretation(
      4 * dim, DataComponentInterpretation::component_is_part_of_vector);

  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    *this->dof_handler_trial_interior,
    *this->present_solution,
    solution_interior_names,
    solution_interior_data_component_interpretation);

  // Skeleton output setup
  // TODO:  it will need its own writer probably because we need to use
  // DataOutFaces object

  return solution_output_structs;
}

template <int dim>
std::vector<double>
TimeHarmonicMaxwell<dim>::calculate_L2_error()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Interior L2 error
  FEValues<dim> fe_values_trial_interior(*this->mapping,
                                         *this->fe_trial_interior,
                                         *this->cell_quadrature,
                                         update_values |
                                           update_quadrature_points |
                                           update_JxW_values);

  const unsigned int n_q_points = this->cell_quadrature->size();

  // The exact solution will be defined by user but will need to be on all
  // possible fields of the ultraweak formulation so we need 4*dim components.
  std::vector<Vector<double>> exact_solution_values(n_q_points,
                                                    Vector<double>(4 * dim));
  auto                       &exact_solution =
    simulation_parameters.analytical_solution->electromagnetics;

  // When looping on each cell we will extract the different field
  // solution obtained numerically. The containers used to store the
  // interpolated solution at the quadrature points are declared below.
  std::vector<Tensor<1, dim>> local_E_values_real(n_q_points);
  std::vector<Tensor<1, dim>> local_E_values_imag(n_q_points);
  std::vector<Tensor<1, dim>> local_H_values_real(n_q_points);
  std::vector<Tensor<1, dim>> local_H_values_imag(n_q_points);

  // We create variables that will store all the integration result we are
  // interested in.
  double L2_error_E_real = 0;
  double L2_error_E_imag = 0;
  double L2_error_H_real = 0;
  double L2_error_H_imag = 0;

  for (const auto &cell : dof_handler_trial_interior->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_trial_interior.reinit(cell);

          // Get the simulated solution at quadrature points
          fe_values_trial_interior[extractor_E_real].get_function_values(
            *present_solution, local_E_values_real);
          fe_values_trial_interior[extractor_E_imag].get_function_values(
            *present_solution, local_E_values_imag);
          fe_values_trial_interior[extractor_H_real].get_function_values(
            *present_solution, local_H_values_real);
          fe_values_trial_interior[extractor_H_imag].get_function_values(
            *present_solution, local_H_values_imag);

          // Get the exact solution at quadrature points
          exact_solution.vector_value_list(
            fe_values_trial_interior.get_quadrature_points(),
            exact_solution_values);

          // Loop on quadrature points to compute the L2 error contributions
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values_trial_interior.JxW(q);

              // Loop on dimensions to compute the squared error
              for (unsigned int d = 0; d < dim; ++d)
                {
                  // E real part
                  L2_error_E_real +=
                    pow(local_E_values_real[q][d] - exact_solution_values[q][d],
                        2) *
                    JxW;

                  // E imag part
                  L2_error_E_imag += pow(local_E_values_imag[q][d] -
                                           exact_solution_values[q][d + dim],
                                         2) *
                                     JxW;
                  // H real part
                  L2_error_H_real +=
                    pow(local_H_values_real[q][d] -
                          exact_solution_values[q][d + 2 * dim],
                        2) *
                    JxW;
                  // H imag part
                  L2_error_H_imag +=
                    pow(local_H_values_imag[q][d] -
                          exact_solution_values[q][d + 3 * dim],
                        2) *
                    JxW;
                }
            }
        }
    }
  // Skeleton L2 error
  // TODO

  L2_error_E_real = Utilities::MPI::sum(L2_error_E_real, mpi_communicator);
  L2_error_E_imag = Utilities::MPI::sum(L2_error_E_imag, mpi_communicator);
  L2_error_H_real = Utilities::MPI::sum(L2_error_H_real, mpi_communicator);
  L2_error_H_imag = Utilities::MPI::sum(L2_error_H_imag, mpi_communicator);

  return {L2_error_E_real, L2_error_E_imag, L2_error_H_real, L2_error_H_imag};
}


template <int dim>
void
TimeHarmonicMaxwell<dim>::finish_simulation()
{
  auto         mpi_communicator = this->triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity !=
        Parameters::Verbosity::quiet)
    {
      this->error_table.evaluate_all_convergence_rates(
        ConvergenceTable::reduction_rate_log2);

      this->error_table.set_scientific("error_E_real", true);
      this->error_table.set_scientific("error_E_imag", true);
      this->error_table.set_scientific("error_H_real", true);
      this->error_table.set_scientific("error_H_imag", true);
      this->error_table.set_precision(
        "error_E_real", this->simulation_control->get_log_precision());
      this->error_table.set_precision(
        "error_E_imag", this->simulation_control->get_log_precision());
      this->error_table.set_precision(
        "error_H_real", this->simulation_control->get_log_precision());
      this->error_table.set_precision(
        "error_H_imag", this->simulation_control->get_log_precision());
      this->error_table.write_text(std::cout);
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::percolate_time_vectors()
{
  // No time-dependent vectors to percolate in time-harmonic Maxwell
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::modify_solution()
{
  // No modification of the solution is required at the moment
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions.time_dependent)
    return;

  AssertThrow(
    false,
    ExcMessage(
      "Time-dependent boundary conditions not yet implemented for TimeHarmonicMaxwell."));
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      std::vector<double> errors       = calculate_L2_error();
      double              E_real_error = errors[0];
      double              E_imag_error = errors[1];
      double              H_real_error = errors[2];
      double              H_imag_error = errors[3];

      this->error_table.add_value("cells",
                                  this->triangulation->n_global_active_cells());
      this->error_table.add_value("error_E_real", E_real_error);
      this->error_table.add_value("error_E_imag", E_imag_error);
      this->error_table.add_value("error_H_real", H_real_error);
      this->error_table.add_value("error_H_imag", H_imag_error);

      if (simulation_parameters.analytical_solution->verbosity !=
          Parameters::Verbosity::quiet)
        {
          this->pcout << "L2 error E real: " << E_real_error << std::endl;
          this->pcout << "L2 error E imag: " << E_imag_error << std::endl;
          this->pcout << "L2 error H real: " << H_real_error << std::endl;
          this->pcout << "L2 error H imag: " << H_imag_error << std::endl;
        }
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "Time Harmonic Electromagnetics");
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::pre_mesh_adaptation()
{
  this->solution_transfer->prepare_for_coarsening_and_refinement(
    *this->present_solution);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Set up the vectors for the transfer
  GlobalVectorType tmp(this->locally_owned_dofs_trial_interior,
                       mpi_communicator);

  // Interpolate the solution at time and previous time
  this->solution_transfer->interpolate(tmp);

  // Distribute constraints
  this->nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  *this->present_solution = tmp;
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::write_checkpoint()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::read_checkpoint()
{
  // TODO
}

template <int dim>
std::vector<OutputStructTableHandler>
TimeHarmonicMaxwell<dim>::gather_tables()
{
  // TODO
  return std::vector<OutputStructTableHandler>();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::compute_energy_norm(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_dofs()
{
  verify_consistency_of_boundary_conditions();

  auto mpi_communicator = triangulation->get_mpi_communicator();

  // Setup each dof handlers
  this->dof_handler_trial_interior->distribute_dofs(*this->fe_trial_interior);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler_trial_interior);
  this->dof_handler_trial_skeleton->distribute_dofs(*this->fe_trial_skeleton);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler_trial_skeleton);
  this->dof_handler_test->distribute_dofs(*this->fe_test);
  DoFRenumbering::Cuthill_McKee(*this->dof_handler_test);

  // Get the locally owned dofs
  this->locally_owned_dofs_trial_interior =
    this->dof_handler_trial_interior->locally_owned_dofs();
  this->locally_owned_dofs_trial_skeleton =
    this->dof_handler_trial_skeleton->locally_owned_dofs();
  this->locally_owned_dofs_test = this->dof_handler_test->locally_owned_dofs();

  // Get the locally relevant dofs
  this->locally_relevant_dofs_trial_interior =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_trial_interior);
  this->locally_relevant_dofs_trial_skeleton =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_trial_skeleton);
  this->locally_relevant_dofs_test =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler_test);

  // Initialize the solution vector
  this->present_solution->reinit(this->locally_owned_dofs_trial_interior,
                                 this->locally_relevant_dofs_trial_interior,
                                 mpi_communicator);
  this->present_solution_skeleton->reinit(
    this->locally_owned_dofs_trial_skeleton,
    this->locally_relevant_dofs_trial_skeleton,
    mpi_communicator);

  // We reinitialize the system rhs with the skeleton dofs because we have
  // performed a static condensation of the interior dofs using the Schur
  // complement.
  this->system_rhs.reinit(this->locally_owned_dofs_trial_skeleton,
                          mpi_communicator);

  // Define constraints
  define_constraints();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(this->locally_relevant_dofs_trial_skeleton);
  DoFTools::make_sparsity_pattern(*this->dof_handler_trial_skeleton,
                                  dsp,
                                  this->nonzero_constraints,
                                  /*keep_constrained_dofs = */ false);
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->locally_owned_dofs_trial_skeleton,
    mpi_communicator,
    this->locally_relevant_dofs_trial_skeleton);

  this->system_matrix.reinit(this->locally_owned_dofs_trial_skeleton,
                             this->locally_owned_dofs_trial_skeleton,
                             dsp,
                             mpi_communicator);

  this->pcout << "   Number of skeleton degrees of freedom: "
              << this->dof_handler_trial_skeleton->n_dofs() << std::endl;
  this->pcout << "   Number of interior degrees of freedom: "
              << this->dof_handler_trial_interior->n_dofs() << std::endl;

  // Provide the TimeHarmonicMaxwell dof_handler and solution to the
  // multiphysics interface. Note that we provide the interior dof_handler and
  // solution has it is the one useful for the other physics.
  multiphysics->set_dof_handler(PhysicsID::electromagnetics,
                                this->dof_handler_trial_interior);
  multiphysics->set_solution(PhysicsID::electromagnetics,
                             this->present_solution);
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::setup_preconditioner()
{
  preconditioner = std::make_shared<TrilinosWrappers::PreconditionIdentity>();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::define_constraints()
{
  // Clear previous constraints
  this->nonzero_constraints.clear();
  this->nonzero_constraints.reinit(this->locally_owned_dofs_trial_skeleton,
                                   this->locally_relevant_dofs_trial_skeleton);

  DoFTools::make_hanging_node_constraints(*this->dof_handler_trial_skeleton,
                                          this->nonzero_constraints);

  // Loop over all boundary conditions defined
  for (const auto &[id, type] :
       this->simulation_parameters
         .boundary_conditions_time_harmonic_electromagnetics.type)
    {
      if (type == BoundaryConditions::BoundaryType::pec)
        {
          // Perfect electric conductor (PEC) boundary condition
          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            0,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            dim,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::pmc)
        {
          // Perfect magnetic conductor (PMC) boundary condition
          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            2 * dim,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            3 * dim,
            Functions::ZeroFunction<dim>(4 * dim),
            id,
            this->nonzero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::imposed_electric_field)
        {
          // Imposed electric field boundary condition

          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            0,
            TimeHarmonicMaxwellElectricFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_imag),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            dim,
            TimeHarmonicMaxwellElectricFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->e_z_imag),
            id,
            this->nonzero_constraints);
        }

      if (type == BoundaryConditions::BoundaryType::imposed_magnetic_field)
        {
          // Imposed magnetic field boundary condition

          // Real
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            0,
            TimeHarmonicMaxwellMagneticFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_imag),
            id,
            this->nonzero_constraints);

          // Imaginary
          VectorTools::project_boundary_values_curl_conforming_l2(
            *this->dof_handler_trial_skeleton,
            dim,
            TimeHarmonicMaxwellMagneticFieldDefined<dim>(
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_real,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_x_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_y_imag,
              &this->simulation_parameters
                 .boundary_conditions_time_harmonic_electromagnetics
                 .imposed_electromagnetic_fields.at(id)
                 ->h_z_imag),
            id,
            this->nonzero_constraints);
        }
    }

  // The DPG method requires the use of skeleton elements and shape functions.
  // Because we want to use high-order Nedelec elements for the face, but dealii
  // does not support a FE_FaceNedelec, we hack our way around by using full
  // cell elements, but freezing the interior dofs using constrain_dofs_to_zero.
  // To do so, we create a container for all dof indices and a container for the
  // face dof indices and we loop on all the faces of each cell to find which
  // dofs are on the face and flag them. The remaining dofs are then constrained
  // to zero.

  std::vector<types::global_dof_index> cell_dof_indices(
    this->fe_trial_skeleton->n_dofs_per_cell());
  std::vector<types::global_dof_index> face_dof_indices(
    this->fe_trial_skeleton->n_dofs_per_face());
  std::vector<bool> is_dof_on_face(this->fe_trial_skeleton->n_dofs_per_cell(),
                                   false);

  // Loop on all skeleton dofs and set the interior constraints to zero
  for (const auto &cell :
       this->dof_handler_trial_skeleton->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Get all dof indices on the cell
          cell->get_dof_indices(cell_dof_indices);

          // Loop on all the faces of the cell
          for (const auto &face : cell->face_iterators())
            {
              face->get_dof_indices(face_dof_indices);

              // Loop on all dofs on the face
              for (const auto &face_dof : face_dof_indices)
                {
                  // Find the first iterator in the cell dof indices that
                  // matches the face dof (find where the face dof is in the
                  // cell dof indices)
                  const auto it = std::find(cell_dof_indices.begin(),
                                            cell_dof_indices.end(),
                                            face_dof);

                  // If the dof is on the face (find returns the second
                  // iterator if no match is found), set the corresponding
                  // flag to true
                  if (it != cell_dof_indices.end())
                    {
                      is_dof_on_face[std::distance(cell_dof_indices.begin(),
                                                   it)] = true;
                    }
                }
            }

          // Loop on all dofs on the cell and constrain the interior ones to
          // zero
          for (unsigned int index = 0; index < cell_dof_indices.size(); ++index)
            {
              // If the dof is not on a face, then it is an interior dof
              if (!is_dof_on_face[index])
                {
                  this->nonzero_constraints.constrain_dof_to_zero(
                    cell_dof_indices[index]);
                }
            }
        }
    }
  this->nonzero_constraints.close();
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::solve_linear_system()
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = this->triangulation->get_mpi_communicator();

  // Define the linear solver tolerance
  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
      .minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
      .relative_residual;

  const double rescale_metric    = this->get_residual_rescale_metric();
  const double rescaled_residual = this->system_rhs.l2_norm() / rescale_metric;
  const double linear_solver_tolerance =
    std::max(relative_residual * rescaled_residual, absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  // Set up the solver control
  setup_preconditioner();

  SolverControl solver_control(this->dof_handler_trial_skeleton->n_dofs(),
                               linear_solver_tolerance);

  // Solve
  GlobalVectorType completely_distributed_solution(
    this->locally_owned_dofs_trial_skeleton, mpi_communicator);
  TrilinosWrappers::SolverCG solver(solver_control);

  solver.solve(this->system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               *this->preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::electromagnetics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -CG iterative solver took : "
                  << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() / rescale_metric << std::endl;
    }

  // Update the solution vector
  nonzero_constraints.distribute(completely_distributed_solution);
  *this->present_solution_skeleton = completely_distributed_solution;
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_matrix()
{
  // TODO
}

template <int dim>
void
TimeHarmonicMaxwell<dim>::assemble_system_rhs()
{
  // TODO
}



template class TimeHarmonicMaxwell<2>;
template class TimeHarmonicMaxwell<3>;
