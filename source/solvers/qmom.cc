// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/boundary_conditions.h>
#include <core/dimensionality.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <solvers/analytical_solutions.h>
#include <solvers/initial_conditions.h>
#include <solvers/qmom.h>
#include <solvers/qmom_scratch_data.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_component_interpretation.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>

template <int dim>
QMOM<dim>::QMOM(
  MultiphysicsInterface<dim>      *multiphysics_interface,
  const SimulationParameters<dim> &p_simulation_parameters,
  const std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                           &p_triangulation,
  const std::shared_ptr<SimulationControl> &p_simulation_control)
  : AuxiliaryPhysics<dim, GlobalVectorType>()
  , n_scalars(
      p_simulation_parameters.multiphysics.qmom_parameters.number_of_scalars)
  , n_moments(
      p_simulation_parameters.multiphysics.qmom_parameters.number_of_moments)
  , n_fields(n_scalars + n_moments)
  , multiphysics(multiphysics_interface)
  , computing_timer(p_triangulation->get_mpi_communicator(),
                    this->pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , simulation_parameters(p_simulation_parameters)
  , triangulation(p_triangulation)
  , simulation_control(p_simulation_control)
  , dof_handler(std::make_shared<DoFHandler<dim>>(*triangulation))
{
  const Parameters::QMOM &qmom_parameters =
    simulation_parameters.multiphysics.qmom_parameters;

  // Relative tolerance used to compare parameter values
  constexpr double parameter_relative_tolerance = 1e-12;
  const auto       values_are_equal             = [](const double first_value,
                                   const double second_value) {
    return std::abs(first_value - second_value) <=
           parameter_relative_tolerance *
             std::max({1., std::abs(first_value), std::abs(second_value)});
  };

  // The QMOM phase properties and operating conditions are given in SI units
  // and are not rescaled by the dimensionality subsection
  const Parameters::Dimensionality &dimensionality =
    simulation_parameters.dimensionality;
  const bool dimensionality_is_si =
    values_are_equal(dimensionality.length, 1.) &&
    values_are_equal(dimensionality.mass, 1.) &&
    values_are_equal(dimensionality.time, 1.) &&
    values_are_equal(dimensionality.temperature, 1.);
  AssertThrow(dimensionality_is_si,
              QMOMDimensionalityError(
                "length = " + std::to_string(dimensionality.length) +
                ", mass = " + std::to_string(dimensionality.mass) +
                ", time = " + std::to_string(dimensionality.time) +
                " and temperature = " +
                std::to_string(dimensionality.temperature)));

  // The transport uses the fluid velocity without mesh velocity correction
  AssertThrow(!simulation_parameters.ale.enabled(),
              ExcMessage("The QMOM physics does not support the arbitrary "
                         "Lagrangian-Eulerian (ALE) formulation yet."));

  // The aqueous and organic phases form a single fluid, whose QMOM properties
  // are used by the phase-space evolution
  AssertThrow(
    simulation_parameters.physical_properties_manager.get_number_of_fluids() ==
      1,
    QMOMInvalidNumberOfFluids(simulation_parameters.physical_properties_manager
                                .get_number_of_fluids()));

  // Features of the QMOM physics that are not implemented yet
  AssertThrow(n_scalars == 0,
              ExcMessage("The transport of the QMOM scalars is not "
                         "implemented yet, but \"number of scalars\" is " +
                         std::to_string(n_scalars) +
                         " in the QMOM subsection. Set it to 0."));
  AssertThrow(
    !time_stepping_is_sdirk(simulation_parameters.simulation_control.method),
    ExcMessage("The SDIRK scheme is not yet supported for this physics"));
  AssertThrow(qmom_parameters.splitting_scheme ==
                Parameters::QMOMSplittingScheme::lie,
              ExcMessage("The strang splitting scheme of the QMOM physics is "
                         "not implemented yet. Set \"splitting scheme\" to lie "
                         "in the QMOM subsection."));
  AssertThrow(!qmom_parameters.enable_nucleation &&
                !qmom_parameters.enable_coalescence &&
                !qmom_parameters.enable_growth,
              ExcMessage("The evolution of the moments in phase space is not "
                         "implemented yet. Set \"enable nucleation\", \"enable "
                         "coalescence\" and \"enable growth\" to false in the "
                         "QMOM subsection."));
  AssertThrow(simulation_parameters.source_term.qmom_source.distribution_type ==
                Parameters::QMOMDistributionType::none,
              ExcMessage("The source terms of the QMOM physics are not "
                         "implemented yet. Set the distribution type of the "
                         "QMOM subsection of the source term subsection to "
                         "none."));
  AssertThrow(
    simulation_parameters.linear_solver.at(PhysicsID::qmom).preconditioner ==
      Parameters::LinearSolver::PreconditionerType::ilu,
    ExcMessage("The QMOM physics only supports the ilu "
               "preconditioner. Set \"preconditioner\" to ilu in the "
               "qmom subsection of the linear solver subsection."));

  // The QMOM boundary conditions are imposed weakly by the upwind fluxes
  for (const auto &[boundary_id, type] :
       simulation_parameters.boundary_conditions_qmom.type)
    AssertThrow(type == BoundaryConditions::BoundaryType::none ||
                  type == BoundaryConditions::BoundaryType::outlet ||
                  type == BoundaryConditions::BoundaryType::qmom_dirichlet,
                ExcMessage("The QMOM boundary condition of boundary id " +
                           std::to_string(boundary_id) +
                           " is not supported. The supported types are none, "
                           "outlet and dirichlet."));

  this->pcout << std::setprecision(simulation_control->get_log_precision())
              << std::scientific;

  // Names of the fields, used for the output and the error table
  field_names = qmom_parameters.scalar_names;
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    field_names.emplace_back("moment_" + std::to_string(i_moment));

  // All the fields share the same scalar discontinuous element
  const unsigned int qmom_degree =
    simulation_parameters.fem_parameters.qmom_degree;
  if (simulation_parameters.mesh.simplex)
    {
      // Simplex meshes (triangles and tetrahedra)
      fe              = std::make_unique<FE_SimplexDGP<dim>>(qmom_degree);
      mapping         = std::make_shared<MappingFE<dim>>(FE_SimplexP<dim>(1));
      cell_quadrature = std::make_unique<QGaussSimplex<dim>>(qmom_degree + 1);
      face_quadrature =
        std::make_unique<QGaussSimplex<dim - 1>>(qmom_degree + 1);
    }
  else
    {
      // Quadrilateral and hexahedral meshes. A degree of 0 still requires a
      // linear mapping.
      fe      = std::make_unique<FE_DGQ<dim>>(qmom_degree);
      mapping = std::make_shared<MappingQ<dim>>(std::max(1U, qmom_degree));
      cell_quadrature = std::make_unique<QGauss<dim>>(qmom_degree + 1);
      face_quadrature = std::make_unique<QGauss<dim - 1>>(qmom_degree + 1);
    }

  // Functions of all the fields. Their constructors verify the number of
  // expressions given for the scalars.
  initial_conditions = std::make_shared<QMOMFieldFunction<dim>>(
    simulation_parameters.initial_condition->qmom,
    n_scalars,
    n_moments,
    "initial conditions");
  for (const auto &[boundary_id, functions] :
       simulation_parameters.boundary_conditions_qmom.dirichlet)
    dirichlet_functions[boundary_id] = std::make_shared<QMOMFieldFunction<dim>>(
      functions,
      n_scalars,
      n_moments,
      "Dirichlet boundary condition of boundary id " +
        std::to_string(boundary_id));
  analytical_solution = std::make_shared<QMOMFieldFunction<dim>>(
    simulation_parameters.analytical_solution->qmom,
    n_scalars,
    n_moments,
    "analytical solution");

  present_solutions.resize(n_fields);
  previous_solutions.assign(n_fields,
                            std::vector<GlobalVectorType>(
                              maximum_number_of_previous_solutions()));
  system_rhs.resize(n_fields);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
QMOM<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;
  solution_output_structs.reserve(n_fields);

  for (unsigned int field = 0; field < n_fields; ++field)
    solution_output_structs.emplace_back(
      std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
      *dof_handler,
      present_solutions[field],
      std::vector<std::string>{field_names[field]},
      std::vector<DataComponentInterpretation::DataComponentInterpretation>{
        DataComponentInterpretation::component_is_scalar});

  return solution_output_structs;
}

template <int dim>
void
QMOM<dim>::finish_simulation()
{
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(triangulation->get_mpi_communicator());

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->calculate_error() &&
      simulation_parameters.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose)
    {
      for (const std::string &field_name : field_names)
        {
          error_table.set_scientific("error_" + field_name, true);
          error_table.set_precision("error_" + field_name,
                                    simulation_control->get_log_precision());
        }
      error_table.write_text(std::cout);
    }

  if (simulation_parameters.timer.type == Parameters::Timer::Type::end)
    {
      announce_string(this->pcout, "QMOM");
      this->pcout << std::defaultfloat;
      computing_timer.print_summary();
      this->pcout << std::scientific;
    }
}

template <int dim>
void
QMOM<dim>::percolate_time_vectors()
{
  for (unsigned int field = 0; field < n_fields; ++field)
    {
      std::vector<GlobalVectorType> &field_previous_solutions =
        previous_solutions[field];
      for (unsigned int i_previous_solution =
             field_previous_solutions.size() - 1;
           i_previous_solution > 0;
           --i_previous_solution)
        field_previous_solutions[i_previous_solution] =
          field_previous_solutions[i_previous_solution - 1];
      field_previous_solutions[0] = present_solutions[field];
    }
}

template <int dim>
void
QMOM<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() &&
      !first_iteration)
    {
      TimerOutput::Scope t(computing_timer, "Calculate L2 errors");

      const std::vector<double> errors = calculate_L2_errors();

      error_table.add_value("time", simulation_control->get_current_time());
      for (unsigned int field = 0; field < n_fields; ++field)
        error_table.add_value("error_" + field_names[field], errors[field]);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        for (unsigned int field = 0; field < n_fields; ++field)
          this->pcout << "L2 error " << field_names[field] << ": "
                      << errors[field] << std::endl;
    }
}

template <int dim>
void
QMOM<dim>::output_per_iteration_timer()
{
  announce_string(this->pcout, "QMOM");
  this->pcout << std::defaultfloat;
  computing_timer.print_summary();
  this->pcout << std::scientific;
  computing_timer.reset();
}

template <int dim>
std::vector<GlobalVectorType *>
QMOM<dim>::gather_all_solutions()
{
  std::vector<GlobalVectorType *> all_solutions;
  all_solutions.reserve(n_fields *
                        (1 + maximum_number_of_previous_solutions()));

  for (GlobalVectorType &present_solution : present_solutions)
    all_solutions.push_back(&present_solution);
  for (std::vector<GlobalVectorType> &field_previous_solutions :
       previous_solutions)
    for (GlobalVectorType &previous_solution : field_previous_solutions)
      all_solutions.push_back(&previous_solution);

  return all_solutions;
}

template <int dim>
template <typename TransferFunction>
void
QMOM<dim>::transfer_all_solutions(const TransferFunction &transfer)
{
  const std::vector<GlobalVectorType *> all_solutions = gather_all_solutions();

  // The solution transfer fills locally owned vectors, which are then copied
  // to the ghosted solution vectors
  std::vector<GlobalVectorType> transferred_solutions(
    all_solutions.size(),
    GlobalVectorType(locally_owned_dofs,
                     triangulation->get_mpi_communicator()));
  std::vector<GlobalVectorType *> transferred_solution_pointers;
  transferred_solution_pointers.reserve(transferred_solutions.size());
  for (GlobalVectorType &transferred_solution : transferred_solutions)
    transferred_solution_pointers.push_back(&transferred_solution);

  transfer(transferred_solution_pointers);

  for (unsigned int i_solution = 0; i_solution < all_solutions.size();
       ++i_solution)
    *all_solutions[i_solution] = transferred_solutions[i_solution];
}

template <int dim>
void
QMOM<dim>::pre_mesh_adaptation()
{
  const std::vector<GlobalVectorType *> all_solutions = gather_all_solutions();
  const std::vector<const GlobalVectorType *> solutions_to_transfer(
    all_solutions.begin(), all_solutions.end());

  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);
  solution_transfer->prepare_for_coarsening_and_refinement(
    solutions_to_transfer);
}

template <int dim>
void
QMOM<dim>::post_mesh_adaptation()
{
  transfer_all_solutions(
    [&](std::vector<GlobalVectorType *> &transferred_solutions) {
      solution_transfer->interpolate(transferred_solutions);
    });
}

template <int dim>
void
QMOM<dim>::compute_error_estimate(
  [[maybe_unused]] const std::pair<const Variable,
                                   Parameters::MultipleAdaptationParameters>
                                         &ivar,
  [[maybe_unused]] dealii::Vector<float> &estimated_error_per_cell)
{
  AssertThrow(false,
              ExcMessage("Mesh adaptation based on QMOM variables is not "
                         "implemented yet. Remove the QMOM variables from the "
                         "mesh adaptation subsection."));
}

template <int dim>
void
QMOM<dim>::write_checkpoint()
{
  const std::vector<GlobalVectorType *> all_solutions = gather_all_solutions();
  const std::vector<const GlobalVectorType *> solutions_to_serialize(
    all_solutions.begin(), all_solutions.end());

  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);
  solution_transfer->prepare_for_serialization(solutions_to_serialize);
}

template <int dim>
void
QMOM<dim>::read_checkpoint()
{
  this->pcout << "Reading QMOM checkpoint" << std::endl;

  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);
  transfer_all_solutions(
    [&](std::vector<GlobalVectorType *> &transferred_solutions) {
      solution_transfer->deserialize(transferred_solutions);
    });
}

template <int dim>
std::vector<OutputStructTableHandler>
QMOM<dim>::gather_tables()
{
  return std::vector<OutputStructTableHandler>();
}

template <int dim>
void
QMOM<dim>::define_constraints()
{
  // The discontinuous Galerkin discretization has no hanging-node constraints
  // and imposes the boundary conditions weakly
  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  constraints.close();
}

template <int dim>
void
QMOM<dim>::setup_dofs()
{
  dof_handler->distribute_dofs(*fe);

  const MPI_Comm mpi_communicator = triangulation->get_mpi_communicator();

  locally_owned_dofs    = dof_handler->locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(*dof_handler);

  for (unsigned int field = 0; field < n_fields; ++field)
    {
      present_solutions[field].reinit(locally_owned_dofs,
                                      locally_relevant_dofs,
                                      mpi_communicator);
      for (GlobalVectorType &previous_solution : previous_solutions[field])
        previous_solution.reinit(locally_owned_dofs,
                                 locally_relevant_dofs,
                                 mpi_communicator);
      system_rhs[field].reinit(locally_owned_dofs, mpi_communicator);
    }

  define_constraints();

  // The upwind fluxes couple the degrees of freedom of neighboring cells
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_flux_sparsity_pattern(*dof_handler,
                                       dsp,
                                       constraints,
                                       /*keep_constrained_dofs = */ true);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of QMOM degrees of freedom: "
              << n_fields * dof_handler->n_dofs() << " (" << n_fields
              << " fields of " << dof_handler->n_dofs() << ")" << std::endl;

  // Only the DoFHandler and the mapping are shared with the other physics:
  // the fields are stored in one vector each, and no other physics uses them
  multiphysics->set_dof_handler(PhysicsID::qmom, dof_handler);
  multiphysics->set_mapping(PhysicsID::qmom, mapping);
}

template <int dim>
void
QMOM<dim>::set_initial_conditions()
{
  initial_conditions->set_time(simulation_control->get_current_time());

  // The fields are interpolated at the support points of the scalar element.
  // All the fields are evaluated at once at each support point, since the
  // moments are computed together from the particle size distribution.
  const Quadrature<dim> support_points_quadrature(
    fe->get_unit_support_points());
  FEValues<dim>                        fe_values(*mapping,
                          *fe,
                          support_points_quadrature,
                          update_quadrature_points);
  const unsigned int                   n_dofs = fe->n_dofs_per_cell();
  std::vector<types::global_dof_index> dof_indices(n_dofs);
  Vector<double>                       field_values(n_fields);

  std::vector<GlobalVectorType> interpolated_solutions(
    n_fields,
    GlobalVectorType(locally_owned_dofs,
                     triangulation->get_mpi_communicator()));

  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          initial_conditions->vector_value(fe_values.quadrature_point(i),
                                           field_values);
          for (unsigned int field = 0; field < n_fields; ++field)
            interpolated_solutions[field](dof_indices[i]) = field_values[field];
        }
    }

  for (unsigned int field = 0; field < n_fields; ++field)
    {
      interpolated_solutions[field].compress(VectorOperation::insert);
      present_solutions[field] = interpolated_solutions[field];
    }

  percolate_time_vectors();
}

template <int dim>
void
QMOM<dim>::update_boundary_conditions()
{
  if (!simulation_parameters.boundary_conditions_qmom.time_dependent)
    return;

  const double time = simulation_control->get_current_time();
  for (const auto &[boundary_id, function] : dirichlet_functions)
    function->set_time(time);
}

template <int dim>
BoundaryConditions::BoundaryType
QMOM<dim>::get_boundary_type(const types::boundary_id boundary_id) const
{
  const auto &boundary_types =
    simulation_parameters.boundary_conditions_qmom.type;
  const auto type_iterator = boundary_types.find(boundary_id);
  return type_iterator == boundary_types.end() ?
           BoundaryConditions::BoundaryType::none :
           type_iterator->second;
}

template <int dim>
void
QMOM<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(computing_timer, "Assemble matrix");

  if (multiphysics->fluid_dynamics_is_block())
    assemble_transport_matrix(
      multiphysics->get_block_solution(PhysicsID::fluid_dynamics));
  else
    assemble_transport_matrix(
      multiphysics->get_solution(PhysicsID::fluid_dynamics));
}

template <int dim>
template <typename FluidVectorType>
void
QMOM<dim>::assemble_transport_matrix(const FluidVectorType &fluid_solution)
{
  system_matrix = 0.;

  const DoFHandler<dim> &fluid_dof_handler =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const auto get_fluid_cell =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell) {
      return typename DoFHandler<dim>::active_cell_iterator(&(*triangulation),
                                                            cell->level(),
                                                            cell->index(),
                                                            &fluid_dof_handler);
    };

  // Coefficient of the present solution in the BDF time derivative, which
  // includes the inverse of the time step. It is zero for steady simulations.
  const bool transient =
    time_stepping_is_bdf(simulation_control->get_assembly_method());
  const double present_bdf_coefficient =
    transient ? simulation_control->get_bdf_coefficients()[0] : 0.;

  QMOMScratchData<dim> scratch_data(*fe,
                                    *cell_quadrature,
                                    *face_quadrature,
                                    *mapping,
                                    fluid_dof_handler.get_fe());
  QMOMCopyData         copy_data(fe->n_dofs_per_cell());

  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        QMOMScratchData<dim>                                 &scratch,
        QMOMCopyData                                         &copy) {
      // The copy data is reset for every cell, including ghost cells, whose
      // faces are assembled when they are finer than an owned neighbor
      copy.local_matrix = 0.;
      copy.face_data.clear();
      copy.cell_is_local = cell->is_locally_owned();
      if (!copy.cell_is_local)
        return;

      cell->get_dof_indices(copy.local_dof_indices);
      scratch.reinit(cell);
      scratch.reinit_velocity(get_fluid_cell(cell), fluid_solution);

      // Weak form of the time derivative and of the advection of a field c:
      // int alpha_0 phi_i c - int (u . grad(phi_i)) c, where alpha_0 is the
      // coefficient of the present solution in the BDF scheme
      for (unsigned int q = 0; q < scratch.n_q_points; ++q)
        {
          const Tensor<1, dim> &velocity = scratch.velocity_values[q];
          for (unsigned int i = 0; i < scratch.n_dofs; ++i)
            {
              const double test_function_term =
                (present_bdf_coefficient * scratch.phi[i][q] -
                 velocity * scratch.grad_phi[i][q]) *
                scratch.JxW[q];
              for (unsigned int j = 0; j < scratch.n_dofs; ++j)
                copy.local_matrix(i, j) +=
                  test_function_term * scratch.phi[j][q];
            }
        }
    };

  const auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        QMOMScratchData<dim>                                 &scratch,
        QMOMCopyData                                         &copy) {
      // No flux crosses a boundary without QMOM boundary condition
      if (get_boundary_type(cell->face(face_no)->boundary_id()) ==
          BoundaryConditions::BoundaryType::none)
        return;

      scratch.reinit_boundary_face(cell, face_no);
      scratch.reinit_face_velocity(get_fluid_cell(cell),
                                   face_no,
                                   numbers::invalid_unsigned_int,
                                   fluid_solution);

      // Outflow: the interior values leave the domain through outlet and
      // Dirichlet boundaries. The inflow through Dirichlet boundaries is
      // known and assembled in the right-hand side.
      const FEFaceValuesBase<dim> &fe_face_values =
        scratch.fe_interface_values.get_fe_face_values(0);
      for (unsigned int q = 0; q < scratch.n_face_q_points; ++q)
        {
          const double normal_velocity =
            scratch.face_velocity_values[q] * fe_face_values.normal_vector(q);
          if (normal_velocity <= 0.)
            continue;

          for (unsigned int i = 0; i < scratch.n_dofs; ++i)
            {
              const double test_function_term =
                normal_velocity * fe_face_values.shape_value(i, q) *
                fe_face_values.JxW(q);
              for (unsigned int j = 0; j < scratch.n_dofs; ++j)
                copy.local_matrix(i, j) +=
                  test_function_term * fe_face_values.shape_value(j, q);
            }
        }
    };

  const auto face_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        const unsigned int                                    sub_face_no,
        const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
        const unsigned int                                    neighbor_face_no,
        const unsigned int    neighbor_sub_face_no,
        QMOMScratchData<dim> &scratch,
        QMOMCopyData         &copy) {
      scratch.reinit_internal_face(cell,
                                   face_no,
                                   sub_face_no,
                                   neighbor_cell,
                                   neighbor_face_no,
                                   neighbor_sub_face_no);
      scratch.reinit_face_velocity(get_fluid_cell(cell),
                                   face_no,
                                   sub_face_no,
                                   fluid_solution);

      const FEInterfaceValues<dim> &fe_interface_values =
        scratch.fe_interface_values;
      const unsigned int n_interface_dofs =
        fe_interface_values.n_current_interface_dofs();

      QMOMCopyData::FaceData &face_data = copy.face_data.emplace_back();
      face_data.joint_dof_indices =
        fe_interface_values.get_interface_dof_indices();
      face_data.face_matrix.reinit(n_interface_dofs, n_interface_dofs);

      // Upwind flux (u . n) c_upwind, with n pointing from the first to the
      // second cell. The jump of the test function, phi_first - phi_second,
      // gives the contributions of opposite signs to the two cells.
      for (unsigned int q = 0; q < scratch.n_face_q_points; ++q)
        {
          const double normal_velocity = scratch.face_velocity_values[q] *
                                         fe_interface_values.normal_vector(q);
          const bool upwind_is_first_cell = normal_velocity > 0.;
          for (unsigned int i = 0; i < n_interface_dofs; ++i)
            {
              const double test_function_term =
                normal_velocity *
                fe_interface_values.jump_in_shape_values(i, q) *
                fe_interface_values.JxW(q);
              for (unsigned int j = 0; j < n_interface_dofs; ++j)
                face_data.face_matrix(i, j) +=
                  test_function_term *
                  fe_interface_values.shape_value(upwind_is_first_cell, j, q);
            }
        }
    };

  // A face between a locally owned and a ghost cell is assembled by the two
  // processes that own its cells. Each process only adds the rows that it
  // owns, so that every row receives the flux of each of its faces exactly
  // once.
  std::vector<types::global_dof_index> row_columns;
  std::vector<double>                  row_values;
  const auto                           copier = [&](const QMOMCopyData &copy) {
    if (copy.cell_is_local)
      constraints.distribute_local_to_global(copy.local_matrix,
                                             copy.local_dof_indices,
                                             system_matrix);

    for (const QMOMCopyData::FaceData &face_data : copy.face_data)
      {
        const unsigned int n_interface_dofs =
          face_data.joint_dof_indices.size();
        for (unsigned int i = 0; i < n_interface_dofs; ++i)
          {
            const types::global_dof_index row = face_data.joint_dof_indices[i];
            if (!locally_owned_dofs.is_element(row))
              continue;

            row_columns.clear();
            row_values.clear();
            for (unsigned int j = 0; j < n_interface_dofs; ++j)
              {
                row_columns.push_back(face_data.joint_dof_indices[j]);
                row_values.push_back(face_data.face_matrix(i, j));
              }
            system_matrix.add(row, row_columns, row_values);
          }
      }
  };

  MeshWorker::mesh_loop(dof_handler->begin_active(),
                        dof_handler->end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_ghost_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_ghost_faces_both,
                        boundary_worker,
                        face_worker);

  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
QMOM<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(computing_timer, "Assemble RHS");

  if (multiphysics->fluid_dynamics_is_block())
    assemble_transport_rhs(
      multiphysics->get_block_solution(PhysicsID::fluid_dynamics));
  else
    assemble_transport_rhs(
      multiphysics->get_solution(PhysicsID::fluid_dynamics));
}

template <int dim>
template <typename FluidVectorType>
void
QMOM<dim>::assemble_transport_rhs(const FluidVectorType &fluid_solution)
{
  for (GlobalVectorType &field_rhs : system_rhs)
    field_rhs = 0.;

  const DoFHandler<dim> &fluid_dof_handler =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  // Number of previous solutions of the BDF time derivative, zero for steady
  // simulations
  const auto         method = simulation_control->get_assembly_method();
  const unsigned int n_previous_solutions =
    time_stepping_is_bdf(method) ? number_of_previous_solutions(method) : 0;

  QMOMScratchData<dim> scratch_data(*fe,
                                    *cell_quadrature,
                                    *face_quadrature,
                                    *mapping,
                                    fluid_dof_handler.get_fe());
  const unsigned int   n_dofs = scratch_data.n_dofs;

  // The mass matrix and the shape functions are the same for all the fields:
  // they are evaluated once per cell and face, then applied to every field
  FullMatrix<double>          mass_matrix(n_dofs, n_dofs);
  Vector<double>              previous_values(n_dofs);
  Vector<double>              previous_terms(n_dofs);
  std::vector<Vector<double>> local_rhs(n_fields, Vector<double>(n_dofs));
  std::vector<types::global_dof_index> dof_indices(n_dofs);
  Vector<double>                       dirichlet_values(n_fields);

  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      cell->get_dof_indices(dof_indices);
      for (Vector<double> &field_rhs : local_rhs)
        field_rhs = 0.;

      // Previous solutions of the BDF time derivative:
      // - int phi_i sum_{p >= 1} alpha_p c_p
      if (n_previous_solutions > 0)
        {
          const Vector<double> &bdf_coefficients =
            simulation_control->get_bdf_coefficients();

          scratch_data.reinit(cell);
          mass_matrix = 0.;
          for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
            for (unsigned int i = 0; i < n_dofs; ++i)
              {
                const double phi_i_JxW =
                  scratch_data.phi[i][q] * scratch_data.JxW[q];
                for (unsigned int j = 0; j < n_dofs; ++j)
                  mass_matrix(i, j) += phi_i_JxW * scratch_data.phi[j][q];
              }

          for (unsigned int field = 0; field < n_fields; ++field)
            {
              previous_terms = 0.;
              for (unsigned int p = 0; p < n_previous_solutions; ++p)
                {
                  cell->get_dof_values(previous_solutions[field][p],
                                       previous_values);
                  previous_terms.add(-bdf_coefficients[p + 1], previous_values);
                }
              mass_matrix.vmult_add(local_rhs[field], previous_terms);
            }
        }

      // Inflow through the Dirichlet boundaries: - int phi_i (u . n) c_D
      for (const unsigned int face_no : cell->face_indices())
        {
          if (!cell->at_boundary(face_no))
            continue;
          const types::boundary_id boundary_id =
            cell->face(face_no)->boundary_id();
          if (get_boundary_type(boundary_id) !=
              BoundaryConditions::BoundaryType::qmom_dirichlet)
            continue;

          scratch_data.reinit_boundary_face(cell, face_no);
          scratch_data.reinit_face_velocity(
            typename DoFHandler<dim>::active_cell_iterator(&(*triangulation),
                                                           cell->level(),
                                                           cell->index(),
                                                           &fluid_dof_handler),
            face_no,
            numbers::invalid_unsigned_int,
            fluid_solution);

          const FEFaceValuesBase<dim> &fe_face_values =
            scratch_data.fe_interface_values.get_fe_face_values(0);
          const QMOMFieldFunction<dim> &dirichlet_function =
            *dirichlet_functions.at(boundary_id);
          for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
            {
              const double normal_velocity =
                scratch_data.face_velocity_values[q] *
                fe_face_values.normal_vector(q);
              if (normal_velocity >= 0.)
                continue;

              dirichlet_function.vector_value(
                fe_face_values.quadrature_point(q), dirichlet_values);
              for (unsigned int i = 0; i < n_dofs; ++i)
                {
                  const double test_function_term =
                    -normal_velocity * fe_face_values.shape_value(i, q) *
                    fe_face_values.JxW(q);
                  for (unsigned int field = 0; field < n_fields; ++field)
                    local_rhs[field](i) +=
                      test_function_term * dirichlet_values[field];
                }
            }
        }

      for (unsigned int field = 0; field < n_fields; ++field)
        constraints.distribute_local_to_global(local_rhs[field],
                                               dof_indices,
                                               system_rhs[field]);
    }

  for (GlobalVectorType &field_rhs : system_rhs)
    field_rhs.compress(VectorOperation::add);
}

template <int dim>
void
QMOM<dim>::setup_preconditioner()
{
  TimerOutput::Scope t(computing_timer, "Setup ILU");

  const Parameters::LinearSolver &linear_solver_parameters =
    simulation_parameters.linear_solver.at(PhysicsID::qmom);
  const TrilinosWrappers::PreconditionILU::AdditionalData
    preconditioner_options(linear_solver_parameters.ilu_precond_fill,
                           linear_solver_parameters.ilu_precond_atol,
                           linear_solver_parameters.ilu_precond_rtol,
                           0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();
  ilu_preconditioner->initialize(system_matrix, preconditioner_options);
}

template <int dim>
void
QMOM<dim>::solve_linear_system()
{
  switch (simulation_parameters.multiphysics.qmom_parameters.splitting_scheme)
    {
      case Parameters::QMOMSplittingScheme::lie:
        advance_physical_space();
        advance_phase_space(simulation_control->get_time_step());
        break;
      case Parameters::QMOMSplittingScheme::strang:
        AssertThrow(false,
                    ExcMessage("The strang splitting scheme of the QMOM "
                               "physics is not implemented yet."));
        break;
    }
}

template <int dim>
void
QMOM<dim>::advance_physical_space()
{
  TimerOutput::Scope t(computing_timer, "Solve linear system");

  const Parameters::LinearSolver &linear_solver_parameters =
    simulation_parameters.linear_solver.at(PhysicsID::qmom);
  const bool verbose =
    linear_solver_parameters.verbosity != Parameters::Verbosity::quiet;
  const double rescale_metric = this->get_residual_rescale_metric();

  GlobalVectorType completely_distributed_solution(
    locally_owned_dofs, triangulation->get_mpi_communicator());

  // All the fields share the matrix and the preconditioner: only the
  // right-hand side changes
  for (unsigned int field = 0; field < n_fields; ++field)
    {
      // The fields have very different magnitudes, so the tolerance is
      // relative to the right-hand side of each field
      const double current_residual =
        system_rhs[field].l2_norm() / rescale_metric;
      const double linear_solver_tolerance =
        std::max(linear_solver_parameters.relative_residual * current_residual,
                 linear_solver_parameters.minimum_residual);
      if (verbose)
        this->pcout << "  -Tolerance of iterative solver of "
                    << field_names[field] << " is : " << linear_solver_tolerance
                    << std::endl;

      SolverControl solver_control(linear_solver_parameters.max_iterations,
                                   linear_solver_tolerance * rescale_metric,
                                   true,
                                   true);
      const TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
        false, linear_solver_parameters.max_krylov_vectors);
      TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

      // The solution at the beginning of the time step is the initial guess
      completely_distributed_solution = present_solutions[field];
      solver.solve(system_matrix,
                   completely_distributed_solution,
                   system_rhs[field],
                   *ilu_preconditioner);

      if (verbose)
        this->pcout << "  -Iterative solver of " << field_names[field]
                    << " took : " << solver_control.last_step()
                    << " steps to reach a residual norm of "
                    << solver_control.last_value() / rescale_metric
                    << std::endl;

      constraints.distribute(completely_distributed_solution);
      present_solutions[field] = completely_distributed_solution;
    }
}

template <int dim>
void
QMOM<dim>::advance_phase_space([[maybe_unused]] const double time_step)
{
  // The evolution of the moments in phase space is not implemented yet: the
  // constructor rejects nucleation, growth and coalescence, so the moments
  // are constant in phase space.
}

template <int dim>
std::vector<double>
QMOM<dim>::calculate_L2_errors()
{
  analytical_solution->set_time(simulation_control->get_current_time());

  FEValues<dim>                    fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int               n_q_points = cell_quadrature->size();
  std::vector<std::vector<double>> values(n_fields,
                                          std::vector<double>(n_q_points));
  Vector<double>                   exact_values(n_fields);
  std::vector<double>              squared_errors(n_fields, 0.);

  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      for (unsigned int field = 0; field < n_fields; ++field)
        fe_values.get_function_values(present_solutions[field], values[field]);

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          analytical_solution->vector_value(fe_values.quadrature_point(q),
                                            exact_values);
          for (unsigned int field = 0; field < n_fields; ++field)
            squared_errors[field] += Utilities::fixed_power<2>(
                                       values[field][q] - exact_values[field]) *
                                     fe_values.JxW(q);
        }
    }

  std::vector<double> errors(n_fields);
  Utilities::MPI::sum(squared_errors,
                      triangulation->get_mpi_communicator(),
                      errors);
  for (double &error : errors)
    error = std::sqrt(error);
  return errors;
}

template <int dim>
double
QMOM<dim>::get_residual_rescale_metric() const
{
  return simulation_parameters.linear_solver.at(PhysicsID::qmom)
             .rescale_residual_by_volume ?
           std::sqrt(GridTools::volume(*triangulation, *mapping)) :
           1.;
}

template class QMOM<2>;
template class QMOM<3>;
