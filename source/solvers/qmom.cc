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
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
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
  , n_components(n_scalars + n_moments)
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

  // Names of the components, used for the output and the error table
  component_names = qmom_parameters.scalar_names;
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    component_names.emplace_back("moment_" + std::to_string(i_moment));

  // All the components share the same scalar discontinuous element
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

  // Functions of all the components. The scalars are given by one expression
  // each, while the moments are computed from the particle size distribution,
  // which always gives n_moments moments: the constructors therefore only
  // verify the number of expressions given for the scalars.
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
  source_term = std::make_shared<QMOMFieldFunction<dim>>(
    simulation_parameters.source_term.qmom_source,
    n_scalars,
    n_moments,
    "source term");
  analytical_solution = std::make_shared<QMOMFieldFunction<dim>>(
    simulation_parameters.analytical_solution->qmom,
    n_scalars,
    n_moments,
    "analytical solution");

  present_solutions.resize(n_components);
  previous_solutions.assign(n_components,
                            std::vector<GlobalVectorType>(
                              maximum_number_of_previous_solutions()));
  system_rhs.resize(n_components);

  // Solution transfer, for mesh adaptation and restarts
  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
QMOM<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;
  solution_output_structs.reserve(n_components);

  // Each component is stored in its own vector, which is output as a scalar
  // named after the component
  for (unsigned int component = 0; component < n_components; ++component)
    solution_output_structs.emplace_back(
      std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
      *dof_handler,
      present_solutions[component],
      std::vector<std::string>{component_names[component]},
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
      for (const std::string &component_name : component_names)
        {
          error_table.set_scientific("error_" + component_name, true);
          error_table.set_precision("error_" + component_name,
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
  for (unsigned int component = 0; component < n_components; ++component)
    {
      std::vector<GlobalVectorType> &component_previous_solutions =
        previous_solutions[component];
      for (unsigned int i_previous_solution =
             component_previous_solutions.size() - 1;
           i_previous_solution > 0;
           --i_previous_solution)
        component_previous_solutions[i_previous_solution] =
          component_previous_solutions[i_previous_solution - 1];
      component_previous_solutions[0] = present_solutions[component];
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
      for (unsigned int component = 0; component < n_components; ++component)
        error_table.add_value("error_" + component_names[component],
                              errors[component]);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        for (unsigned int component = 0; component < n_components; ++component)
          this->pcout << "L2 error " << component_names[component] << ": "
                      << errors[component] << std::endl;
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
void
QMOM<dim>::pre_mesh_adaptation()
{
  // The present and previous solutions of all the components are transferred
  // together: present solution of component 0, previous solutions of
  // component 0, present solution of component 1, etc.
  std::vector<const GlobalVectorType *> sol_set_transfer;
  for (unsigned int component = 0; component < n_components; ++component)
    {
      sol_set_transfer.emplace_back(&present_solutions[component]);
      for (const GlobalVectorType &previous_solution :
           previous_solutions[component])
        sol_set_transfer.emplace_back(&previous_solution);
    }
  solution_transfer->prepare_for_coarsening_and_refinement(sol_set_transfer);
}

template <int dim>
void
QMOM<dim>::post_mesh_adaptation()
{
  const MPI_Comm mpi_communicator = triangulation->get_mpi_communicator();

  // Set up the vectors for the transfer, in the order of pre_mesh_adaptation()
  std::vector<GlobalVectorType> distributed_solutions(
    n_components * (1 + maximum_number_of_previous_solutions()),
    GlobalVectorType(locally_owned_dofs, mpi_communicator));
  std::vector<GlobalVectorType *> output_vectors;
  for (GlobalVectorType &distributed_solution : distributed_solutions)
    output_vectors.emplace_back(&distributed_solution);

  // Interpolate the present and previous solutions
  solution_transfer->interpolate(output_vectors);

  // Fix on the new mesh
  unsigned int i_vector = 0;
  for (unsigned int component = 0; component < n_components; ++component)
    {
      present_solutions[component] = distributed_solutions[i_vector++];
      for (GlobalVectorType &previous_solution : previous_solutions[component])
        previous_solution = distributed_solutions[i_vector++];
    }
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
  // Same order as in pre_mesh_adaptation()
  std::vector<const GlobalVectorType *> sol_set_transfer;
  for (unsigned int component = 0; component < n_components; ++component)
    {
      sol_set_transfer.emplace_back(&present_solutions[component]);
      for (const GlobalVectorType &previous_solution :
           previous_solutions[component])
        sol_set_transfer.emplace_back(&previous_solution);
    }

  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);
  solution_transfer->prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
QMOM<dim>::read_checkpoint()
{
  const MPI_Comm mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading QMOM checkpoint" << std::endl;

  // Set up the vectors for the deserialization, in the order of
  // write_checkpoint()
  std::vector<GlobalVectorType> distributed_solutions(
    n_components * (1 + maximum_number_of_previous_solutions()),
    GlobalVectorType(locally_owned_dofs, mpi_communicator));
  std::vector<GlobalVectorType *> input_vectors;
  for (GlobalVectorType &distributed_solution : distributed_solutions)
    input_vectors.emplace_back(&distributed_solution);

  solution_transfer->deserialize(input_vectors);

  unsigned int i_vector = 0;
  for (unsigned int component = 0; component < n_components; ++component)
    {
      present_solutions[component] = distributed_solutions[i_vector++];
      for (GlobalVectorType &previous_solution : previous_solutions[component])
        previous_solution = distributed_solutions[i_vector++];
    }
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

  for (unsigned int component = 0; component < n_components; ++component)
    {
      present_solutions[component].reinit(locally_owned_dofs,
                                          locally_relevant_dofs,
                                          mpi_communicator);
      for (GlobalVectorType &previous_solution : previous_solutions[component])
        previous_solution.reinit(locally_owned_dofs,
                                 locally_relevant_dofs,
                                 mpi_communicator);
      system_rhs[component].reinit(locally_owned_dofs, mpi_communicator);
    }

  define_constraints();

  // Sparse matrix initialization. The upwind fluxes couple the degrees of
  // freedom of neighboring cells, hence the flux sparsity pattern.
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
              << n_components * dof_handler->n_dofs() << " (" << n_components
              << " components of " << dof_handler->n_dofs() << ")" << std::endl;

  // Provide the QMOM dof_handler and mapping to the multiphysics interface.
  // The solution vectors, one per component, are not provided.
  multiphysics->set_dof_handler(PhysicsID::qmom, dof_handler);
  multiphysics->set_mapping(PhysicsID::qmom, mapping);
}

template <int dim>
void
QMOM<dim>::set_initial_conditions()
{
  initial_conditions->set_time(simulation_control->get_current_time());

  GlobalVectorType interpolated_solution(locally_owned_dofs,
                                         triangulation->get_mpi_communicator());
  for (unsigned int component = 0; component < n_components; ++component)
    {
      // Each component is stored in its own vector, with the scalar element:
      // the component of the initial conditions is interpolated as a scalar
      // function
      const ScalarFunctionFromFunctionObject<dim> component_initial_condition(
        [&](const Point<dim> &p) {
          return initial_conditions->value(p, component);
        });
      VectorTools::interpolate(*mapping,
                               *dof_handler,
                               component_initial_condition,
                               interpolated_solution);
      constraints.distribute(interpolated_solution);
      present_solutions[component] = interpolated_solution;
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
void
QMOM<dim>::reinit_velocity_with_adequate_solution(
  const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
  QMOMScratchData<dim>                                 &scratch_data)
{
  if (multiphysics->fluid_dynamics_is_block())
    scratch_data.reinit_velocity(velocity_cell,
                                 multiphysics->get_block_solution(
                                   PhysicsID::fluid_dynamics));
  else
    scratch_data.reinit_velocity(
      velocity_cell, multiphysics->get_solution(PhysicsID::fluid_dynamics));
}

template <int dim>
void
QMOM<dim>::reinit_face_velocity_with_adequate_solution(
  const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
  const unsigned int                                    face_no,
  const unsigned int                                    sub_face_no,
  QMOMScratchData<dim>                                 &scratch_data)
{
  if (multiphysics->fluid_dynamics_is_block())
    scratch_data.reinit_face_velocity(velocity_cell,
                                      face_no,
                                      sub_face_no,
                                      multiphysics->get_block_solution(
                                        PhysicsID::fluid_dynamics));
  else
    scratch_data.reinit_face_velocity(velocity_cell,
                                      face_no,
                                      sub_face_no,
                                      multiphysics->get_solution(
                                        PhysicsID::fluid_dynamics));
}

template <int dim>
void
QMOM<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(computing_timer, "Assemble matrix");

  system_matrix = 0.;

  const DoFHandler<dim> &dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  // Get the cell of the fluid dynamics that corresponds to a QMOM cell
  const auto get_velocity_cell =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell) {
      return typename DoFHandler<dim>::active_cell_iterator(&(*triangulation),
                                                            cell->level(),
                                                            cell->index(),
                                                            &dof_handler_fluid);
    };

  // Coefficient of the present solution in the BDF time derivative, which
  // includes the inverse of the time step. It is zero for steady simulations.
  const bool transient =
    time_stepping_is_bdf(simulation_control->get_assembly_method());
  const double bdf_coefficient_present =
    transient ? simulation_control->get_bdf_coefficients()[0] : 0.;

  QMOMScratchData<dim> scratch_data(*fe,
                                    *cell_quadrature,
                                    *face_quadrature,
                                    *mapping,
                                    dof_handler_fluid.get_fe(),
                                    n_components);
  QMOMCopyData         copy_data(fe->n_dofs_per_cell(), n_components);

  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        QMOMScratchData<dim>                                 &scratch_data,
        QMOMCopyData                                         &copy_data) {
      // The copy data is reset for every cell, including ghost cells, whose
      // faces are assembled when they are finer than an owned neighbor
      copy_data.local_matrix = 0.;
      copy_data.face_data.clear();
      copy_data.cell_is_local = cell->is_locally_owned();
      if (!copy_data.cell_is_local)
        return;

      cell->get_dof_indices(copy_data.local_dof_indices);
      scratch_data.reinit(cell);
      reinit_velocity_with_adequate_solution(get_velocity_cell(cell),
                                             scratch_data);

      for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
        {
          const Tensor<1, dim> &velocity = scratch_data.velocity_values[q];
          const double          JxW      = scratch_data.JxW[q];

          for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
            {
              const double          phi_i      = scratch_data.phi[q][i];
              const Tensor<1, dim> &grad_phi_i = scratch_data.grad_phi[q][i];

              for (unsigned int j = 0; j < scratch_data.n_dofs; ++j)
                {
                  const double phi_j = scratch_data.phi[q][j];

                  // Time derivative with the BDF scheme:
                  // \phi_i * \alpha_0 * \phi_j
                  copy_data.local_matrix(i, j) +=
                    phi_i * bdf_coefficient_present * phi_j * JxW;

                  // Advection in conservative form, integrated by parts:
                  // - (u . \nabla \phi_i) * \phi_j
                  copy_data.local_matrix(i, j) -=
                    velocity * grad_phi_i * phi_j * JxW;
                }
            }
        }
    };

  const auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        QMOMScratchData<dim>                                 &scratch_data,
        QMOMCopyData                                         &copy_data) {
      // Only the outlet and Dirichlet boundaries contribute to the matrix: no
      // flux crosses a boundary of type none
      if (simulation_parameters.boundary_conditions_qmom.type.at(
            cell->face(face_no)->boundary_id()) ==
          BoundaryConditions::BoundaryType::none)
        return;

      scratch_data.reinit_boundary_face(cell, face_no);
      reinit_face_velocity_with_adequate_solution(get_velocity_cell(cell),
                                                  face_no,
                                                  numbers::invalid_unsigned_int,
                                                  scratch_data);

      const FEFaceValuesBase<dim> &fe_face =
        scratch_data.fe_interface_values.get_fe_face_values(0);

      for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
        {
          const double velocity_dot_n =
            scratch_data.face_velocity_values[q] * fe_face.normal_vector(q);

          // The advection term is only applied to the outflow part of the
          // boundary. The inflow through the Dirichlet boundaries is known and
          // assembled in the right-hand side.
          if (velocity_dot_n <= 0.)
            continue;

          const double advection_factor = velocity_dot_n * fe_face.JxW(q);

          for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
            {
              const double phi_i = fe_face.shape_value(i, q);
              for (unsigned int j = 0; j < scratch_data.n_dofs; ++j)
                {
                  // Outflow of the interior values:
                  // \phi_i * \phi_j * (u . n)
                  copy_data.local_matrix(i, j) +=
                    phi_i * fe_face.shape_value(j, q) * advection_factor;
                }
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
        QMOMScratchData<dim> &scratch_data,
        QMOMCopyData         &copy_data) {
      scratch_data.reinit_internal_face(cell,
                                        face_no,
                                        sub_face_no,
                                        neighbor_cell,
                                        neighbor_face_no,
                                        neighbor_sub_face_no);
      reinit_face_velocity_with_adequate_solution(get_velocity_cell(cell),
                                                  face_no,
                                                  sub_face_no,
                                                  scratch_data);

      const unsigned int n_interface_dofs = scratch_data.n_interface_dofs;

      // Pad copy_data memory for the internal faces elementary matrices
      QMOMCopyData::FaceData &copy_data_face =
        copy_data.face_data.emplace_back();
      copy_data_face.joint_dof_indices =
        scratch_data.fe_interface_values.get_interface_dof_indices();
      copy_data_face.face_matrix.reinit(n_interface_dofs, n_interface_dofs);

      const std::vector<double> &JxW =
        scratch_data.fe_interface_values.get_JxW_values();
      const std::vector<Tensor<1, dim>> &normals =
        scratch_data.fe_interface_values.get_normal_vectors();

      for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
        {
          // The normal points from the "here" cell to the "there" cell
          const double velocity_dot_n =
            scratch_data.face_velocity_values[q] * normals[q];

          // Shape functions of the upwind cell
          const auto &phi_upwind = (velocity_dot_n > 0.) ?
                                     scratch_data.face_phi_here[q] :
                                     scratch_data.face_phi_there[q];
          const auto &jump_phi   = scratch_data.face_jump_phi[q];

          const double advection_factor = velocity_dot_n * JxW[q];

          for (unsigned int i = 0; i < n_interface_dofs; ++i)
            {
              const double jump_phi_i = jump_phi[i];
              for (unsigned int j = 0; j < n_interface_dofs; ++j)
                {
                  // Advection term with upwinding:
                  // [\phi_i] * \phi_j^{upwind} * (u . n)
                  copy_data_face.face_matrix(i, j) +=
                    jump_phi_i * phi_upwind[j] * advection_factor;
                }
            }
        }
    };

  // A face between a locally owned and a ghost cell is assembled by the two
  // processes that own its cells. Each process only adds the rows that it
  // owns, so that every row receives the flux of each of its faces exactly
  // once.
  std::vector<types::global_dof_index> row_columns;
  std::vector<double>                  row_values;
  const auto copier = [&](const QMOMCopyData &copy_data) {
    if (copy_data.cell_is_local)
      constraints.distribute_local_to_global(copy_data.local_matrix,
                                             copy_data.local_dof_indices,
                                             system_matrix);

    for (const QMOMCopyData::FaceData &copy_data_face : copy_data.face_data)
      {
        const unsigned int n_interface_dofs =
          copy_data_face.joint_dof_indices.size();
        for (unsigned int i = 0; i < n_interface_dofs; ++i)
          {
            const types::global_dof_index row =
              copy_data_face.joint_dof_indices[i];
            if (!locally_owned_dofs.is_element(row))
              continue;

            row_columns.clear();
            row_values.clear();
            for (unsigned int j = 0; j < n_interface_dofs; ++j)
              {
                row_columns.push_back(copy_data_face.joint_dof_indices[j]);
                row_values.push_back(copy_data_face.face_matrix(i, j));
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

  for (GlobalVectorType &component_rhs : system_rhs)
    component_rhs = 0.;

  // Update the source term time
  source_term->set_time(simulation_control->get_current_time());

  const DoFHandler<dim> &dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  // Get the cell of the fluid dynamics that corresponds to a QMOM cell
  const auto get_velocity_cell =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell) {
      return typename DoFHandler<dim>::active_cell_iterator(&(*triangulation),
                                                            cell->level(),
                                                            cell->index(),
                                                            &dof_handler_fluid);
    };

  // Number of previous solutions of the BDF time derivative, zero for steady
  // simulations
  const auto         method = simulation_control->get_assembly_method();
  const unsigned int n_previous_solutions =
    time_stepping_is_bdf(method) ? number_of_previous_solutions(method) : 0;
  const Vector<double> &bdf_coefficients =
    simulation_control->get_bdf_coefficients();

  QMOMScratchData<dim> scratch_data(*fe,
                                    *cell_quadrature,
                                    *face_quadrature,
                                    *mapping,
                                    dof_handler_fluid.get_fe(),
                                    n_components);
  QMOMCopyData         copy_data(fe->n_dofs_per_cell(), n_components);

  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        QMOMScratchData<dim>                                 &scratch_data,
        QMOMCopyData                                         &copy_data) {
      for (Vector<double> &component_rhs : copy_data.local_rhs)
        component_rhs = 0.;
      copy_data.cell_is_local = cell->is_locally_owned();
      if (!copy_data.cell_is_local)
        return;

      cell->get_dof_indices(copy_data.local_dof_indices);
      scratch_data.reinit(cell);
      if (n_previous_solutions > 0)
        scratch_data.reinit_previous_solutions(previous_solutions);
      scratch_data.reinit_source(*source_term);

      for (unsigned int q = 0; q < scratch_data.n_q_points; ++q)
        {
          const double JxW = scratch_data.JxW[q];

          for (unsigned int component = 0; component < n_components;
               ++component)
            {
              // Previous solutions of the BDF scheme:
              // sum_{p >= 1} \alpha_p * \zeta_c^{n+1-p}
              double previous_terms = 0.;
              for (unsigned int p = 0; p < n_previous_solutions; ++p)
                previous_terms += bdf_coefficients[p + 1] *
                                  scratch_data.previous_values[component][p][q];

              const double source = scratch_data.source_values[q][component];

              Vector<double> &local_rhs = copy_data.local_rhs[component];
              for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                {
                  const double phi_i = scratch_data.phi[q][i];

                  // Time derivative with the BDF scheme, previous solutions:
                  // - \phi_i * sum_{p >= 1} \alpha_p * \zeta_c^{n+1-p}
                  local_rhs(i) -= phi_i * previous_terms * JxW;

                  // Source term: \phi_i * S_c
                  local_rhs(i) += phi_i * source * JxW;
                }
            }
        }
    };

  const auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        QMOMScratchData<dim>                                 &scratch_data,
        QMOMCopyData                                         &copy_data) {
      // Only the inflow through the Dirichlet boundaries contributes to the
      // right-hand side
      const types::boundary_id boundary_id = cell->face(face_no)->boundary_id();
      if (simulation_parameters.boundary_conditions_qmom.type.at(boundary_id) !=
          BoundaryConditions::BoundaryType::qmom_dirichlet)
        return;

      scratch_data.reinit_boundary_face(cell, face_no);
      reinit_face_velocity_with_adequate_solution(get_velocity_cell(cell),
                                                  face_no,
                                                  numbers::invalid_unsigned_int,
                                                  scratch_data);

      const FEFaceValuesBase<dim> &fe_face =
        scratch_data.fe_interface_values.get_fe_face_values(0);
      const QMOMFieldFunction<dim> &dirichlet_function =
        *dirichlet_functions.at(boundary_id);
      Vector<double> dirichlet_values(n_components);

      for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
        {
          const double velocity_dot_n =
            scratch_data.face_velocity_values[q] * fe_face.normal_vector(q);

          // The Dirichlet values only enter the domain through the inflow part
          // of the boundary. The outflow is assembled in the matrix.
          if (velocity_dot_n >= 0.)
            continue;

          // All the components are evaluated at once, since the moments are
          // computed together from the particle size distribution
          dirichlet_function.vector_value(fe_face.quadrature_point(q),
                                          dirichlet_values);

          const double advection_factor = velocity_dot_n * fe_face.JxW(q);

          for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
            {
              const double phi_i = fe_face.shape_value(i, q);
              for (unsigned int component = 0; component < n_components;
                   ++component)
                {
                  // Inflow of the Dirichlet values:
                  // - \phi_i * \zeta_{c,D} * (u . n)
                  copy_data.local_rhs[component](i) -=
                    phi_i * dirichlet_values[component] * advection_factor;
                }
            }
        }
    };

  const auto copier = [&](const QMOMCopyData &copy_data) {
    if (!copy_data.cell_is_local)
      return;

    for (unsigned int component = 0; component < n_components; ++component)
      constraints.distribute_local_to_global(copy_data.local_rhs[component],
                                             copy_data.local_dof_indices,
                                             system_rhs[component]);
  };

  // The right-hand side has no interior face terms: the upwind fluxes between
  // cells only involve the unknown solution, and are assembled in the matrix
  MeshWorker::mesh_loop(dof_handler->begin_active(),
                        dof_handler->end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces,
                        boundary_worker);

  for (GlobalVectorType &component_rhs : system_rhs)
    component_rhs.compress(VectorOperation::add);
}

template <int dim>
void
QMOM<dim>::setup_preconditioner()
{
  const auto preconditioner_type =
    simulation_parameters.linear_solver.at(PhysicsID::qmom).preconditioner;

  if (preconditioner_type == Parameters::LinearSolver::PreconditionerType::ilu)
    setup_ilu();
  else if (preconditioner_type ==
           Parameters::LinearSolver::PreconditionerType::amg)
    setup_amg();
  else
    AssertThrow(false,
                ExcMessage("This linear solver does not support this "
                           "preconditioner. Only <ilu> and <amg> "
                           "preconditioners are supported."));
}

template <int dim>
void
QMOM<dim>::setup_ilu()
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
QMOM<dim>::setup_amg()
{
  TimerOutput::Scope t(computing_timer, "Setup AMG");

  const Parameters::LinearSolver &linear_solver_parameters =
    simulation_parameters.linear_solver.at(PhysicsID::qmom);

  // Constant modes of the scalar element shared by all the components
  const ComponentMask            components(1, true);
  std::vector<std::vector<bool>> constant_modes =
    DoFTools::extract_constant_modes(*dof_handler, components);

  // The transport operator is non-symmetric since the components are
  // advected by the fluid velocity
  const bool  elliptic              = false;
  const bool  higher_order_elements = fe->degree > 1;
  const bool  output_details        = false;
  const char *smoother_type         = "ILU";
  const char *coarse_type           = "ILU";
  TrilinosWrappers::PreconditionAMG::AdditionalData preconditioner_options(
    elliptic,
    higher_order_elements,
    linear_solver_parameters.amg_n_cycles,
    linear_solver_parameters.amg_w_cycles,
    linear_solver_parameters.amg_aggregation_threshold,
    constant_modes,
    linear_solver_parameters.amg_smoother_sweeps,
    linear_solver_parameters.amg_smoother_overlap,
    output_details,
    smoother_type,
    coarse_type);

  Teuchos::ParameterList              parameter_ml;
  std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
  preconditioner_options.set_parameters(parameter_ml,
                                        distributed_constant_modes,
                                        system_matrix);
  const double ilu_fill = linear_solver_parameters.amg_precond_ilu_fill;
  const double ilu_atol = linear_solver_parameters.amg_precond_ilu_atol;
  const double ilu_rtol = linear_solver_parameters.amg_precond_ilu_rtol;
  parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

  parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);

  amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
  amg_preconditioner->initialize(system_matrix, parameter_ml);
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

  const auto preconditioner_type = linear_solver_parameters.preconditioner;

  // All the components share the matrix and the preconditioner: only the
  // right-hand side changes
  for (unsigned int component = 0; component < n_components; ++component)
    {
      // The components have very different magnitudes, so the tolerance is
      // relative to the right-hand side of each component
      const double current_residual =
        system_rhs[component].l2_norm() / rescale_metric;
      const double linear_solver_tolerance =
        std::max(linear_solver_parameters.relative_residual * current_residual,
                 linear_solver_parameters.minimum_residual);
      if (verbose)
        this->pcout << "  -Tolerance of iterative solver of "
                    << component_names[component]
                    << " is : " << linear_solver_tolerance << std::endl;

      SolverControl solver_control(linear_solver_parameters.max_iterations,
                                   linear_solver_tolerance * rescale_metric,
                                   true,
                                   true);
      const TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
        false, linear_solver_parameters.max_krylov_vectors);
      TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);

      // The solution at the beginning of the time step is the initial guess
      completely_distributed_solution = present_solutions[component];
      if (preconditioner_type ==
          Parameters::LinearSolver::PreconditionerType::ilu)
        solver.solve(system_matrix,
                     completely_distributed_solution,
                     system_rhs[component],
                     *ilu_preconditioner);
      else
        solver.solve(system_matrix,
                     completely_distributed_solution,
                     system_rhs[component],
                     *amg_preconditioner);

      if (verbose)
        this->pcout << "  -Iterative solver of " << component_names[component]
                    << " took : " << solver_control.last_step()
                    << " steps to reach a residual norm of "
                    << solver_control.last_value() / rescale_metric
                    << std::endl;

      constraints.distribute(completely_distributed_solution);
      present_solutions[component] = completely_distributed_solution;
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
  std::vector<std::vector<double>> values(n_components,
                                          std::vector<double>(n_q_points));
  Vector<double>                   exact_values(n_components);
  std::vector<double>              squared_errors(n_components, 0.);

  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      for (unsigned int component = 0; component < n_components; ++component)
        fe_values.get_function_values(present_solutions[component],
                                      values[component]);

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          analytical_solution->vector_value(fe_values.quadrature_point(q),
                                            exact_values);
          for (unsigned int component = 0; component < n_components;
               ++component)
            squared_errors[component] +=
              Utilities::fixed_power<2>(values[component][q] -
                                        exact_values[component]) *
              fe_values.JxW(q);
        }
    }

  std::vector<double> errors(n_components);
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
