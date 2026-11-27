// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/dimensionality.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
#include <core/qmom_boundary_conditions.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <solvers/analytical_solutions.h>
#include <solvers/copy_data.h>
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

#ifdef LETHE_WITH_QMOMKIT
#  include <qmomkit/moment_set.hpp>
#  include <qmomkit/moment_source.hpp>
#  include <qmomkit/node_advection.hpp>
#  include <qmomkit/phase_space_evolution.hpp>
#  include <qmomkit/qmom.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
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
  , has_phase_space_mechanism(
      p_simulation_parameters.multiphysics.qmom_parameters.enable_growth ||
      p_simulation_parameters.multiphysics.qmom_parameters.enable_nucleation ||
      p_simulation_parameters.multiphysics.qmom_parameters.enable_coalescence)
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
  // The steady state of the Lie splitting depends on the time step, and the
  // phase-space evolution needs a physical time step
  using TimeSteppingMethod = Parameters::SimulationControl::TimeSteppingMethod;
  const TimeSteppingMethod method =
    simulation_parameters.simulation_control.method;
  AssertThrow(method == TimeSteppingMethod::bdf1 ||
                method == TimeSteppingMethod::bdf2 ||
                method == TimeSteppingMethod::bdf3,
              ExcMessage("The QMOM physics only supports the bdf1, bdf2 and "
                         "bdf3 time-stepping methods, but \"method\" is " +
                         enum_to_string(method) +
                         " in the simulation control subsection. Set it to "
                         "bdf1, bdf2 or bdf3."));
  // With the Lie splitting, the previous solutions of the BDF scheme include
  // the changes of the previous phase-space steps. A BDF scheme of order 2 or
  // 3 extrapolates them, so that the transport repeats part of these changes,
  // even without velocity. The Lie splitting is first-order accurate, so the
  // bdf1 scheme loses no accuracy.
  AssertThrow(!has_phase_space_mechanism || method == TimeSteppingMethod::bdf1,
              ExcMessage("The growth, nucleation and coalescence of the QMOM "
                         "physics only support the bdf1 time-stepping method, "
                         "but \"method\" is " +
                         enum_to_string(method) +
                         " in the simulation control subsection. Set it to "
                         "bdf1."));
  AssertThrow(qmom_parameters.splitting_scheme ==
                Parameters::QMOMSplittingScheme::lie,
              ExcMessage("The strang splitting scheme of the QMOM physics is "
                         "not implemented yet. Set \"splitting scheme\" to lie "
                         "in the QMOM subsection."));

  // Growth, nucleation and coalescence of the particles in phase space
  phase_space_source_terms = std::make_unique<QMOMSourceTerms>(qmom_parameters);

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
       simulation_parameters.qmom_boundary_conditions.dirichlet)
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
        {
          error_table.add_value("relative_error_" + component_names[component],
                                errors[component]);
          error_table.set_scientific("relative_error_" +
                                       component_names[component],
                                     true);
        }

      if (Utilities::MPI::this_mpi_process(
            triangulation->get_mpi_communicator()) == 0)
        {
          const std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.analytical_solution->get_filename() +
            "_qmom.dat";
          std::ofstream output(filename.c_str());
          error_table.write_text(output);
        }

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        for (unsigned int component = 0; component < n_components; ++component)
          this->pcout << "Relative L2 error " << component_names[component]
                      << ": " << errors[component] << std::endl;
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
  output_vectors.reserve(distributed_solutions.size());
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
QMOM<dim>::write_checkpoint(const std::string & /*checkpoint_slot_label*/)
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
QMOM<dim>::read_checkpoint(const std::string & /*checkpoint_slot_label*/)
{
  const MPI_Comm mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading QMOM checkpoint" << std::endl;

  // Set up the vectors for the deserialization, in the order of
  // write_checkpoint()
  std::vector<GlobalVectorType> distributed_solutions(
    n_components * (1 + maximum_number_of_previous_solutions()),
    GlobalVectorType(locally_owned_dofs, mpi_communicator));
  std::vector<GlobalVectorType *> input_vectors;
  input_vectors.reserve(distributed_solutions.size());
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
  if (!simulation_parameters.qmom_boundary_conditions.time_dependent)
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
  // includes the inverse of the time step
  const double bdf_coefficient_present =
    simulation_control->get_bdf_coefficients()[0];

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
      // flux crosses a boundary of type noflux
      if (simulation_parameters.qmom_boundary_conditions.type.at(
            cell->face(face_no)->boundary_id()) ==
          BoundaryConditions::QMOMBoundaryType::noflux)
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

  // Number of previous solutions of the BDF time derivative
  const unsigned int n_previous_solutions =
    number_of_previous_solutions(simulation_control->get_assembly_method());
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
      if (simulation_parameters.qmom_boundary_conditions.type.at(boundary_id) !=
          BoundaryConditions::QMOMBoundaryType::dirichlet)
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
        if (has_phase_space_mechanism)
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

  // present_solutions are ghosted, hence read-only: the solution is updated
  // in a locally owned copy
  GlobalVectorType locally_owned_present_solution(
    locally_owned_dofs, triangulation->get_mpi_communicator());
  GlobalVectorType residual(locally_owned_dofs,
                            triangulation->get_mpi_communicator());
  GlobalVectorType completely_distributed_solution(
    locally_owned_dofs, triangulation->get_mpi_communicator());

  const auto preconditioner_type = linear_solver_parameters.preconditioner;

  // All the components share the matrix and the preconditioner: only the
  // right-hand side changes
  for (unsigned int component = 0; component < n_components; ++component)
    {
      // The system is solved for the correction of the present solution: the
      // right-hand side is the residual of that solution, and the tolerance
      // is relative to it.
      locally_owned_present_solution = present_solutions[component];
      const double current_residual =
        system_matrix.residual(residual,
                               locally_owned_present_solution,
                               system_rhs[component]) /
        rescale_metric;
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

      completely_distributed_solution = 0.;
      if (preconditioner_type ==
          Parameters::LinearSolver::PreconditionerType::ilu)
        solver.solve(system_matrix,
                     completely_distributed_solution,
                     residual,
                     *ilu_preconditioner);
      else
        solver.solve(system_matrix,
                     completely_distributed_solution,
                     residual,
                     *amg_preconditioner);

      if (verbose)
        this->pcout << "  -Iterative solver of " << component_names[component]
                    << " took : " << solver_control.last_step()
                    << " steps to reach a residual norm of "
                    << solver_control.last_value() / rescale_metric
                    << std::endl;

      locally_owned_present_solution += completely_distributed_solution;
      constraints.distribute(locally_owned_present_solution);
      present_solutions[component] = locally_owned_present_solution;
    }
}

template <int dim>
void
QMOM<dim>::advance_phase_space([[maybe_unused]] const double time_step)
{
#ifdef LETHE_WITH_QMOMKIT
  TimerOutput::Scope t(computing_timer, "Phase-space evolution");

  const Parameters::QMOM &qmom_parameters =
    simulation_parameters.multiphysics.qmom_parameters;

  // The phase-space step starts where the physical-space transport started,
  // at the beginning of the time step
  const double phase_space_start_time =
    simulation_control->get_current_time() - time_step;

  // Callbacks of QMOMKit, built once by the source terms: the growth moves
  // the quadrature nodes, the coalescence changes the moments through the
  // quadrature, and the nucleation adds particles. A disabled mechanism
  // contributes zero.
  const qmomkit::PhaseSpaceVelocity &phase_space_velocity =
    phase_space_source_terms->get_phase_space_velocity();
  const qmomkit::QmomMomentSource &moment_source =
    phase_space_source_terms->get_moment_source();
  const qmomkit::QmomMomentBoundaryContribution &boundary_contribution =
    phase_space_source_terms->get_boundary_contribution();

  // Node count and inversion options shared by all the degrees of freedom.
  // The moments of the case are not used: each degree of freedom builds its
  // own state from its moments.
  qmomkit::QmomCase qmom_case;
  qmom_case.node_count = qmom_parameters.number_of_nodes;
  qmom_case.qmom_options.adaptive =
    qmom_parameters.phase_evolution.enable_adaptive_node_count;

  // Tolerances of the adaptive Runge-Kutta integration of the trajectories
  // of the nodes
  qmomkit::QmomPhaseSpaceEvolutionOptions evolution_options;
  evolution_options.node_options.adaptive_options.relative_tolerance =
    qmom_parameters.phase_evolution.relative_tolerance;
  evolution_options.node_options.adaptive_options.absolute_tolerance =
    qmom_parameters.phase_evolution.absolute_tolerance;

  // The moments are read from the ghosted present solutions and written into
  // locally owned (non-ghosted) vectors, which are copied back at the end
  std::vector<GlobalVectorType> evolved_moments(n_moments);
  for (GlobalVectorType &evolved_moment : evolved_moments)
    evolved_moment.reinit(locally_owned_dofs,
                          triangulation->get_mpi_communicator());

  // The discontinuous elements are nodal and no degree of freedom is shared
  // between cells, so the moments of each degree of freedom are the moments
  // at its support point and are evolved independently of the others
  std::vector<double> local_scalars(n_scalars);
  std::vector<double> local_moments(n_moments);
  std::vector<double> nucleation_moment_rates(n_moments);
  bool                evolution_failed = false;
  for (const types::global_dof_index dof : locally_owned_dofs)
    {
      // The models of the mechanisms may depend on the scalars at this degree
      // of freedom (e.g. through the supersaturation of a process). The
      // scalars are not evolved in phase space: they keep their values at the
      // end of the transport in physical space during the whole step.
      for (unsigned int i_scalar = 0; i_scalar < n_scalars; ++i_scalar)
        local_scalars[i_scalar] = present_solutions[i_scalar][dof];
      phase_space_source_terms->set_local_scalars(local_scalars);

      // The moments are left unchanged unless they are evolved below. They
      // are only inserted: switching to addition is a collective operation.
      for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
        {
          local_moments[i_moment] =
            present_solutions[n_scalars + i_moment][dof];
          evolved_moments[i_moment][dof] = local_moments[i_moment];
        }

      const qmomkit::QmomStateCandidate candidate =
        qmom_case.state_candidate(qmomkit::MomentSet(local_moments),
                                  evolution_options.realizability_policy);

      // Without particles, only the nucleation changes the moments. Its rate
      // is integrated with the midpoint rule, exact for the constant model.
      if (candidate.status == qmomkit::StateCandidateStatus::empty)
        {
          if (phase_space_source_terms->nucleation_is_enabled())
            {
              phase_space_source_terms->compute_nucleation_moment_rates(
                phase_space_start_time + 0.5 * time_step,
                nucleation_moment_rates);
              for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
                evolved_moments[i_moment][dof] =
                  local_moments[i_moment] +
                  time_step * nucleation_moment_rates[i_moment];
            }
          continue;
        }

      // The moments cannot be inverted into a quadrature
      if (!candidate.ok())
        {
          evolution_failed = true;
          break;
        }

      const qmomkit::QmomPhaseSpaceEvolutionResult result =
        qmomkit::advance_qmom_phase_space(*candidate.state,
                                          phase_space_start_time,
                                          time_step,
                                          phase_space_velocity,
                                          moment_source,
                                          boundary_contribution,
                                          evolution_options);

      // The evolved quadrature cannot be rebuilt into moments
      if (!result.ok())
        {
          evolution_failed = true;
          break;
        }

      const std::vector<double> &moment_values = result.state.moment_values();
      for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
        evolved_moments[i_moment][dof] = moment_values[i_moment];
    }

  // A rank that throws alone would leave the others waiting in the next
  // collective call, so the failure is shared by all the ranks
  AssertThrow(
    !Utilities::MPI::logical_or(evolution_failed,
                                triangulation->get_mpi_communicator()),
    ExcMessage(
      "The QMOM phase-space evolution failed: the moments of at least one "
      "degree of freedom are not realizable and cannot be inverted or "
      "evolved by QMOMKit. Refine the mesh or reduce the time step to limit "
      "the oscillations of the transport in physical space, or set \"enable "
      "adaptive node count\" to true in the \"phase evolution\" subsection "
      "of the QMOM subsection."));

  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    {
      evolved_moments[i_moment].compress(VectorOperation::insert);
      present_solutions[n_scalars + i_moment] = evolved_moments[i_moment];
    }
#else
  // The constructor of the multiphysics interface already stops the
  // simulation in this case
  AssertThrow(false, QMOMWithoutQMOMKitError());
#endif
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
  std::vector<double>              squared_exact_norms(n_components, 0.);

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
            {
              squared_errors[component] +=
                Utilities::fixed_power<2>(values[component][q] -
                                          exact_values[component]) *
                fe_values.JxW(q);
              squared_exact_norms[component] +=
                Utilities::fixed_power<2>(exact_values[component]) *
                fe_values.JxW(q);
            }
        }
    }

  std::vector<double> global_squared_errors(n_components);
  std::vector<double> global_squared_exact_norms(n_components);
  Utilities::MPI::sum(squared_errors,
                      triangulation->get_mpi_communicator(),
                      global_squared_errors);
  Utilities::MPI::sum(squared_exact_norms,
                      triangulation->get_mpi_communicator(),
                      global_squared_exact_norms);

  // The moments differ by orders of magnitude: each error is relative to the
  // norm of its analytical solution, unless that solution vanishes
  std::vector<double> errors(n_components);
  for (unsigned int component = 0; component < n_components; ++component)
    errors[component] = global_squared_exact_norms[component] > 0. ?
                          std::sqrt(global_squared_errors[component] /
                                    global_squared_exact_norms[component]) :
                          std::sqrt(global_squared_errors[component]);
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
