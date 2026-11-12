// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/boundary_conditions.h>
#include <core/dimensionality.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
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
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_component_interpretation.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <span>
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
  , n_scalars(p_simulation_parameters.multiphysics.qmom_parameters
                .get_number_of_nodes())
  , n_moments(
      p_simulation_parameters.multiphysics.qmom_parameters.number_of_moments)
  , n_components(n_moments + n_scalars)
  , first_moment_component(n_scalars)
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
  // Relative tolerance used to compare parameter values
  constexpr double parameter_relative_tolerance = 1e-12;
  const auto       values_are_equal = [](const double first_value,
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
  // are used by the phase-space step
  AssertThrow(
    simulation_parameters.physical_properties_manager.get_number_of_fluids() ==
      1,
    InvalidNumberOfFluid(simulation_parameters.physical_properties_manager
                           .get_number_of_fluids()));

  lipid_nanoparticle_production =
    simulation_parameters.multiphysics.qmom_parameters
      .lipid_nanoparticle_production;

  this->pcout << std::setprecision(simulation_control->get_log_precision())
              << std::scientific;

  // Names of the fields, used for the output and the error table
  field_names = simulation_parameters.multiphysics.qmom_parameters.scalar_names;
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    field_names.emplace_back("moment_" + std::to_string(i_moment));

  // All the components share the same discontinuous scalar element, so that
  // they have the same support points.
  const unsigned int qmom_degree =
    simulation_parameters.fem_parameters.qmom_degree;
  if (simulation_parameters.mesh.simplex)
    {
      // Simplex meshes (triangles and tetrahedra)
      const FE_SimplexDGP<dim> scalar_fe(qmom_degree);
      fe      = std::make_unique<FESystem<dim>>(scalar_fe, n_components);
      mapping = std::make_shared<MappingFE<dim>>(FE_SimplexP<dim>(1));
      cell_quadrature = std::make_unique<QGaussSimplex<dim>>(qmom_degree + 1);
      face_quadrature =
        std::make_unique<QGaussSimplex<dim - 1>>(qmom_degree + 1);
    }
  else
    {
      // Quadrilateral and hexahedral meshes. A degree of 0 still requires a
      // linear mapping.
      fe =
        std::make_unique<FESystem<dim>>(FE_DGQ<dim>(qmom_degree), n_components);
      mapping = std::make_shared<MappingQ<dim>>(std::max(1U, qmom_degree));
      cell_quadrature = std::make_unique<QGauss<dim>>(qmom_degree + 1);
      face_quadrature = std::make_unique<QGauss<dim - 1>>(qmom_degree + 1);
    }

  // Functions of all the fields. Their constructors verify the number of
  // expressions given for the moments.
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

  present_solution = std::make_shared<GlobalVectorType>();
  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);

  previous_solutions = std::make_shared<std::vector<GlobalVectorType>>(
    maximum_number_of_previous_solutions());
  previous_solutions_transfer.reserve(previous_solutions->size());
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    previous_solutions_transfer.emplace_back(
      SolutionTransfer<dim, GlobalVectorType>(*dof_handler));
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
QMOM<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;

  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation(n_components,
                             DataComponentInterpretation::component_is_scalar);

  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    *dof_handler,
    *present_solution,
    field_names,
    component_interpretation);
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
  for (unsigned int i_previous_solution = previous_solutions->size() - 1;
       i_previous_solution > 0;
       --i_previous_solution)
    (*previous_solutions)[i_previous_solution] =
      (*previous_solutions)[i_previous_solution - 1];
  (*previous_solutions)[0] = *present_solution;
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
        error_table.add_value("error_" + field_names[component],
                              errors[component]);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        for (unsigned int component = 0; component < n_components; ++component)
          this->pcout << "L2 error " << field_names[component] << ": "
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
  solution_transfer->prepare_for_coarsening_and_refinement(*present_solution);
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    previous_solutions_transfer[i_previous_solution]
      .prepare_for_coarsening_and_refinement(
        (*previous_solutions)[i_previous_solution]);
}

template <int dim>
void
QMOM<dim>::post_mesh_adaptation()
{
  const auto mpi_communicator = triangulation->get_mpi_communicator();

  GlobalVectorType interpolated_solution(locally_owned_dofs, mpi_communicator);
  solution_transfer->interpolate(interpolated_solution);
  *present_solution = interpolated_solution;

  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    {
      GlobalVectorType interpolated_previous_solution(locally_owned_dofs,
                                                      mpi_communicator);
      previous_solutions_transfer[i_previous_solution].interpolate(
        interpolated_previous_solution);
      (*previous_solutions)[i_previous_solution] =
        interpolated_previous_solution;
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
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer =
    std::make_unique<SolutionTransfer<dim, GlobalVectorType>>(*dof_handler);

  sol_set_transfer.emplace_back(&(*present_solution));
  for (const auto &previous_solution : *previous_solutions)
    sol_set_transfer.emplace_back(&previous_solution);

  solution_transfer->prepare_for_serialization(sol_set_transfer);
}

template <int dim>
void
QMOM<dim>::read_checkpoint()
{
  const auto mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading QMOM checkpoint" << std::endl;

  std::vector<GlobalVectorType *> input_vectors(1 + previous_solutions->size());
  GlobalVectorType distributed_system(locally_owned_dofs, mpi_communicator);
  input_vectors[0] = &distributed_system;

  std::vector<GlobalVectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions->size());
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(locally_owned_dofs, mpi_communicator));
      input_vectors[i_previous_solution + 1] =
        &distributed_previous_solutions[i_previous_solution];
    }

  solution_transfer->deserialize(input_vectors);

  *present_solution = distributed_system;
  for (unsigned int i_previous_solution = 0;
       i_previous_solution < previous_solutions->size();
       ++i_previous_solution)
    (*previous_solutions)[i_previous_solution] =
      distributed_previous_solutions[i_previous_solution];
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
  nonzero_constraints.clear();
  nonzero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  nonzero_constraints.close();
}

template <int dim>
void
QMOM<dim>::setup_dofs()
{
  dof_handler->distribute_dofs(*fe);

  const auto mpi_communicator = triangulation->get_mpi_communicator();

  locally_owned_dofs    = dof_handler->locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(*dof_handler);

  present_solution->reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           mpi_communicator);
  for (auto &solution : *previous_solutions)
    solution.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    mpi_communicator);

  evaluation_point.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  newton_update.reinit(locally_owned_dofs, mpi_communicator);
  local_evaluation_point.reinit(locally_owned_dofs, mpi_communicator);

  ghosted_stage_solution.reinit(locally_owned_dofs,
                                locally_relevant_dofs,
                                mpi_communicator);
  step_initial_solution.reinit(locally_owned_dofs, mpi_communicator);
  stage_solution.reinit(locally_owned_dofs, mpi_communicator);
  transport_residual.reinit(locally_owned_dofs, mpi_communicator);
  transport_time_derivative.reinit(locally_owned_dofs, mpi_communicator);

  define_constraints();

  this->pcout << "   Number of QMOM degrees of freedom: "
              << dof_handler->n_dofs() << std::endl;

  // Share the QMOM fields with the other physics
  multiphysics->set_dof_handler(PhysicsID::qmom, dof_handler);
  multiphysics->set_mapping(PhysicsID::qmom, mapping);
  multiphysics->set_solution(PhysicsID::qmom, present_solution);
  multiphysics->set_previous_solutions(PhysicsID::qmom, previous_solutions);
}

template <int dim>
void
QMOM<dim>::set_initial_conditions()
{
  initial_conditions->set_time(simulation_control->get_current_time());

  GlobalVectorType interpolated_solution(locally_owned_dofs,
                                         triangulation->get_mpi_communicator());
  VectorTools::interpolate(*mapping,
                           *dof_handler,
                           *initial_conditions,
                           interpolated_solution);
  *present_solution = interpolated_solution;

  percolate_time_vectors();
}

template <int dim>
void
QMOM<dim>::assemble_system_matrix()
{}

template <int dim>
void
QMOM<dim>::assemble_system_rhs()
{}

template <int dim>
void
QMOM<dim>::solve_linear_system()
{
  advance_time_step();
}

template <int dim>
void
QMOM<dim>::advance_time_step()
{
  // The time of the simulation control is already at the end of the step
  const double time_step = simulation_control->get_time_step();
  const double time_step_start =
    simulation_control->get_current_time() - time_step;

  switch (simulation_parameters.multiphysics.qmom_parameters.splitting_scheme)
    {
      case Parameters::QMOMSplittingScheme::lie:
        advance_transport(time_step_start, time_step);
        apply_phase_space_step(time_step);
        break;
      case Parameters::QMOMSplittingScheme::strang:
        apply_phase_space_step(0.5 * time_step);
        advance_transport(time_step_start, time_step);
        apply_phase_space_step(0.5 * time_step);
        break;
    }
}

template <int dim>
void
QMOM<dim>::advance_transport(const double time_step_start,
                             const double time_step)
{
  TimerOutput::Scope t(computing_timer, "Transport");
}

template <int dim>
void
QMOM<dim>::apply_phase_space_step(const double time_step)
{
  TimerOutput::Scope t(computing_timer, "Phase-space step");
}

template <int dim>
template <typename FluidVectorType>
void
QMOM<dim>::evaluate_transport_operator(
  const GlobalVectorType &solution,
  const double            time,
  const double            interpolation_factor,
  const FluidVectorType  &previous_fluid_solution,
  const FluidVectorType  &present_fluid_solution,
  GlobalVectorType       &time_derivative)
{
  for (const auto &[boundary_id, function] : dirichlet_functions)
    function->set_time(time);

  const DoFHandler<dim> &fluid_dof_handler =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const auto get_fluid_cell =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell) {
      return typename DoFHandler<dim>::active_cell_iterator(&(*triangulation),
                                                            cell->level(),
                                                            cell->index(),
                                                            &fluid_dof_handler);
    };

  // 1. Residual R(u) of the upwind discontinuous Galerkin discretization of
  // du/dt + div(u_fluid u) = 0 for every field: cell terms, interior face
  // fluxes and boundary fluxes
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
      // faces may be assembled when they are finer than an owned neighbor
      copy.local_rhs = 0.;
      copy.face_data.clear();
      copy.cell_is_local = cell->is_locally_owned();
      if (!copy.cell_is_local)
        return;

      cell->get_dof_indices(copy.local_dof_indices);
      scratch.reinit(cell, solution);
      scratch.reinit_velocity(get_fluid_cell(cell),
                              previous_fluid_solution,
                              present_fluid_solution,
                              interpolation_factor);

      // Weak form of the advective flux: int_K c u . grad(phi)
      const FEValues<dim> &fe_values = scratch.fe_values;
      for (unsigned int q = 0; q < scratch.n_q_points; ++q)
        {
          const double JxW = fe_values.JxW(q);
          for (unsigned int i = 0; i < scratch.n_dofs; ++i)
            copy.local_rhs(i) +=
              scratch.values[q](scratch.component_of_dof[i]) *
              (scratch.velocity_values[q] * fe_values.shape_grad(i, q)) * JxW;
        }
    };

  const auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        QMOMScratchData<dim>                                 &scratch,
        QMOMCopyData                                         &copy) {
      scratch.reinit_boundary_face(cell, face_no, solution);
      scratch.reinit_face_velocity(get_fluid_cell(cell),
                                   face_no,
                                   numbers::invalid_unsigned_int,
                                   previous_fluid_solution,
                                   present_fluid_solution,
                                   interpolation_factor);
      assemble_boundary_face_flux(scratch,
                                  cell->face(face_no)->boundary_id(),
                                  copy.local_rhs);
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
                                   neighbor_sub_face_no,
                                   solution);
      scratch.reinit_face_velocity(get_fluid_cell(cell),
                                   face_no,
                                   sub_face_no,
                                   previous_fluid_solution,
                                   present_fluid_solution,
                                   interpolation_factor);

      QMOMCopyData::FaceData &face_data = copy.face_data.emplace_back();
      face_data.joint_dof_indices =
        scratch.fe_interface_values.get_interface_dof_indices();
      face_data.face_rhs.reinit(scratch.n_interface_dofs);
      assemble_interior_face_flux(scratch, face_data.face_rhs);
    };

  // Each face is assembled once per process that owns one of its cells, with
  // the contributions to both cells. Only the locally owned entries are
  // added, so that every degree of freedom receives the flux of each of its
  // faces exactly once and no off-process entry is written.
  std::vector<types::global_dof_index> owned_face_indices;
  std::vector<double>                  owned_face_values;
  const auto                           copier = [&](const QMOMCopyData &copy) {
    if (copy.cell_is_local)
      transport_residual.add(copy.local_dof_indices, copy.local_rhs);

    for (const QMOMCopyData::FaceData &face_data : copy.face_data)
      {
        owned_face_indices.clear();
        owned_face_values.clear();
        for (unsigned int i = 0; i < face_data.joint_dof_indices.size(); ++i)
          if (locally_owned_dofs.is_element(face_data.joint_dof_indices[i]))
            {
              owned_face_indices.push_back(face_data.joint_dof_indices[i]);
              owned_face_values.push_back(face_data.face_rhs(i));
            }
        transport_residual.add(owned_face_indices, owned_face_values);
      }
  };

  transport_residual = 0.;
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
  transport_residual.compress(VectorOperation::add);

  // 2. L(u) = M^{-1} R(u). The mass matrix is block diagonal: one block per
  // cell, itself made of identical blocks for every component, since all the
  // components share the same scalar element.
  const unsigned int n_base_dofs = fe->base_element(0).n_dofs_per_cell();
  FEValues<dim>      fe_values(*mapping,
                               *fe,
                               *cell_quadrature,
                               update_values | update_JxW_values);
  FullMatrix<double> inverse_mass_matrix(n_base_dofs, n_base_dofs);
  Vector<double>     local_residual(n_base_dofs);
  Vector<double>     local_time_derivative(n_base_dofs);
  std::vector<types::global_dof_index> dof_indices(fe->n_dofs_per_cell());

  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);

      inverse_mass_matrix = 0.;
      for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
        for (unsigned int i = 0; i < n_base_dofs; ++i)
          {
            const double phi_i =
              fe_values.shape_value(fe->component_to_system_index(0, i), q);
            for (unsigned int j = 0; j < n_base_dofs; ++j)
              inverse_mass_matrix(i, j) +=
                phi_i *
                fe_values.shape_value(fe->component_to_system_index(0, j), q) *
                fe_values.JxW(q);
          }
      inverse_mass_matrix.gauss_jordan();

      for (unsigned int component = 0; component < n_components; ++component)
        {
          for (unsigned int i = 0; i < n_base_dofs; ++i)
            local_residual(i) = transport_residual(
              dof_indices[fe->component_to_system_index(component, i)]);
          inverse_mass_matrix.vmult(local_time_derivative, local_residual);
          for (unsigned int i = 0; i < n_base_dofs; ++i)
            time_derivative(
              dof_indices[fe->component_to_system_index(component, i)]) =
              local_time_derivative(i);
        }
    }
  time_derivative.compress(VectorOperation::insert);
}

template <int dim>
void
QMOM<dim>::assemble_interior_face_flux(const QMOMScratchData<dim> &scratch_data,
                                       Vector<double> &face_rhs) const
{
  const FEInterfaceValues<dim> &fe_interface_values =
    scratch_data.fe_interface_values;
  const std::vector<double>         &JxW = fe_interface_values.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals =
    fe_interface_values.get_normal_vectors();

  // Upwind flux (u . n) c^upwind, with n pointing from the first to the
  // second cell. The jump of the test function, phi_here - phi_there, gives
  // the contributions of opposite signs to the two cells.
  for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
    {
      const double velocity_dot_n =
        scratch_data.face_velocity_values[q] * normals[q];
      const Vector<double> &upwind_values = velocity_dot_n > 0. ?
                                              scratch_data.values_here[q] :
                                              scratch_data.values_there[q];

      for (unsigned int i = 0; i < scratch_data.n_interface_dofs; ++i)
        {
          const unsigned int component =
            scratch_data.component_of_interface_dof[i];
          face_rhs(i) -=
            velocity_dot_n * upwind_values(component) *
            fe_interface_values.jump_in_shape_values(i, q, component) * JxW[q];
        }
    }
}

template <int dim>
void
QMOM<dim>::assemble_boundary_face_flux(QMOMScratchData<dim>    &scratch_data,
                                       const types::boundary_id boundary_id,
                                       Vector<double>          &local_rhs) const
{
  // Boundaries without QMOM boundary condition are treated as none
  const auto &boundary_types =
    simulation_parameters.boundary_conditions_qmom.type;
  const auto type_iterator = boundary_types.find(boundary_id);
  if (type_iterator == boundary_types.end() ||
      type_iterator->second == BoundaryConditions::BoundaryType::none)
    return;
  const BoundaryConditions::BoundaryType type = type_iterator->second;
  AssertThrow(type == BoundaryConditions::BoundaryType::qmom_dirichlet ||
                type == BoundaryConditions::BoundaryType::outlet,
              ExcInternalError());

  const FEFaceValuesBase<dim> &fe_face =
    scratch_data.fe_interface_values.get_fe_face_values(0);
  const std::vector<double>         &JxW      = fe_face.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals  = fe_face.get_normal_vectors();
  const std::vector<Point<dim>>     &q_points = fe_face.get_quadrature_points();
  Vector<double>                    &upwind_values = scratch_data.upwind_values;

  for (unsigned int q = 0; q < scratch_data.n_face_q_points; ++q)
    {
      const double velocity_dot_n =
        scratch_data.face_velocity_values[q] * normals[q];

      // Outflow: the interior values leave the domain. Inflow: the Dirichlet
      // values enter it, and nothing enters through an outlet.
      if (velocity_dot_n > 0.)
        upwind_values = scratch_data.values_here[q];
      else if (type == BoundaryConditions::BoundaryType::qmom_dirichlet)
        dirichlet_functions.at(boundary_id)
          ->vector_value(q_points[q], upwind_values);
      else
        upwind_values = 0.;

      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
        local_rhs(i) -= velocity_dot_n *
                        upwind_values(scratch_data.component_of_dof[i]) *
                        fe_face.shape_value(i, q) * JxW[q];
    }
}

template <int dim>
std::vector<double>
QMOM<dim>::calculate_L2_errors()
{
  analytical_solution->set_time(simulation_control->get_current_time());

  FEValues<dim>      fe_values(*mapping,
                               *fe,
                               *cell_quadrature,
                               update_values | update_quadrature_points |
                                 update_JxW_values);
  const unsigned int n_q_points = cell_quadrature->size();
  std::vector<Vector<double>> values(n_q_points, Vector<double>(n_components));
  Vector<double>              exact_values(n_components);
  std::vector<double>         squared_errors(n_components, 0.);

  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      fe_values.get_function_values(*present_solution, values);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          analytical_solution->vector_value(fe_values.quadrature_point(q),
                                            exact_values);
          for (unsigned int component = 0; component < n_components;
               ++component)
            squared_errors[component] +=
              Utilities::fixed_power<2>(values[q](component) -
                                        exact_values(component)) *
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
void
QMOM<dim>::output_newton_update_norms()
{
  this->pcout << "\t  ||dx||_L2 = " << std::setw(6) << newton_update.l2_norm()
              << std::setw(6)
              << "\t  ||dx||_Linfty = " << newton_update.linfty_norm()
              << std::endl;
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
