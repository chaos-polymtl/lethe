// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/lethe_grid_tools.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/tracer.h>
#include <solvers/tracer_assemblers.h>
#include <solvers/tracer_scratch_data.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <map>

template <int dim>
void
Tracer<dim>::setup_assemblers()
{
  AssertThrow(
    is_sdirk(this->simulation_control->get_assembly_method()) == false,
    ExcMessage("The SDIRK scheme is not yet supported for this physics"));

  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.emplace_back(
        std::make_shared<TracerAssemblerBDF<dim>>(this->simulation_control));
    }
  // Core assemblers are different between DG and CG versions.
  if (simulation_parameters.fem_parameters.tracer_uses_dg)
    {
      this->assemblers.emplace_back(
        std::make_shared<TracerAssemblerDGCore<dim>>());
      this->inner_face_assembler = std::make_shared<TracerAssemblerSIPG<dim>>();
      this->boundary_face_assembler =
        std::make_shared<TracerAssemblerBoundaryNitsche<dim>>(
          simulation_parameters.boundary_conditions_tracer);
    }
  else
    {
      this->assemblers.emplace_back(
        std::make_shared<TracerAssemblerCore<dim>>(this->simulation_control));
    }

  auto reaction_constant_model =
    simulation_parameters.physical_properties_manager
      .get_tracer_reaction_prefactor();
  if (!std::dynamic_pointer_cast<NoneTracerReactionPrefactor>(
        reaction_constant_model))
    this->assemblers.emplace_back(
      std::make_shared<TracerAssemblerReaction<dim>>(this->simulation_control));
}

template <int dim>
void
Tracer<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");

  this->system_matrix = 0;
  setup_assemblers();

  // Update the source term time
  simulation_parameters.source_term.tracer_source->set_time(
    simulation_control->get_current_time());

  if (simulation_parameters.fem_parameters.tracer_uses_dg)
    assemble_system_matrix_dg();
  else
    assemble_system_matrix_cg();
}

template <int dim>
void
Tracer<dim>::assemble_system_matrix_cg()
{
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = TracerScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->face_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &Tracer::assemble_local_system_matrix,
                  &Tracer::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  system_matrix.compress(VectorOperation::add);
}


template <int dim>
void
Tracer<dim>::assemble_system_matrix_dg()
{
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = TracerScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->face_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe());

  StabilizedDGMethodsCopyData copy_data(this->fe->n_dofs_per_cell(),
                                        this->cell_quadrature->size());

  // We first wrap the assembly of the matrix within a cell_worker lambda
  // function. This is only done for compatibility reasons with the MeshWorker
  // paradigm and does not have any functional purpose.
  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        TracerScratchData<dim>                               &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      this->assemble_local_system_matrix(cell, scratch_data, copy_data);
    };

  const auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                   &face_no,
        TracerScratchData<dim>                               &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      const auto boundary_index = cell->face(face_no)->boundary_id();

      scratch_data.reinit_boundary_face(cell,
                                        face_no,
                                        boundary_index,
                                        this->evaluation_point);

      scratch_data.reinit_signed_distance_at_face(
        cell, face_no, this->multiphysics->get_immersed_solid_shape());

      // Gather velocity information at the face to properly advect
      const DoFHandler<dim> *dof_handler_fluid =
        multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
      // Get the cell that corresponds to the fluid dynamics
      typename DoFHandler<dim>::active_cell_iterator velocity_cell(
        &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

      // Reinit the internal face velocity within the scratch data
      reinit_face_velocity_with_adequate_solution(velocity_cell,
                                                  face_no,
                                                  scratch_data);

      scratch_data.calculate_face_physical_properties();

      this->boundary_face_assembler->assemble_matrix(scratch_data, copy_data);
    };

  const auto face_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                   &face_no,
        const unsigned int                                   &sub_face_no,
        const typename DoFHandler<dim>::active_cell_iterator &neigh_cell,
        const unsigned int                                   &neigh_face_no,
        const unsigned int                                   &neigh_sub_face_no,
        TracerScratchData<dim>                               &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      scratch_data.reinit_internal_face(cell,
                                        face_no,
                                        sub_face_no,
                                        neigh_cell,
                                        neigh_face_no,
                                        neigh_sub_face_no,
                                        this->evaluation_point);

      scratch_data.reinit_signed_distance_at_face(
        cell, face_no, this->multiphysics->get_immersed_solid_shape());

      // Pad copy_data memory for the internal faces elementary matrices
      // BB note : Array could be pre-allocated
      copy_data.face_data.emplace_back();
      auto &copy_data_face = copy_data.face_data.back();
      copy_data_face.face_matrix.reinit(scratch_data.n_interface_dofs,
                                        scratch_data.n_interface_dofs);
      copy_data_face.joint_dof_indices =
        scratch_data.fe_interface_values_tracer.get_interface_dof_indices();

      // Gather velocity information at the face to properly advect
      const DoFHandler<dim> *dof_handler_fluid =
        multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
      // Get the cell that corresponds to the fluid dynamics
      typename DoFHandler<dim>::active_cell_iterator velocity_cell(
        &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

      // Reinit the internal face velocity within the scratch data
      reinit_face_velocity_with_adequate_solution(velocity_cell,
                                                  face_no,
                                                  scratch_data);

      scratch_data.calculate_face_physical_properties();

      this->inner_face_assembler->assemble_matrix(scratch_data, copy_data);
    };

  const auto copier = [&](const StabilizedDGMethodsCopyData &copy_data) {
    this->copy_local_matrix_to_global_matrix(copy_data);

    const AffineConstraints<double> &constraints_used = this->zero_constraints;

    for (const auto &cdf : copy_data.face_data)
      {
        constraints_used.distribute_local_to_global(cdf.face_matrix,
                                                    cdf.joint_dof_indices,
                                                    system_matrix);
      }
  };

  MeshWorker::mesh_loop(this->dof_handler.begin_active(),
                        this->dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_ghost_faces_both,
                        boundary_worker,
                        face_worker);
}



template <int dim>
void
Tracer<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  TracerScratchData<dim>                               &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*simulation_parameters.source_term.tracer_source));

  scratch_data.reinit_signed_distance(
    cell, &(*this->multiphysics->get_immersed_solid_shape()));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
    }
  else
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
Tracer<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              system_matrix);
}

template <int dim>
void
Tracer<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  // TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  this->system_rhs = 0;
  setup_assemblers();

  // Update the source term time
  simulation_parameters.source_term.tracer_source->set_time(
    simulation_control->get_current_time());


  if (simulation_parameters.fem_parameters.tracer_uses_dg)
    assemble_system_rhs_dg();

  else
    assemble_system_rhs_cg();
}


template <int dim>
void
Tracer<dim>::assemble_system_rhs_cg()
{
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = TracerScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->face_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &Tracer::assemble_local_system_rhs,
                  &Tracer::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
Tracer<dim>::assemble_system_rhs_dg()
{
  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = TracerScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->face_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe());

  StabilizedDGMethodsCopyData copy_data(this->fe->n_dofs_per_cell(),
                                        this->cell_quadrature->size());

  const auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        TracerScratchData<dim>                               &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      this->assemble_local_system_rhs(cell, scratch_data, copy_data);
    };

  const auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                   &face_no,
        TracerScratchData<dim>                               &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data) {
      // Identify which boundary condition corresponds to the boundary id.
      const auto boundary_index = cell->face(face_no)->boundary_id();

      scratch_data.reinit_boundary_face(cell,
                                        face_no,
                                        boundary_index,
                                        this->evaluation_point);

      scratch_data.reinit_signed_distance_at_face(
        cell, face_no, this->multiphysics->get_immersed_solid_shape());

      // Gather velocity information at the face to properly advect
      // First gather the dof handler for the fluid dynamics
      const DoFHandler<dim> *dof_handler_fluid =
        multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
      // Get the cell that corresponds to the fluid dynamics
      typename DoFHandler<dim>::active_cell_iterator velocity_cell(
        &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

      reinit_face_velocity_with_adequate_solution(velocity_cell,
                                                  face_no,
                                                  scratch_data);

      scratch_data.calculate_face_physical_properties();

      this->boundary_face_assembler->assemble_rhs(scratch_data, copy_data);
    };

  const auto face_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                   &face_no,
        const unsigned int                                   &sub_face_no,
        const typename DoFHandler<dim>::active_cell_iterator &neigh_cell,
        const unsigned int                                   &neigh_face_no,
        const unsigned int                                   &neigh_sub_face_no,
        TracerScratchData<dim>                               &scratch_data,
        StabilizedDGMethodsCopyData                          &copy_data)

  {
    scratch_data.reinit_internal_face(cell,
                                      face_no,
                                      sub_face_no,
                                      neigh_cell,
                                      neigh_face_no,
                                      neigh_sub_face_no,
                                      this->evaluation_point);

    scratch_data.reinit_signed_distance_at_face(
      cell, face_no, this->multiphysics->get_immersed_solid_shape());

    copy_data.face_data.emplace_back();
    auto &copy_data_face = copy_data.face_data.back();
    copy_data_face.joint_dof_indices =
      scratch_data.fe_interface_values_tracer.get_interface_dof_indices();
    copy_data_face.face_rhs.reinit(scratch_data.n_interface_dofs);

    // Gather velocity information at the face to properly advect
    // First gather the dof handler for the fluid dynamics
    const DoFHandler<dim> *dof_handler_fluid =
      multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
    // Get the cell that corresponds to the fluid dynamics
    typename DoFHandler<dim>::active_cell_iterator velocity_cell(
      &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

    reinit_face_velocity_with_adequate_solution(velocity_cell,
                                                face_no,
                                                scratch_data);

    scratch_data.calculate_face_physical_properties();

    this->inner_face_assembler->assemble_rhs(scratch_data, copy_data);
  };


  const auto copier = [&](const StabilizedDGMethodsCopyData &c) {
    this->copy_local_rhs_to_global_rhs(c);
    const AffineConstraints<double> &constraints_used = this->zero_constraints;
    for (const auto &cdf : c.face_data)
      {
        constraints_used.distribute_local_to_global(cdf.face_rhs,
                                                    cdf.joint_dof_indices,
                                                    system_rhs);
      }
  };

  MeshWorker::mesh_loop(this->dof_handler.begin_active(),
                        this->dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once |
                          MeshWorker::assemble_ghost_faces_both,
                        boundary_worker,
                        face_worker);
}

template <int dim>
void
Tracer<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  TracerScratchData<dim>                               &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  auto source_term = simulation_parameters.source_term.tracer_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  scratch_data.reinit_signed_distance(
    cell, &(*this->multiphysics->get_immersed_solid_shape()));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
    }
  else
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::FluidDynamicsInitialConditionType::
              average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
      else
        {
          scratch_data.reinit_velocity(
            velocity_cell,
            *multiphysics->get_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale,
            this->simulation_parameters.tracer_drift_velocity.drift_velocity);
        }
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
Tracer<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsCopyData &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              system_rhs);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
Tracer<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;
  std::vector<std::string>                         solution_names(1, "tracer");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation(1,
                             DataComponentInterpretation::component_is_scalar);
  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    dof_handler,
    present_solution,
    solution_names,
    component_interpretation);
  return solution_output_structs;
}

template <int dim>
double
Tracer<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = cell_quadrature->size();

  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->tracer;
  exact_solution.set_time(simulation_control->get_current_time());

  double l2error = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_scalar_values);

          // Get the exact solution at all gauss points
          exact_solution.value_list(fe_values.get_quadrature_points(),
                                    q_exact_solution);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              double sim   = q_scalar_values[q];
              double exact = q_exact_solution[q];
              l2error += (sim - exact) * (sim - exact) * fe_values.JxW(q);
            }
        }
    }
  l2error = Utilities::MPI::sum(l2error, mpi_communicator);
  return l2error;
}

template <int dim>
void
Tracer<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_mpi_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity ==
        Parameters::Verbosity::verbose)
    {
      error_table.omit_column_from_convergence_rate_evaluation("cells");

      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          error_table.evaluate_all_convergence_rates(
            ConvergenceTable::reduction_rate_log2);
        }
      error_table.set_scientific("error_tracer", true);
      error_table.set_precision("error_tracer",
                                simulation_control->get_log_precision());
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
Tracer<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
Tracer<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      double tracer_error = calculate_L2_error();

      error_table.add_value("cells",
                            this->triangulation->n_global_active_cells());
      error_table.add_value("error_tracer", tracer_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error tracer: " << tracer_error << std::endl;
        }
    }

  if (simulation_parameters.post_processing.calculate_tracer_statistics)
    {
      calculate_tracer_statistics();
      if (simulation_control->get_step_number() %
            this->simulation_parameters.post_processing.output_frequency ==
          0)
        this->write_tracer_statistics();
    }

  // Calculate tracer flow rate at every boundary
  if (this->simulation_parameters.post_processing.calculate_tracer_flow_rate)
    {
      TimerOutput::Scope t(this->computing_timer, "Calculate tracer flow rate");

      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          announce_string(this->pcout, "Tracer flow rates");
        }

      std::map<types::boundary_id, double> tracer_flow_rates;
      if (multiphysics->fluid_dynamics_is_block())
        {
          tracer_flow_rates = postprocess_tracer_flow_rate(
            *multiphysics->get_block_solution(PhysicsID::fluid_dynamics));
        }
      else
        {
          tracer_flow_rates = postprocess_tracer_flow_rate(
            *multiphysics->get_solution(PhysicsID::fluid_dynamics));
        }
      this->write_tracer_flow_rates(tracer_flow_rates);
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "Tracer");
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }
}

template <int dim>
void
Tracer<dim>::calculate_tracer_statistics()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  FEValues<dim> fe_values(*mapping,
                          *fe,
                          *cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int  n_q_points = cell_quadrature->size();
  std::vector<double> q_tracer_values(n_q_points);

  double volume_integral  = 0;
  double max_tracer_value = DBL_MIN;
  double min_tracer_value = DBL_MAX;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_tracer_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              volume_integral += q_tracer_values[q] * fe_values.JxW(q);
              max_tracer_value = std::max(q_tracer_values[q], max_tracer_value);
              min_tracer_value = std::min(q_tracer_values[q], min_tracer_value);
            }
        }
    }
  volume_integral      = Utilities::MPI::sum(volume_integral, mpi_communicator);
  double global_volume = GridTools::volume(*triangulation, *mapping);
  double tracer_average = volume_integral / global_volume;

  double variance_integral = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values.get_function_values(present_solution, q_tracer_values);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              variance_integral += (q_tracer_values[q] - tracer_average) *
                                   (q_tracer_values[q] - tracer_average) *
                                   fe_values.JxW(q);
            }
        }
    }

  variance_integral = Utilities::MPI::sum(variance_integral, mpi_communicator);
  double tracer_variance      = variance_integral / global_volume;
  double tracer_std_deviation = std::sqrt(tracer_variance);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << std::setprecision(
        this->simulation_control->get_log_precision());
      this->pcout << "Tracer statistics: " << std::endl;
      this->pcout << "\t     Min: " << min_tracer_value << std::endl;
      this->pcout << "\t     Max: " << max_tracer_value << std::endl;
      this->pcout << "\t Average: " << tracer_average << std::endl;
      this->pcout << "\t Std-Dev: " << tracer_std_deviation << std::endl;
    }

  statistics_table.add_value("time", simulation_control->get_current_time());
  statistics_table.set_scientific("time", true);
  statistics_table.set_precision("time", 12);

  statistics_table.add_value("min", min_tracer_value);
  statistics_table.set_scientific("min", true);
  statistics_table.set_precision("min", 12);

  statistics_table.add_value("max", max_tracer_value);
  statistics_table.set_scientific("max", true);
  statistics_table.set_precision("max", 12);

  statistics_table.add_value("average", tracer_average);
  statistics_table.set_scientific("average", true);
  statistics_table.set_precision("average", 12);

  statistics_table.add_value("std-dev", tracer_std_deviation);
  statistics_table.set_scientific("std-dev", true);
  statistics_table.set_precision("std-dev", 12);
}


template <int dim>
template <typename VectorType>
std::map<types::boundary_id, double>
Tracer<dim>::postprocess_tracer_flow_rate(const VectorType &current_solution_fd)
{
  const unsigned int n_q_points_face = this->face_quadrature->size();
  const MPI_Comm mpi_communicator    = this->dof_handler.get_mpi_communicator();

  // Initialize tracer information
  std::vector<double>         tracer_values(n_q_points_face);
  std::vector<Tensor<1, dim>> tracer_gradient(n_q_points_face);
  FEFaceValues<dim>           fe_face_values_tracer(*this->mapping,
                                          this->dof_handler.get_fe(),
                                          *this->face_quadrature,
                                          update_values | update_gradients |
                                            update_JxW_values |
                                            update_normal_vectors |
                                            update_quadrature_points);

  // Initialize fluid dynamics information
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      velocity_values(n_q_points_face);
  FEFaceValues<dim>                fe_face_values_fd(*this->mapping,
                                      dof_handler_fd->get_fe(),
                                      *this->face_quadrature,
                                      update_values);

  // Initialize fluid properties
  auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::map<field, std::vector<double>> fields;

  const auto diffusivity_model = properties_manager.get_tracer_diffusivity();

  std::vector<double>     tracer_diffusivity(n_q_points_face);
  std::vector<Point<dim>> face_quadrature_points;

  std::map<types::boundary_id, double> tracer_flow_rate;

  // Get vector of all boundary conditions
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          for (const auto face : cell->face_indices())
            {
              if (cell->face(face)->at_boundary())
                {
                  const auto boundary_id = cell->face(face)->boundary_id();

                  // Gather tracer information
                  fe_face_values_tracer.reinit(cell, face);
                  fe_face_values_tracer.get_function_values(
                    this->present_solution, tracer_values);
                  fe_face_values_tracer.get_function_gradients(
                    this->present_solution, tracer_gradient);

                  // We update the fields required by the diffusivity
                  // model
                  fields.clear();
                  if (diffusivity_model->depends_on(field::levelset))
                    {
                      std::vector<double> levelset_values(n_q_points_face);
                      fields.insert(
                        std::pair<field, std::vector<double>>(field::levelset,
                                                              n_q_points_face));
                      face_quadrature_points =
                        fe_face_values_tracer.get_quadrature_points();
                      this->multiphysics->get_immersed_solid_shape()
                        ->value_list(face_quadrature_points, levelset_values);
                      set_field_vector(field::levelset,
                                       levelset_values,
                                       fields);
                    }

                  diffusivity_model->vector_value(fields, tracer_diffusivity);

                  // Get fluid dynamics active cell iterator
                  typename DoFHandler<dim>::active_cell_iterator cell_fd(
                    &(*(this->triangulation)),
                    cell->level(),
                    cell->index(),
                    dof_handler_fd);

                  // Gather fluid dynamics information
                  fe_face_values_fd.reinit(cell_fd, face);
                  fe_face_values_fd[velocities].get_function_values(
                    current_solution_fd, velocity_values);

                  // Loop on the quadrature points
                  for (unsigned int q = 0; q < n_q_points_face; q++)
                    {
                      Tensor<1, dim> normal_vector_tracer =
                        -fe_face_values_tracer.normal_vector(q);

                      tracer_flow_rate[boundary_id] +=
                        (-tracer_diffusivity[q] * tracer_gradient[q] *
                           normal_vector_tracer +
                         tracer_values[q] * velocity_values[q] *
                           normal_vector_tracer) *
                        fe_face_values_tracer.JxW(q);
                    } // end loop on quadrature points
                }     // end face is a boundary face
            }         // end loop on faces
        }             // end condition cell at boundary
    }                 // end loop on cells


  // Sum across all cores
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    tracer_flow_rate[id] =
      Utilities::MPI::sum(tracer_flow_rate[id], mpi_communicator);

  return tracer_flow_rate;
}

template <int dim>
void
Tracer<dim>::write_tracer_flow_rates(
  const std::map<types::boundary_id, double> &tracer_flow_rate_map)
{
  // Fill table
  this->tracer_flow_rate_table.add_value(
    "time", this->simulation_control->get_current_time());
  this->tracer_flow_rate_table.set_precision("time", 12);
  for (auto const &[id, value] : tracer_flow_rate_map)
    {
      this->tracer_flow_rate_table.add_value("flow-rate-" +
                                               Utilities::int_to_string(id, 2),
                                             value);
      this->tracer_flow_rate_table.set_scientific(
        "flow-rate-" + Utilities::int_to_string(id, 2), true);
    }

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Tracer flow rate at the boundaries: " << std::endl;
      for (auto const &[id, value] : tracer_flow_rate_map)
        this->pcout << "\t boundary " << id << ": " << value << std::endl;
    }

  auto mpi_communicator = triangulation->get_mpi_communicator();
  if ((simulation_control->get_step_number() %
         this->simulation_parameters.post_processing.output_frequency ==
       0) &&
      Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.tracer_flow_rate_output_name +
        ".dat";
      std::ofstream output(filename.c_str());
      this->tracer_flow_rate_table.write_text(output);
    }
}

template <int dim>
void
Tracer<dim>::write_tracer_statistics()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.tracer_output_name + ".dat";
      std::ofstream output(filename.c_str());


      statistics_table.write_text(output);
    }
}

template <int dim>
void
Tracer<dim>::pre_mesh_adaptation()
{
  solution_transfer->prepare_for_coarsening_and_refinement(present_solution);

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }
}

template <int dim>
void
Tracer<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();

  // Set up the vectors for the transfer
  GlobalVectorType tmp(locally_owned_dofs, mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer->interpolate(tmp);

  // Distribute constraints
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;

  // Transfer previous solutions
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      GlobalVectorType tmp_previous_solution(locally_owned_dofs,
                                             mpi_communicator);
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }
}

template <int dim>
std::vector<OutputStructTableHandler>
Tracer<dim>::gather_tables()
{
  std::vector<OutputStructTableHandler> table_output_structs;

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";

  if (this->simulation_parameters.analytical_solution->calculate_error())
    table_output_structs.emplace_back(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_tracer" + suffix);

  if (this->simulation_parameters.post_processing.calculate_tracer_statistics)
    table_output_structs.emplace_back(
      this->statistics_table,
      prefix + this->simulation_parameters.post_processing.tracer_output_name +
        suffix);

  if (this->simulation_parameters.post_processing.calculate_tracer_flow_rate)
    table_output_structs.emplace_back(
      this->tracer_flow_rate_table,
      prefix +
        this->simulation_parameters.post_processing
          .tracer_flow_rate_output_name +
        suffix);

  return table_output_structs;
}

template <int dim>
void
Tracer<dim>::write_checkpoint()
{
  auto mpi_communicator = this->triangulation->get_mpi_communicator();
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer =
    std::make_shared<SolutionTransfer<dim, GlobalVectorType>>(dof_handler);

  sol_set_transfer.emplace_back(&present_solution);
  for (const auto &previous_solution : previous_solutions)
    {
      sol_set_transfer.emplace_back(&previous_solution);
    }
  solution_transfer->prepare_for_serialization(sol_set_transfer);

  // Serialize all post-processing tables that are currently used with the
  // Tracer solver
  const std::vector<OutputStructTableHandler> &table_output_structs =
    this->gather_tables();
  serialize_tables_vector(table_output_structs, mpi_communicator);
}

template <int dim>
void
Tracer<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_mpi_communicator();
  this->pcout << "Reading tracer checkpoint" << std::endl;

  std::vector<GlobalVectorType *> input_vectors(1 + previous_solutions.size());
  GlobalVectorType distributed_system(locally_owned_dofs, mpi_communicator);
  input_vectors[0] = &distributed_system;


  std::vector<GlobalVectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        GlobalVectorType(locally_owned_dofs, mpi_communicator));
      input_vectors[i + 1] = &distributed_previous_solutions[i];
    }

  solution_transfer->deserialize(input_vectors);

  present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }

  // Deserialize all post-processing tables that are currently used with the
  // Tracer solver
  std::vector<OutputStructTableHandler> table_output_structs =
    this->gather_tables();
  deserialize_tables_vector(table_output_structs, mpi_communicator);
}


template <int dim>
void
Tracer<dim>::setup_dofs()
{
  verify_consistency_of_boundary_conditions();

  dof_handler.distribute_dofs(*fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  auto mpi_communicator = triangulation->get_mpi_communicator();


  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  present_solution.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);

  // Previous solutions for transient schemes
  for (auto &solution : this->previous_solutions)
    {
      solution.reinit(locally_owned_dofs,
                      locally_relevant_dofs,
                      mpi_communicator);
    }

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

  newton_update.reinit(locally_owned_dofs, mpi_communicator);

  local_evaluation_point.reinit(this->locally_owned_dofs, mpi_communicator);

  {
    nonzero_constraints.clear();
    nonzero_constraints.reinit(locally_owned_dofs, this->locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_tracer.type)
      {
        // Dirichlet condition : imposed temperature at i_bc
        if (type == BoundaryConditions::BoundaryType::tracer_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              id,
              *this->simulation_parameters.boundary_conditions_tracer.tracer.at(
                id),
              nonzero_constraints);
          }
        if (type == BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              id,
              this->simulation_parameters.boundary_conditions_tracer
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_tracer
                .periodic_direction.at(id),
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();

  // Boundary conditions for Newton correction
  {
    zero_constraints.clear();
    zero_constraints.reinit(this->locally_owned_dofs,
                            this->locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            zero_constraints);

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_tracer.type)
      {
        if (type == BoundaryConditions::BoundaryType::tracer_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              id,
              Functions::ZeroFunction<dim>(),
              zero_constraints);
          }
        if (type == BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              id,
              this->simulation_parameters.boundary_conditions_tracer
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_tracer
                .periodic_direction.at(id),
              zero_constraints);
          }
      }
  }
  zero_constraints.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(locally_relevant_dofs);


  if (simulation_parameters.fem_parameters.tracer_uses_dg)
    {
      DoFTools::make_flux_sparsity_pattern(this->dof_handler,
                                           dsp,
                                           nonzero_constraints,
                                           /*keep_constrained_dofs = */ true);
    }
  else
    {
      DoFTools::make_sparsity_pattern(this->dof_handler,
                                      dsp,
                                      nonzero_constraints,
                                      /*keep_constrained_dofs = */ true);
    }



  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of tracer degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the tracer dof_handler and present solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::tracer, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::tracer, &this->present_solution);
  multiphysics->set_previous_solutions(PhysicsID::tracer,
                                       &this->previous_solutions);
}

template <int dim>
void
Tracer<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions_tracer.time_dependent)
    return;

  double time = this->simulation_control->get_current_time();
  nonzero_constraints.clear();
  nonzero_constraints.reinit(this->locally_owned_dofs,
                             this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          nonzero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_tracer.type)
    {
      // Dirichlet condition : imposed temperature at i_bc
      if (type == BoundaryConditions::BoundaryType::tracer_dirichlet)
        {
          this->simulation_parameters.boundary_conditions_tracer.tracer.at(id)
            ->set_time(time);
          VectorTools::interpolate_boundary_values(
            this->dof_handler,
            id,
            *this->simulation_parameters.boundary_conditions_tracer.tracer.at(
              id),
            nonzero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions_tracer
              .periodic_neighbor_id.at(id),
            this->simulation_parameters.boundary_conditions_tracer
              .periodic_direction.at(id),
            nonzero_constraints);
        }
    }
  nonzero_constraints.close();
  nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim>
void
Tracer<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*mapping,
                           dof_handler,
                           simulation_parameters.initial_condition->tracer,
                           newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  percolate_time_vectors();
}


template <int dim>
void
Tracer<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (ivar.first == Variable::tracer)
    {
      const FEValuesExtractors::Scalar tracer(0);

      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(tracer));
    }
}

template <int dim>
void
Tracer<dim>::solve_linear_system(const bool initial_step,
                                 const bool /*renewed_matrix*/)
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = triangulation->get_mpi_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::tracer).minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::tracer).relative_residual;

  const double rescale_metric   = this->get_residual_rescale_metric();
  const double current_residual = system_rhs.l2_norm() / rescale_metric;
  const double linear_solver_tolerance =
    std::max(relative_residual * current_residual, absolute_residual);
  if (this->simulation_parameters.linear_solver.at(PhysicsID::tracer)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  const double non_rescaled_linear_solver_tolerance =
    linear_solver_tolerance * rescale_metric;

  const unsigned int ilu_fill =
    simulation_parameters.linear_solver.at(PhysicsID::tracer).ilu_precond_fill;
  const double ilu_atol =
    simulation_parameters.linear_solver.at(PhysicsID::tracer).ilu_precond_atol;
  const double ilu_rtol =
    simulation_parameters.linear_solver.at(PhysicsID::tracer).ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(
    simulation_parameters.linear_solver.at(PhysicsID::tracer).max_iterations,
    non_rescaled_linear_solver_tolerance,
    true,
    true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    simulation_parameters.linear_solver.at(PhysicsID::tracer)
      .max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::tracer).verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() / rescale_metric << std::endl;
    }

  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}



template class Tracer<2>;
template class Tracer<3>;
