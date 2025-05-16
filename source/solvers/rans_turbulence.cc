// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/rans_turbulence.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

template <int dim>
void
RANSTurbulence<dim>::assemble_matrix_and_rhs()
{
  assemble_system_matrix();
  assemble_system_rhs();
}


template <int dim>
void
RANSTurbulence<dim>::assemble_rhs()
{
  assemble_system_rhs();
}

template <int dim>
void
RANSTurbulence<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.emplace_back(
        std::make_shared<RANSTurbulenceAssemblerBDF<dim>>(
          this->simulation_control));
    }

  // Core assembler
  this->assemblers.emplace_back(
    std::make_shared<RANSTurbulenceAssemblerCore<dim>>(
      this->simulation_control));
}

template <int dim>
void
RANSTurbulence<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");

  this->system_matrix = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = RANSTurbulenceScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->rans_mapping,
    dof_handler_fluid->get_fe(),
    *this->face_quadrature);

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &RANSTurbulence::assemble_local_system_matrix,
                  &RANSTurbulence::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
RANSTurbulence<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  RANSTurbulenceScratchData<dim>                       &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  auto source_term = simulation_parameters.source_term.rans_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator fd_cell(&(*triangulation),
                                                         cell->level(),
                                                         cell->index(),
                                                         dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::InitialConditionType::average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_fluid_dynamics(
            fd_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          if (!this->simulation_parameters.ale.enabled())
            scratch_data.reinit_fluid_dynamics(
              fd_cell,
              *multiphysics->get_block_solution(PhysicsID::fluid_dynamics),
              this->simulation_parameters.ale);
        }
    }
  else
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::InitialConditionType::average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_fluid_dynamics(
            fd_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_fluid_dynamics(fd_cell,
                                             *multiphysics->get_solution(
                                               PhysicsID::fluid_dynamics),
                                             this->simulation_parameters.ale);
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
RANSTurbulence<dim>::copy_local_matrix_to_global_matrix(
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
RANSTurbulence<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_rhs = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = RANSTurbulenceScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->turbulence_mapping,
    dof_handler_fluid->get_fe(),
    *this->face_quadrature);

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &RANSTurbulence::assemble_local_system_rhs,
                  &RANSTurbulence::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
RANSTurbulence<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  RANSTurbulenceScratchData<dim>                       &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  auto source_term = simulation_parameters.source_term.rans_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator fd_cell(&(*triangulation),
                                                         cell->level(),
                                                         cell->index(),
                                                         dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::InitialConditionType::average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_fluid_dynamics(
            fd_cell,
            *multiphysics->get_block_time_average_solution(
              PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_fluid_dynamics(fd_cell,
                                             *multiphysics->get_block_solution(
                                               PhysicsID::fluid_dynamics),
                                             this->simulation_parameters.ale);
        }
      scratch_data.reinit_velocity_gradient(
        *multiphysics->get_block_solution(PhysicsID::fluid_dynamics));
    }
  else
    {
      // Check if the post processed variable needs to be calculated with the
      // average velocity profile or the fluid solution.
      if (this->simulation_parameters.initial_condition->type ==
            Parameters::InitialConditionType::average_velocity_profile &&
          !this->simulation_parameters.multiphysics.fluid_dynamics &&
          simulation_control->get_current_time() >
            this->simulation_parameters.post_processing
              .initial_time_for_average_velocities)
        {
          scratch_data.reinit_fluid_dynamics(
            fd_cell,
            *multiphysics->get_time_average_solution(PhysicsID::fluid_dynamics),
            this->simulation_parameters.ale);
        }
      else
        {
          scratch_data.reinit_fluid_dynamics(fd_cell,
                                             *multiphysics->get_solution(
                                               PhysicsID::fluid_dynamics),
                                             this->simulation_parameters.ale);
        }

      scratch_data.reinit_velocity_gradient(
        *multiphysics->get_solution(PhysicsID::fluid_dynamics));
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
RANSTurbulence<dim>::copy_local_rhs_to_global_rhs(
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
void
RANSTurbulence<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  // TODO compute turbulent viscosity from k-e model
  data_out.add_data_vector(dof_handler,
                           present_solution,
                           "turbulent viscosity");
}

template <int dim>
void
RANSTurbulence<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_communicator();
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  if (this_mpi_process == 0 &&
      simulation_parameters.analytical_solution->verbosity !=
        Parameters::Verbosity::quiet)
    {
      error_table.omit_column_from_convergence_rate_evaluation("cells");


      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          error_table.evaluate_all_convergence_rates(
            ConvergenceTable::reduction_rate_log2);
        }
      error_table.set_scientific("error_turbulent_viscosity", true);
      error_table.set_precision("error_turbulent_viscosity",
                                simulation_control->get_log_precision());
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
RANSTurbulence<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
RANSTurbulence<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      double turbulent_viscosity_error = calculate_L2_error();

      error_table.add_value("cells",
                            this->triangulation->n_global_active_cells());
      error_table.add_value("error_turbulent_viscosity",
                            turbulent_viscosity_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error turbulent viscosity : "
                      << turbulent_viscosity_error << std::endl;
        }
    }

  // Set-up domain name for output files
  Parameters::FluidIndicator monitored_fluid =
    this->simulation_parameters.post_processing.postprocessed_fluid;
  // default: monophase simulations
  std::string domain_name("fluid");
  bool        gather_vof(false);

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "RANS Turbulence model");
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }
}

template <int dim>
void
RANSTurbulence<dim>::pre_mesh_adaptation()
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
RANSTurbulence<dim>::post_mesh_adaptation()
{
  auto mpi_communicator = triangulation->get_communicator();


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
void
RANSTurbulence<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  if (ivar.first == Variable::turbulent_viscosity)
    {
      const FEValuesExtractors::Scalar turbulent_viscosity(0);

      KellyErrorEstimator<dim>::estimate(
        *this->rans_mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(turbulent_viscosity));
    }
}

// TODO Continue checking from here

template <int dim>
void
RANSTurbulence<dim>::write_checkpoint()
{
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer = std::make_shared<
    parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>(
    dof_handler);

  sol_set_transfer.emplace_back(&present_solution);
  for (const auto &previous_solution : previous_solutions)
    {
      sol_set_transfer.emplace_back(&previous_solution);
    }

  std::string checkpoint_file_prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  if (simulation_parameters.post_processing.calculate_average_temp_and_hf)
    {
      std::vector<const GlobalVectorType *> avg_scalar_set_transfer =
        this->average_temperature->save(checkpoint_file_prefix);

      sol_set_transfer.insert(sol_set_transfer.end(),
                              avg_scalar_set_transfer.begin(),
                              avg_scalar_set_transfer.end());
    }

  solution_transfer->prepare_for_serialization(sol_set_transfer);

  // Serialize error table
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (this->simulation_parameters.analytical_solution->calculate_error())
    serialize_table(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_RANS" + suffix);
  // if (this->simulation_parameters.post_processing.calculate_heat_flux)
  //   serialize_table(
  //     this->heat_flux_table,
  //     prefix +
  //       this->simulation_parameters.post_processing.heat_flux_output_name +
  //       suffix);
  // if (this->simulation_parameters.post_processing
  //       .calculate_temperature_statistics)
  //   serialize_table(
  //     this->statistics_table,
  //     prefix +
  //       this->simulation_parameters.post_processing.temperature_output_name +
  //       suffix);

  if (this->simulation_parameters.post_processing.calculate_liquid_fraction)
    serialize_table(this->liquid_fraction_table,
                    prefix +
                      this->simulation_parameters.post_processing
                        .liquid_fraction_output_name +
                      suffix);
}

template <int dim>
void
RANSTurbulence<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_communicator();
  this->pcout << "Reading RANS turbulence checkpoint" << std::endl;

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

  parallel::distributed::SolutionTransfer<dim, GlobalVectorType>
    system_trans_vectors(this->dof_handler);

  std::string checkpoint_file_prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;

  if (simulation_parameters.post_processing.calculate_average_temp_and_hf)
    {
      std::vector<GlobalVectorType *> sum_vectors =
        this->average_temperature->read(checkpoint_file_prefix);

      input_vectors.insert(input_vectors.end(),
                           sum_vectors.begin(),
                           sum_vectors.end());
    }

  solution_transfer->deserialize(input_vectors);

  if (simulation_parameters.post_processing.calculate_average_temp_and_hf)
    {
      // Reset the time-averaged turbulence variables if the initial time
      // for averaging has not been reached
      if ((this->simulation_parameters.post_processing
             .initial_time_for_average_temp_and_hf +
           1e-6 * simulation_control->get_time_step()) >
          this->simulation_control->get_current_time())
        {
          this->pcout
            << "Warning: The checkpointed time-averaged turbulence have been reinitialized because the initial averaging time has not yet been reached."
            << std::endl;
          this->average_temperature->zero_average_after_restart();
        }
      this->average_temperature->sanitize_after_restart();
    }

  present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }

  // Deserialize error table
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (this->simulation_parameters.analytical_solution->calculate_error())
    deserialize_table(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_RANS" + suffix);
  // if (this->simulation_parameters.post_processing.calculate_heat_flux)
  //   deserialize_table(
  //     this->heat_flux_table,
  //     prefix +
  //       this->simulation_parameters.post_processing.heat_flux_output_name +
  //       suffix);
  // if (this->simulation_parameters.post_processing
  //       .calculate_temperature_statistics)
  //   deserialize_table(
  //     this->statistics_table,
  //     prefix +
  //       this->simulation_parameters.post_processing.temperature_output_name +
  //       suffix);
  // if (this->simulation_parameters.post_processing.calculate_liquid_fraction)
  //   deserialize_table(this->liquid_fraction_table,
  //                     prefix +
  //                       this->simulation_parameters.post_processing
  //                         .liquid_fraction_output_name +
  //                       suffix);
}


template <int dim>
void
RANSTurbulence<dim>::setup_dofs()
{
  verify_consistency_of_boundary_conditions();

  // Proceed with setting up the DoFs
  dof_handler.distribute_dofs(*fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  auto mpi_communicator = triangulation->get_communicator();


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
    nonzero_constraints.reinit(this->locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_ht.type)
      {
        // Dirichlet condition : imposed temperature at i_bc
        if (type == BoundaryConditions::BoundaryType::temperature)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              id,
              *this->simulation_parameters.boundary_conditions_ht
                 .dirichlet_value.at(id),
              nonzero_constraints);
          }
        if (type == BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              id,
              this->simulation_parameters.boundary_conditions_ht
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_ht
                .periodic_direction.at(id),
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();

  // Boundary conditions for Newton correction
  {
    zero_constraints.clear();
    zero_constraints.reinit(this->locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            zero_constraints);

    for (auto const &[id, type] :
         this->simulation_parameters.boundary_conditions_rans.type)
      {
        // if (type == BoundaryConditions::BoundaryType::temperature)
        //   {
        //     VectorTools::interpolate_boundary_values(
        //       this->dof_handler,
        //       id,
        //       Functions::ZeroFunction<dim>(),
        //       zero_constraints);
        //   }
        // if (type == BoundaryConditions::BoundaryType::periodic)
        //   {
        //     DoFTools::make_periodicity_constraints(
        //       this->dof_handler,
        //       id,
        //       this->simulation_parameters.boundary_conditions_ht
        //         .periodic_neighbor_id.at(id),
        //       this->simulation_parameters.boundary_conditions_ht
        //         .periodic_direction.at(id),
        //       zero_constraints);
        //   }
      }
  }
  zero_constraints.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  /*keep_constrained_dofs = */ true);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of thermal degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the heat transfer dof_handler and present solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::rans_turbulence, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::rans_turbulence, &this->present_solution);
  multiphysics->set_previous_solutions(PhysicsID::rans_turbulence,
                                       &this->previous_solutions);
}

template <int dim>
void
RANSTurbulence<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions_ht.time_dependent)
    return;

  double time = this->simulation_control->get_current_time();
  // We begin by setting the new time for all expressions, although the change
  // for the convection-radiation-flux boundary conditions won't be applied in
  // this function
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_ht.type)
    {
      this->simulation_parameters.boundary_conditions_ht.dirichlet_value.at(id)
        ->set_time(time);
      this->simulation_parameters.boundary_conditions_ht.h.at(id)->set_time(
        time);
      this->simulation_parameters.boundary_conditions_ht.Tinf.at(id)->set_time(
        time);
      this->simulation_parameters.boundary_conditions_ht.emissivity.at(id)
        ->set_time(time);
      this->simulation_parameters.boundary_conditions_ht.heat_flux_bc.at(id)
        ->set_time(time);
    }

  nonzero_constraints.clear();
  nonzero_constraints.reinit(this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          nonzero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_ht.type)
    {
      // Dirichlet condition : imposed temperature at i_bc
      if (type == BoundaryConditions::BoundaryType::temperature)
        {
          VectorTools::interpolate_boundary_values(
            this->dof_handler,
            id,
            *this->simulation_parameters.boundary_conditions_ht.dirichlet_value
               .at(id),
            nonzero_constraints);
        }
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions_ht
              .periodic_neighbor_id.at(id),
            this->simulation_parameters.boundary_conditions_ht
              .periodic_direction.at(id),
            nonzero_constraints);
        }
    }
  nonzero_constraints.close();
}

template <int dim>
void
RANSTurbulence<dim>::set_initial_conditions()
{
  VectorTools::interpolate(*this->rans_mapping,
                           dof_handler,
                           simulation_parameters.initial_condition->temperature,
                           newton_update);
  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  percolate_time_vectors();
}

template <int dim>
void
RANSTurbulence<dim>::solve_linear_system(const bool initial_step,
                                         const bool /*renewed_matrix*/)
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = triangulation->get_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
      .minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
      .relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const double ilu_fill =
    simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
      .ilu_precond_fill;
  const double ilu_atol =
    simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
      .ilu_precond_atol;
  const double ilu_rtol =
    simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(simulation_parameters.linear_solver
                                 .at(PhysicsID::rans_turbulence)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
      .max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::rans_turbulence)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

template class RANSTurbulence<2>;
template class RANSTurbulence<3>;
