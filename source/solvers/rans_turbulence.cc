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
    std::make_shared<RANSTurbulenceAssemblerCore<dim>>(this->simulation_control));
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
  RANSTurbulenceScratchData<dim>                         &scratch_data,
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
    *this->temperature_mapping,
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
  RANSTurbulenceScratchData<dim>                         &scratch_data,
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
  data_out.add_data_vector(dof_handler, present_solution, "turbulent viscosity");
}

template <int dim>
double
RANSTurbulence<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values(*this->temperature_mapping,
                          *fe,
                          *this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = this->cell_quadrature->size();

  std::vector<double> q_exact_solution(n_q_points);
  std::vector<double> q_scalar_values(n_q_points);

  auto &exact_solution = simulation_parameters.analytical_solution->temperature;
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
      error_table.add_value("error_turbulent_viscosity", turbulent_viscosity_error);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "L2 error turbulent viscosity : " << turbulent_viscosity_error
                      << std::endl;
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
      announce_string(this->pcout, "Heat Transfer");
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
        "_HT" + suffix);
  if (this->simulation_parameters.post_processing.calculate_heat_flux)
    serialize_table(
      this->heat_flux_table,
      prefix +
        this->simulation_parameters.post_processing.heat_flux_output_name +
        suffix);
  if (this->simulation_parameters.post_processing
        .calculate_temperature_statistics)
    serialize_table(
      this->statistics_table,
      prefix +
        this->simulation_parameters.post_processing.temperature_output_name +
        suffix);

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
  this->pcout << "Reading heat transfer checkpoint" << std::endl;

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
      // Reset the time-averaged temperature and heat flux if the initial time
      // for averaging has not been reached
      if ((this->simulation_parameters.post_processing
             .initial_time_for_average_temp_and_hf +
           1e-6 * simulation_control->get_time_step()) >
          this->simulation_control->get_current_time())
        {
          this->pcout
            << "Warning: The checkpointed time-averaged temperature and heat flux have been reinitialized because the initial averaging time has not yet been reached."
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
        "_HT" + suffix);
  if (this->simulation_parameters.post_processing.calculate_heat_flux)
    deserialize_table(
      this->heat_flux_table,
      prefix +
        this->simulation_parameters.post_processing.heat_flux_output_name +
        suffix);
  if (this->simulation_parameters.post_processing
        .calculate_temperature_statistics)
    deserialize_table(
      this->statistics_table,
      prefix +
        this->simulation_parameters.post_processing.temperature_output_name +
        suffix);
  if (this->simulation_parameters.post_processing.calculate_liquid_fraction)
    deserialize_table(this->liquid_fraction_table,
                      prefix +
                        this->simulation_parameters.post_processing
                          .liquid_fraction_output_name +
                        suffix);
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

    // Initialize the vectors used in the time average temperature calculation
    if (this->simulation_parameters.post_processing
          .calculate_average_temp_and_hf)
      {
        this->average_temperature->initialize_vectors(
          this->locally_owned_dofs,
          this->locally_relevant_dofs,
          mpi_communicator);
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
         this->simulation_parameters.boundary_conditions_ht.type)
      {
        if (type == BoundaryConditions::BoundaryType::temperature)
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
              this->simulation_parameters.boundary_conditions_ht
                .periodic_neighbor_id.at(id),
              this->simulation_parameters.boundary_conditions_ht
                .periodic_direction.at(id),
              zero_constraints);
          }
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
  multiphysics->set_dof_handler(PhysicsID::heat_transfer, &this->dof_handler);
  multiphysics->set_solution(PhysicsID::heat_transfer, &this->present_solution);
  multiphysics->set_previous_solutions(PhysicsID::heat_transfer,
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
  VectorTools::interpolate(*this->temperature_mapping,
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
    simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
      .minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
      .relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const double ilu_fill =
    simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
      .ilu_precond_fill;
  const double ilu_atol =
    simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
      .ilu_precond_atol;
  const double ilu_rtol =
    simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(simulation_parameters.linear_solver
                                 .at(PhysicsID::heat_transfer)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
      .max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::heat_transfer)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  constraints_used.distribute(completely_distributed_solution);
  newton_update = completely_distributed_solution;
}

template <int dim>
void
RANSTurbulence<dim>::postprocess_temperature_statistics(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const bool                       time_average)
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double> local_temperature_values(n_q_points);
  FEValues<dim>       fe_values_ht(*this->temperature_mapping,
                             *this->fe,
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  // Initialize VOF information
  const DoFHandler<dim>         *dof_handler_vof = nullptr;
  std::shared_ptr<FEValues<dim>> fe_values_vof;
  std::vector<double>            filtered_phase_values(n_q_points);

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_values_vof =
        std::make_shared<FEValues<dim>>(*this->temperature_mapping,
                                        dof_handler_vof->get_fe(),
                                        *this->cell_quadrature,
                                        update_values);
    }

  // Other initializations
  double phase_coefficient(0.);
  double volume_integral(0.);
  double temperature_integral(0.);
  double minimum_temperature(std::numeric_limits<double>::max());
  double maximum_temperature(std::numeric_limits<double>::lowest());
  // boolean used to state if the quadrature point is in the
  // postprocessed_fluid, i.e. its temperature should be considered
  // in the minimum and maximum temperature calculation
  bool point_is_in_postprocessed_fluid(false);

  // Calculate min, max and average
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Gather heat transfer information
          fe_values_ht.reinit(cell);

          if (!time_average)
            {
              fe_values_ht.get_function_values(this->present_solution,
                                               local_temperature_values);
            }
          else
            {
              // calculate the average using the time-averaged temperature
              fe_values_ht.get_function_values(
                this->average_temperature->get_average_scalar(),
                local_temperature_values);
            }

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              std::tie(phase_coefficient, point_is_in_postprocessed_fluid) =
                set_phase_coefficient(gather_vof,
                                      monitored_fluid,
                                      filtered_phase_values[q]);

              if (point_is_in_postprocessed_fluid)
                {
                  maximum_temperature =
                    std::max(local_temperature_values[q], maximum_temperature);
                  minimum_temperature =
                    std::min(local_temperature_values[q], minimum_temperature);
                }

              volume_integral += fe_values_ht.JxW(q) * phase_coefficient;
              temperature_integral += local_temperature_values[q] *
                                      phase_coefficient * fe_values_ht.JxW(q);
            } // end loop on quadrature points
        }
    } // end loop on cell

  minimum_temperature =
    Utilities::MPI::min(minimum_temperature, mpi_communicator);
  maximum_temperature =
    Utilities::MPI::max(maximum_temperature, mpi_communicator);

  volume_integral = Utilities::MPI::sum(volume_integral, mpi_communicator);
  temperature_integral =
    Utilities::MPI::sum(temperature_integral, mpi_communicator);
  double temperature_average = temperature_integral / volume_integral;

  // Calculate standard deviation
  double temperature_variance_integral = 0;
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values_ht.reinit(cell);

          if (!time_average)
            {
              fe_values_ht.get_function_values(this->present_solution,
                                               local_temperature_values);
            }
          else
            {
              // calculate the average using the time-averaged temperature
              fe_values_ht.get_function_values(
                this->average_temperature->get_average_scalar(),
                local_temperature_values);
            }

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              phase_coefficient =
                set_phase_coefficient(gather_vof,
                                      monitored_fluid,
                                      filtered_phase_values[q])
                  .first;

              temperature_variance_integral +=
                (local_temperature_values[q] - temperature_average) *
                (local_temperature_values[q] - temperature_average) *
                phase_coefficient * fe_values_ht.JxW(q);
            } // end loop on quadrature points
        }
    } // end loop on cell

  temperature_variance_integral =
    Utilities::MPI::sum(temperature_variance_integral, mpi_communicator);
  double temperature_variance = temperature_variance_integral / volume_integral;
  double temperature_std_deviation = std::sqrt(temperature_variance);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      std::string temperature_label =
        time_average ? "Time-averaged temperature statistics on " :
                       "Temperature statistics on ";
      this->pcout << temperature_label << domain_name << ": " << std::endl;
      this->pcout << "\t     Min: " << minimum_temperature << std::endl;
      this->pcout << "\t     Max: " << maximum_temperature << std::endl;
      this->pcout << "\t Average: " << temperature_average << std::endl;
      this->pcout << "\t Std-Dev: " << temperature_std_deviation << std::endl;
    }

  // Fill table
  if (!time_average)
    {
      this->statistics_table.add_value(
        "time", this->simulation_control->get_current_time());
      this->statistics_table.add_value("min", minimum_temperature);
      this->statistics_table.add_value("max", maximum_temperature);
      this->statistics_table.add_value("average", temperature_average);
      this->statistics_table.add_value("std-dev", temperature_std_deviation);
    }
  else
    {
      this->statistics_table.add_value("min_time_average", minimum_temperature);
      this->statistics_table.add_value("max_time_average", maximum_temperature);
      this->statistics_table.add_value("time_average_average",
                                       temperature_average);
      this->statistics_table.add_value("std-dev_time_average",
                                       temperature_std_deviation);
    }
}

template <int dim>
void
RANSTurbulence<dim>::write_temperature_statistics(const std::string domain_name)
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.temperature_output_name + "_" +
        domain_name + ".dat";
      std::ofstream output(filename.c_str());

      this->statistics_table.write_text(output);
    }
}

template <int dim>
void
RANSTurbulence<dim>::postprocess_liquid_fraction(const bool gather_vof)
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double> local_temperature_values(n_q_points);
  FEValues<dim>       fe_values_ht(*this->temperature_mapping,
                             *this->fe,
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  // Initialize VOF information
  const DoFHandler<dim>         *dof_handler_vof = nullptr;
  std::shared_ptr<FEValues<dim>> fe_values_vof;
  std::vector<double>            filtered_phase_values(n_q_points);

  // Get the raw physical properties parameters to calculate the liquid fraction
  // in-situ
  const auto &physical_properties_parameters =
    this->simulation_parameters.physical_properties_manager
      .get_physical_properties_parameters();

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_values_vof =
        std::make_shared<FEValues<dim>>(*this->temperature_mapping,
                                        dof_handler_vof->get_fe(),
                                        *this->cell_quadrature,
                                        update_values);
    }

  // Variables for the integration
  double volume_integral(0.);
  double liquid_volume_integral(0.);

  // Calculate min, max and average
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Gather heat transfer information
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // If VOF is enabled, gather the liquid fraction if the fluid has
              // a phase fraction
              if (!gather_vof)
                {
                  liquid_volume_integral +=
                    calculate_liquid_fraction(local_temperature_values[q],
                                              physical_properties_parameters
                                                .fluids[0]
                                                .phase_change_parameters) *
                    fe_values_ht.JxW(q);
                  volume_integral += fe_values_ht.JxW(q);
                }
              else
                {
                  // Case of fluid 0 being a phase change
                  if (physical_properties_parameters.fluids[0]
                        .specific_heat_model ==
                      Parameters::Material::SpecificHeatModel::phase_change)
                    {
                      liquid_volume_integral +=
                        (1. - filtered_phase_values[q]) *
                        calculate_liquid_fraction(local_temperature_values[q],
                                                  physical_properties_parameters
                                                    .fluids[0]
                                                    .phase_change_parameters) *
                        fe_values_ht.JxW(q);
                      volume_integral +=
                        (1. - filtered_phase_values[q]) * fe_values_ht.JxW(q);
                    }

                  // Case of fluid 1 being a phase change
                  if (physical_properties_parameters.fluids[1]
                        .specific_heat_model ==
                      Parameters::Material::SpecificHeatModel::phase_change)
                    {
                      liquid_volume_integral +=
                        (filtered_phase_values[q]) *
                        calculate_liquid_fraction(local_temperature_values[q],
                                                  physical_properties_parameters
                                                    .fluids[1]
                                                    .phase_change_parameters) *
                        fe_values_ht.JxW(q);
                      volume_integral +=
                        (filtered_phase_values[q]) * fe_values_ht.JxW(q);
                    }
                }
            } // end loop on quadrature points
        }
    } // end loop on cell

  volume_integral = Utilities::MPI::sum(volume_integral, mpi_communicator);
  liquid_volume_integral =
    Utilities::MPI::sum(liquid_volume_integral, mpi_communicator);
  const double liquid_fraction = liquid_volume_integral / volume_integral;

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Liquid fraction"
                  << ": " << liquid_fraction << std::endl;
    }

  // Fill table
  this->liquid_fraction_table.add_value(
    "time", this->simulation_control->get_current_time());
  this->liquid_fraction_table.add_value("liquid fraction", liquid_fraction);
}

template <int dim>
void
RANSTurbulence<dim>::write_liquid_fraction()
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.liquid_fraction_output_name +
        ".dat";
      std::ofstream output(filename.c_str());

      this->liquid_fraction_table.write_text(output);
    }
}



template <int dim>
template <typename VectorType>
void
RANSTurbulence<dim>::postprocess_heat_flux_on_bc(
  const bool        gather_vof,
  const VectorType &current_solution_fd)
{
  const unsigned int n_q_points_face  = this->face_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double>         local_temperature_values(n_q_points_face);
  std::vector<Tensor<1, dim>> temperature_gradient(n_q_points_face);
  FEFaceValues<dim>           fe_face_values_ht(*this->temperature_mapping,
                                      this->dof_handler.get_fe(),
                                      *this->face_quadrature,
                                      update_values | update_gradients |
                                        update_JxW_values |
                                        update_normal_vectors);

  // Initialize fluid dynamics information
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      local_velocity_values(n_q_points_face);
  FEFaceValues<dim>                fe_face_values_fd(*this->temperature_mapping,
                                      dof_handler_fd->get_fe(),
                                      *this->face_quadrature,
                                      update_values);

  // Initialize VOF information
  DoFHandler<dim>                   *dof_handler_vof = nullptr;
  std::shared_ptr<FEFaceValues<dim>> fe_face_values_vof;
  std::vector<double>                filtered_phase_values(n_q_points_face);

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_face_values_vof =
        std::make_shared<FEFaceValues<dim>>(*this->temperature_mapping,
                                            dof_handler_vof->get_fe(),
                                            *this->face_quadrature,
                                            update_values);
    }

  // Initialize fluid properties
  auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::map<field, double>              field_values;
  std::map<field, std::vector<double>> fields;

  // monophase flow
  double density(0.);
  double specific_heat(0.);
  double thermal_conductivity(0.);

  // multiphase flow
  std::vector<double> density_0(n_q_points_face);
  std::vector<double> specific_heat_0(n_q_points_face);
  std::vector<double> thermal_conductivity_0(n_q_points_face);
  std::vector<double> density_1(n_q_points_face);
  std::vector<double> specific_heat_1(n_q_points_face);
  std::vector<double> thermal_conductivity_1(n_q_points_face);

  switch (properties_manager.get_number_of_fluids())
    {
      default:
        {
          // Get values for monophase flow
          const auto density_model = properties_manager.get_density();
          const auto specific_heat_model =
            properties_manager.get_specific_heat();
          const auto conductivity_model =
            properties_manager.get_thermal_conductivity();

          density              = density_model->value(field_values);
          specific_heat        = specific_heat_model->value(field_values);
          thermal_conductivity = conductivity_model->value(field_values);

          break;
        }
      case 2:
        {
          // Get prm values for multiphase flow - will be blended in the
          // integration loop
          const auto density_models = properties_manager.get_density_vector();
          const auto specific_heat_models =
            properties_manager.get_specific_heat_vector();
          const auto thermal_conductivity_models =
            properties_manager.get_thermal_conductivity_vector();

          density_models[0]->vector_value(fields, density_0);
          specific_heat_models[0]->vector_value(fields, specific_heat_0);
          thermal_conductivity_models[0]->vector_value(fields,
                                                       thermal_conductivity_0);

          density_models[1]->vector_value(fields, density_1);
          specific_heat_models[1]->vector_value(fields, specific_heat_1);
          thermal_conductivity_models[1]->vector_value(fields,
                                                       thermal_conductivity_1);

          break;
        }
    } // end switch on number of fluids

  std::map<types::boundary_id, double> heat_flux;
  std::map<types::boundary_id, double> convective_flux;

  // Get vector of heat transfer boundary conditions
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          for (const auto face : cell->face_indices())
            {
              if (cell->face(face)->at_boundary())
                {
                  types::boundary_id boundary_index =
                    cell->face(face)->boundary_id();

                  // Gather h coefficient and T_inf
                  const double h_coefficient =
                    this->simulation_parameters.boundary_conditions_ht.h
                      .at(boundary_index)
                      ->value(Point<dim>());
                  const double T_inf =
                    this->simulation_parameters.boundary_conditions_ht.Tinf
                      .at(boundary_index)
                      ->value(Point<dim>());


                  // Gather heat transfer information
                  fe_face_values_ht.reinit(cell, face);
                  fe_face_values_ht.get_function_values(
                    this->present_solution, local_temperature_values);
                  fe_face_values_ht.get_function_gradients(
                    this->present_solution, temperature_gradient);

                  // Get fluid dynamics active cell iterator
                  typename DoFHandler<dim>::active_cell_iterator cell_fd(
                    &(*(this->triangulation)),
                    cell->level(),
                    cell->index(),
                    dof_handler_fd);

                  // Gather fluid dynamics information
                  fe_face_values_fd.reinit(cell_fd, face);
                  fe_face_values_fd[velocities].get_function_values(
                    current_solution_fd, local_velocity_values);

                  if (gather_vof)
                    {
                      // Get VOF active cell iterator
                      typename DoFHandler<dim>::active_cell_iterator cell_vof(
                        &(*(this->triangulation)),
                        cell->level(),
                        cell->index(),
                        dof_handler_vof);

                      // Gather VOF information
                      fe_face_values_vof->reinit(cell_vof, face);
                      fe_face_values_vof->get_function_values(
                        *this->multiphysics->get_filtered_solution(
                          PhysicsID::VOF),
                        filtered_phase_values);
                    }

                  // Loop on the quadrature points
                  for (unsigned int q = 0; q < n_q_points_face; q++)
                    {
                      if (properties_manager.get_number_of_fluids() == 2)
                        {
                          // Blend the physical properties using the VOF
                          // field
                          thermal_conductivity =
                            calculate_point_property(filtered_phase_values[q],
                                                     thermal_conductivity_0[q],
                                                     thermal_conductivity_1[q]);

                          density =
                            calculate_point_property(filtered_phase_values[q],
                                                     density_0[q],
                                                     density_1[q]);

                          specific_heat =
                            calculate_point_property(filtered_phase_values[q],
                                                     specific_heat_0[q],
                                                     specific_heat_1[q]);
                        }

                      Tensor<1, dim> normal_vector_ht =
                        -fe_face_values_ht.normal_vector(q);

                      heat_flux[boundary_index] +=
                        (-thermal_conductivity * temperature_gradient[q] *
                           normal_vector_ht +
                         local_temperature_values[q] * density * specific_heat *
                           local_velocity_values[q] * normal_vector_ht) *
                        fe_face_values_ht.JxW(q);

                      convective_flux[boundary_index] +=
                        h_coefficient * (local_temperature_values[q] - T_inf) *
                        fe_face_values_ht.JxW(q);

                    } // end loop on quadrature points
                }     // end loop on faces
            }         // end face is a boundary face
        }             // end condition cell at boundary
    }                 // end loop on cells


  // Sum across all cores
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_ht.type)
    {
      heat_flux[id] = Utilities::MPI::sum(heat_flux[id], mpi_communicator);

      convective_flux[id] =
        Utilities::MPI::sum(convective_flux[id], mpi_communicator);
    }

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout
        << "Total heat flux at the heat transfer boundary conditions : "
        << std::endl;
      for (auto const &[id, type] :
           this->simulation_parameters.boundary_conditions_ht.type)
        this->pcout << "\t boundary " << id << " : " << heat_flux.at(id)
                    << std::endl;


      this->pcout
        << "Convective heat flux at the heat transfer boundary conditions : "
        << std::endl;
      for (auto const &[id, type] :
           this->simulation_parameters.boundary_conditions_ht.type)
        this->pcout << "\t boundary " << id << " : " << convective_flux.at(id)
                    << std::endl;
    }

  // Filling table
  this->heat_flux_table.add_value("time",
                                  this->simulation_control->get_current_time());
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions_ht.type)
    {
      this->heat_flux_table.add_value("total_flux_bc_" +
                                        Utilities::int_to_string(id, 1),
                                      heat_flux.at(id));
      this->heat_flux_table.add_value("convective_flux_bc_" +
                                        Utilities::int_to_string(id, 1),
                                      convective_flux.at(id));
    }
}

template void
RANSTurbulence<2>::postprocess_heat_flux_on_bc<GlobalVectorType>(
  const bool              gather_vof,
  const GlobalVectorType &current_solution_fd);

template void
RANSTurbulence<3>::postprocess_heat_flux_on_bc<GlobalVectorType>(
  const bool              gather_vof,
  const GlobalVectorType &current_solution_fd);

template void
RANSTurbulence<2>::postprocess_heat_flux_on_bc<GlobalBlockVectorType>(
  const bool                   gather_vof,
  const GlobalBlockVectorType &current_solution_fd);

template void
RANSTurbulence<3>::postprocess_heat_flux_on_bc<GlobalBlockVectorType>(
  const bool                   gather_vof,
  const GlobalBlockVectorType &current_solution_fd);

template <int dim>
template <typename VectorType>
void
RANSTurbulence<dim>::postprocess_thermal_energy_in_fluid(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const VectorType                &current_solution_fd)
{
  const unsigned int n_q_points       = this->cell_quadrature->size();
  const MPI_Comm     mpi_communicator = this->dof_handler.get_communicator();

  // Initialize heat transfer information
  std::vector<double>         local_temperature_values(n_q_points);
  std::vector<Tensor<1, dim>> temperature_gradient(n_q_points);
  FEValues<dim>               fe_values_ht(*this->temperature_mapping,
                             this->dof_handler.get_fe(),
                             *this->cell_quadrature,
                             update_values | update_JxW_values);

  // Initialize fluid dynamics information
  const DoFHandler<dim> *dof_handler_fd =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      local_velocity_values(n_q_points);
  FEValues<dim>                    fe_values_fd(*this->temperature_mapping,
                             dof_handler_fd->get_fe(),
                             *this->cell_quadrature,
                             update_values);

  // Initialize VOF information
  const DoFHandler<dim>         *dof_handler_vof = nullptr;
  std::shared_ptr<FEValues<dim>> fe_values_vof;
  std::vector<double>            filtered_phase_values(n_q_points);

  if (gather_vof)
    {
      dof_handler_vof = this->multiphysics->get_dof_handler(PhysicsID::VOF);
      fe_values_vof =
        std::make_shared<FEValues<dim>>(*this->temperature_mapping,
                                        dof_handler_vof->get_fe(),
                                        *this->cell_quadrature,
                                        update_values);
    }

  // Initialize fluid properties
  auto &properties_manager =
    this->simulation_parameters.physical_properties_manager;
  std::map<field, double>              field_values;
  std::map<field, std::vector<double>> fields;

  // monophase flow
  double density(0.);
  double specific_heat(0.);

  // multiphase flow
  std::vector<double> density_0(n_q_points);
  std::vector<double> specific_heat_0(n_q_points);
  std::vector<double> density_1(n_q_points);
  std::vector<double> specific_heat_1(n_q_points);

  switch (properties_manager.get_number_of_fluids())
    {
      default:
        {
          // Get values for monophase flow
          const auto density_model = properties_manager.get_density();
          const auto specific_heat_model =
            properties_manager.get_specific_heat();

          density       = density_model->value(field_values);
          specific_heat = specific_heat_model->value(field_values);

          break;
        }
      case 2:
        {
          // Get prm values for multiphase flow - will be blended in the
          // integration loop
          const auto density_models = properties_manager.get_density_vector();
          const auto specific_heat_models =
            properties_manager.get_specific_heat_vector();

          density_models[0]->vector_value(fields, density_0);
          specific_heat_models[0]->vector_value(fields, specific_heat_0);

          density_models[1]->vector_value(fields, density_1);
          specific_heat_models[1]->vector_value(fields, specific_heat_1);

          break;
        }
    } // end switch on number of fluids

  // Other initializations
  double phase_coefficient(0.);
  double heat_in_domain(0.);

  // Integrate on all domain
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Gather heat transfer information
          fe_values_ht.reinit(cell);
          fe_values_ht.get_function_values(this->present_solution,
                                           local_temperature_values);

          // Get fluid dynamics active cell iterator
          typename DoFHandler<dim>::active_cell_iterator cell_fd(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            dof_handler_fd);

          // Gather fluid dynamics information
          fe_values_fd.reinit(cell_fd);
          fe_values_fd[velocities].get_function_values(current_solution_fd,
                                                       local_velocity_values);

          if (gather_vof)
            {
              // Get VOF active cell iterator
              typename DoFHandler<dim>::active_cell_iterator cell_vof(
                &(*(this->triangulation)),
                cell->level(),
                cell->index(),
                dof_handler_vof);

              // Gather VOF information
              fe_values_vof->reinit(cell_vof);
              fe_values_vof->get_function_values(
                *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
                filtered_phase_values);
            }

          // Loop on the quadrature points
          for (unsigned int q = 0; q < n_q_points; q++)
            {
              phase_coefficient =
                set_phase_coefficient(gather_vof,
                                      monitored_fluid,
                                      filtered_phase_values[q])
                  .first;

              if (properties_manager.get_number_of_fluids() == 2)
                {
                  // Blend the physical properties using the VOF
                  // field
                  density = calculate_point_property(filtered_phase_values[q],
                                                     density_0[q],
                                                     density_1[q]);

                  specific_heat =
                    calculate_point_property(filtered_phase_values[q],
                                             specific_heat_0[q],
                                             specific_heat_1[q]);
                }

              heat_in_domain +=
                (density * specific_heat * local_temperature_values[q]) *
                phase_coefficient * fe_values_ht.JxW(q);

            } // end loop on quadrature points
        }
    } // end loop on cells

  heat_in_domain = Utilities::MPI::sum(heat_in_domain, mpi_communicator);

  // Console output
  if (simulation_parameters.post_processing.verbosity ==
      Parameters::Verbosity::verbose)
    {
      this->pcout << "Thermal energy in " << domain_name << " : "
                  << heat_in_domain << std::endl;
    }

  // Filling table
  this->heat_flux_table.add_value("thermal_energy_" + domain_name,
                                  heat_in_domain);
}

template void
RANSTurbulence<2>::postprocess_thermal_energy_in_fluid<GlobalVectorType>(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const GlobalVectorType          &current_solution_fd);

template void
RANSTurbulence<3>::postprocess_thermal_energy_in_fluid<GlobalVectorType>(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const GlobalVectorType          &current_solution_fd);

template void
RANSTurbulence<2>::postprocess_thermal_energy_in_fluid<GlobalBlockVectorType>(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const GlobalBlockVectorType     &current_solution_fd);

template void
RANSTurbulence<3>::postprocess_thermal_energy_in_fluid<GlobalBlockVectorType>(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const std::string                domain_name,
  const GlobalBlockVectorType     &current_solution_fd);

template <int dim>
void
RANSTurbulence<dim>::write_heat_flux(const std::string domain_name)
{
  auto mpi_communicator = triangulation->get_communicator();

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.post_processing.heat_flux_output_name + "_" +
        domain_name + ".dat";
      std::ofstream output(filename.c_str());

      this->heat_flux_table.write_text(output);
    }
}

template <int dim>
std::pair<double, bool>
RANSTurbulence<dim>::set_phase_coefficient(
  const bool                       gather_vof,
  const Parameters::FluidIndicator monitored_fluid,
  const double                     phase_value_q)
{
  double phase_coefficient(0.);
  bool   point_is_in_postprocessed_fluid(false);

  switch (monitored_fluid)
    {
      default:
        {
          // Parameters::FluidIndicator::both
          phase_coefficient               = 1.;
          point_is_in_postprocessed_fluid = true;
          break;
        }
      case Parameters::FluidIndicator::fluid0:
        {
          if (gather_vof)
            {
              phase_coefficient = 1. - phase_value_q;
              if (phase_value_q < 0.5)
                point_is_in_postprocessed_fluid = true;
            }
          else
            {
              // In the case of monophase simulations, "fluid0" is
              // equivalent to "both" (all domain) calculation
              phase_coefficient               = 1.;
              point_is_in_postprocessed_fluid = true;
            }
          break;
        }
      case Parameters::FluidIndicator::fluid1:
        {
          if (gather_vof)
            {
              phase_coefficient = phase_value_q;
              if (phase_value_q > 0.5)
                point_is_in_postprocessed_fluid = true;
            }
          else
            {
              // Defensive only here, this has been checked already
              // in simulation_parameters.h
              throw std::logic_error(
                "Inconsistency in .prm!\n when VOF = false"
                "\n use (default value): set postprocessed fluid = both"
                "\n or: set postprocessed fluid = fluid 0");
            }
          break;
        }
    } // end switch on monitored fluid

  return std::make_pair(phase_coefficient, point_is_in_postprocessed_fluid);
}

template class RANSTurbulence<2>;
template class RANSTurbulence<3>;
