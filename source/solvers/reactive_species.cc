/*---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <core/bdf.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/reactive_species.h>
#include <solvers/reactive_species_assemblers.h>
#include <solvers/reactive_species_scratch_data.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>

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
ReactiveSpecies<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Time-stepping schemes
  if (is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.push_back(
        std::make_shared<ReactiveSpeciesAssemblerBDF<dim>>(
          this->simulation_control));
    }

  // For all the assemblers below, the parameter epsilon is passed to the
  // constructor explicitly because it might be calculated directly using
  // properties from the triangulation (e.g. minimal cell-size). Consequently,
  // it cannot be used directly from the simulation parameters and it is passed
  // to the constructor separately.

  // Core assembler
  this->assemblers.push_back(
    std::make_shared<ReactiveSpeciesAssemblerCore<dim>>(
      this->simulation_control,
      this->simulation_parameters.multiphysics.reactive_species_parameters));
}

template <int dim>
void
ReactiveSpecies<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");

  this->system_matrix = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = ReactiveSpeciesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &ReactiveSpecies::assemble_local_system_matrix,
                  &ReactiveSpecies::copy_local_matrix_to_global_matrix,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
ReactiveSpecies<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  ReactiveSpeciesScratchData<dim>                      &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!copy_data.cell_is_local)
    return;

  auto source_term = simulation_parameters.source_term.reactive_species_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_block_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
    }
  else
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
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
ReactiveSpecies<dim>::copy_local_matrix_to_global_matrix(
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
ReactiveSpecies<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_rhs = 0;
  setup_assemblers();

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  auto scratch_data = ReactiveSpeciesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    dof_handler_fluid->get_fe());

  WorkStream::run(this->dof_handler.begin_active(),
                  this->dof_handler.end(),
                  *this,
                  &ReactiveSpecies::assemble_local_system_rhs,
                  &ReactiveSpecies::copy_local_rhs_to_global_rhs,
                  scratch_data,
                  StabilizedMethodsCopyData(this->fe->n_dofs_per_cell(),
                                            this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
ReactiveSpecies<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  ReactiveSpeciesScratchData<dim>                      &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!copy_data.cell_is_local)
    return;

  auto source_term = simulation_parameters.source_term.reactive_species_source;
  source_term->set_time(simulation_control->get_current_time());

  scratch_data.reinit(cell,
                      this->evaluation_point,
                      this->previous_solutions,
                      &(*source_term));

  const DoFHandler<dim> *dof_handler_fluid =
    multiphysics->get_dof_handler(PhysicsID::fluid_dynamics);

  typename DoFHandler<dim>::active_cell_iterator velocity_cell(
    &(*triangulation), cell->level(), cell->index(), dof_handler_fluid);

  if (multiphysics->fluid_dynamics_is_block())
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_block_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
    }
  else
    {
      scratch_data.reinit_velocity(velocity_cell,
                                   *multiphysics->get_solution(
                                     PhysicsID::fluid_dynamics),
                                   this->simulation_parameters.ale);
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
ReactiveSpecies<dim>::copy_local_rhs_to_global_rhs(
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
ReactiveSpecies<dim>::attach_solution_to_output(DataOut<dim> &data_out)
{
  // Add the interpretation of the solution. The first component is the
  // phase order (Phi) and the following one is the chemical potential (eta)

  std::vector<std::string> solution_names;

  // TODO Change to flexible number of species
  unsigned int number_of_reactive_species = 4;
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      solution_names.push_back("species_" + std::to_string(i));
    }

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      // TODO Change to flexible number of species
      4,
      DataComponentInterpretation::component_is_scalar);

  data_out.add_data_vector(dof_handler,
                           present_solution,
                           solution_names,
                           data_component_interpretation);
}

template <int dim>
std::vector<double>
ReactiveSpecies<dim>::calculate_L2_error()
{
  auto mpi_communicator = triangulation->get_communicator();

  FEValues<dim> fe_values(*this->mapping,
                          *fe,
                          *this->cell_quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int n_q_points = this->cell_quadrature->size();

  // TODO Change to flexible number of species
  unsigned int                            number_of_reactive_species = 4;
  std::vector<FEValuesExtractors::Scalar> fe_values_extractors;
  std::vector<std::vector<double>>        local_species_values;
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      fe_values_extractors.emplace_back(i);
      local_species_values.emplace_back(std::vector<double>(n_q_points));
    }

  // TODO Change to flexible number of species
  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(4));

  auto &exact_solution =
    simulation_parameters.analytical_solution->reactive_species;
  exact_solution.set_time(simulation_control->get_current_time());

  // TODO Change to flexible number of species
  std::vector<double> l2_error_species(4);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          // TODO Change to flexible number of species
          for (unsigned int i = 0; i < number_of_reactive_species; i++)
            {
              fe_values[fe_values_extractors[i]].get_function_values(
                evaluation_point, local_species_values[i]);
            }

          // Get the exact solution at all Gauss points
          exact_solution.vector_value_list(fe_values.get_quadrature_points(),
                                           q_exactSol);

          // TODO Change to flexible number of species
          for (unsigned int i = 0; i < number_of_reactive_species; i++)
            {
              for (unsigned int q = 0; q < n_q_points; q++)
                {
                  // Find the values of x and u_h (the finite element solution)
                  // at the quadrature points
                  double species_sim   = local_species_values[i][q];
                  double species_exact = q_exactSol[q][i];

                  l2_error_species[i] += (species_sim - species_exact) *
                                         (species_sim - species_exact) *
                                         fe_values.JxW(q);
                }
            }
        }
    }

  // TODO Change to flexible number of species
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      l2_error_species[i] =
        Utilities::MPI::sum(l2_error_species[i], mpi_communicator);
      l2_error_species[i] = std::sqrt(l2_error_species[i]);
    }

  return l2_error_species;
}

template <int dim>
void
ReactiveSpecies<dim>::finish_simulation()
{
  auto         mpi_communicator = triangulation->get_communicator();
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

      // TODO Change to flexible number of species
      unsigned int number_of_reactive_species = 4;
      for (unsigned int i = 0; i < number_of_reactive_species; i++)
        {
          error_table.set_scientific("error_species_" + std::to_string(i),
                                     true);
          error_table.set_precision("error_species_" + std::to_string(i),
                                    simulation_control->get_log_precision());
        }
      error_table.write_text(std::cout);
    }
}

template <int dim>
void
ReactiveSpecies<dim>::percolate_time_vectors()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim>
void
ReactiveSpecies<dim>::modify_solution()
{
  // TODO Delete if unused
}

template <int dim>
void
ReactiveSpecies<dim>::postprocess(bool first_iteration)
{
  if (simulation_parameters.analytical_solution->calculate_error() == true &&
      !first_iteration)
    {
      std::vector<double> reactive_species_error = calculate_L2_error();

      error_table.add_value("cells",
                            this->triangulation->n_global_active_cells());

      // TODO Change to flexible number of species
      unsigned int number_of_reactive_species = 4;
      for (unsigned int i = 0; i < number_of_reactive_species; i++)
        {
          error_table.add_value("error_species_" + std::to_string(i),
                                reactive_species_error[i]);
          error_table.set_scientific("error_species_" + std::to_string(i),
                                     true);
        }

      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.analytical_solution->get_filename() +
        "_reactive_species.dat";
      std::ofstream output(filename.c_str());
      error_table.write_text(output);

      if (simulation_parameters.analytical_solution->verbosity ==
          Parameters::Verbosity::verbose)
        {
          for (unsigned int i = 0; i < number_of_reactive_species; i++)
            {
              this->pcout << "L2 error reactive species : "
                          << reactive_species_error[i] << std::endl;
            }
        }
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "Reactive Species");
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }
}


template <int dim>
void
ReactiveSpecies<dim>::pre_mesh_adaptation()
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
ReactiveSpecies<dim>::post_mesh_adaptation()
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
ReactiveSpecies<dim>::compute_kelly(
  const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                        &ivar,
  dealii::Vector<float> &estimated_error_per_cell)
{
  // TODO Change to flexible number of species
  unsigned int                            number_of_reactive_species = 4;
  std::vector<FEValuesExtractors::Scalar> fe_values_extractors;
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      fe_values_extractors.emplace_back(i);
      if (ivar.first == Variable::concentration_reactive_species)
        {
          KellyErrorEstimator<dim>::estimate(
            *this->mapping,
            this->dof_handler,
            *this->face_quadrature,
            typename std::map<types::boundary_id,
                              const Function<dim, double> *>(),
            present_solution,
            estimated_error_per_cell,
            this->fe->component_mask(fe_values_extractors[i]));
        }
    }
}


template <int dim>
void
ReactiveSpecies<dim>::write_checkpoint()
{
  std::vector<const GlobalVectorType *> sol_set_transfer;

  solution_transfer = std::make_shared<
    parallel::distributed::SolutionTransfer<dim, GlobalVectorType>>(
    dof_handler);

  sol_set_transfer.push_back(&present_solution);
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&previous_solutions[i]);
    }
  solution_transfer->prepare_for_serialization(sol_set_transfer);

  // Serialize tables
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (this->simulation_parameters.analytical_solution->calculate_error())
    serialize_table(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_RS" + suffix);
}

template <int dim>
void
ReactiveSpecies<dim>::read_checkpoint()
{
  auto mpi_communicator = triangulation->get_communicator();
  this->pcout << "Reading Reactive species checkpoint" << std::endl;

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

  // Deserialize tables
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (this->simulation_parameters.analytical_solution->calculate_error())
    deserialize_table(
      this->error_table,
      prefix + this->simulation_parameters.analytical_solution->get_filename() +
        "_RS" + suffix);
}


template <int dim>
void
ReactiveSpecies<dim>::setup_dofs()
{
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
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);

    for (unsigned int i_bc = 0;
         i_bc <
         this->simulation_parameters.boundary_conditions_reactive_species.size;
         ++i_bc)
      {
        // Dirichlet condition: imposed species value at i_bc
        if (this->simulation_parameters.boundary_conditions_reactive_species
              .type[i_bc] ==
            BoundaryConditions::BoundaryType::reactive_species_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_reactive_species
                .id[i_bc],
              ReactiveSpeciesFunctionDefined<dim>(
                &this->simulation_parameters
                   .boundary_conditions_reactive_species.bcFunctions[i_bc]
                   .dirichlet),
              nonzero_constraints);
          }
      }
  }
  nonzero_constraints.close();

  // Boundary conditions for Newton correction
  {
    zero_constraints.clear();
    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            zero_constraints);

    for (unsigned int i_bc = 0;
         i_bc <
         this->simulation_parameters.boundary_conditions_reactive_species.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions_reactive_species
              .type[i_bc] ==
            BoundaryConditions::BoundaryType::reactive_species_dirichlet)
          {
            VectorTools::interpolate_boundary_values(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions_reactive_species
                .id[i_bc],
              Functions::ZeroFunction<dim>(4),
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

  this->pcout << "   Number of Reactive physics degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the reactive_species dof_handler and present solution pointers to
  // the multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::reactive_species,
                                &this->dof_handler);
  multiphysics->set_solution(PhysicsID::reactive_species,
                             &this->present_solution);
  multiphysics->set_previous_solutions(PhysicsID::reactive_species,
                                       &this->previous_solutions);
}

template <int dim>
void
ReactiveSpecies<dim>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions_reactive_species
         .time_dependent)
    return;

  double time = this->simulation_control->get_current_time();
  nonzero_constraints.clear();
  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          nonzero_constraints);

  for (unsigned int i_bc = 0;
       i_bc <
       this->simulation_parameters.boundary_conditions_reactive_species.size;
       ++i_bc)
    {
      // Dirichlet condition: imposed species value at i_bc
      if (this->simulation_parameters.boundary_conditions_reactive_species
            .type[i_bc] ==
          BoundaryConditions::BoundaryType::reactive_species_dirichlet)
        {
          this->simulation_parameters.boundary_conditions_reactive_species
            .bcFunctions[i_bc]
            .dirichlet.set_time(time);
          VectorTools::interpolate_boundary_values(
            this->dof_handler,
            this->simulation_parameters.boundary_conditions_reactive_species
              .id[i_bc],
            ReactiveSpeciesFunctionDefined<dim>(
              &this->simulation_parameters.boundary_conditions_reactive_species
                 .bcFunctions[i_bc]
                 .dirichlet),
            nonzero_constraints);
        }
    }

  nonzero_constraints.close();
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim>
void
ReactiveSpecies<dim>::set_initial_conditions()
{
  // TODO Change to flexible number of species
  unsigned int                            number_of_reactive_species = 4;
  std::vector<FEValuesExtractors::Scalar> fe_values_extractors;
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      fe_values_extractors.emplace_back(i);

      VectorTools::interpolate(
        *mapping,
        dof_handler,
        simulation_parameters.initial_condition->reactive_species,
        newton_update,
        fe->component_mask(fe_values_extractors[i]));
    }

  nonzero_constraints.distribute(newton_update);
  present_solution = newton_update;
  percolate_time_vectors();
}

template <int dim>
void
ReactiveSpecies<dim>::solve_linear_system(const bool initial_step,
                                          const bool /*renewed_matrix*/)
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  auto mpi_communicator = triangulation->get_communicator();

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double absolute_residual =
    simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
      .minimum_residual;
  const double relative_residual =
    simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
      .relative_residual;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  const double ilu_fill =
    simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
      .ilu_precond_fill;
  const double ilu_atol =
    simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
      .ilu_precond_atol;
  const double ilu_rtol =
    simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);

  TrilinosWrappers::PreconditionILU ilu_preconditioner;

  ilu_preconditioner.initialize(system_matrix, preconditionerOptions);

  GlobalVectorType completely_distributed_solution(locally_owned_dofs,
                                                   mpi_communicator);

  SolverControl solver_control(simulation_parameters.linear_solver
                                 .at(PhysicsID::reactive_species)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    false,
    simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
      .max_krylov_vectors);


  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);


  solver.solve(system_matrix,
               completely_distributed_solution,
               system_rhs,
               ilu_preconditioner);

  if (simulation_parameters.linear_solver.at(PhysicsID::reactive_species)
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
ReactiveSpecies<dim>::output_newton_update_norms(
  const unsigned int display_precision)
{
  auto mpi_communicator = triangulation->get_communicator();

  // TODO Change to flexible number of species
  unsigned int                            number_of_reactive_species = 4;
  std::vector<FEValuesExtractors::Scalar> fe_values_extractors;
  std::vector<ComponentMask>              component_mask_species;
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      fe_values_extractors.emplace_back(i);
      component_mask_species.emplace_back(
        fe->component_mask(fe_values_extractors[i]));
    }

  // TODO Change to flexible number of species
  std::vector<double> local_sum(4);
  std::vector<double> local_max(4);
  std::vector<double> global_l2_norm(4);
  std::vector<double> global_linfty_norm(4);
  for (unsigned int i = 0; i < number_of_reactive_species; i++)
    {
      local_max[i] = std::numeric_limits<double>::lowest();
      const std::vector<IndexSet> index_set_species =
        DoFTools::locally_owned_dofs_per_component(dof_handler,
                                                   component_mask_species[i]);

      for (auto j = index_set_species[i].begin();
           j != index_set_species[i].end();
           j++)
        {
          double dof_newton_update = newton_update[*j];

          local_sum[i] += dof_newton_update * dof_newton_update;

          local_max[i] = std::max(local_max[i], std::abs(dof_newton_update));
        }

      global_l2_norm[i] =
        std::sqrt(Utilities::MPI::sum(local_sum[i], mpi_communicator));
      global_linfty_norm[i] =
        Utilities::MPI::max(local_max[i], mpi_communicator);

      this->pcout << std::setprecision(display_precision)
                  << "\n\t||dphi||_L2 = " << std::setw(6) << global_l2_norm[i]
                  << std::setw(6) << "\t||dphi||_Linfty = "
                  << std::setprecision(display_precision)
                  << global_linfty_norm[i] << std::endl;
    }
}



template class ReactiveSpecies<2>;
template class ReactiveSpecies<3>;
