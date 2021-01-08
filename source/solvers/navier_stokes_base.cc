/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out_faces.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <core/grids.h>
#include <core/solutions_output.h>
#include <core/utilities.h>
#include <solvers/flow_control.h>
#include <solvers/navier_stokes_base.h>
#include <solvers/post_processors.h>
#include <solvers/postprocessing_cfd.h>
#include <solvers/postprocessing_velocities.h>

#include "core/time_integration_utilities.h"


/*
 * Constructor for the Navier-Stokes base class
 */
template <int dim, typename VectorType, typename DofsType>
NavierStokesBase<dim, VectorType, DofsType>::NavierStokesBase(
  NavierStokesSolverParameters<dim> &p_nsparam)
  : PhysicsSolver<VectorType>(p_nsparam.non_linear_solver,
                              p_nsparam.multiphysics.number_physic_total)
  , mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , triangulation(dynamic_cast<parallel::DistributedTriangulationBase<dim> *>(
      new parallel::distributed::Triangulation<dim>(
        this->mpi_communicator,
        typename Triangulation<dim>::MeshSmoothing(
          Triangulation<dim>::smoothing_on_refinement |
          Triangulation<dim>::smoothing_on_coarsening))))
  , dof_handler(*this->triangulation)
  , fe(FE_Q<dim>(p_nsparam.fem_parameters.velocity_order),
       dim,
       FE_Q<dim>(p_nsparam.fem_parameters.pressure_order),
       1)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , nsparam(p_nsparam)
  , flow_control(nsparam.flow_control)
  , velocity_fem_degree(p_nsparam.fem_parameters.velocity_order)
  , pressure_fem_degree(p_nsparam.fem_parameters.pressure_order)
  , number_quadrature_points(p_nsparam.fem_parameters.velocity_order + 1)
{
  this->pcout.set_condition(
    Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0);

  if (nsparam.simulation_control.method ==
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    simulation_control =
      std::make_shared<SimulationControlSteady>(nsparam.simulation_control);
  else if (nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    {
      simulation_control = std::make_shared<SimulationControlAdjointSteady>(
        nsparam.simulation_control);
    }
  else
    {
      if (nsparam.simulation_control.output_control ==
          Parameters::SimulationControl::OutputControl::time)
        simulation_control =
          std::make_shared<SimulationControlTransientDynamicOutput>(
            nsparam.simulation_control);
      else
        simulation_control = std::make_shared<SimulationControlTransient>(
          nsparam.simulation_control);
    }

  // Change the behavior of the timer for situations when you don't want outputs
  if (nsparam.timer.type == Parameters::Timer::Type::none)
    this->computing_timer.disable_output();

  // Pre-allocate the force tables to match the number of boundary conditions
  forces_on_boundaries.resize(nsparam.boundary_conditions.size);
  torques_on_boundaries.resize(nsparam.boundary_conditions.size);
  forces_tables.resize(nsparam.boundary_conditions.size);
  torques_tables.resize(nsparam.boundary_conditions.size);

  // Get the exact solution from the parser
  exact_solution = &nsparam.analytical_solution->velocity;

  // If there is a forcing function, get it from the parser
  if (nsparam.sourceTerm->source_term())
    forcing_function = &nsparam.sourceTerm->navier_stokes_source;
  else
    forcing_function = new NoForce<dim>;

  this->pcout << "Running on "
              << Utilities::MPI::n_mpi_processes(this->mpi_communicator)
              << " MPI rank(s)..." << std::endl;

  this->pcout << std::setprecision(nsparam.simulation_control.log_precision);
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_flow_rate(
  const VectorType &evaluation_point)
{
  this->flow_rate = calculate_flow_rate(this->dof_handler,
                                        evaluation_point,
                                        nsparam.flow_control.boundary_flow_id,
                                        nsparam.fem_parameters,
                                        mpi_communicator);

  // Showing results (area and flow rate)
  if (nsparam.flow_control.verbosity == Parameters::Verbosity::verbose &&
      simulation_control->get_step_number() > 0 && this->this_mpi_process == 0)
    {
      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Flow control summary                    |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      this->pcout << "Inlet area : " << flow_rate.second << std::endl;
      this->pcout << "Flow rate : " << flow_rate.first << std::endl;
      this->pcout << "Beta applied : "
                  << beta[nsparam.flow_control.flow_direction] << std::endl;
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::dynamic_flow_control()
{
  if (nsparam.flow_control.enable_flow_control &&
      nsparam.simulation_control.method !=
        Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      this->flow_control.calculate_beta(flow_rate,
                                        simulation_control->get_time_step(),
                                        simulation_control->get_step_number());

      this->beta = flow_control.get_beta();
    }
}


// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_forces(
  const VectorType &evaluation_point)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_forces");

  this->forces_on_boundaries = calculate_forces(this->dof_handler,
                                                evaluation_point,
                                                nsparam.physical_properties,
                                                nsparam.boundary_conditions,
                                                mpi_communicator);

  if (nsparam.forces_parameters.verbosity == Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      std::cout << std::endl;
      std::string independent_column_names = "Boundary ID";

      std::vector<std::string> dependent_column_names;
      dependent_column_names.push_back("f_x");
      dependent_column_names.push_back("f_y");
      if (dim == 3)
        dependent_column_names.push_back("f_z");

      TableHandler table = make_table_scalars_tensors(
        nsparam.boundary_conditions.id,
        independent_column_names,
        this->forces_on_boundaries,
        dependent_column_names,
        this->nsparam.simulation_control.log_precision);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Force  summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (unsigned int i_boundary = 0;
       i_boundary < nsparam.boundary_conditions.size;
       ++i_boundary)
    {
      this->forces_tables[i_boundary].add_value(
        "time", simulation_control->get_current_time());
      this->forces_tables[i_boundary].add_value(
        "f_x", this->forces_on_boundaries[i_boundary][0]);
      this->forces_tables[i_boundary].add_value(
        "f_y", this->forces_on_boundaries[i_boundary][1]);
      if (dim == 3)
        this->forces_tables[i_boundary].add_value(
          "f_z", this->forces_on_boundaries[i_boundary][2]);
      else
        this->forces_tables[i_boundary].add_value("f_z", 0.);

      // Precision
      this->forces_tables[i_boundary].set_precision(
        "f_x", nsparam.forces_parameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "f_y", nsparam.forces_parameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "f_z", nsparam.forces_parameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "time", nsparam.forces_parameters.output_precision);
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_torques(
  const VectorType &evaluation_point)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_torques");

  this->torques_on_boundaries = calculate_torques(this->dof_handler,
                                                  evaluation_point,
                                                  nsparam.physical_properties,
                                                  nsparam.fem_parameters,
                                                  nsparam.boundary_conditions,
                                                  mpi_communicator);

  if (nsparam.forces_parameters.verbosity == Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      this->pcout << std::endl;
      std::string independent_column_names = "Boundary ID";

      std::vector<std::string> dependent_column_names;
      dependent_column_names.push_back("T_x");
      dependent_column_names.push_back("T_y");
      dependent_column_names.push_back("T_z");

      TableHandler table = make_table_scalars_tensors(
        nsparam.boundary_conditions.id,
        independent_column_names,
        this->torques_on_boundaries,
        dependent_column_names,
        this->nsparam.simulation_control.log_precision);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Torque summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (unsigned int boundary_id = 0;
       boundary_id < nsparam.boundary_conditions.size;
       ++boundary_id)
    {
      this->torques_tables[boundary_id].add_value(
        "time", simulation_control->get_current_time());
      this->torques_tables[boundary_id].add_value(
        "T_x", this->torques_on_boundaries[boundary_id][0]);
      this->torques_tables[boundary_id].add_value(
        "T_y", this->torques_on_boundaries[boundary_id][1]);
      this->torques_tables[boundary_id].add_value(
        "T_z", this->torques_on_boundaries[boundary_id][2]);

      // Precision
      this->torques_tables[boundary_id].set_precision(
        "T_x", nsparam.forces_parameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "T_y", nsparam.forces_parameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "T_z", nsparam.forces_parameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "time", nsparam.forces_parameters.output_precision);
    }
}



template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::finish_simulation()
{
  if (nsparam.forces_parameters.calculate_force)
    this->write_output_forces();

  if (nsparam.forces_parameters.calculate_torque)
    this->write_output_torques();

  if (nsparam.analytical_solution->calculate_error())
    {
      if (nsparam.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          error_table.set_scientific("error_pressure", true);
          error_table.omit_column_from_convergence_rate_evaluation("cells");
          error_table.omit_column_from_convergence_rate_evaluation(
            "total_time");
          error_table.evaluate_all_convergence_rates(
            ConvergenceTable::reduction_rate_log2);
        }
      error_table.set_scientific("error_velocity", true);

      if (this->this_mpi_process == 0)
        {
          std::string filename =
            nsparam.analytical_solution->get_filename() + ".dat";
          std::ofstream output(filename.c_str());
          error_table.write_text(output);
          std::vector<std::string> sub_columns;
          if (nsparam.simulation_control.method ==
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              sub_columns.push_back("cells");
              sub_columns.push_back("error_velocity");
              sub_columns.push_back("error_pressure");
              error_table.set_column_order(sub_columns);
            }
          error_table.write_text(std::cout);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::finish_time_step()
{
  if (nsparam.simulation_control.method !=
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      this->solution_m3      = this->solution_m2;
      this->solution_m2      = this->solution_m1;
      auto &present_solution = this->get_present_solution();
      this->solution_m1      = present_solution;
      const double CFL       = calculate_CFL(this->dof_handler,
                                       present_solution,
                                       simulation_control->get_time_step(),
                                       mpi_communicator);
      this->simulation_control->set_CFL(CFL);
    }
  if (this->nsparam.restart_parameters.checkpoint &&
      simulation_control->get_step_number() %
          this->nsparam.restart_parameters.frequency ==
        0)
    {
      this->write_checkpoint();
    }

  if (this->nsparam.timer.type == Parameters::Timer::Type::iteration)
    {
      this->computing_timer.print_summary();
      this->computing_timer.reset();
    }
}

// Do an iteration with the NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial condition
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::iterate()
{
  auto &present_solution = this->get_present_solution();
  if (nsparam.simulation_control.method ==
      Parameters::SimulationControl::TimeSteppingMethod::sdirk22)
    {
      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1,
        false,
        false);
      this->solution_m2 = present_solution;

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2,
        false,
        false);
    }

  else if (nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::sdirk33)
    {
      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1,
        false,
        false);

      this->solution_m2 = present_solution;

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2,
        false,
        false);

      this->solution_m3 = present_solution;

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3,
        false,
        false);
    }
  else
    {
      PhysicsSolver<VectorType>::solve_non_linear_system(
        nsparam.simulation_control.method, false, false);
    }
}

// Do an iteration with the NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial condition
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::first_iteration()
{
  auto &present_solution = this->get_present_solution();

  // First step if the method is not a multi-step method
  if (!is_bdf_high_order(nsparam.simulation_control.method))
    {
      iterate();
    }

  // Taking care of the multi-step methods
  else if (nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    {
      Parameters::SimulationControl timeParameters = nsparam.simulation_control;

      // Start the BDF2 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;
      simulation_control->set_current_time_step(time_step);
      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = present_solution;

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps to complete the full step

      time_step =
        timeParameters.dt * (1. - timeParameters.startup_timestep_scaling);

      simulation_control->set_current_time_step(time_step);

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf2, false, true);

      simulation_control->set_suggested_time_step(timeParameters.dt);
    }

  else if (nsparam.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    {
      Parameters::SimulationControl timeParameters = nsparam.simulation_control;

      // Start the BDF3 with a single Euler time step with a lower time step
      double time_step =
        timeParameters.dt * timeParameters.startup_timestep_scaling;

      simulation_control->set_current_time_step(time_step);

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = present_solution;

      // Reset the time step and do a bdf 2 newton iteration using the two
      // steps

      simulation_control->set_current_time_step(time_step);

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf1, false, true);
      this->solution_m3 = this->solution_m2;
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = present_solution;

      // Reset the time step and do a bdf 3 newton iteration using the two
      // steps to complete the full step
      time_step =
        timeParameters.dt * (1. - 2. * timeParameters.startup_timestep_scaling);
      simulation_control->set_current_time_step(time_step);

      PhysicsSolver<VectorType>::solve_non_linear_system(
        Parameters::SimulationControl::TimeSteppingMethod::bdf3, false, true);
      simulation_control->set_suggested_time_step(timeParameters.dt);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh()
{
  if (simulation_control->get_step_number() %
        this->nsparam.mesh_adaptation.frequency ==
      0)
    {
      if (this->nsparam.mesh_adaptation.type ==
          Parameters::MeshAdaptation::Type::kelly)
        refine_mesh_kelly();

      else if (this->nsparam.mesh_adaptation.type ==
               Parameters::MeshAdaptation::Type::uniform)
        refine_mesh_uniform();
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh_kelly()
{
  if (dynamic_cast<parallel::distributed::Triangulation<dim> *>(
        this->triangulation.get()) == nullptr)
    return;

  auto &tria = *dynamic_cast<parallel::distributed::Triangulation<dim> *>(
    this->triangulation.get());

  // Time monitoring
  TimerOutput::Scope t(this->computing_timer, "refine");

  Vector<float>       estimated_error_per_cell(tria.n_active_cells());
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);
  auto &present_solution = this->get_present_solution();
  if (this->nsparam.mesh_adaptation.variable ==
      Parameters::MeshAdaptation::Variable::pressure)
    {
      KellyErrorEstimator<dim>::estimate(
        mapping,
        this->dof_handler,
        QGauss<dim - 1>(this->number_quadrature_points + 1),
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        this->fe.component_mask(pressure));
    }
  else if (this->nsparam.mesh_adaptation.variable ==
           Parameters::MeshAdaptation::Variable::velocity)
    {
      KellyErrorEstimator<dim>::estimate(
        mapping,
        this->dof_handler,
        QGauss<dim - 1>(this->number_quadrature_points + 1),
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        this->fe.component_mask(velocity));
    }

  if (this->nsparam.mesh_adaptation.fractionType ==
      Parameters::MeshAdaptation::FractionType::number)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      tria,
      estimated_error_per_cell,
      this->nsparam.mesh_adaptation.refinement_fraction,
      this->nsparam.mesh_adaptation.coarsening_fraction,
      this->nsparam.mesh_adaptation.maximum_number_elements);

  else if (this->nsparam.mesh_adaptation.fractionType ==
           Parameters::MeshAdaptation::FractionType::fraction)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      tria,
      estimated_error_per_cell,
      this->nsparam.mesh_adaptation.refinement_fraction,
      this->nsparam.mesh_adaptation.coarsening_fraction);

  if (tria.n_levels() > this->nsparam.mesh_adaptation.maximum_refinement_level)
    for (typename Triangulation<dim>::active_cell_iterator cell =
           tria.begin_active(
             this->nsparam.mesh_adaptation.maximum_refinement_level);
         cell != tria.end();
         ++cell)
      cell->clear_refine_flag();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active(
           this->nsparam.mesh_adaptation.minimum_refinement_level);
       cell !=
       tria.end_active(this->nsparam.mesh_adaptation.minimum_refinement_level);
       ++cell)
    cell->clear_coarsen_flag();

  tria.prepare_coarsening_and_refinement();

  // Solution transfer objects for all the solutions
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m1(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m2(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m3(
    this->dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(this->solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(this->solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(this->solution_m3);

  tria.execute_coarsening_and_refinement();
  this->setup_dofs();

  // Set up the vectors for the transfer
  VectorType tmp(locally_owned_dofs, this->mpi_communicator);
  VectorType tmp_m1(locally_owned_dofs, this->mpi_communicator);
  VectorType tmp_m2(locally_owned_dofs, this->mpi_communicator);
  VectorType tmp_m3(locally_owned_dofs, this->mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);
  solution_transfer_m1.interpolate(tmp_m1);
  solution_transfer_m2.interpolate(tmp_m2);
  solution_transfer_m3.interpolate(tmp_m3);

  // Distribute constraints
  auto &nonzero_constraints = this->get_nonzero_constraints();
  nonzero_constraints.distribute(tmp);
  nonzero_constraints.distribute(tmp_m1);
  nonzero_constraints.distribute(tmp_m2);
  nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  present_solution  = tmp;
  this->solution_m1 = tmp_m1;
  this->solution_m2 = tmp_m2;
  this->solution_m3 = tmp_m3;
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh_uniform()
{
  TimerOutput::Scope t(this->computing_timer, "refine");

  // Solution transfer objects for all the solutions
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m1(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m2(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m3(
    this->dof_handler);
  auto &present_solution = this->get_present_solution();
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(this->solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(this->solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(this->solution_m3);

  // Refine
  this->triangulation->refine_global(1);

  setup_dofs();

  // Set up the vectors for the transfer
  VectorType tmp(locally_owned_dofs, this->mpi_communicator);
  VectorType tmp_m1(locally_owned_dofs, this->mpi_communicator);
  VectorType tmp_m2(locally_owned_dofs, this->mpi_communicator);
  VectorType tmp_m3(locally_owned_dofs, this->mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);
  solution_transfer_m1.interpolate(tmp_m1);
  solution_transfer_m2.interpolate(tmp_m2);
  solution_transfer_m3.interpolate(tmp_m3);

  // Distribute constraints
  auto &nonzero_constraints = this->get_nonzero_constraints();
  nonzero_constraints.distribute(tmp);
  nonzero_constraints.distribute(tmp_m1);
  nonzero_constraints.distribute(tmp_m2);
  nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  present_solution  = tmp;
  this->solution_m1 = tmp_m1;
  this->solution_m2 = tmp_m2;
  this->solution_m3 = tmp_m3;
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocess(bool firstIter)
{
  auto &present_solution = this->get_present_solution();

  if (this->nsparam.post_processing.calculate_enstrophy)
    {
      double enstrophy = calculate_enstrophy(this->dof_handler,
                                             present_solution,
                                             nsparam.fem_parameters,
                                             mpi_communicator);

      this->enstrophy_table.add_value("time",
                                      simulation_control->get_current_time());
      this->enstrophy_table.add_value("enstrophy", enstrophy);

      // Display Enstrophy to screen if verbosity is enabled
      if (this->nsparam.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Enstrophy  : " << enstrophy << std::endl;
        }

      // Output Enstrophy to a text file from processor 0
      if (simulation_control->get_step_number() %
              this->nsparam.post_processing.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            nsparam.post_processing.enstrophy_output_name + ".dat";
          std::ofstream output(filename.c_str());
          enstrophy_table.set_precision("time", 12);
          enstrophy_table.set_precision("enstrophy", 12);
          this->enstrophy_table.write_text(output);
        }
    }

  // The average velocities and reynolds stresses are calculated when the
  // time reaches the initial time. (time >= initial time) with 1e-6 as
  // tolerance.
  if (this->nsparam.post_processing.calculate_average_velocities)
    {
      // Calculate average velocities when the time reaches the initial time.
      // time >= initial time with the epsilon as tolerance.
      const double dt = simulation_control->get_time_step();
      if (simulation_control->get_current_time() >
          (nsparam.post_processing.initial_time - 1e-6 * dt))
        {
          this->average_velocities.calculate_average_velocities(
            this->get_local_evaluation_point(),
            nsparam.post_processing,
            simulation_control->get_current_time(),
            dt);
        }
    }

  if (this->nsparam.post_processing.calculate_kinetic_energy)
    {
      TimerOutput::Scope t(this->computing_timer, "kinetic_energy_calculation");
      double             kE = calculate_kinetic_energy(this->dof_handler,
                                           present_solution,
                                           nsparam.fem_parameters,
                                           mpi_communicator);
      this->kinetic_energy_table.add_value(
        "time", simulation_control->get_current_time());
      this->kinetic_energy_table.add_value("kinetic-energy", kE);
      if (this->nsparam.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Kinetic energy : " << kE << std::endl;
        }

      // Output Kinetic Energy to a text file from processor 0
      if ((simulation_control->get_step_number() %
             this->nsparam.post_processing.output_frequency ==
           0) &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            nsparam.post_processing.kinetic_energy_output_name + ".dat";
          std::ofstream output(filename.c_str());
          kinetic_energy_table.set_precision("time", 12);
          kinetic_energy_table.set_precision("kinetic-energy", 12);
          this->kinetic_energy_table.write_text(output);
        }
    }

  // Calculate inlet flow rate and area
  if (this->nsparam.flow_control.enable_flow_control)
    {
      this->postprocessing_flow_rate(this->get_present_solution());
    }

  if (!firstIter)
    {
      // Calculate forces on the boundary conditions
      if (this->nsparam.forces_parameters.calculate_force)
        {
          if (simulation_control->get_step_number() %
                this->nsparam.forces_parameters.calculation_frequency ==
              0)
            this->postprocessing_forces(present_solution);
          if (simulation_control->get_step_number() %
                this->nsparam.forces_parameters.output_frequency ==
              0)
            this->write_output_forces();
        }

      // Calculate torques on the boundary conditions
      if (this->nsparam.forces_parameters.calculate_torque)
        {
          if (simulation_control->get_step_number() %
                this->nsparam.forces_parameters.calculation_frequency ==
              0)
            this->postprocessing_torques(present_solution);
          if (simulation_control->get_step_number() %
                this->nsparam.forces_parameters.output_frequency ==
              0)
            this->write_output_torques();
        }

      // Calculate error with respect to analytical solution
      if (this->nsparam.analytical_solution->calculate_error())
        {
          // Update the time of the exact solution to the actual time
          this->exact_solution->set_time(
            simulation_control->get_current_time());
          const std::pair<double, double> errors =
            calculate_L2_error(dof_handler,
                               present_solution,
                               exact_solution,
                               nsparam.fem_parameters,
                               mpi_communicator);
          const double error_velocity = errors.first;
          const double error_pressure = errors.second;
          if (nsparam.simulation_control.method ==
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              this->error_table.add_value(
                "cells", this->triangulation->n_global_active_cells());
              this->error_table.add_value("error_velocity", error_velocity);
              this->error_table.add_value("error_pressure", error_pressure);
              auto summary = computing_timer.get_summary_data(
                computing_timer.total_wall_time);
              double total_time = 0;
              for (auto it = summary.begin(); it != summary.end(); ++it)
                {
                  total_time += summary[it->first];
                }

              this->error_table.add_value("total_time", total_time);
            }
          else
            {
              this->error_table.add_value(
                "time", simulation_control->get_current_time());
              this->error_table.add_value("error_velocity", error_velocity);
            }
          if (this->nsparam.analytical_solution->verbosity ==
              Parameters::Verbosity::verbose)
            {
              this->pcout << "L2 error velocity : " << error_velocity
                          << std::endl;
            }
        }
    }
  if (this->simulation_control->is_output_iteration())
    this->write_output_results(present_solution);
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::set_nodal_values()
{
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  const MappingQ<dim>              mapping(this->velocity_fem_degree,
                              this->nsparam.fem_parameters.qmapping_all);
  auto &                           newton_update = this->get_newton_update();
  VectorTools::interpolate(mapping,
                           this->dof_handler,
                           this->nsparam.initial_condition->uvwp,
                           newton_update,
                           this->fe.component_mask(velocities));
  VectorTools::interpolate(mapping,
                           this->dof_handler,
                           this->nsparam.initial_condition->uvwp,
                           newton_update,
                           this->fe.component_mask(pressure));
  auto &nonzero_constraints = this->get_nonzero_constraints();
  nonzero_constraints.distribute(newton_update);
  auto &present_solution = this->get_present_solution();
  present_solution       = newton_update;
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "read_checkpoint");
  std::string        prefix = this->nsparam.restart_parameters.filename;
  this->simulation_control->read(prefix);
  this->pvdhandler.read(prefix);

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  try
    {
      if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim> *>(
            this->triangulation.get()))
        tria->load(filename.c_str());
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }
  setup_dofs();
  std::vector<VectorType *> x_system(4);



  /**
   * The newton_update vector is used to initialize the x_system to which
   * the serialized restart data is read. This is a necessary step since the
   * serialized data must be read into a locally_owned vector. The
   * present_solution, solution_m1, solution_m2 and solution_m3 vectors are
   * locally_relevant vector. Consequently, they cannot be written directly to
   * and these intermediary vectors must be used to store the data.
   */
  VectorType distributed_system(locally_owned_dofs, this->mpi_communicator);
  VectorType distributed_system_m1(locally_owned_dofs, this->mpi_communicator);
  VectorType distributed_system_m2(locally_owned_dofs, this->mpi_communicator);
  VectorType distributed_system_m3(locally_owned_dofs, this->mpi_communicator);
  x_system[0] = &(distributed_system);
  x_system[1] = &(distributed_system_m1);
  x_system[2] = &(distributed_system_m2);
  x_system[3] = &(distributed_system_m3);
  parallel::distributed::SolutionTransfer<dim, VectorType> system_trans_vectors(
    this->dof_handler);

  if (nsparam.post_processing.calculate_average_velocities)
    {
      std::vector<VectorType *> sum_vectors =
        this->average_velocities.read(prefix);

      x_system.insert(x_system.end(), sum_vectors.begin(), sum_vectors.end());
    }

  system_trans_vectors.deserialize(x_system);
  auto &present_solution = this->get_present_solution();
  present_solution       = distributed_system;
  this->solution_m1      = distributed_system_m1;
  this->solution_m2      = distributed_system_m2;
  this->solution_m3      = distributed_system_m3;

  if (nsparam.flow_control.enable_flow_control)
    {
      this->flow_control.read(prefix);

      this->flow_rate =
        calculate_flow_rate(this->dof_handler,
                            present_solution,
                            nsparam.flow_control.boundary_flow_id,
                            nsparam.fem_parameters,
                            mpi_communicator);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::output_field_hook(
  DataOut<dim> & /*data_out*/)
{}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_results(
  const VectorType &solution)
{
  TimerOutput::Scope  t(this->computing_timer, "output");
  const MappingQ<dim> mapping(this->velocity_fem_degree,
                              nsparam.fem_parameters.qmapping_all);

  const std::string  folder        = simulation_control->get_output_path();
  const std::string  solution_name = simulation_control->get_output_name();
  const unsigned int iter          = simulation_control->get_step_number();
  const double       time          = simulation_control->get_current_time();
  const unsigned int subdivision = simulation_control->get_number_subdivision();
  const unsigned int group_files = simulation_control->get_group_files();

  // Add the interpretation of the solution. The dim first components are the
  // velocity vectors and the following one is the pressure.
  std::vector<std::string> solution_names(dim, "velocity");
  solution_names.push_back("pressure");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);


  DataOut<dim> data_out;

  // Additional flag to enable the output of high-order elements
  DataOutBase::VtkFlags flags;
  if (this->velocity_fem_degree > 1)
    flags.write_higher_order_cells = true;
  data_out.set_flags(flags);

  // Attach the solution data to data_out object
  data_out.attach_dof_handler(this->dof_handler);
  data_out.add_data_vector(solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);

  if (this->nsparam.post_processing.calculate_average_velocities)
    {
      // Add the interpretation of the average solution. The dim first
      // components are the average velocity vectors and the following one is
      // the average pressure. (<u>, <v>, <w>, <p>)
      std::vector<std::string> average_solution_names(dim, "average_velocity");
      average_solution_names.push_back("average_pressure");

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        average_data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
      average_data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector(
        this->average_velocities.get_average_velocities(),
        average_solution_names,
        DataOut<dim>::type_dof_data,
        average_data_component_interpretation);

      // Add the interpretation of the reynolds stresses of solution.
      // The dim first components are the normal reynolds stress vectors and
      // the following ones are others resolved reynolds stresses.
      std::vector<std::string> reynolds_normal_stress_names(
        dim, "reynolds_normal_stress");
      reynolds_normal_stress_names.push_back("turbulent_kinetic_energy");
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        reynolds_normal_stress_data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
      reynolds_normal_stress_data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);


      std::vector<std::string> reynolds_shear_stress_names = {
        "reynolds_shear_stress_uv"};
      if (dim == 2)
        {
          reynolds_shear_stress_names.push_back("dummy_rss_2d");
        }
      if (dim == 3)
        {
          reynolds_shear_stress_names.push_back("reynolds_shear_stress_vw");
          reynolds_shear_stress_names.push_back("reynolds_shear_stress_uw");
        }
      reynolds_shear_stress_names.push_back("dummy_rss");

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        reynolds_shear_stress_data_component_interpretation(
          dim, DataComponentInterpretation::component_is_scalar);
      reynolds_shear_stress_data_component_interpretation.push_back(
        DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector(
        this->average_velocities.get_reynolds_normal_stresses(),
        reynolds_normal_stress_names,
        DataOut<dim>::type_dof_data,
        reynolds_normal_stress_data_component_interpretation);

      data_out.add_data_vector(
        this->average_velocities.get_reynolds_shear_stresses(),
        reynolds_shear_stress_names,
        DataOut<dim>::type_dof_data,
        reynolds_shear_stress_data_component_interpretation);
    }

  Vector<float> subdomain(this->triangulation->n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = this->triangulation->locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");


  // Create additional post-processor that derives information from the solution
  VorticityPostprocessor<dim> vorticity;
  data_out.add_data_vector(solution, vorticity);

  QCriterionPostprocessor<dim> qcriterion;
  data_out.add_data_vector(solution, qcriterion);

  SRFPostprocessor<dim> srf(nsparam.velocitySource.omega_x,
                            nsparam.velocitySource.omega_y,
                            nsparam.velocitySource.omega_z);

  if (nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::srf)
    data_out.add_data_vector(solution, srf);

  output_field_hook(data_out);

  // Build the patches and write the output

  data_out.build_patches(mapping,
                         subdivision,
                         DataOut<dim>::curved_inner_cells);

  write_vtu_and_pvd<dim>(this->pvdhandler,
                         data_out,
                         folder,
                         solution_name,
                         time,
                         iter,
                         group_files,
                         this->mpi_communicator);

  if (simulation_control->get_output_boundaries())
    {
      DataOutFaces<dim>          data_out_faces;
      BoundaryPostprocessor<dim> boundary_id;
      data_out_faces.attach_dof_handler(this->dof_handler);
      data_out_faces.add_data_vector(solution, boundary_id);
      data_out_faces.build_patches();

      write_boundaries_vtu<dim>(
        data_out_faces, folder, time, iter, this->mpi_communicator);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_forces()
{
  TimerOutput::Scope t(this->computing_timer, "output_forces");
  for (unsigned int boundary_id = 0;
       boundary_id < nsparam.boundary_conditions.size;
       ++boundary_id)
    {
      std::string filename = nsparam.forces_parameters.force_output_name + "." +
                             Utilities::int_to_string(boundary_id, 2) + ".dat";
      std::ofstream output(filename.c_str());

      forces_tables[boundary_id].write_text(output);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_torques()
{
  TimerOutput::Scope t(this->computing_timer, "output_torques");
  for (unsigned int boundary_id = 0;
       boundary_id < nsparam.boundary_conditions.size;
       ++boundary_id)
    {
      std::string filename = nsparam.forces_parameters.torque_output_name +
                             "." + Utilities::int_to_string(boundary_id, 2) +
                             ".dat";
      std::ofstream output(filename.c_str());

      this->torques_tables[boundary_id].write_text(output);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "write_checkpoint");
  std::string        prefix = this->nsparam.restart_parameters.filename;
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      simulation_control->save(prefix);
      this->pvdhandler.save(prefix);

      if (nsparam.flow_control.enable_flow_control)
        this->flow_control.save(prefix);
    }

  std::vector<const VectorType *> sol_set_transfer;
  auto &present_solution = this->get_present_solution();
  sol_set_transfer.push_back(&present_solution);
  sol_set_transfer.push_back(&this->solution_m1);
  sol_set_transfer.push_back(&this->solution_m2);
  sol_set_transfer.push_back(&this->solution_m3);

  if (nsparam.post_processing.calculate_average_velocities)
    {
      std::vector<const VectorType *> av_set_transfer =
        this->average_velocities.save(prefix);

      // Insert average velocities vectors into the set transfer vector
      sol_set_transfer.insert(sol_set_transfer.end(),
                              av_set_transfer.begin(),
                              av_set_transfer.end());
    }

  parallel::distributed::SolutionTransfer<dim, VectorType> system_trans_vectors(
    this->dof_handler);
  system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim> *>(
        this->triangulation.get()))
    {
      std::string triangulationName = prefix + ".triangulation";
      tria->save(prefix + ".triangulation");
    }
}

// Pre-compile the 2D and 3D version with the types that can occur
template class NavierStokesBase<2, TrilinosWrappers::MPI::Vector, IndexSet>;
template class NavierStokesBase<3, TrilinosWrappers::MPI::Vector, IndexSet>;
template class NavierStokesBase<2,
                                TrilinosWrappers::MPI::BlockVector,
                                std::vector<IndexSet>>;
template class NavierStokesBase<3,
                                TrilinosWrappers::MPI::BlockVector,
                                std::vector<IndexSet>>;
