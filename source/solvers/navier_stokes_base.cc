// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/grids.h>
#include <core/lethe_grid_tools.h>
#include <core/mesh_controller.h>
#include <core/mortar_coupling_manager.h>
#include <core/solutions_output.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/flow_control.h>
#include <solvers/navier_stokes_base.h>
#include <solvers/postprocessing_cfd.h>
#include <solvers/postprocessing_velocities.h>
#include <solvers/postprocessors.h>
#include <solvers/postprocessors_smoothing.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_bubbles.h>
#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/data_out_resample.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <sys/stat.h>

/*
 * Constructor for the Navier-Stokes base class
 */
template <int dim, typename VectorType, typename DofsType>
NavierStokesBase<dim, VectorType, DofsType>::NavierStokesBase(
  SimulationParameters<dim> &p_nsparam)
  : PhysicsSolver<VectorType>(
      p_nsparam.non_linear_solver.at(PhysicsID::fluid_dynamics))
  , mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , dof_handler()
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , simulation_parameters(p_nsparam)
  , flow_control(simulation_parameters.flow_control)
  , velocity_fem_degree(p_nsparam.fem_parameters.velocity_order)
  , pressure_fem_degree(p_nsparam.fem_parameters.pressure_order)
  , number_quadrature_points(p_nsparam.fem_parameters.velocity_order + 1)
  , mesh_controller(p_nsparam.mesh_adaptation.maximum_number_elements)
{
  if (simulation_parameters.mesh.simplex)
    {
      // for simplex meshes
      const FE_SimplexP<dim> velocity_fe(
        p_nsparam.fem_parameters.velocity_order);
      const FE_SimplexP<dim> pressure_fe(
        p_nsparam.fem_parameters.pressure_order);
      fe = std::make_shared<FESystem<dim>>(velocity_fe, dim, pressure_fe, 1);

      AssertThrow(
        !p_nsparam.fem_parameters.enable_bubble_function_velocity ||
          !p_nsparam.fem_parameters.enable_bubble_function_pressure,
        ExcMessage(
          "Bubble enrichment functions are not compatible with simplex meshes."));

      mapping = std::make_shared<MappingFE<dim>>(velocity_fe);
      cell_quadrature =
        std::make_shared<QGaussSimplex<dim>>(number_quadrature_points);
      face_quadrature =
        std::make_shared<QGaussSimplex<dim - 1>>(number_quadrature_points);
      triangulation =
        std::make_shared<parallel::fullydistributed::Triangulation<dim>>(
          this->mpi_communicator);
      dof_handler.clear();
      dof_handler.reinit(*this->triangulation);
    }
  else
    {
      // Usual case, for quad/hex meshes
      std::shared_ptr<FiniteElement<dim>> fe_u;
      std::shared_ptr<FiniteElement<dim>> fe_p;

      if (p_nsparam.fem_parameters.enable_bubble_function_velocity)
        fe_u = std::make_shared<FE_Q_Bubbles<dim>>(
          p_nsparam.fem_parameters.velocity_order);
      else
        fe_u =
          std::make_shared<FE_Q<dim>>(p_nsparam.fem_parameters.velocity_order);

      if (p_nsparam.fem_parameters.enable_bubble_function_pressure)
        fe_p = std::make_shared<FE_Q_Bubbles<dim>>(
          p_nsparam.fem_parameters.pressure_order);
      else
        fe_p =
          std::make_shared<FE_Q<dim>>(p_nsparam.fem_parameters.pressure_order);

      fe              = std::make_shared<FESystem<dim>>(*fe_u, dim, *fe_p, 1);
      mapping         = std::make_shared<MappingQ<dim>>(velocity_fem_degree);
      cell_quadrature = std::make_shared<QGauss<dim>>(number_quadrature_points);
      face_quadrature =
        std::make_shared<QGauss<dim - 1>>(number_quadrature_points);
      triangulation =
        std::make_shared<parallel::distributed::Triangulation<dim>>(
          this->mpi_communicator,
          Triangulation<dim>::limit_level_difference_at_vertices,
          (p_nsparam.linear_solver.at(fluid_dynamics).preconditioner ==
           Parameters::LinearSolver::PreconditionerType::lsmg) ?
            parallel::distributed::Triangulation<
              dim>::construct_multigrid_hierarchy :
            parallel::distributed::Triangulation<dim>::default_setting);
      dof_handler.clear();
      dof_handler.reinit(*this->triangulation);
    }

  this->pcout.set_condition(
    Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0);

  // Check if the output directory exists
  std::string output_dir_name =
    simulation_parameters.simulation_control.output_folder;
  struct stat buffer;

  // If output directory does not exist, create it
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      if (stat(output_dir_name.c_str(), &buffer) != 0)
        {
          create_output_folder(output_dir_name);
        }
    }

  if (simulation_parameters.simulation_control.method ==
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    simulation_control = std::make_shared<SimulationControlSteady>(
      simulation_parameters.simulation_control);
  else if (simulation_parameters.simulation_control.method ==
           Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    {
      simulation_control = std::make_shared<SimulationControlAdjointSteady>(
        simulation_parameters.simulation_control);
    }
  else
    {
      simulation_control = std::make_shared<SimulationControlTransient>(
        simulation_parameters.simulation_control);
    }

  // Provide the simulation control object to the physical property manager
  simulation_parameters.physical_properties_manager.provide_simulation_control(
    simulation_control);

  multiphysics = std::make_shared<MultiphysicsInterface<dim>>(
    simulation_parameters, triangulation, simulation_control, this->pcout);

  // Pre-allocate memory for the previous solutions using the information
  // of the BDF schemes
  previous_solutions.resize(maximum_number_of_previous_solutions());

  // Change the behavior of the timer for situations when you don't want
  // outputs
  if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
    this->computing_timer.disable_output();

  // Get the exact solution from the parser
  exact_solution = &simulation_parameters.analytical_solution->uvwp;

  // If there is a forcing function, get it from the parser
  forcing_function = simulation_parameters.source_term.navier_stokes_source;

  if (this->simulation_parameters.post_processing
        .calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    average_velocities =
      std::make_shared<AverageVelocities<dim, VectorType, DofsType>>(
        dof_handler);

  this->pcout << "Running on "
              << Utilities::MPI::n_mpi_processes(this->mpi_communicator)
              << " MPI rank(s)..." << std::endl;

  this->pcout << std::setprecision(
    simulation_parameters.simulation_control.log_precision);
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::dynamic_flow_control()
{
  if (simulation_parameters.flow_control.enable_flow_control &&
      simulation_parameters.simulation_control.method !=
        Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      // Calculate the average velocity
      double average_velocity = calculate_average_velocity(
        this->dof_handler,
        this->present_solution,
        simulation_parameters.flow_control.boundary_flow_id,
        *this->face_quadrature,
        *this->get_mapping());

      this->flow_control.calculate_beta(average_velocity,
                                        simulation_control->get_time_step(),
                                        simulation_control->get_step_number());

      // Showing results
      if (simulation_parameters.flow_control.verbosity ==
            Parameters::Verbosity::verbose &&
          simulation_control->get_step_number() > 0 &&
          this->this_mpi_process == 0)
        {
          announce_string(this->pcout, "Flow control summary");
          this->pcout << "Space-average velocity: " << average_velocity
                      << std::endl;
          this->pcout
            << "Beta force:  "
            << flow_control
                 .get_beta()[simulation_parameters.flow_control.flow_direction]
            << std::endl;
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
bool
NavierStokesBase<dim, VectorType, DofsType>::check_existance_of_bc(
  BoundaryConditions::BoundaryType bc)
{
  bool bc_exist = false;
  // Loop over the boundary to check if they need assembler on their face.
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == bc)
        bc_exist = true;
    }
  return bc_exist;
}


// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_forces(
  const VectorType &evaluation_point)
{
  TimerOutput::Scope t(this->computing_timer, "Calculate forces");

  this->forces_on_boundaries =
    calculate_forces(this->dof_handler,
                     evaluation_point,
                     simulation_parameters.physical_properties_manager,
                     simulation_parameters.boundary_conditions,
                     *this->face_quadrature,
                     *this->get_mapping());

  if (simulation_parameters.forces_parameters.verbosity ==
        Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      std::cout << std::endl;
      std::string independent_column_names = "Boundary ID";

      std::vector<std::string> dependent_column_names;
      dependent_column_names.emplace_back("f_x");
      dependent_column_names.emplace_back("f_y");
      if (dim == 3)
        dependent_column_names.emplace_back("f_z");
      dependent_column_names.emplace_back("f_xv");
      dependent_column_names.emplace_back("f_yv");
      if (dim == 3)
        dependent_column_names.emplace_back("f_zv");
      dependent_column_names.emplace_back("f_xp");
      dependent_column_names.emplace_back("f_yp");
      if (dim == 3)
        dependent_column_names.emplace_back("f_zp");

      std::vector<unsigned int> boundary_ids =
        extract_keys_from_map(simulation_parameters.boundary_conditions.type);

      std::vector<Tensor<1, dim>> total_forces =
        extract_values_from_map(this->forces_on_boundaries[0]);

      std::vector<Tensor<1, dim>> viscous_forces =
        extract_values_from_map(this->forces_on_boundaries[1]);

      std::vector<Tensor<1, dim>> pressure_forces =
        extract_values_from_map(this->forces_on_boundaries[2]);

      std::vector<std::vector<Tensor<1, dim>>> all_forces{total_forces,
                                                          viscous_forces,
                                                          pressure_forces};

      TableHandler table = make_table_scalars_tensors(
        boundary_ids,
        independent_column_names,
        all_forces,
        dependent_column_names,
        this->simulation_parameters.simulation_control.log_precision,
        true);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Force summary                           |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (auto const &[id, type] : simulation_parameters.boundary_conditions.type)
    {
      if (simulation_control->is_steady())
        {
          this->forces_tables[id].add_value(
            "cells", this->triangulation->n_global_active_cells());
        }
      else
        {
          this->forces_tables[id].add_value(
            "time", simulation_control->get_current_time());
          this->forces_tables[id].set_precision(
            "time", simulation_parameters.forces_parameters.output_precision);
        }

      this->forces_tables[id].add_value("f_x",
                                        this->forces_on_boundaries[0][id][0]);
      this->forces_tables[id].add_value("f_y",
                                        this->forces_on_boundaries[0][id][1]);
      if (dim == 3)
        this->forces_tables[id].add_value("f_z",
                                          this->forces_on_boundaries[0][id][2]);
      else
        this->forces_tables[id].add_value("f_z", 0.);

      this->forces_tables[id].add_value("f_xv",
                                        this->forces_on_boundaries[1][id][0]);
      this->forces_tables[id].add_value("f_yv",
                                        this->forces_on_boundaries[1][id][1]);
      if (dim == 3)
        this->forces_tables[id].add_value("f_zv",
                                          this->forces_on_boundaries[1][id][2]);
      else
        this->forces_tables[id].add_value("f_zv", 0.);

      this->forces_tables[id].add_value("f_xp",
                                        this->forces_on_boundaries[2][id][0]);
      this->forces_tables[id].add_value("f_yp",
                                        this->forces_on_boundaries[2][id][1]);
      if (dim == 3)
        this->forces_tables[id].add_value("f_zp",
                                          this->forces_on_boundaries[2][id][2]);
      else
        this->forces_tables[id].add_value("f_zp", 0.);

      // Precision
      this->forces_tables[id].set_precision(
        "f_x", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[id].set_precision(
        "f_y", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[id].set_precision(
        "f_z", simulation_parameters.forces_parameters.output_precision);

      this->forces_tables[id].set_precision(
        "f_xv", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[id].set_precision(
        "f_yv", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[id].set_precision(
        "f_zv", simulation_parameters.forces_parameters.output_precision);

      this->forces_tables[id].set_precision(
        "f_xp", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[id].set_precision(
        "f_yp", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[id].set_precision(
        "f_zp", simulation_parameters.forces_parameters.output_precision);
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_torques(
  const VectorType &evaluation_point)
{
  TimerOutput::Scope t(this->computing_timer, "Calculate torques");

  std::map<types::boundary_id, Tensor<1, 3>> torques_on_boundaries =
    calculate_torques(this->dof_handler,
                      evaluation_point,
                      simulation_parameters.physical_properties_manager,
                      simulation_parameters.boundary_conditions,
                      *this->face_quadrature,
                      *this->get_mapping());

  if (simulation_parameters.forces_parameters.verbosity ==
        Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      this->pcout << std::endl;
      std::string independent_column_names = "Boundary ID";

      std::vector<std::string> dependent_column_names;
      dependent_column_names.emplace_back("T_x");
      dependent_column_names.emplace_back("T_y");
      dependent_column_names.emplace_back("T_z");

      std::vector<unsigned int> boundary_ids =
        extract_keys_from_map(simulation_parameters.boundary_conditions.type);
      std::vector<Tensor<1, 3>> torques =
        extract_values_from_map(torques_on_boundaries);

      TableHandler table = make_table_scalars_tensors(
        boundary_ids,
        independent_column_names,
        torques,
        dependent_column_names,
        this->simulation_parameters.simulation_control.log_precision,
        true);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Torque summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (auto const &[boundary_id, type] :
       simulation_parameters.boundary_conditions.type)
    {
      if (simulation_control->is_steady())
        {
          this->torques_tables[boundary_id].add_value(
            "cells", this->triangulation->n_global_active_cells());
        }
      else
        {
          this->torques_tables[boundary_id].add_value(
            "time", simulation_control->get_current_time());
          this->torques_tables[boundary_id].set_precision(
            "time", simulation_parameters.forces_parameters.output_precision);
        }
      this->torques_tables[boundary_id].add_value(
        "T_x", torques_on_boundaries[boundary_id][0]);
      this->torques_tables[boundary_id].add_value(
        "T_y", torques_on_boundaries[boundary_id][1]);
      this->torques_tables[boundary_id].add_value(
        "T_z", torques_on_boundaries[boundary_id][2]);

      // Precision
      this->torques_tables[boundary_id].set_precision(
        "T_x", simulation_parameters.forces_parameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "T_y", simulation_parameters.forces_parameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "T_z", simulation_parameters.forces_parameters.output_precision);
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::finish_simulation_fd()
{
  if (simulation_parameters.forces_parameters.calculate_force)
    this->write_output_forces();

  if (simulation_parameters.forces_parameters.calculate_torque)
    this->write_output_torques();

  if (simulation_parameters.analytical_solution->calculate_error())
    {
      if (simulation_parameters.simulation_control.method ==
          Parameters::SimulationControl::TimeSteppingMethod::steady)
        {
          error_table.set_scientific("error_pressure", true);
          error_table.omit_column_from_convergence_rate_evaluation("cells");
          error_table.evaluate_all_convergence_rates(
            ConvergenceTable::reduction_rate_log2);
        }
      error_table.set_scientific("error_velocity", true);

      if (this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.analytical_solution->get_filename() + ".dat";
          std::ofstream output(filename.c_str());
          error_table.write_text(output);
          std::vector<std::string> sub_columns;
          if (simulation_parameters.simulation_control.method ==
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              sub_columns.emplace_back("cells");
              sub_columns.emplace_back("error_velocity");
              sub_columns.emplace_back("error_pressure");
              error_table.set_precision(
                "error_pressure", simulation_control->get_log_precision());
              error_table.set_column_order(sub_columns);
            }
          error_table.set_precision("error_velocity",
                                    simulation_control->get_log_precision());

          error_table.write_text(std::cout);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::percolate_time_vectors_fd()
{
  for (unsigned int i = previous_solutions.size() - 1; i > 0; --i)
    {
      previous_solutions[i] = previous_solutions[i - 1];
    }
  previous_solutions[0] = this->present_solution;
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::finish_time_step()
{
  if (simulation_parameters.simulation_control.method !=
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      TimerOutput::Scope t(this->computing_timer,
                           "Calculate CFL and percolate time vectors");

      percolate_time_vectors_fd();
      const double CFL = calculate_CFL(this->dof_handler,
                                       this->present_solution,
                                       simulation_control->get_time_step(),
                                       *this->cell_quadrature,
                                       *this->get_mapping());
      this->simulation_control->set_CFL(CFL);
    }
  if (this->simulation_parameters.restart_parameters.checkpoint &&
      simulation_control->get_step_number() != 0 &&
      (simulation_control->get_step_number() %
           this->simulation_parameters.restart_parameters.frequency ==
         0 ||
       simulation_control->is_at_end()))
    {
      this->write_checkpoint();
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
    {
      announce_string(this->pcout, "Fluid Dynamics");
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
  auto &present_solution = this->present_solution;

  if (simulation_parameters.multiphysics.fluid_dynamics)
    {
      // Solve and percolate the auxiliary physics that should be treated BEFORE
      // the fluid dynamics
      multiphysics->solve(false,
                          simulation_parameters.simulation_control.method);
      multiphysics->percolate_time_vectors(false);

      if (simulation_parameters.non_linear_solver.at(PhysicsID::fluid_dynamics)
              .verbosity != Parameters::Verbosity::quiet ||
          simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
              .verbosity != Parameters::Verbosity::quiet)
        announce_string(this->pcout, "Fluid Dynamics");
      PhysicsSolver<VectorType>::solve_non_linear_system(false);

      // If the physics need to be solved after the physics, the matrix free
      // solver requires to update the value here. This is due to the different
      // type of vectors.
      if (this->multiphysics->get_active_physics().size() > 1)
        this->update_solutions_for_multiphysics();

      // Solve and percolate the auxiliary physics that should be treated AFTER
      // the fluid dynamics
      multiphysics->solve(true,
                          simulation_parameters.simulation_control.method);
      // Dear future Bruno, percolating auxiliary physics before fluid dynamics
      // is necessary because of the checkpointing mechanism. You spent an
      // evening debugging this, trust me.
      multiphysics->percolate_time_vectors(true);
    }
  else
    {
      // Fluid dynamics is not to be solved, but rather specified. Update
      // condition and move on.

      // Solve and percolate the auxiliary physics that should be treated
      // BEFORE the fluid dynamics
      multiphysics->solve(false,
                          simulation_parameters.simulation_control.method);
      multiphysics->percolate_time_vectors(false);

      if (this->simulation_parameters.initial_condition->type ==
          Parameters::InitialConditionType::average_velocity_profile)
        {
          // We get the solution via the average solution
          this->local_evaluation_point =
            this->average_velocities->get_average_velocities();
          present_solution = this->local_evaluation_point;
        }
      else
        {
          // We get the solution via an initial condition
          this->simulation_parameters.initial_condition->uvwp.set_time(
            this->simulation_control->get_current_time());
          set_initial_condition_fd(
            this->simulation_parameters.initial_condition->type);
        }

      // Solve and percolate the auxiliary physics that should be treated
      // AFTER the fluid dynamics
      multiphysics->solve(true,
                          simulation_parameters.simulation_control.method);
      multiphysics->percolate_time_vectors(true);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::
  enable_dynamic_zero_constraints_fd()
{
  if (!this->simulation_parameters.constrain_solid_domain.enable)
    return;

  // Initialize StasisConstraintWithTemperature structs for each constraint
  for (unsigned int c_id = 0;
       c_id <
       this->simulation_parameters.constrain_solid_domain.number_of_constraints;
       c_id++)
    {
      StasisConstraintWithTemperature stasis_constraint_struct(
        this->simulation_parameters.constrain_solid_domain.fluid_ids[c_id],
        this->simulation_parameters.constrain_solid_domain
          .temperature_min_values[c_id],
        this->simulation_parameters.constrain_solid_domain
          .temperature_max_values[c_id],
        this->simulation_parameters.constrain_solid_domain
          .filtered_phase_fraction_tolerance[c_id]);
      this->stasis_constraint_structs.emplace_back(stasis_constraint_struct);
    }

  // For temperature-dependent constraints
  const DoFHandler<dim> *dof_handler_ht =
    this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

  this->fe_values_temperature =
    std::make_shared<FEValues<dim>>(*this->get_mapping(),
                                    dof_handler_ht->get_fe(),
                                    *this->cell_quadrature,
                                    update_values);

  // For VOF simulations
  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      this->fe_values_vof =
        std::make_shared<FEValues<dim>>(*this->get_mapping(),
                                        dof_handler_vof->get_fe(),
                                        *this->cell_quadrature,
                                        update_values);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh()
{
  bool refinement_step;
  if (this->simulation_parameters.mesh_adaptation.refinement_at_frequency)
    refinement_step = this->simulation_control->get_step_number() %
                        this->simulation_parameters.mesh_adaptation.frequency ==
                      0;
  else
    refinement_step = this->simulation_control->get_step_number() == 0;
  if (refinement_step)
    {
      if (this->simulation_parameters.mesh_adaptation.type ==
          Parameters::MeshAdaptation::Type::kelly)
        refine_mesh_kelly();

      else if (this->simulation_parameters.mesh_adaptation.type ==
               Parameters::MeshAdaptation::Type::uniform)
        refine_mesh_uniform();
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::box_refine_mesh(const bool restart)
{
  if (restart)
    {
      return;
    }
  // Read the mesh that define the box use in this function
  Triangulation<dim> box_to_refine;
  if (this->simulation_parameters.mesh_box_refinement->box_mesh->type ==
      Parameters::Mesh::Type::gmsh)
    {
      if (this->simulation_parameters.mesh_box_refinement->box_mesh->simplex)
        {
          Triangulation<dim> basetria(
            Triangulation<dim>::limit_level_difference_at_vertices);

          GridIn<dim> grid_in;
          grid_in.attach_triangulation(basetria);
          std::ifstream input_file(this->simulation_parameters
                                     .mesh_box_refinement->box_mesh->file_name);

          grid_in.read_msh(input_file);

          // By default uses the METIS partitioner.
          // A user parameter option could be made to choose a partitionner.
          GridTools::partition_triangulation(0, basetria);


          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(basetria, nullptr);

          triangulation->create_triangulation(construction_data);
        }
      else
        {
          GridIn<dim> grid_in;
          grid_in.attach_triangulation(box_to_refine);
          std::ifstream input_file(this->simulation_parameters
                                     .mesh_box_refinement->box_mesh->file_name);
          grid_in.read_msh(input_file);
        }
    }
  // Dealii grids
  else if (this->simulation_parameters.mesh_box_refinement->box_mesh->type ==
           Parameters::Mesh::Type::dealii)
    {
      if (this->simulation_parameters.mesh_box_refinement->box_mesh->simplex)
        {
          Triangulation<dim> temporary_quad_triangulation;
          GridGenerator::generate_from_name_and_arguments(
            temporary_quad_triangulation,
            this->simulation_parameters.mesh_box_refinement->box_mesh
              ->grid_type,
            this->simulation_parameters.mesh_box_refinement->box_mesh
              ->grid_arguments);

          // initial refinement
          const int initial_refinement =
            this->simulation_parameters.mesh_box_refinement->box_mesh
              ->initial_refinement;
          temporary_quad_triangulation.refine_global(initial_refinement);
          // flatten the triangulation
          Triangulation<dim> flat_temp_quad_triangulation;
          GridGenerator::flatten_triangulation(temporary_quad_triangulation,
                                               flat_temp_quad_triangulation);

          Triangulation<dim> temporary_tri_triangulation(
            Triangulation<dim>::limit_level_difference_at_vertices);
          GridGenerator::convert_hypercube_to_simplex_mesh(
            flat_temp_quad_triangulation, temporary_tri_triangulation);

          GridTools::partition_triangulation_zorder(
            0, temporary_tri_triangulation);
          GridTools::partition_multigrid_levels(temporary_tri_triangulation);

          // extract relevant information from distributed triangulation
          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(
              temporary_tri_triangulation,
              nullptr,
              TriangulationDescription::Settings::
                construct_multigrid_hierarchy);
          box_to_refine.create_triangulation(construction_data);
        }
      else
        {
          GridGenerator::generate_from_name_and_arguments(
            box_to_refine,
            this->simulation_parameters.mesh_box_refinement->box_mesh
              ->grid_type,
            this->simulation_parameters.mesh_box_refinement->box_mesh
              ->grid_arguments);
        }
    }

  // Define a local dofhandler of this mesh. This won't be needed in later
  // version of LetheGridTools

  box_to_refine.refine_global(this->simulation_parameters.mesh_box_refinement
                                ->box_mesh->initial_refinement);
  DoFHandler<dim> box_to_refine_dof_handler(box_to_refine);
  // Refine the number of time needed
  for (unsigned int i = 0;
       i < this->simulation_parameters.mesh_box_refinement->initial_refinement;
       ++i)
    {
      if (dynamic_cast<parallel::distributed::Triangulation<dim> *>(
            this->triangulation.get()) == nullptr)
        return;

      auto &tria = *dynamic_cast<parallel::distributed::Triangulation<dim> *>(
        this->triangulation.get());

      // Time monitoring
      TimerOutput::Scope t(this->computing_timer, "Box refine");
      this->pcout
        << "Initial refinement in box - Step  " << i + 1 << " of "
        << this->simulation_parameters.mesh_box_refinement->initial_refinement
        << std::endl;


      Vector<float> estimated_error_per_cell(tria.n_active_cells());
      const FEValuesExtractors::Vector velocity(0);
      const FEValuesExtractors::Scalar pressure(dim);
      auto &present_solution = this->present_solution;

      const auto &cell_iterator =
        box_to_refine_dof_handler.active_cell_iterators();

      // Find all the cells of the principal mesh that are partially contained
      // inside the box_mesh and set them up for refinement.
      for (const auto &cell : cell_iterator)
        {
          std::vector<typename DoFHandler<dim>::active_cell_iterator>
            cell_to_refine;
          cell_to_refine =
            (LetheGridTools::find_cells_in_cells(this->dof_handler, cell));
          for (unsigned int j = 0; j < cell_to_refine.size(); ++j)
            {
              cell_to_refine[j]->set_refine_flag();
            }
        }

      tria.prepare_coarsening_and_refinement();

      // Solution transfer objects for all the solutions
      SolutionTransfer<dim, VectorType> solution_transfer(this->dof_handler,
                                                          true);
      std::vector<SolutionTransfer<dim, VectorType>>
        previous_solutions_transfer;
      // Important to reserve to prevent pointer dangling
      previous_solutions_transfer.reserve(previous_solutions.size());
      for (unsigned int i = 0; i < previous_solutions.size(); ++i)
        {
          previous_solutions_transfer.emplace_back(
            SolutionTransfer<dim, VectorType>(this->dof_handler, true));
          if constexpr (std::is_same_v<
                          VectorType,
                          LinearAlgebra::distributed::Vector<double>>)
            previous_solutions[i].update_ghost_values();
          previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
            previous_solutions[i]);
        }

      SolutionTransfer<dim, VectorType> solution_transfer_m1(this->dof_handler,
                                                             true);
      SolutionTransfer<dim, VectorType> solution_transfer_m2(this->dof_handler,
                                                             true);
      SolutionTransfer<dim, VectorType> solution_transfer_m3(this->dof_handler,
                                                             true);

      if constexpr (std::is_same_v<VectorType,
                                   LinearAlgebra::distributed::Vector<double>>)
        present_solution.update_ghost_values();
      solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

      multiphysics->prepare_for_mesh_adaptation();
      if (this->simulation_parameters.post_processing
            .calculate_average_velocities)
        average_velocities->prepare_for_mesh_adaptation();

      tria.execute_coarsening_and_refinement();
      this->setup_dofs();

      // Set up the vectors for the transfer
      VectorType tmp = init_temporary_vector();

      // Interpolate the solution at time and previous time
      solution_transfer.interpolate(tmp);

      // Distribute constraints
      auto &nonzero_constraints = this->nonzero_constraints;
      nonzero_constraints.distribute(tmp);

      // Fix on the new mesh
      present_solution = tmp;

      for (unsigned int i = 0; i < previous_solutions.size(); ++i)
        {
          VectorType tmp_previous_solution = init_temporary_vector();
          previous_solutions_transfer[i].interpolate(tmp_previous_solution);
          nonzero_constraints.distribute(tmp_previous_solution);
          previous_solutions[i] = tmp_previous_solution;
        }

      multiphysics->post_mesh_adaptation();
      if (this->simulation_parameters.post_processing
            .calculate_average_velocities)
        average_velocities->post_mesh_adaptation();
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
  TimerOutput::Scope t(this->computing_timer, "Refine");

  Vector<float> estimated_error_per_cell(tria.n_active_cells());
  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);
  auto                            &present_solution = this->present_solution;
  VectorType                       locally_relevant_solution;
  locally_relevant_solution.reinit(this->locally_owned_dofs,
                                   this->locally_relevant_dofs,
                                   this->mpi_communicator);
  locally_relevant_solution = this->present_solution;
  locally_relevant_solution.update_ghost_values();

  // Global flags
  // Their dimension is consistent with the dimension returned by
  // save_refine_flags(), in order to be able to use load_refine_flags()
  std::vector<bool> global_refine_flags(dim * tria.n_active_cells(), false);
  std::vector<bool> global_coarsen_flags(dim * tria.n_active_cells(), false);

  bool         first_variable(true);
  const double coarsening_factor = mesh_controller.calculate_coarsening_factor(
    this->triangulation->n_global_active_cells());

  unsigned int maximal_number_of_elements =
    this->simulation_parameters.mesh_adaptation.maximum_number_elements;
  // Override the maximal number of elements if the controller is used. The
  // controller will find a coarsening_factor that respects the user-defined
  // maximal_number_of_elements.
  if (this->simulation_parameters.mesh_adaptation.mesh_controller_is_enabled)
    maximal_number_of_elements = INT_MAX;

  for (const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
         &ivar : this->simulation_parameters.mesh_adaptation.variables)
    {
      double ivar_coarsening_factor = ivar.second.coarsening_fraction;
      if (this->simulation_parameters.mesh_adaptation
            .mesh_controller_is_enabled)
        ivar_coarsening_factor = coarsening_factor;

      if (ivar.first == Variable::pressure)
        {
          KellyErrorEstimator<dim>::estimate(
            *this->get_mapping(),
            this->dof_handler,
            *this->face_quadrature,
            typename std::map<types::boundary_id,
                              const Function<dim, double> *>(),
            locally_relevant_solution,
            estimated_error_per_cell,
            this->fe->component_mask(pressure));
        }
      else if (ivar.first == Variable::velocity)
        {
          KellyErrorEstimator<dim>::estimate(
            *this->get_mapping(),
            this->dof_handler,
            *this->face_quadrature,
            typename std::map<types::boundary_id,
                              const Function<dim, double> *>(),
            locally_relevant_solution,
            estimated_error_per_cell,
            this->fe->component_mask(velocity));
        }
      else
        {
          // refine_mesh on an auxiliary physic parameter
          multiphysics->compute_kelly(ivar, estimated_error_per_cell);
        }

      if (this->simulation_parameters.mesh_adaptation.fractionType ==
          Parameters::MeshAdaptation::FractionType::number)
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
          tria,
          estimated_error_per_cell,
          ivar.second.refinement_fraction,
          ivar_coarsening_factor,
          maximal_number_of_elements);

      else if (this->simulation_parameters.mesh_adaptation.fractionType ==
               Parameters::MeshAdaptation::FractionType::fraction)
        parallel::distributed::GridRefinement::
          refine_and_coarsen_fixed_fraction(tria,
                                            estimated_error_per_cell,
                                            ivar.second.refinement_fraction,
                                            ivar_coarsening_factor);

      // Remove the flags if the cell is at the boundary and is set as do not
      // touch in the parameter file
      if (this->simulation_parameters.mesh_adaptation
            .is_boundary_refinement_fixed)
        for (const auto &cell : tria.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              if (cell->at_boundary())
                for (const auto &face : cell->face_iterators())
                  {
                    if (std::find(this->simulation_parameters.mesh_adaptation
                                    .boundaries_to_fix.begin(),
                                  this->simulation_parameters.mesh_adaptation
                                    .boundaries_to_fix.end(),
                                  face->boundary_id()) !=
                        this->simulation_parameters.mesh_adaptation
                          .boundaries_to_fix.end())
                      {
                        cell->clear_refine_flag();
                        cell->clear_coarsen_flag();
                      }
                  }
            }

      std::vector<bool> current_refine_flags;
      std::vector<bool> current_coarsen_flags;

      tria.save_refine_flags(current_refine_flags);
      tria.save_coarsen_flags(current_coarsen_flags);

      // Fill global flags
      if (first_variable)
        {
          // special case of the first refinement variable
          global_refine_flags  = current_refine_flags;
          global_coarsen_flags = current_coarsen_flags;

          first_variable = false;
        }
      else
        {
          // for subsequent refinement variables
          for (std::vector<bool>::size_type i = 0;
               i != global_refine_flags.size();
               ++i)
            {
              // refine if at least refinement flag on one variable
              global_refine_flags[i] =
                global_refine_flags[i] || current_refine_flags[i];
            }

          // for subsequent coarsen variables
          for (std::vector<bool>::size_type i = 0;
               i != global_coarsen_flags.size();
               ++i)
            {
              // coarsen if refinement flag on all variables
              global_coarsen_flags[i] =
                global_coarsen_flags[i] && current_coarsen_flags[i];
            }
        }

      // Clear flags in the triangulation
      for (const auto &cell : tria.active_cell_iterators())
        {
          cell->clear_coarsen_flag();
          cell->clear_refine_flag();
        }
    }

  // Load global flags
  tria.load_refine_flags(global_refine_flags);
  tria.load_coarsen_flags(global_coarsen_flags);

  if (tria.n_levels() >
      this->simulation_parameters.mesh_adaptation.maximum_refinement_level)
    for (typename Triangulation<dim>::active_cell_iterator cell =
           tria.begin_active(this->simulation_parameters.mesh_adaptation
                               .maximum_refinement_level);
         cell != tria.end();
         ++cell)
      cell->clear_refine_flag();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active(this->simulation_parameters.mesh_adaptation
                             .minimum_refinement_level);
       cell !=
       tria.end_active(
         this->simulation_parameters.mesh_adaptation.minimum_refinement_level);
       ++cell)
    cell->clear_coarsen_flag();

  tria.prepare_coarsening_and_refinement();

  // Solution transfer objects for all the solutions
  SolutionTransfer<dim, VectorType> solution_transfer(this->dof_handler, true);
  std::vector<SolutionTransfer<dim, VectorType>> previous_solutions_transfer;
  // Important to reserve to prevent pointer dangling
  previous_solutions_transfer.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer.emplace_back(
        SolutionTransfer<dim, VectorType>(this->dof_handler, true));
      if constexpr (std::is_same_v<VectorType,
                                   LinearAlgebra::distributed::Vector<double>>)
        previous_solutions[i].update_ghost_values();
      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }

  if constexpr (std::is_same_v<VectorType,
                               LinearAlgebra::distributed::Vector<double>>)
    present_solution.update_ghost_values();
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

  multiphysics->prepare_for_mesh_adaptation();
  if (this->simulation_parameters.post_processing
        .calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    average_velocities->prepare_for_mesh_adaptation();

  tria.execute_coarsening_and_refinement();
  this->setup_dofs();

  // Set up the vectors for the transfer
  VectorType tmp = init_temporary_vector();

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);

  // Distribute constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      VectorType tmp_previous_solution = init_temporary_vector();
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }

  multiphysics->post_mesh_adaptation();
  if (this->simulation_parameters.post_processing
        .calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    average_velocities->post_mesh_adaptation();

  // Only needed if other physics apart from fluid dynamics are enabled.
  if (this->multiphysics->get_active_physics().size() > 1)
    this->update_multiphysics_time_average_solution();
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh_uniform()
{
  if (this->triangulation->n_global_levels() >
      this->simulation_parameters.mesh_adaptation.maximum_refinement_level)
    return;
  AssertThrow(this->triangulation->all_reference_cells_are_hyper_cube(),
              ExcMessage("Uniform refinement is not supported for "
                         "simplex meshes."));
  TimerOutput::Scope t(this->computing_timer, "Refine");

  // Solution transfer objects for all the solutions
  SolutionTransfer<dim, VectorType> solution_transfer(this->dof_handler, true);
  SolutionTransfer<dim, VectorType> solution_transfer_m2(this->dof_handler,
                                                         true);
  SolutionTransfer<dim, VectorType> solution_transfer_m3(this->dof_handler,
                                                         true);

  if constexpr (std::is_same_v<VectorType,
                               LinearAlgebra::distributed::Vector<double>>)
    present_solution.update_ghost_values();

  solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);

  std::vector<SolutionTransfer<dim, VectorType>> previous_solutions_transfer;
  // Important to reserve to prevent pointer dangling
  previous_solutions_transfer.reserve(previous_solutions.size());

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer.emplace_back(
        SolutionTransfer<dim, VectorType>(this->dof_handler, true));

      if constexpr (std::is_same_v<VectorType,
                                   LinearAlgebra::distributed::Vector<double>>)
        previous_solutions[i].update_ghost_values();

      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }

  multiphysics->prepare_for_mesh_adaptation();

  // Refine
  this->triangulation->refine_global(1);

  // If mortar is enabled, update mapping cache with refined triangulation
  if (this->simulation_parameters.mortar_parameters.enable)
    this->mapping_cache->initialize(*this->mapping, *this->triangulation);

  setup_dofs();

  // Set up the vectors for the transfer
  VectorType tmp = init_temporary_vector();

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);

  // Distribute constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;

  // Set up the vectors for the transfer
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      VectorType tmp_previous_solution = init_temporary_vector();
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }

  multiphysics->post_mesh_adaptation();
  if (this->simulation_parameters.post_processing
        .calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    average_velocities->post_mesh_adaptation();

  // Only needed if other physics apart from fluid dynamics are enabled.
  if (this->multiphysics->get_active_physics().size() > 1)
    this->update_multiphysics_time_average_solution();
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocess_fd(bool firstIter)
{
  auto &present_solution = this->present_solution;

  // Enstrophy
  if (this->simulation_parameters.post_processing.calculate_enstrophy)
    {
      TimerOutput::Scope t(this->computing_timer, "Calculate enstrophy");

      double enstrophy = calculate_enstrophy(this->dof_handler,
                                             present_solution,
                                             *this->cell_quadrature,
                                             *this->get_mapping());

      this->enstrophy_table.add_value("time",
                                      simulation_control->get_current_time());
      this->enstrophy_table.add_value("enstrophy", enstrophy);

      // Display Enstrophy to screen if verbosity is enabled
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Enstrophy: " << enstrophy << std::endl;
        }

      // Output Enstrophy to a text file from processor 0
      if (simulation_control->get_step_number() %
              this->simulation_parameters.post_processing.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing.enstrophy_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          enstrophy_table.set_precision("time", 12);
          enstrophy_table.set_precision("enstrophy", 12);
          this->enstrophy_table.write_text(output);
        }
    }

  // Pressure power
  if (this->simulation_parameters.post_processing.calculate_pressure_power)
    {
      TimerOutput::Scope t(this->computing_timer, "Calculate pressure power");

      const double pressure_power =
        calculate_pressure_power(this->dof_handler,
                                 present_solution,
                                 *this->cell_quadrature,
                                 *this->get_mapping());

      this->pressure_power_table.add_value(
        "time", simulation_control->get_current_time());
      this->pressure_power_table.add_value("pressure_power", pressure_power);

      // Display pressure power to screen if verbosity is enabled
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Pressure power: " << pressure_power << std::endl;
        }

      // Output pressure power to a text file from processor 0
      if (simulation_control->get_step_number() %
              this->simulation_parameters.post_processing.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing.pressure_power_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          pressure_power_table.set_precision("time", 12);
          pressure_power_table.set_precision("pressure_power", 12);
          this->pressure_power_table.write_text(output);
        }
    }


  // Viscous dissipation
  if (this->simulation_parameters.post_processing.calculate_viscous_dissipation)
    {
      TimerOutput::Scope t(this->computing_timer,
                           "Calculate viscous dissipation");

      const double viscous_dissipation = calculate_viscous_dissipation(
        this->dof_handler,
        present_solution,
        *this->cell_quadrature,
        *this->get_mapping(),
        simulation_parameters.physical_properties_manager);

      this->viscous_dissipation_table.add_value(
        "time", simulation_control->get_current_time());
      this->viscous_dissipation_table.add_value("viscous_dissipation",
                                                viscous_dissipation);

      // Display pressure power to screen if verbosity is enabled
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Viscous dissipation: " << viscous_dissipation
                      << std::endl;
        }

      // Output pressure power to a text file from processor 0
      if (simulation_control->get_step_number() %
              this->simulation_parameters.post_processing.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing
              .viscous_dissipation_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          viscous_dissipation_table.set_precision("time", 12);
          viscous_dissipation_table.set_precision("viscous_dissipation", 12);
          this->viscous_dissipation_table.write_text(output);
        }
    }

  // The average velocities and reynolds stresses are calculated when the
  // time reaches the initial time. (time >= initial time) with 1e-6 as
  // tolerance.
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      TimerOutput::Scope t(this->computing_timer,
                           "Calculate average velocities");
      this->average_velocities->calculate_average_velocities(
        this->local_evaluation_point,
        simulation_parameters.post_processing,
        simulation_control->get_current_time(),
        simulation_control->get_time_step());
    }

  if (this->simulation_parameters.post_processing.calculate_kinetic_energy)
    {
      TimerOutput::Scope t(this->computing_timer, "Calculate kinetic energy");

      double kE = calculate_kinetic_energy(this->dof_handler,
                                           present_solution,
                                           *this->cell_quadrature,
                                           *this->get_mapping());
      this->kinetic_energy_table.add_value(
        "time", simulation_control->get_current_time());
      this->kinetic_energy_table.add_value("kinetic-energy", kE);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Kinetic energy: " << kE << std::endl;
        }

      // Output Kinetic Energy to a text file from processor 0
      if ((simulation_control->get_step_number() %
             this->simulation_parameters.post_processing.output_frequency ==
           0) &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing.kinetic_energy_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          kinetic_energy_table.set_precision("time", 12);
          kinetic_energy_table.set_precision("kinetic-energy", 12);
          this->kinetic_energy_table.write_text(output);
        }
    }

  // Calculate apparent viscosity
  if (this->simulation_parameters.post_processing.calculate_apparent_viscosity)
    {
      TimerOutput::Scope t(this->computing_timer,
                           "Calculate apparent viscosity");

      double apparent_viscosity = calculate_apparent_viscosity(
        this->dof_handler,
        this->present_solution,
        *this->cell_quadrature,
        *this->get_mapping(),
        this->simulation_parameters.physical_properties_manager);

      this->apparent_viscosity_table.add_value(
        "time", simulation_control->get_current_time());
      this->apparent_viscosity_table.add_value("apparent_viscosity",
                                               apparent_viscosity);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Apparent viscosity: " << apparent_viscosity
                      << " m^2/s" << std::endl;
        }

      // Output apparent viscosity to a text file from processor 0
      if (simulation_control->get_step_number() %
              this->simulation_parameters.post_processing.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing
              .apparent_viscosity_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          apparent_viscosity_table.set_precision("time", 12);
          apparent_viscosity_table.set_precision("apparent_viscosity", 12);
          this->apparent_viscosity_table.write_text(output);
        }
    }

  // Calculate pressure drop between two boundaries
  if (this->simulation_parameters.post_processing.calculate_pressure_drop)
    {
      TimerOutput::Scope t(this->computing_timer, "Calculate pressure drop");

      double pressure_drop, total_pressure_drop;
      std::tie(pressure_drop, total_pressure_drop) = calculate_pressure_drop(
        this->dof_handler,
        *this->get_mapping(),
        this->evaluation_point,
        *this->cell_quadrature,
        *this->face_quadrature,
        this->simulation_parameters.post_processing.inlet_boundary_id,
        this->simulation_parameters.post_processing.outlet_boundary_id);
      this->pressure_drop_table.add_value(
        "time", simulation_control->get_current_time());
      this->pressure_drop_table.add_value("pressure-drop", pressure_drop);
      this->pressure_drop_table.add_value("total-pressure-drop",
                                          total_pressure_drop);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Pressure drop: "
                      << std::setprecision(
                           simulation_control->get_log_precision())
                      << this->simulation_parameters.physical_properties_manager
                             .get_density_scale() *
                           pressure_drop
                      << " Pa" << std::endl;
          this->pcout << "Total pressure drop: "
                      << std::setprecision(
                           simulation_control->get_log_precision())
                      << this->simulation_parameters.physical_properties_manager
                             .get_density_scale() *
                           total_pressure_drop
                      << " Pa" << std::endl;
        }

      // Output pressure drop to a text file from processor 0
      if ((simulation_control->get_step_number() %
             this->simulation_parameters.post_processing.output_frequency ==
           0) &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing.pressure_drop_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          pressure_drop_table.set_precision("time", 12);
          pressure_drop_table.set_precision("pressure-drop", 12);
          pressure_drop_table.set_precision("total-pressure-drop", 12);
          this->pressure_drop_table.write_text(output);
        }
    }

  // Calculate flow rate at every boundary
  if (this->simulation_parameters.post_processing.calculate_flow_rate)
    {
      TimerOutput::Scope t(this->computing_timer, "Calculate flow rate");

      this->flow_rate_table.add_value("time",
                                      simulation_control->get_current_time());
      this->flow_rate_table.set_scientific("time", true);

      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          announce_string(this->pcout, "Flow rates");
        }

      for (auto const &[boundary_id, type] :
           simulation_parameters.boundary_conditions.type)
        {
          std::pair<double, double> boundary_flow_rate =
            calculate_flow_rate(this->dof_handler,
                                this->present_solution,
                                boundary_id,
                                *this->face_quadrature,
                                *this->get_mapping());

          this->flow_rate_table.add_value(
            "flow-rate-" + Utilities::int_to_string(boundary_id, 2),
            boundary_flow_rate.first);
          this->flow_rate_table.set_scientific(
            "flow-rate-" + Utilities::int_to_string(boundary_id, 2), true);
          if (this->simulation_parameters.post_processing.verbosity ==
              Parameters::Verbosity::verbose)
            {
              this->pcout << "Flow rate at boundary " +
                               std::to_string(boundary_id) + ": "
                          << std::setprecision(
                               simulation_control->get_log_precision())
                          << boundary_flow_rate.first << " m^3/s" << std::endl;
            }
        }

      // Output flow rate to a text file from processor 0
      if ((simulation_control->get_step_number() %
             this->simulation_parameters.post_processing.output_frequency ==
           0) &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.post_processing.flow_rate_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          flow_rate_table.set_precision("time", 12);
          for (auto const &[boundary_id, type] :
               simulation_parameters.boundary_conditions.type)
            flow_rate_table.set_precision(
              "flow-rate-" + Utilities::int_to_string(boundary_id, 2), 12);
          this->flow_rate_table.write_text(output);
        }
    }

  if (!firstIter)
    {
      // Calculate forces on the boundary conditions
      if (this->simulation_parameters.forces_parameters.calculate_force)
        {
          if (simulation_control->get_step_number() %
                this->simulation_parameters.forces_parameters
                  .calculation_frequency ==
              0)
            this->postprocessing_forces(present_solution);
          if (simulation_control->get_step_number() %
                this->simulation_parameters.forces_parameters
                  .output_frequency ==
              0)
            this->write_output_forces();
        }

      // Calculate torques on the boundary conditions
      if (this->simulation_parameters.forces_parameters.calculate_torque)
        {
          if (simulation_control->get_step_number() %
                this->simulation_parameters.forces_parameters
                  .calculation_frequency ==
              0)
            this->postprocessing_torques(present_solution);
          if (simulation_control->get_step_number() %
                this->simulation_parameters.forces_parameters
                  .output_frequency ==
              0)
            this->write_output_torques();
        }

      // Calculate error with respect to analytical solution
      if (this->simulation_parameters.analytical_solution->calculate_error())
        {
          TimerOutput::Scope t(this->computing_timer,
                               "Calculate error w.r.t. analytical solution");

          // Update the time of the exact solution to the actual time
          this->exact_solution->set_time(
            simulation_control->get_current_time());

          present_solution.update_ghost_values();

          const std::pair<double, double> errors =
            calculate_L2_error(dof_handler,
                               present_solution,
                               exact_solution,
                               *this->cell_quadrature,
                               *this->get_mapping());
          const double error_velocity = errors.first;
          const double error_pressure = errors.second;
          if (simulation_parameters.simulation_control.method ==
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

              if (this->simulation_parameters.timer.write_time_in_error_table)
                {
                  auto summary = computing_timer.get_summary_data(
                    computing_timer.total_wall_time);
                  double total_time = 0;
                  for (auto it = summary.begin(); it != summary.end(); ++it)
                    {
                      total_time += summary[it->first];
                    }
                  this->error_table.add_value("total_time", total_time);
                }

              // Calculate error on pressure for VOF or Cahn-Hilliard
              // simulations
              if (this->simulation_parameters.multiphysics.VOF ||
                  this->simulation_parameters.multiphysics.cahn_hilliard)
                this->error_table.add_value("error_pressure", error_pressure);
            }
          if (this->simulation_parameters.analytical_solution->verbosity ==
              Parameters::Verbosity::verbose)
            {
              this->pcout << "L2 error velocity: "
                          << std::setprecision(
                               simulation_control->get_log_precision())
                          << error_velocity << std::endl;
              if (this->simulation_parameters.multiphysics.cahn_hilliard)
                {
                  this->pcout << "L2 error pressure: "
                              << std::setprecision(
                                   simulation_control->get_log_precision())
                              << error_pressure << std::endl;
                }
            }
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::set_nodal_values()
{
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  VectorTools::interpolate(*this->get_mapping(),
                           this->dof_handler,
                           this->simulation_parameters.initial_condition->uvwp,
                           this->newton_update,
                           this->fe->component_mask(velocities));
  VectorTools::interpolate(*this->get_mapping(),
                           this->dof_handler,
                           this->simulation_parameters.initial_condition->uvwp,
                           this->newton_update,
                           this->fe->component_mask(pressure));
  this->nonzero_constraints.distribute(this->newton_update);
  this->present_solution = this->newton_update;
  if (this->simulation_parameters.simulation_control.bdf_startup_method ==
      Parameters::SimulationControl::BDFStartupMethods::initial_solution)
    {
      for (unsigned int i = 1; i < this->previous_solutions.size(); ++i)
        {
          double previous_solution_time =
            -this->simulation_parameters.simulation_control.dt * i;
          this->simulation_parameters.initial_condition->uvwp.set_time(
            previous_solution_time);
          const FEValuesExtractors::Vector velocities(0);
          const FEValuesExtractors::Scalar pressure(dim);
          VectorTools::interpolate(
            *this->get_mapping(),
            this->dof_handler,
            this->simulation_parameters.initial_condition->uvwp,
            this->newton_update,
            this->fe->component_mask(velocities));
          VectorTools::interpolate(
            *this->get_mapping(),
            this->dof_handler,
            this->simulation_parameters.initial_condition->uvwp,
            this->newton_update,
            this->fe->component_mask(pressure));
          this->previous_solutions[i - 1] = this->newton_update;
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::define_non_zero_constraints()
{
  double time = this->simulation_control->get_current_time();
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);

  // Non-zero constraints
  auto &nonzero_constraints = this->get_nonzero_constraints();
  nonzero_constraints.clear();

  if constexpr (std::is_same_v<VectorType, GlobalBlockVectorType>)
    {
      std::vector<unsigned int> block_component(dim + 1, 0);
      block_component[dim] = 1;
      DoFRenumbering::component_wise(this->dof_handler, block_component);
      std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(this->dof_handler, block_component);

      unsigned int dof_u = dofs_per_block[0];
      unsigned int dof_p = dofs_per_block[1];

      IndexSet locally_relevant_dofs_acquisition;
      locally_relevant_dofs_acquisition =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);
      this->locally_relevant_dofs.resize(2);
      this->locally_relevant_dofs[0] =
        locally_relevant_dofs_acquisition.get_view(0, dof_u);
      this->locally_relevant_dofs[1] =
        locally_relevant_dofs_acquisition.get_view(dof_u, dof_u + dof_p);
      nonzero_constraints.reinit(this->dof_handler.locally_owned_dofs(),
                                 locally_relevant_dofs_acquisition);
    }
  else
    {
      nonzero_constraints.reinit(this->dof_handler.locally_owned_dofs(),
                                 this->locally_relevant_dofs);
    }

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          nonzero_constraints);
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::noslip)
        {
          VectorTools::interpolate_boundary_values(
            *this->get_mapping(),
            this->dof_handler,
            id,
            dealii::Functions::ZeroFunction<dim>(dim + 1),
            nonzero_constraints,
            this->fe->component_mask(velocities));
        }
      else if (type == BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> no_normal_flux_boundaries;
          no_normal_flux_boundaries.insert(id);
          VectorTools::compute_no_normal_flux_constraints(
            this->dof_handler,
            0,
            no_normal_flux_boundaries,
            nonzero_constraints,
            *this->get_mapping());
        }
      else if (type == BoundaryConditions::BoundaryType::function)
        {
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->u.set_time(time);
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->v.set_time(time);
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->w.set_time(time);
          VectorTools::interpolate_boundary_values(
            *this->get_mapping(),
            this->dof_handler,
            id,
            NavierStokesFunctionDefined<dim>(
              &this->simulation_parameters.boundary_conditions
                 .navier_stokes_functions.at(id)
                 ->u,
              &this->simulation_parameters.boundary_conditions
                 .navier_stokes_functions.at(id)
                 ->v,
              &this->simulation_parameters.boundary_conditions
                 .navier_stokes_functions.at(id)
                 ->w),
            nonzero_constraints,
            this->fe->component_mask(velocities));
        }
      else if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions.periodic_neighbor_id
              .at(id),
            this->simulation_parameters.boundary_conditions.periodic_direction
              .at(id),
            nonzero_constraints);
        }
    }

  this->establish_solid_domain(true);

  nonzero_constraints.close();
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::define_zero_constraints()
{
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  this->zero_constraints.clear();

  if constexpr (std::is_same_v<VectorType, GlobalBlockVectorType>)
    {
      std::vector<unsigned int> block_component(dim + 1, 0);
      block_component[dim] = 1;
      DoFRenumbering::component_wise(this->dof_handler, block_component);
      std::vector<types::global_dof_index> dofs_per_block =
        DoFTools::count_dofs_per_fe_block(this->dof_handler, block_component);

      unsigned int dof_u = dofs_per_block[0];
      unsigned int dof_p = dofs_per_block[1];

      IndexSet locally_relevant_dofs_acquisition;
      locally_relevant_dofs_acquisition =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);
      this->locally_relevant_dofs.resize(2);
      this->locally_relevant_dofs[0] =
        locally_relevant_dofs_acquisition.get_view(0, dof_u);
      this->locally_relevant_dofs[1] =
        locally_relevant_dofs_acquisition.get_view(dof_u, dof_u + dof_p);
      this->zero_constraints.reinit(this->dof_handler.locally_owned_dofs(),
                                    locally_relevant_dofs_acquisition);
    }
  else
    {
      this->locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(this->dof_handler);
      this->zero_constraints.reinit(this->dof_handler.locally_owned_dofs(),
                                    this->locally_relevant_dofs);
    }

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::noslip ||
          type == BoundaryConditions::BoundaryType::function)
        {
          VectorTools::interpolate_boundary_values(
            *this->get_mapping(),
            this->dof_handler,
            id,
            dealii::Functions::ZeroFunction<dim>(dim + 1),
            this->zero_constraints,
            this->fe->component_mask(velocities));
        }
      else if (type == BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> no_normal_flux_boundaries;
          no_normal_flux_boundaries.insert(id);
          VectorTools::compute_no_normal_flux_constraints(
            this->dof_handler,
            0,
            no_normal_flux_boundaries,
            this->zero_constraints,
            *this->get_mapping());
        }
      else if (type == BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            id,
            this->simulation_parameters.boundary_conditions.periodic_neighbor_id
              .at(id),
            this->simulation_parameters.boundary_conditions.periodic_direction
              .at(id),
            this->zero_constraints);
        }
      else if (type == BoundaryConditions::BoundaryType::pressure)
        {
          /*The pressure boundary condition is implemented in the matrix-based
           * assemblers*/
        }
      else if (type == BoundaryConditions::BoundaryType::function_weak)
        {
          /*The function weak boundary condition is implemented in the
           * matrix-based assemblers and the matrix-free operators*/
        }
      else if (type == BoundaryConditions::BoundaryType::partial_slip)
        {
          /*The partial slip boundary condition is implemented in the
           * matrix-based assemblers*/
        }
      else if (type == BoundaryConditions::BoundaryType::outlet)
        {
          /*The directional do-nothing boundary condition is implemented
           * in the matrix-based assemblers and the matrix-free operators*/
        }
      else if (type == BoundaryConditions::BoundaryType::none)
        {
          /*Default boundary condition*/
        }
    }

  this->establish_solid_domain(false);

  this->zero_constraints.close();
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::update_mortar_configuration()
{
  if (!this->simulation_parameters.mortar_parameters.enable)
    return;

  bool refinement_step;
  if (this->simulation_parameters.mesh_adaptation.refinement_at_frequency)
    refinement_step = this->simulation_control->get_step_number() %
                        this->simulation_parameters.mesh_adaptation.frequency ==
                      0;
  else
    refinement_step = this->simulation_control->get_step_number() == 0;

  if (this->simulation_control->is_at_start() || !refinement_step ||
      this->simulation_parameters.mesh_adaptation.type ==
        Parameters::MeshAdaptation::Type::none)
    {
      setup_dofs();

      // Set up the vectors for the transfer
      VectorType tmp = init_temporary_vector();
      tmp            = this->present_solution;

      if constexpr (std::is_same_v<VectorType,
                                   LinearAlgebra::distributed::Vector<double>>)
        tmp.update_ghost_values();

      // Distribute constraints
      auto &nonzero_constraints = this->nonzero_constraints;
      nonzero_constraints.distribute(tmp);

      // Fix on the new mesh
      this->present_solution = tmp;

      for (unsigned int i = 0; i < previous_solutions.size(); ++i)
        {
          VectorType tmp_previous_solution = init_temporary_vector();
          if constexpr (std::is_same_v<
                          VectorType,
                          LinearAlgebra::distributed::Vector<double>>)
            tmp_previous_solution.update_ghost_values();

          nonzero_constraints.distribute(tmp_previous_solution);
          previous_solutions[i] = tmp_previous_solution;
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::reinit_mortar()
{
  if (!this->simulation_parameters.mortar_parameters.enable)
    return;

  rotate_rotor_mapping(true);

  TimerOutput::Scope t(this->computing_timer, "Reinit mortar");

  // Create mortar manager
  this->mortar_manager = std::make_shared<MortarManagerCircle<dim>>(
    *this->cell_quadrature,
    this->dof_handler,
    this->simulation_parameters.mortar_parameters);

  // Create mortar coupling evaluator
  this->mortar_coupling_evaluator =
    std::make_shared<NavierStokesCouplingEvaluation<dim, double>>(
      *this->get_mapping(),
      this->dof_handler,
      this->simulation_parameters.physical_properties_manager
        .get_kinematic_viscosity_scale());

  this->mortar_coupling_operator =
    std::make_shared<CouplingOperator<dim, double>>(
      *this->get_mapping(),
      this->dof_handler,
      this->zero_constraints,
      this->mortar_coupling_evaluator,
      this->mortar_manager,
      this->simulation_parameters.mortar_parameters.rotor_boundary_id,
      this->simulation_parameters.mortar_parameters.stator_boundary_id,
      this->simulation_parameters.mortar_parameters.sip_factor);
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::rotate_rotor_mapping(bool is_first)
{
  if (!this->simulation_parameters.mortar_parameters.enable)
    return;

  TimerOutput::Scope t(this->computing_timer, "Rotate mortar");

  // Get updated rotation angle (radians)
  simulation_parameters.mortar_parameters.rotor_rotation_angle->set_time(
    this->simulation_control->get_current_time());
  const double rotation_angle =
    simulation_parameters.mortar_parameters.rotor_rotation_angle->value(
      Point<dim>());

  // Get updated angular velocity (radians/time)
  simulation_parameters.mortar_parameters.rotor_angular_velocity->set_time(
    this->simulation_control->get_current_time());
  const double angular_velocity =
    simulation_parameters.mortar_parameters.rotor_angular_velocity->value(
      Point<dim>());

  if (simulation_parameters.mortar_parameters.verbosity ==
        Parameters::Verbosity::verbose &&
      is_first)
    {
      this->pcout << "Mortar - Rotor grid angle is: " << rotation_angle
                  << " rad \n"
                  << "         Rotor grid velocity is: " << angular_velocity
                  << " rad/time \n"
                  << std::endl;
    }

  // Create new mapping cache
  this->mapping_cache =
    std::make_shared<MappingQCache<dim>>(this->velocity_fem_degree);

  LetheGridTools::rotate_mapping(
    this->dof_handler,
    *this->mapping_cache,
    *this->mapping,
    compute_n_subdivisions_and_radius(
      *this->triangulation, this->simulation_parameters.mortar_parameters)
      .second,
    rotation_angle);
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::update_boundary_conditions()
{
  if (!this->simulation_parameters.boundary_conditions.time_dependent)
    return;

  // We can never assume in the code anywhere that the local_evaluation_point is
  // at the right value its value must always be reinitialized from the present
  // solution. This may appear trivial, but this is extremely important when we
  // are checkpointing. Trust me future Bruno.
  this->local_evaluation_point = this->present_solution;

  double time = this->simulation_control->get_current_time();

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      // Only set the time if it actually makes sense to set the time.
      if (type == BoundaryConditions::BoundaryType::function ||
          type == BoundaryConditions::BoundaryType::function_weak ||
          type == BoundaryConditions::BoundaryType::pressure ||
          type == BoundaryConditions::BoundaryType::partial_slip)
        {
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->u.set_time(time);
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->v.set_time(time);
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->w.set_time(time);
          this->simulation_parameters.boundary_conditions
            .navier_stokes_functions.at(id)
            ->p.set_time(time);
        }
    }
  this->define_non_zero_constraints();
  // Distribute constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "Read checkpoint");
  std::string        prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  this->simulation_control->read(prefix);
  this->pvdhandler.read(prefix);

  this->set_solution_from_checkpoint(prefix);

  // Calculate the initial condition for the average velocity profile
  if (simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->average_velocities->calculate_average_velocities(
        this->local_evaluation_point,
        simulation_parameters.post_processing,
        simulation_control->get_current_time(),
        simulation_control->get_time_step());
    }

  if (simulation_parameters.flow_control.enable_flow_control)
    {
      this->flow_control.read(prefix);
    }

  multiphysics->read_checkpoint();

  // Deserialize all post-processing tables that are currently used
  {
    const Parameters::PostProcessing post_processing =
      this->simulation_parameters.post_processing;
    std::string prefix =
      this->simulation_parameters.simulation_control.output_folder;
    std::string suffix = ".checkpoint";
    if (post_processing.calculate_enstrophy)
      deserialize_table(this->enstrophy_table,
                        prefix + post_processing.enstrophy_output_name +
                          suffix);
    if (post_processing.calculate_pressure_power)
      deserialize_table(this->pressure_power_table,
                        prefix + post_processing.pressure_power_output_name +
                          suffix);
    if (post_processing.calculate_viscous_dissipation)
      deserialize_table(this->viscous_dissipation_table,
                        prefix +
                          post_processing.viscous_dissipation_output_name +
                          suffix);
    if (post_processing.calculate_kinetic_energy)
      deserialize_table(this->kinetic_energy_table,
                        prefix + post_processing.kinetic_energy_output_name +
                          suffix);
    if (post_processing.calculate_apparent_viscosity)
      deserialize_table(this->apparent_viscosity_table,
                        prefix +
                          post_processing.apparent_viscosity_output_name +
                          suffix);
    if (post_processing.calculate_flow_rate)
      deserialize_table(this->flow_rate_table,
                        prefix + post_processing.flow_rate_output_name +
                          suffix);
    if (post_processing.calculate_pressure_drop)
      deserialize_table(this->pressure_drop_table,
                        prefix + post_processing.pressure_drop_output_name +
                          suffix);
    if (this->simulation_parameters.forces_parameters.calculate_force)
      for (auto const &[id, type] :
           this->simulation_parameters.boundary_conditions.type)
        {
          deserialize_table(
            this->forces_tables[id],
            prefix +
              this->simulation_parameters.forces_parameters.force_output_name +
              "_" + Utilities::int_to_string(id, 2) + suffix);
        }
    if (this->simulation_parameters.forces_parameters.calculate_torque)
      for (auto const &[id, type] :
           this->simulation_parameters.boundary_conditions.type)
        {
          deserialize_table(
            this->torques_tables[id],
            prefix +
              this->simulation_parameters.forces_parameters.torque_output_name +
              "_" + Utilities::int_to_string(id, 2) + suffix);
        }
    if (this->simulation_parameters.analytical_solution->calculate_error())
      deserialize_table(
        this->error_table,
        prefix +
          this->simulation_parameters.analytical_solution->get_filename() +
          "_FD" + suffix);
  }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::set_solution_from_checkpoint(
  std::string checkpoint_file_prefix)
{
  const std::string filename = checkpoint_file_prefix + ".triangulation";
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
        {
          tria->load(filename.c_str());
        }
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }


  setup_dofs();
  enable_dynamic_zero_constraints_fd();

  // BB note: There is an issue right now that will prevent this code from
  // running in debug mode with Trilinos vectors Deal.II vectors require that
  // the vectors used in the checkpointing mechanism have their relevant dofs
  // whereas Trilinos vectors do not allow for this. Right now this code works
  // well in release mode for both vector types, but will not work in debug mode
  // for Trilinos vectors because of an assertion. A workaround will be
  // implemented in a near future

  std::vector<VectorType *> x_system(1 + previous_solutions.size());

  VectorType distributed_system(locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);
  x_system[0] = &(distributed_system);

  std::vector<VectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        VectorType(locally_owned_dofs,
                   this->locally_relevant_dofs,
                   this->mpi_communicator));
      x_system[i + 1] = &distributed_previous_solutions[i];
    }
  SolutionTransfer<dim, VectorType> system_trans_vectors(this->dof_handler);

  if (simulation_parameters.post_processing.calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    {
      std::vector<VectorType *> sum_vectors =
        this->average_velocities->read(checkpoint_file_prefix);

      x_system.insert(x_system.end(), sum_vectors.begin(), sum_vectors.end());
    }

  system_trans_vectors.deserialize(x_system);
  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }

  // Reset the average velocity profile if the initial time to average the
  // velocities has not been reached. Disabled if the initial condition is an
  // average velocity profile.
  if (simulation_parameters.post_processing.calculate_average_velocities &
      (this->simulation_parameters.initial_condition->type !=
       Parameters::InitialConditionType::average_velocity_profile))
    {
      if ((this->simulation_parameters.post_processing
             .initial_time_for_average_velocities +
           1e-6 * simulation_control->get_time_step()) >
          this->simulation_control->get_current_time())
        {
          this->pcout
            << "Warning: The checkpointed time-averaged velocity has been reinitialized because the initial averaging time has not yet been reached."
            << std::endl;
          this->average_velocities->zero_average_after_restart();
        }
    }

  if (simulation_parameters.post_processing.calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    {
      this->average_velocities->sanitize_after_restart();
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::establish_solid_domain(
  const bool non_zero_constraints)
{
  // If there are no solid regions, there is no work to be done and we can
  // return.
  if (simulation_parameters.physical_properties_manager
        .get_number_of_solids() == 0)
    return;

  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // We will need to identify which pressure degrees of freedom are connected to
  // fluid region. For these, we won't establish a zero pressure constraint.
  std::unordered_set<types::global_dof_index> dofs_are_connected_to_fluid;

  // Loop through all cells to identify which cells are solid. This first step
  // is used to 1) constraint the velocity degree of freedom to be zero in the
  // solid region and 2) to identify which pressure degrees of freedom are
  // connected to fluid cells
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          // If the material_id is higher than 0, the region is a solid region.
          // Constrain the velocity DOFs to be zero.
          if (cell->material_id() > 0)
            {
              constrain_solid_cell_velocity_dofs(non_zero_constraints,
                                                 local_dof_indices,
                                                 this->zero_constraints);
            }
          else
            {
              // Cell is a fluid cell and as such all the pressure DOFs of that
              // cell are connected to the fluid. This will be used later on to
              // identify which pressure cells to constrain.
              flag_dofs_connected_to_fluid(local_dof_indices,
                                           dofs_are_connected_to_fluid);
            }
        }
    }

  // All pressure DOFs that are not connected to a fluid cell are constrained
  // to ensure that the system matrix has adequate conditioning.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          // If the material_id is > 0, the region is a solid region
          if (cell->material_id() > 0)
            {
              // First check if the cell is connected to a fluid cell by
              // checking if one of the DOF of the cell is connected to a fluid
              // cell.
              bool connected_to_fluid =
                check_cell_is_connected_to_fluid(dofs_are_connected_to_fluid,
                                                 local_dof_indices);

              // All the pressure DOFs with the cell are not connected to the
              // fluid. Consequently, we fix a constraint on these pressure
              // DOFs.
              if (!connected_to_fluid)
                {
                  constrain_pressure(non_zero_constraints,
                                     local_dof_indices,
                                     this->zero_constraints);
                }
            }
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::constrain_stasis_with_temperature(
  const DoFHandler<dim> *dof_handler_ht)
{
  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Get struct containing temperature range information and flag
  // containers for DOFs.
  StasisConstraintWithTemperature &stasis_constraint_struct =
    this->stasis_constraint_structs[0];

  // Get domain restriction with plane information
  bool restrain_domain_with_plane =
    this->simulation_parameters.constrain_solid_domain
      .enable_domain_restriction_with_plane;
  Point<dim> plane_point(
    this->simulation_parameters.constrain_solid_domain.restriction_plane_point);
  Tensor<1, dim> plane_normal_vector(
    this->simulation_parameters.constrain_solid_domain
      .restriction_plane_normal_vector);

  // Get temperature solution
  const auto temperature_solution =
    *this->multiphysics->get_solution(PhysicsID::heat_transfer);
  std::vector<double> local_temperature_values(this->cell_quadrature->size());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);

          // If a restriction plane is defined, check if the cell is in the
          // valid domain.
          if (!restrain_domain_with_plane ||
              cell_in_constraining_domain(cell,
                                          plane_point,
                                          plane_normal_vector))
            {
              get_cell_temperature_values(cell,
                                          dof_handler_ht,
                                          temperature_solution,
                                          local_temperature_values);
              identify_cell_and_constrain_velocity(local_dof_indices,
                                                   local_temperature_values,
                                                   stasis_constraint_struct);
            }
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::
  constrain_stasis_with_temperature_vof(const DoFHandler<dim> *dof_handler_vof,
                                        const DoFHandler<dim> *dof_handler_ht)
{
  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Get domain restriction with plane information
  bool restrain_domain_with_plane =
    this->simulation_parameters.constrain_solid_domain
      .enable_domain_restriction_with_plane;
  Point<dim> plane_point(
    this->simulation_parameters.constrain_solid_domain.restriction_plane_point);
  Tensor<1, dim> plane_normal_vector(
    this->simulation_parameters.constrain_solid_domain
      .restriction_plane_normal_vector);

  // Get filtered phase fraction solution
  const auto filtered_phase_fraction_solution =
    *this->multiphysics->get_filtered_solution(PhysicsID::VOF);
  std::vector<double> local_filtered_phase_fraction_values(
    this->cell_quadrature->size());

  // Get temperature solution
  const auto temperature_solution =
    *this->multiphysics->get_solution(PhysicsID::heat_transfer);
  std::vector<double> local_temperature_values(this->cell_quadrature->size());

  // Loop over structs containing fluid id, temperature and phase fraction range
  // information, and flag containers for DOFs.
  for (StasisConstraintWithTemperature &stasis_constraint_struct :
       this->stasis_constraint_structs)
    {
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned() || cell->is_ghost())
            {
              cell->get_dof_indices(local_dof_indices);

              // If a restriction plane is defined, check if the cell is in the
              // valid domain.
              if (!restrain_domain_with_plane ||
                  cell_in_constraining_domain(cell,
                                              plane_point,
                                              plane_normal_vector))
                {
                  bool cell_is_in_right_fluid = true;
                  get_cell_filtered_phase_fraction_values(
                    cell,
                    dof_handler_vof,
                    filtered_phase_fraction_solution,
                    local_filtered_phase_fraction_values);

                  // Check if cell is only in the fluid of interest. As soon as
                  // one filtered phase fraction value is outside the tolerated
                  // range, the cell is perceived as being in the wrong fluid.
                  for (const double &filtered_phase_fraction :
                       local_filtered_phase_fraction_values)
                    {
                      if (abs(stasis_constraint_struct.fluid_id -
                              filtered_phase_fraction) >=
                          stasis_constraint_struct
                            .filtered_phase_fraction_tolerance)
                        {
                          cell_is_in_right_fluid = false;
                          break;
                        }
                    }

                  // If the cell is not in the right fluid, no solid constraint
                  // will be applied on the cell's DOFs; we skip to the next
                  // cell.
                  if (!cell_is_in_right_fluid)
                    continue;

                  get_cell_temperature_values(cell,
                                              dof_handler_ht,
                                              temperature_solution,
                                              local_temperature_values);
                  identify_cell_and_constrain_velocity(
                    local_dof_indices,
                    local_temperature_values,
                    stasis_constraint_struct);
                }
            }
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::
  identify_cell_and_constrain_velocity(
    const std::vector<types::global_dof_index> &local_dof_indices,
    const std::vector<double>                  &local_temperature_values,
    StasisConstraintWithTemperature            &stasis_constraint_struct)
{
  for (const double &temperature : local_temperature_values)
    {
      // Skip cells with at least 1 DOF that is outbound the temperature limits.
      if (temperature < stasis_constraint_struct.min_solid_temperature ||
          temperature > stasis_constraint_struct.max_solid_temperature)
        return;
    }
  constrain_solid_cell_velocity_dofs(false,
                                     local_dof_indices,
                                     this->dynamic_zero_constraints);
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
  TimerOutput::Scope t(this->computing_timer, "Output VTU");

  const std::string  folder        = simulation_control->get_output_path();
  const std::string  solution_name = simulation_control->get_output_name();
  const unsigned int iter          = simulation_control->get_step_number();
  const double       time          = simulation_control->get_current_time();
  const unsigned int subdivision = simulation_control->get_number_subdivision();
  const unsigned int group_files = simulation_control->get_group_files();

  // Add the interpretation of the solution. The dim first components are the
  // velocity vectors and the following one is the pressure.
  std::vector<std::string> solution_names(dim, "velocity");
  solution_names.emplace_back("pressure");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.emplace_back(
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

  if (this->simulation_parameters.post_processing
        .calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    {
      // Add the interpretation of the average solution. The dim first
      // components are the average velocity vectors and the following one is
      // the average pressure. (<u>, <v>, <w>, <p>)
      std::vector<std::string> average_solution_names(dim, "average_velocity");
      average_solution_names.emplace_back("average_pressure");

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        average_data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
      average_data_component_interpretation.emplace_back(
        DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector(
        this->average_velocities->get_average_velocities(),
        average_solution_names,
        DataOut<dim>::type_dof_data,
        average_data_component_interpretation);

      // Add the interpretation of the reynolds stresses of solution.
      // The dim first components are the normal reynolds stress vectors and
      // the following ones are others resolved reynolds stresses.
      std::vector<std::string> reynolds_normal_stress_names(
        dim, "reynolds_normal_stress");
      reynolds_normal_stress_names.emplace_back("turbulent_kinetic_energy");
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        reynolds_normal_stress_data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
      reynolds_normal_stress_data_component_interpretation.emplace_back(
        DataComponentInterpretation::component_is_scalar);


      std::vector<std::string> reynolds_shear_stress_names = {
        "reynolds_shear_stress_uv"};
      if (dim == 2)
        {
          reynolds_shear_stress_names.emplace_back("dummy_rss_2d");
        }
      if (dim == 3)
        {
          reynolds_shear_stress_names.emplace_back("reynolds_shear_stress_vw");
          reynolds_shear_stress_names.emplace_back("reynolds_shear_stress_uw");
        }
      reynolds_shear_stress_names.emplace_back("dummy_rss");

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        reynolds_shear_stress_data_component_interpretation(
          dim, DataComponentInterpretation::component_is_scalar);
      reynolds_shear_stress_data_component_interpretation.emplace_back(
        DataComponentInterpretation::component_is_scalar);

      data_out.add_data_vector(
        this->average_velocities->get_reynolds_normal_stresses(),
        reynolds_normal_stress_names,
        DataOut<dim>::type_dof_data,
        reynolds_normal_stress_data_component_interpretation);

      data_out.add_data_vector(
        this->average_velocities->get_reynolds_shear_stresses(),
        reynolds_shear_stress_names,
        DataOut<dim>::type_dof_data,
        reynolds_shear_stress_data_component_interpretation);
    }

  Vector<float> subdomain(this->triangulation->n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = this->triangulation->locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");


  // Create the post-processors to have derived information about the velocity
  // They are generated outside the if condition for smoothing to ensure
  // that the objects still exist when the write output of DataOut is called
  // Regular discontinuous postprocessors
  QCriterionPostprocessor<dim> qcriterion;
  DivergencePostprocessor<dim> divergence;
  VorticityPostprocessor<dim>  vorticity;
  data_out.add_data_vector(solution, vorticity);

  // Get physical properties models
  std::vector<std::shared_ptr<DensityModel>> density_models =
    this->simulation_parameters.physical_properties_manager
      .get_density_vector();
  std::vector<std::shared_ptr<RheologicalModel>> rheological_models =
    this->simulation_parameters.physical_properties_manager
      .get_rheology_vector();

  const double n_fluids = this->simulation_parameters
                            .physical_properties_manager.get_number_of_fluids();

  std::vector<DensityPostprocessor<dim>> density_postprocessors;
  density_postprocessors.reserve(n_fluids);
  std::vector<KinematicViscosityPostprocessor<dim>>
    kinematic_viscosity_postprocessors;
  kinematic_viscosity_postprocessors.reserve(n_fluids);
  std::vector<DynamicViscosityPostprocessor<dim>>
    dynamic_viscosity_postprocessors;
  dynamic_viscosity_postprocessors.reserve(n_fluids);

  for (unsigned int f_id = 0; f_id < n_fluids; ++f_id)
    {
      density_postprocessors.emplace_back(
        DensityPostprocessor<dim>(density_models[f_id], f_id));
      kinematic_viscosity_postprocessors.emplace_back(
        KinematicViscosityPostprocessor<dim>(rheological_models[f_id], f_id));
      dynamic_viscosity_postprocessors.emplace_back(
        DynamicViscosityPostprocessor<dim>(
          rheological_models[f_id],
          density_models[f_id]->get_density_ref(),
          f_id));

      // Only output when density is not constant or if it is a multiphase
      // flow
      if (!density_models[f_id]->is_constant_density_model() ||
          this->simulation_parameters.multiphysics.VOF ||
          this->simulation_parameters.multiphysics.cahn_hilliard)
        data_out.add_data_vector(solution, density_postprocessors[f_id]);

      // Only output the kinematic viscosity for non-newtonian rheology
      if (rheological_models[f_id]->is_non_newtonian_rheological_model())
        {
          data_out.add_data_vector(solution,
                                   kinematic_viscosity_postprocessors[f_id]);

          // Only output the dynamic viscosity for multiphase flows
          if (this->simulation_parameters.multiphysics.VOF ||
              this->simulation_parameters.multiphysics.cahn_hilliard)
            data_out.add_data_vector(solution,
                                     dynamic_viscosity_postprocessors[f_id]);
        }
    }

  ShearRatePostprocessor<dim> shear_rate_processor;
  if (this->simulation_parameters.physical_properties_manager
        .is_non_newtonian())
    data_out.add_data_vector(solution, shear_rate_processor);

  // Trilinos vector for the smoothed output fields
  QcriterionPostProcessorSmoothing<dim, VectorType> qcriterion_smoothing(
    *this->triangulation,
    this->simulation_parameters,
    number_quadrature_points);

  ContinuityPostProcessorSmoothing<dim, VectorType> continuity_smoothing(
    *this->triangulation,
    this->simulation_parameters,
    number_quadrature_points);

  // TODO: generalize this to VectorType
  GlobalVectorType qcriterion_field;
  GlobalVectorType continuity_field;

  if (this->simulation_parameters.post_processing.smoothed_output_fields)
    {
      // Qcriterion smoothing
      {
        qcriterion_field =
          qcriterion_smoothing.calculate_smoothed_field(solution,
                                                        this->dof_handler,
                                                        this->get_mapping());

        std::vector<DataComponentInterpretation::DataComponentInterpretation>
          data_component_interpretation(
            1, DataComponentInterpretation::component_is_scalar);

        std::vector<std::string> qcriterion_name = {"qcriterion"};
        const DoFHandler<dim>   &dof_handler_qcriterion =
          qcriterion_smoothing.get_dof_handler();
        data_out.add_data_vector(dof_handler_qcriterion,
                                 qcriterion_field,
                                 qcriterion_name,
                                 data_component_interpretation);
      }
      // Continuity smoothing
      {
        continuity_field =
          continuity_smoothing.calculate_smoothed_field(solution,
                                                        this->dof_handler,
                                                        this->get_mapping());

        std::vector<DataComponentInterpretation::DataComponentInterpretation>
          data_component_interpretation(
            1, DataComponentInterpretation::component_is_scalar);

        std::vector<std::string> continuity_name = {"velocity_divergence"};
        const DoFHandler<dim>   &dof_handler_qcriterion =
          continuity_smoothing.get_dof_handler();
        data_out.add_data_vector(dof_handler_qcriterion,
                                 continuity_field,
                                 continuity_name,
                                 data_component_interpretation);
      }
    }
  else
    {
      // Use the non-smoothed version of the post-processors
      data_out.add_data_vector(solution, divergence);
      data_out.add_data_vector(solution, qcriterion);
    }


  SRFPostprocessor<dim> srf(simulation_parameters.velocity_sources.omega_x,
                            simulation_parameters.velocity_sources.omega_y,
                            simulation_parameters.velocity_sources.omega_z);

  if (simulation_parameters.velocity_sources.rotating_frame_type ==
      Parameters::VelocitySource::RotatingFrameType::srf)
    data_out.add_data_vector(solution, srf);

  output_field_hook(data_out);

  multiphysics->attach_solution_to_output(data_out);

  // Build the patches and write the output

  data_out.build_patches(*this->get_mapping(),
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

  if (simulation_control->get_output_boundaries() &&
      simulation_control->get_step_number() == 0)
    {
      DataOutFaces<dim> data_out_faces;

      // Add the additional flag to enable high-order cells output when the
      // velocity interpolation order is larger than 1
      DataOutBase::VtkFlags flags;
      if (this->velocity_fem_degree > 1)
        flags.write_higher_order_cells = true;
      data_out_faces.set_flags(flags);

      BoundaryPostprocessor<dim> boundary_id;
      data_out_faces.attach_dof_handler(this->dof_handler);
      data_out_faces.add_data_vector(solution, boundary_id);
      data_out_faces.build_patches(*this->get_mapping(), subdivision);

      write_boundaries_vtu<dim>(
        data_out_faces, folder, time, iter, this->mpi_communicator);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_forces()
{
  TimerOutput::Scope t(this->computing_timer, "Output forces");
  if (this->this_mpi_process == 0)
    {
      for (auto const &[boundary_id, type] :
           simulation_parameters.boundary_conditions.type)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.forces_parameters.force_output_name + "." +
            Utilities::int_to_string(boundary_id, 2) + ".dat";
          std::ofstream output(filename.c_str());

          forces_tables[boundary_id].write_text(output);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_torques()
{
  TimerOutput::Scope t(this->computing_timer, "Output torques");
  if (this->this_mpi_process == 0)
    {
      for (auto const &[boundary_id, type] :
           simulation_parameters.boundary_conditions.type)
        {
          std::string filename =
            simulation_parameters.simulation_control.output_folder +
            simulation_parameters.forces_parameters.torque_output_name + "." +
            Utilities::int_to_string(boundary_id, 2) + ".dat";
          std::ofstream output(filename.c_str());

          this->torques_tables[boundary_id].write_text(output);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "Write checkpoint");
  std::string        prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      simulation_control->save(prefix);
      this->pvdhandler.save(prefix);

      if (simulation_parameters.flow_control.enable_flow_control)
        this->flow_control.save(prefix);
    }

  std::vector<const VectorType *> sol_set_transfer;
  sol_set_transfer.emplace_back(&this->present_solution);
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      sol_set_transfer.emplace_back(&previous_solutions[i]);
    }

  if (simulation_parameters.post_processing.calculate_average_velocities ||
      this->simulation_parameters.initial_condition->type ==
        Parameters::InitialConditionType::average_velocity_profile)
    {
      std::vector<const VectorType *> av_set_transfer =
        this->average_velocities->save(prefix);

      // Insert average velocities vectors into the set transfer vector
      sol_set_transfer.insert(sol_set_transfer.end(),
                              av_set_transfer.begin(),
                              av_set_transfer.end());
    }

  SolutionTransfer<dim, VectorType> system_trans_vectors(this->dof_handler);
  system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  multiphysics->write_checkpoint();

  if (auto tria = dynamic_cast<parallel::distributed::Triangulation<dim> *>(
        this->triangulation.get()))
    {
      std::string triangulationName = prefix + ".triangulation";
      tria->save(prefix + ".triangulation");
    }

  // Serialize all post-processing tables that are currently used
  {
    const Parameters::PostProcessing post_processing =
      this->simulation_parameters.post_processing;
    std::string prefix =
      this->simulation_parameters.simulation_control.output_folder;
    std::string suffix = ".checkpoint";
    if (post_processing.calculate_enstrophy)
      serialize_table(this->enstrophy_table,
                      prefix + post_processing.enstrophy_output_name + suffix);
    if (post_processing.calculate_kinetic_energy)
      serialize_table(this->kinetic_energy_table,
                      prefix + post_processing.kinetic_energy_output_name +
                        suffix);
    if (post_processing.calculate_pressure_power)
      serialize_table(this->pressure_power_table,
                      prefix + post_processing.pressure_power_output_name +
                        suffix);
    if (post_processing.calculate_viscous_dissipation)
      serialize_table(this->viscous_dissipation_table,
                      prefix + post_processing.viscous_dissipation_output_name +
                        suffix);
    if (post_processing.calculate_apparent_viscosity)
      serialize_table(this->apparent_viscosity_table,
                      prefix + post_processing.apparent_viscosity_output_name +
                        suffix);
    if (post_processing.calculate_flow_rate)
      serialize_table(this->flow_rate_table,
                      prefix + post_processing.flow_rate_output_name + suffix);
    if (post_processing.calculate_pressure_drop)
      serialize_table(this->pressure_drop_table,
                      prefix + post_processing.pressure_drop_output_name +
                        suffix);
    if (this->simulation_parameters.forces_parameters.calculate_force)
      for (auto const &[boundary_id, type] :
           this->simulation_parameters.boundary_conditions.type)
        {
          serialize_table(
            this->forces_tables[boundary_id],
            prefix +
              this->simulation_parameters.forces_parameters.force_output_name +
              "_" + Utilities::int_to_string(boundary_id, 2) + suffix);
        }
    if (this->simulation_parameters.forces_parameters.calculate_torque)
      for (auto const &[boundary_id, type] :
           this->simulation_parameters.boundary_conditions.type)
        {
          serialize_table(
            this->torques_tables[boundary_id],
            prefix +
              this->simulation_parameters.forces_parameters.torque_output_name +
              "_" + Utilities::int_to_string(boundary_id, 2) + suffix);
        }
    if (this->simulation_parameters.analytical_solution->calculate_error())
      serialize_table(
        this->error_table,
        prefix +
          this->simulation_parameters.analytical_solution->get_filename() +
          "_FD" + suffix);
  }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::
  rescale_pressure_dofs_in_newton_update()
{
  const double pressure_scaling_factor =
    simulation_parameters.stabilization.pressure_scaling_factor;

  // We skip the function if the factor has a value of 1
  if (abs(pressure_scaling_factor - 1) < 1e-8)
    return;

  TimerOutput::Scope t(this->computing_timer, "Rescale pressure");

  const unsigned int                   dofs_per_cell = this->fe->dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  // Map used to keep track of which DOFs have been looped over
  std::unordered_set<unsigned int> rescaled_dofs_set;
  rescaled_dofs_set.clear();
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          for (unsigned int j = 0; j < local_dof_indices.size(); ++j)
            {
              const unsigned int global_id = local_dof_indices[j];
              // We check if we have already checked this DOF
              auto iterator = rescaled_dofs_set.find(global_id);
              if (iterator == rescaled_dofs_set.end())
                {
                  const unsigned int component_index =
                    this->fe->system_to_component_index(j).first;
                  if (this->dof_handler.locally_owned_dofs().is_element(
                        global_id) &&
                      component_index == dim)
                    {
                      this->newton_update(global_id) =
                        this->newton_update(global_id) *
                        pressure_scaling_factor;
                    }
                  rescaled_dofs_set.insert(global_id);
                }
            }
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::output_newton_update_norms(
  const unsigned int display_precision)
{
  TimerOutput::Scope t(this->computing_timer,
                       "Calculate and output norms after Newton its");

  if constexpr (std::is_same_v<VectorType, GlobalVectorType> ||
                std::is_same_v<VectorType,
                               LinearAlgebra::distributed::Vector<double>>)
    {
      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar pressure(dim);

      ComponentMask velocity_mask = fe->component_mask(velocities);
      ComponentMask pressure_mask = fe->component_mask(pressure);

      const std::vector<IndexSet> index_set_velocity =
        DoFTools::locally_owned_dofs_per_component(dof_handler, velocity_mask);
      const std::vector<IndexSet> index_set_pressure =
        DoFTools::locally_owned_dofs_per_component(dof_handler, pressure_mask);

      double local_sum = 0.0;
      double local_max = std::numeric_limits<double>::lowest();

      for (unsigned int d = 0; d < dim; ++d)
        {
          for (const auto &j : index_set_velocity[d])
            {
              double dof_newton_update = newton_update[j];
              local_sum += dof_newton_update * dof_newton_update;
              local_max = std::max(local_max, std::abs(dof_newton_update));
            }
        }

      double global_velocity_l2_norm =
        std::sqrt(Utilities::MPI::sum(local_sum, this->mpi_communicator));
      double global_velocity_linfty_norm =
        Utilities::MPI::max(local_max, this->mpi_communicator);

      local_sum = 0.0;
      local_max = std::numeric_limits<double>::lowest();

      for (const auto &j : index_set_pressure[dim])
        {
          double dof_newton_update = newton_update[j];
          local_sum += dof_newton_update * dof_newton_update;
          local_max = std::max(local_max, std::abs(dof_newton_update));
        }

      double global_pressure_l2_norm =
        std::sqrt(Utilities::MPI::sum(local_sum, this->mpi_communicator));
      double global_pressure_linfty_norm =
        Utilities::MPI::max(local_max, this->mpi_communicator);

      this->pcout << std::setprecision(display_precision)
                  << "\n\t||du||_L2 = " << std::setw(6)
                  << global_velocity_l2_norm << std::setw(6)
                  << "\t||du||_Linfty = "
                  << std::setprecision(display_precision)
                  << global_velocity_linfty_norm << std::endl;
      this->pcout << std::setprecision(display_precision)
                  << "\t||dp||_L2 = " << std::setw(6) << global_pressure_l2_norm
                  << std::setw(6) << "\t||dp||_Linfty = "
                  << std::setprecision(display_precision)
                  << global_pressure_linfty_norm << std::endl;
    }
  if constexpr (std::is_same_v<VectorType, GlobalBlockVectorType>)
    {
      this->pcout << std::setprecision(display_precision)
                  << "\t||du||_L2 = " << std::setw(6)
                  << newton_update.block(0).l2_norm() << std::setw(6)
                  << "\t||du||_Linfty = "
                  << std::setprecision(display_precision)
                  << newton_update.block(0).linfty_norm() << std::endl;
      this->pcout << std::setprecision(display_precision)
                  << "\t||dp||_L2 = " << std::setw(6)
                  << newton_update.block(1).l2_norm() << std::setw(6)
                  << "\t||dp||_Linfty = "
                  << std::setprecision(display_precision)
                  << newton_update.block(1).linfty_norm() << std::endl;
    }
}

template <int dim, typename VectorType, typename DofsType>
inline VectorType
NavierStokesBase<dim, VectorType, DofsType>::init_temporary_vector()
{
  VectorType tmp;

  if constexpr (std::is_same_v<VectorType, GlobalVectorType> ||
                std::is_same_v<VectorType, GlobalBlockVectorType>)
    tmp.reinit(locally_owned_dofs, this->mpi_communicator);

  else if constexpr (std::is_same_v<VectorType,
                                    LinearAlgebra::distributed::Vector<double>>)
    tmp.reinit(locally_owned_dofs,
               locally_relevant_dofs,
               this->mpi_communicator);

  return tmp;
}

// Pre-compile the 2D and 3D version with the types that can occur
template class NavierStokesBase<2, GlobalVectorType, IndexSet>;
template class NavierStokesBase<3, GlobalVectorType, IndexSet>;
template class NavierStokesBase<2,
                                GlobalBlockVectorType,
                                std::vector<IndexSet>>;
template class NavierStokesBase<3,
                                GlobalBlockVectorType,
                                std::vector<IndexSet>>;

#ifndef LETHE_USE_LDV
template class NavierStokesBase<2,
                                LinearAlgebra::distributed::Vector<double>,
                                IndexSet>;
template class NavierStokesBase<3,
                                LinearAlgebra::distributed::Vector<double>,
                                IndexSet>;
#endif
