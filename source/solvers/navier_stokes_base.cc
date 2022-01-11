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

#include <core/bdf.h>
#include <core/grids.h>
#include <core/lethegridtools.h>
#include <core/sdirk.h>
#include <core/solutions_output.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/flow_control.h>
#include <solvers/navier_stokes_base.h>
#include <solvers/post_processors.h>
#include <solvers/postprocessing_cfd.h>
#include <solvers/postprocessing_velocities.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>


/*
 * Constructor for the Navier-Stokes base class
 */
template <int dim, typename VectorType, typename DofsType>
NavierStokesBase<dim, VectorType, DofsType>::NavierStokesBase(
  SimulationParameters<dim> &p_nsparam)
  : PhysicsSolver<VectorType>(p_nsparam.non_linear_solver)
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
{
  if (simulation_parameters.mesh.simplex)
    {
      // for simplex meshes
      const FE_SimplexP<dim> velocity_fe(
        p_nsparam.fem_parameters.velocity_order);
      const FE_SimplexP<dim> pressure_fe(
        p_nsparam.fem_parameters.pressure_order);
      fe = std::make_shared<FESystem<dim>>(velocity_fe, dim, pressure_fe, 1);
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
      fe = std::make_shared<FESystem<dim>>(
        FE_Q<dim>(p_nsparam.fem_parameters.velocity_order),
        dim,
        FE_Q<dim>(p_nsparam.fem_parameters.pressure_order),
        1);
      mapping = std::make_shared<MappingQ<dim>>(
        velocity_fem_degree, simulation_parameters.fem_parameters.qmapping_all);
      cell_quadrature = std::make_shared<QGauss<dim>>(number_quadrature_points);
      face_quadrature =
        std::make_shared<QGauss<dim - 1>>(number_quadrature_points);
      triangulation =
        std::make_shared<parallel::distributed::Triangulation<dim>>(
          this->mpi_communicator,
          typename Triangulation<dim>::MeshSmoothing(
            Triangulation<dim>::smoothing_on_refinement |
            Triangulation<dim>::smoothing_on_coarsening));
      dof_handler.clear();
      dof_handler.reinit(*this->triangulation);
    }

  this->pcout.set_condition(
    Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0);

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
      if (simulation_parameters.simulation_control.output_control ==
          Parameters::SimulationControl::OutputControl::time)
        simulation_control =
          std::make_shared<SimulationControlTransientDynamicOutput>(
            simulation_parameters.simulation_control);
      else
        simulation_control = std::make_shared<SimulationControlTransient>(
          simulation_parameters.simulation_control);
    }

  multiphysics = std::make_shared<MultiphysicsInterface<dim>>(
    simulation_parameters, triangulation, simulation_control, this->pcout);

  // Pre-allocate memory for the previous solutions using the information
  // of the BDF schemes
  previous_solutions.resize(maximum_number_of_previous_solutions());

  // Pre-allocate memory for intermediary stages if there are any
  if (this->simulation_parameters.simulation_control.bdf_startup_method ==
        Parameters::SimulationControl::BDFStartupMethods::sdirk_step &&
      this->simulation_parameters.simulation_control.method ==
        Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    solution_stages.resize(1);
  else if (this->simulation_parameters.simulation_control.bdf_startup_method ==
             Parameters::SimulationControl::BDFStartupMethods::sdirk_step &&
           this->simulation_parameters.simulation_control.method ==
             Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    solution_stages.resize(2);
  else
    solution_stages.resize(number_of_intermediary_stages(
      simulation_parameters.simulation_control.method));


  // Change the behavior of the timer for situations when you don't want
  // outputs
  if (simulation_parameters.timer.type == Parameters::Timer::Type::none)
    this->computing_timer.disable_output();

  // Pre-allocate the force tables to match the number of boundary conditions
  forces_on_boundaries.resize(simulation_parameters.boundary_conditions.size);
  torques_on_boundaries.resize(simulation_parameters.boundary_conditions.size);
  forces_tables.resize(simulation_parameters.boundary_conditions.size);
  torques_tables.resize(simulation_parameters.boundary_conditions.size);

  // Get the exact solution from the parser
  exact_solution = &simulation_parameters.analytical_solution->uvwp;

  // If there is a forcing function, get it from the parser
  if (simulation_parameters.source_term->source_term())
    forcing_function = &simulation_parameters.source_term->navier_stokes_source;
  else
    forcing_function = new NoForce<dim>;

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
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
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_flow_rate(
  const VectorType &evaluation_point)
{
  this->flow_rate =
    calculate_flow_rate(this->dof_handler,
                        evaluation_point,
                        simulation_parameters.flow_control.boundary_flow_id,
                        *this->face_quadrature,
                        *this->mapping);

  // Showing results (area and flow rate)
  if (simulation_parameters.flow_control.verbosity ==
        Parameters::Verbosity::verbose &&
      simulation_control->get_step_number() > 0 && this->this_mpi_process == 0)
    {
      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Flow control summary                    |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      this->pcout << "Inlet area : " << flow_rate.second << std::endl;
      this->pcout << "Flow rate : " << flow_rate.first << std::endl;
      this->pcout << "Beta applied : "
                  << beta[simulation_parameters.flow_control.flow_direction]
                  << std::endl;
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::dynamic_flow_control()
{
  if (simulation_parameters.flow_control.enable_flow_control &&
      simulation_parameters.simulation_control.method !=
        Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      this->flow_control.calculate_beta(flow_rate,
                                        simulation_control->get_time_step(),
                                        simulation_control->get_step_number());

      this->beta = flow_control.get_beta();
    }
}

template <int dim, typename VectorType, typename DofsType>
bool
NavierStokesBase<dim, VectorType, DofsType>::check_existance_of_bc(
  BoundaryConditions::BoundaryType bc)
{
  bool bc_exist = false;
  // Loop over the boundary to check if they need assembler on their face.
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] == bc)
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
  TimerOutput::Scope t(this->computing_timer, "calculate_forces");

  this->forces_on_boundaries =
    calculate_forces(this->dof_handler,
                     evaluation_point,
                     simulation_parameters.physical_properties,
                     simulation_parameters.boundary_conditions,
                     *this->face_quadrature,
                     *this->mapping);

  if (simulation_parameters.forces_parameters.verbosity ==
        Parameters::Verbosity::verbose &&
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
        simulation_parameters.boundary_conditions.id,
        independent_column_names,
        this->forces_on_boundaries,
        dependent_column_names,
        this->simulation_parameters.simulation_control.log_precision);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Force  summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (unsigned int i_boundary = 0;
       i_boundary < simulation_parameters.boundary_conditions.size;
       ++i_boundary)
    {
      if (simulation_control->is_steady())
        {
          this->forces_tables[i_boundary].add_value(
            "cells", this->triangulation->n_global_active_cells());
        }
      else
        {
          this->forces_tables[i_boundary].add_value(
            "time", simulation_control->get_current_time());
          this->forces_tables[i_boundary].set_precision(
            "time", simulation_parameters.forces_parameters.output_precision);
        }

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
        "f_x", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "f_y", simulation_parameters.forces_parameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "f_z", simulation_parameters.forces_parameters.output_precision);
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_torques(
  const VectorType &evaluation_point)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_torques");

  this->torques_on_boundaries =
    calculate_torques(this->dof_handler,
                      evaluation_point,
                      simulation_parameters.physical_properties,
                      simulation_parameters.boundary_conditions,
                      *this->face_quadrature,
                      *this->mapping);

  if (simulation_parameters.forces_parameters.verbosity ==
        Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      this->pcout << std::endl;
      std::string independent_column_names = "Boundary ID";

      std::vector<std::string> dependent_column_names;
      dependent_column_names.push_back("T_x");
      dependent_column_names.push_back("T_y");
      dependent_column_names.push_back("T_z");

      TableHandler table = make_table_scalars_tensors(
        simulation_parameters.boundary_conditions.id,
        independent_column_names,
        this->torques_on_boundaries,
        dependent_column_names,
        this->simulation_parameters.simulation_control.log_precision);

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Torque summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (unsigned int boundary_id = 0;
       boundary_id < simulation_parameters.boundary_conditions.size;
       ++boundary_id)
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
        "T_x", this->torques_on_boundaries[boundary_id][0]);
      this->torques_tables[boundary_id].add_value(
        "T_y", this->torques_on_boundaries[boundary_id][1]);
      this->torques_tables[boundary_id].add_value(
        "T_z", this->torques_on_boundaries[boundary_id][2]);

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
            simulation_parameters.analytical_solution->get_filename() + ".dat";
          std::ofstream output(filename.c_str());
          error_table.write_text(output);
          std::vector<std::string> sub_columns;
          if (simulation_parameters.simulation_control.method ==
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              sub_columns.push_back("cells");
              sub_columns.push_back("error_velocity");
              sub_columns.push_back("error_pressure");
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
NavierStokesBase<dim, VectorType, DofsType>::finish_time_step_fd()
{
  if (simulation_parameters.simulation_control.method !=
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      percolate_time_vectors_fd();
      const double CFL = calculate_CFL(this->dof_handler,
                                       this->present_solution,
                                       simulation_control->get_time_step(),
                                       *this->cell_quadrature,
                                       *this->mapping);
      this->simulation_control->set_CFL(CFL);
    }
  if (this->simulation_parameters.restart_parameters.checkpoint &&
      simulation_control->get_step_number() != 0 &&
      simulation_control->get_step_number() %
          this->simulation_parameters.restart_parameters.frequency ==
        0)
    {
      this->write_checkpoint();
    }

  if (this->simulation_parameters.timer.type ==
      Parameters::Timer::Type::iteration)
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
  auto &present_solution = this->present_solution;
  if (simulation_control->get_assembly_method() ==
      Parameters::SimulationControl::TimeSteppingMethod::sdirk22)
    {
      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk22_1);
      PhysicsSolver<VectorType>::solve_non_linear_system(false);
      this->solution_stages[0] = present_solution;

      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk22_2);
      PhysicsSolver<VectorType>::solve_non_linear_system(false);
    }
  else if (simulation_control->get_assembly_method() ==
           Parameters::SimulationControl::TimeSteppingMethod::sdirk33)
    {
      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_1);
      PhysicsSolver<VectorType>::solve_non_linear_system(false);

      this->solution_stages[0] = present_solution;

      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_2);
      PhysicsSolver<VectorType>::solve_non_linear_system(false);

      this->solution_stages[1] = present_solution;

      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33_3);
      PhysicsSolver<VectorType>::solve_non_linear_system(false);
    }
  else
    {
      // sdirk schemes are not implemented for multiphysics simulations

      // Solve the auxiliary physics that should be treated BEFORE the fluid
      // dynamics
      multiphysics->pre_solve(simulation_parameters.simulation_control.method);

      PhysicsSolver<VectorType>::solve_non_linear_system(false);

      // Solve the auxiliary physics that should be treated AFTER the fluid
      // dynamics
      multiphysics->post_solve(simulation_parameters.simulation_control.method);
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh()
{
  if (simulation_control->get_step_number() %
        this->simulation_parameters.mesh_adaptation.frequency ==
      0)
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
NavierStokesBase<dim, VectorType, DofsType>::box_refine_mesh()
{
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
          // A user parameter option could be made to chose a partitionner.
          GridTools::partition_triangulation(0, basetria);


          auto construction_data = TriangulationDescription::Utilities::
            create_description_from_triangulation(basetria, 0);

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
              0,
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
      TimerOutput::Scope t(this->computing_timer, "box refine");

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
      parallel::distributed::SolutionTransfer<dim, VectorType>
        solution_transfer(this->dof_handler);
      std::vector<parallel::distributed::SolutionTransfer<dim, VectorType>>
        previous_solutions_transfer;
      // Important to reserve to prevent pointer dangling
      previous_solutions_transfer.reserve(previous_solutions.size());
      for (unsigned int i = 0; i < previous_solutions.size(); ++i)
        {
          previous_solutions_transfer.push_back(
            parallel::distributed::SolutionTransfer<dim, VectorType>(
              this->dof_handler));
          previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
            previous_solutions[i]);
        }

      parallel::distributed::SolutionTransfer<dim, VectorType>
        solution_transfer_m1(this->dof_handler);
      parallel::distributed::SolutionTransfer<dim, VectorType>
        solution_transfer_m2(this->dof_handler);
      parallel::distributed::SolutionTransfer<dim, VectorType>
        solution_transfer_m3(this->dof_handler);
      solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

      multiphysics->prepare_for_mesh_adaptation();
      if (this->simulation_parameters.post_processing
            .calculate_average_velocities)
        average_velocities->prepare_for_mesh_adaptation();


      tria.execute_coarsening_and_refinement();
      this->setup_dofs();

      // Set up the vectors for the transfer
      VectorType tmp(locally_owned_dofs, this->mpi_communicator);

      // Interpolate the solution at time and previous time
      solution_transfer.interpolate(tmp);

      // Distribute constraints
      auto &nonzero_constraints = this->nonzero_constraints;
      nonzero_constraints.distribute(tmp);

      // Fix on the new mesh
      present_solution = tmp;

      for (unsigned int i = 0; i < previous_solutions.size(); ++i)
        {
          VectorType tmp_previous_solution(locally_owned_dofs,
                                           this->mpi_communicator);
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
  TimerOutput::Scope t(this->computing_timer, "refine");

  Vector<float> estimated_error_per_cell(tria.n_active_cells());
  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);
  auto &                           present_solution = this->present_solution;

  if (this->simulation_parameters.mesh_adaptation.variable ==
      Parameters::MeshAdaptation::Variable::pressure)
    {
      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(pressure));
    }
  else if (this->simulation_parameters.mesh_adaptation.variable ==
           Parameters::MeshAdaptation::Variable::velocity)
    {
      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        this->dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        this->fe->component_mask(velocity));
    }
  else
    {
      // refine_mesh on an auxiliary physic parameter
      multiphysics->compute_kelly(estimated_error_per_cell);
    }

  if (this->simulation_parameters.mesh_adaptation.fractionType ==
      Parameters::MeshAdaptation::FractionType::number)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      tria,
      estimated_error_per_cell,
      this->simulation_parameters.mesh_adaptation.refinement_fraction,
      this->simulation_parameters.mesh_adaptation.coarsening_fraction,
      this->simulation_parameters.mesh_adaptation.maximum_number_elements);

  else if (this->simulation_parameters.mesh_adaptation.fractionType ==
           Parameters::MeshAdaptation::FractionType::fraction)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      tria,
      estimated_error_per_cell,
      this->simulation_parameters.mesh_adaptation.refinement_fraction,
      this->simulation_parameters.mesh_adaptation.coarsening_fraction);

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
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
    this->dof_handler);
  std::vector<parallel::distributed::SolutionTransfer<dim, VectorType>>
    previous_solutions_transfer;
  // Important to reserve to prevent pointer dangling
  previous_solutions_transfer.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer.push_back(
        parallel::distributed::SolutionTransfer<dim, VectorType>(
          this->dof_handler));
      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }

  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);

  multiphysics->prepare_for_mesh_adaptation();
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    average_velocities->prepare_for_mesh_adaptation();

  tria.execute_coarsening_and_refinement();
  this->setup_dofs();

  // Set up the vectors for the transfer
  VectorType tmp(locally_owned_dofs, this->mpi_communicator);

  // Interpolate the solution at time and previous time
  solution_transfer.interpolate(tmp);

  // Distribute constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(tmp);

  // Fix on the new mesh
  present_solution = tmp;

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      VectorType tmp_previous_solution(locally_owned_dofs,
                                       this->mpi_communicator);
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }

  multiphysics->post_mesh_adaptation();
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    average_velocities->post_mesh_adaptation();
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh_uniform()
{
  TimerOutput::Scope t(this->computing_timer, "refine");

  // Solution transfer objects for all the solutions
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m2(
    this->dof_handler);
  parallel::distributed::SolutionTransfer<dim, VectorType> solution_transfer_m3(
    this->dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);

  std::vector<parallel::distributed::SolutionTransfer<dim, VectorType>>
    previous_solutions_transfer;
  // Important to reserve to prevent pointer dangling
  previous_solutions_transfer.reserve(previous_solutions.size());

  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions_transfer.emplace_back(
        parallel::distributed::SolutionTransfer<dim, VectorType>(
          this->dof_handler));
      previous_solutions_transfer[i].prepare_for_coarsening_and_refinement(
        previous_solutions[i]);
    }

  multiphysics->prepare_for_mesh_adaptation();

  // Refine
  this->triangulation->refine_global(1);

  setup_dofs();

  // Set up the vectors for the transfer
  VectorType tmp(locally_owned_dofs, this->mpi_communicator);

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
      VectorType tmp_previous_solution(locally_owned_dofs,
                                       this->mpi_communicator);
      previous_solutions_transfer[i].interpolate(tmp_previous_solution);
      nonzero_constraints.distribute(tmp_previous_solution);
      previous_solutions[i] = tmp_previous_solution;
    }

  multiphysics->post_mesh_adaptation();
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocess_fd(bool firstIter)
{
  auto &present_solution = this->present_solution;

  if (this->simulation_parameters.post_processing.calculate_enstrophy)
    {
      double enstrophy = calculate_enstrophy(this->dof_handler,
                                             present_solution,
                                             *this->cell_quadrature,
                                             *this->mapping);

      this->enstrophy_table.add_value("time",
                                      simulation_control->get_current_time());
      this->enstrophy_table.add_value("enstrophy", enstrophy);

      // Display Enstrophy to screen if verbosity is enabled
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Enstrophy  : " << enstrophy << std::endl;
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

  // The average velocities and reynolds stresses are calculated when the
  // time reaches the initial time. (time >= initial time) with 1e-6 as
  // tolerance.
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      // Calculate average velocities when the time reaches the initial time.
      // time >= initial time with the epsilon as tolerance.
      const double dt = simulation_control->get_time_step();
      if (simulation_control->get_current_time() >
          (simulation_parameters.post_processing.initial_time - 1e-6 * dt))
        {
          this->average_velocities->calculate_average_velocities(
            this->local_evaluation_point,
            simulation_parameters.post_processing,
            simulation_control->get_current_time(),
            dt);
        }
    }

  if (this->simulation_parameters.post_processing.calculate_kinetic_energy)
    {
      TimerOutput::Scope t(this->computing_timer, "kinetic_energy_calculation");
      double             kE = calculate_kinetic_energy(this->dof_handler,
                                           present_solution,
                                           *this->cell_quadrature,
                                           *this->mapping);
      this->kinetic_energy_table.add_value(
        "time", simulation_control->get_current_time());
      this->kinetic_energy_table.add_value("kinetic-energy", kE);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Kinetic energy : " << kE << std::endl;
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

  // Calculte apparent viscosity
  if (this->simulation_parameters.post_processing.calculate_apparent_viscosity)
    {
      TimerOutput::Scope t(this->computing_timer,
                           "apparent_viscosity_calculation");
      double             apparent_viscosity = calculate_apparent_viscosity(
        this->dof_handler,
        this->present_solution,
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.physical_properties);

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
      TimerOutput::Scope t(this->computing_timer, "pressure_drop_calculation");
      double             pressure_drop = calculate_pressure_drop(
        this->dof_handler,
        this->mapping,
        this->evaluation_point,
        *this->cell_quadrature,
        *this->face_quadrature,
        this->simulation_parameters.post_processing.inlet_boundary_id,
        this->simulation_parameters.post_processing.outlet_boundary_id);
      this->pressure_drop_table.add_value(
        "time", simulation_control->get_current_time());
      this->pressure_drop_table.add_value("pressure-drop", pressure_drop);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Pressure drop: "
                      << this->simulation_parameters.physical_properties
                             .fluids[0]
                             .density *
                           pressure_drop
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
          this->pressure_drop_table.write_text(output);
        }
    }


  // Calculate inlet flow rate and area
  if (this->simulation_parameters.flow_control.enable_flow_control)
    {
      this->postprocessing_flow_rate(this->present_solution);
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
          // Update the time of the exact solution to the actual time
          this->exact_solution->set_time(
            simulation_control->get_current_time());
          const std::pair<double, double> errors =
            calculate_L2_error(dof_handler,
                               present_solution,
                               exact_solution,
                               *this->cell_quadrature,
                               *this->mapping);
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

              // Calculate error on pressure for VOF simulations
              if (this->simulation_parameters.multiphysics.VOF)
                this->error_table.add_value("error_pressure", error_pressure);
            }
          if (this->simulation_parameters.analytical_solution->verbosity ==
              Parameters::Verbosity::verbose)
            {
              this->pcout << "L2 error velocity : "
                          << std::setprecision(
                               simulation_control->get_log_precision())
                          << error_velocity << std::endl;
            }
        }
    }
  if (this->simulation_control->is_output_iteration())
    {
      this->write_output_results(present_solution);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::set_nodal_values()
{
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  VectorTools::interpolate(*this->mapping,
                           this->dof_handler,
                           this->simulation_parameters.initial_condition->uvwp,
                           this->newton_update,
                           this->fe->component_mask(velocities));
  VectorTools::interpolate(*this->mapping,
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
            *this->mapping,
            this->dof_handler,
            this->simulation_parameters.initial_condition->uvwp,
            this->newton_update,
            this->fe->component_mask(velocities));
          VectorTools::interpolate(
            *this->mapping,
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
NavierStokesBase<dim, VectorType, DofsType>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "read_checkpoint");
  std::string        prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
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
  std::vector<VectorType *> x_system(1 + previous_solutions.size());

  VectorType distributed_system(locally_owned_dofs, this->mpi_communicator);
  x_system[0] = &(distributed_system);

  std::vector<VectorType> distributed_previous_solutions;
  distributed_previous_solutions.reserve(previous_solutions.size());
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        VectorType(locally_owned_dofs, this->mpi_communicator));
      x_system[i + 1] = &distributed_previous_solutions[i];
    }
  parallel::distributed::SolutionTransfer<dim, VectorType> system_trans_vectors(
    this->dof_handler);

  if (simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<VectorType *> sum_vectors =
        this->average_velocities->read(prefix);

      x_system.insert(x_system.end(), sum_vectors.begin(), sum_vectors.end());
    }

  system_trans_vectors.deserialize(x_system);
  this->present_solution = distributed_system;
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      previous_solutions[i] = distributed_previous_solutions[i];
    }

  if (simulation_parameters.flow_control.enable_flow_control)
    {
      this->flow_control.read(prefix);

      this->flow_rate =
        calculate_flow_rate(this->dof_handler,
                            present_solution,
                            simulation_parameters.flow_control.boundary_flow_id,
                            *this->face_quadrature,
                            *this->mapping);
    }


  multiphysics->read_checkpoint();
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
  TimerOutput::Scope t(this->computing_timer, "output");

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

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
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
        this->average_velocities->get_average_velocities(),
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


  // Create additional post-processor that derives information from the
  // solution
  VorticityPostprocessor<dim> vorticity;
  data_out.add_data_vector(solution, vorticity);

  QCriterionPostprocessor<dim> qcriterion;
  data_out.add_data_vector(solution, qcriterion);

  SRFPostprocessor<dim> srf(simulation_parameters.velocity_sources.omega_x,
                            simulation_parameters.velocity_sources.omega_y,
                            simulation_parameters.velocity_sources.omega_z);

  if (simulation_parameters.velocity_sources.type ==
      Parameters::VelocitySource::VelocitySourceType::srf)
    data_out.add_data_vector(solution, srf);

  NonNewtonianViscosityPostprocessor<dim> non_newtonian_viscosity(
    simulation_parameters.physical_properties);
  if (simulation_parameters.physical_properties.non_newtonian_flow)
    data_out.add_data_vector(solution, non_newtonian_viscosity);


  output_field_hook(data_out);

  multiphysics->attach_solution_to_output(data_out);

  // Build the patches and write the output

  data_out.build_patches(*this->mapping,
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
      data_out_faces.build_patches(*this->mapping);

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
       boundary_id < simulation_parameters.boundary_conditions.size;
       ++boundary_id)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.forces_parameters.force_output_name + "." +
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
       boundary_id < simulation_parameters.boundary_conditions.size;
       ++boundary_id)
    {
      std::string filename =
        simulation_parameters.simulation_control.output_folder +
        simulation_parameters.forces_parameters.torque_output_name + "." +
        Utilities::int_to_string(boundary_id, 2) + ".dat";
      std::ofstream output(filename.c_str());

      this->torques_tables[boundary_id].write_text(output);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "write_checkpoint");
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
  sol_set_transfer.push_back(&this->present_solution);
  for (unsigned int i = 0; i < previous_solutions.size(); ++i)
    {
      sol_set_transfer.push_back(&previous_solutions[i]);
    }

  if (simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<const VectorType *> av_set_transfer =
        this->average_velocities->save(prefix);

      // Insert average velocities vectors into the set transfer vector
      sol_set_transfer.insert(sol_set_transfer.end(),
                              av_set_transfer.begin(),
                              av_set_transfer.end());
    }



  parallel::distributed::SolutionTransfer<dim, VectorType> system_trans_vectors(
    this->dof_handler);
  system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  multiphysics->write_checkpoint();

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
