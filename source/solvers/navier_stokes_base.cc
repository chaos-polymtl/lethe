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

#include <solvers/navier_stokes_base.h>
#include <solvers/postprocessing_enstrophy.h>
#include <solvers/postprocessing_force.h>
#include <solvers/postprocessing_kinetic_energy.h>
#include <solvers/postprocessing_torque.h>


/*
 * Constructor for the Navier-Stokes base class
 */
template <int dim, typename VectorType, typename DofsType>
NavierStokesBase<dim, VectorType, DofsType>::NavierStokesBase(
  NavierStokesSolverParameters<dim> &p_nsparam,
  const unsigned int                 p_degreeVelocity,
  const unsigned int                 p_degreePressure)
  : PhysicsSolver<VectorType>(p_nsparam.nonLinearSolver)
  //      new NewtonNonLinearSolver<VectorType>(this,
  //      p_nsparam.nonLinearSolver))
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
  , fe(FE_Q<dim>(p_degreeVelocity), dim, FE_Q<dim>(p_degreePressure), 1)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , nsparam(p_nsparam)
  , degreeVelocity_(p_degreeVelocity)
  , degreePressure_(p_degreePressure)
  , degreeQuadrature_(p_degreeVelocity + 1)
{
  this->pcout.set_condition(
    Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0);
  this->simulationControl = nsparam.simulationControl;

  // Overide default value of quadrature point if they are specified
  if (nsparam.femParameters.quadraturePoints > 0)
    degreeQuadrature_ = nsparam.femParameters.quadraturePoints;

  // Change the behavior of the timer for situations when you don't want outputs
  if (nsparam.timer.type == Parameters::Timer::Type::none)
    this->computing_timer.disable_output();

  // Pre-allocate the force tables to match the number of boundary conditions
  forces_.resize(nsparam.boundaryConditions.size);
  torques_.resize(nsparam.boundaryConditions.size);
  forces_tables.resize(nsparam.boundaryConditions.size);
  torques_tables.resize(nsparam.boundaryConditions.size);

  // Get the exact solution from the parser
  exact_solution = &nsparam.analyticalSolution->velocity;

  // If there is a forcing function, get it from the parser
  if (nsparam.sourceTerm->source_term())
    forcing_function = &nsparam.sourceTerm->source;
  else
    forcing_function = new NoForce<dim>;

  this->pcout << "Running on "
              << Utilities::MPI::n_mpi_processes(this->mpi_communicator)
              << " MPI rank(s)..." << std::endl;
}

template <int dim, typename VectorType, typename DofsType>
double
NavierStokesBase<dim, VectorType, DofsType>::calculate_CFL(
  const VectorType &evaluation_point)
{
  QGauss<dim>                          quadrature_formula(1);
  const MappingQ<dim>                  mapping(this->degreeVelocity_,
                              nsparam.femParameters.qmapping_all);
  FEValues<dim>                        fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points);
  const unsigned int                   dofs_per_cell = this->fe.dofs_per_cell;
  const unsigned int                   n_q_points = quadrature_formula.size();
  std::vector<Vector<double>>          initial_velocity(n_q_points,
                                                        Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);

  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);

  // Element size
  double h;

  // Current time step
  double timeStep = this->simulationControl.getCurrentTimeStep();

  // CFL
  double CFL = 0;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          if (dim == 2)
            h = std::sqrt(4. * cell->measure() / M_PI) / this->degreeVelocity_;
          else if (dim == 3)
            h =
              pow(6 * cell->measure() / M_PI, 1. / 3.) / this->degreeVelocity_;
          fe_values.reinit(cell);
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double localCFL =
                present_velocity_values[q].norm() / h * timeStep;
              CFL = std::max(CFL, localCFL);
            }
        }
    }
  CFL = Utilities::MPI::max(CFL, this->mpi_communicator);

  return CFL;
}

// This is a primitive first implementation that could be greatly improved by
// doing a single pass instead of N boundary passes
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_forces(
  const VectorType &       evaluation_point,
  const SimulationControl &simulationControl)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_forces");

  this->forces_ = calculate_forces(this->fe,
                                   this->dof_handler,
                                   evaluation_point,
                                   nsparam.physicalProperties,
                                   nsparam.femParameters,
                                   nsparam.boundaryConditions,
                                   mpi_communicator);

  if (nsparam.forcesParameters.verbosity == Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      std::cout << std::endl;
      TableHandler table;

      for (unsigned int i_boundary = 0;
           i_boundary < nsparam.boundaryConditions.size;
           ++i_boundary)
        {
          table.add_value("Boundary ID", i_boundary);
          table.add_value("f_x", this->forces_[i_boundary][0]);
          table.add_value("f_y", this->forces_[i_boundary][1]);
          table.set_precision("f_x",
                              nsparam.forcesParameters.display_precision);
          table.set_precision("f_y",
                              nsparam.forcesParameters.display_precision);
          if (dim == 3)
            {
              table.add_value("f_z", this->forces_[i_boundary][2]);
              table.set_precision("f_z",
                                  nsparam.forcesParameters.display_precision);
            }
        }
      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Force  summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (unsigned int i_boundary = 0;
       i_boundary < nsparam.boundaryConditions.size;
       ++i_boundary)
    {
      this->forces_tables[i_boundary].add_value("time",
                                                simulationControl.getTime());
      this->forces_tables[i_boundary].add_value("f_x",
                                                this->forces_[i_boundary][0]);
      this->forces_tables[i_boundary].add_value("f_y",
                                                this->forces_[i_boundary][1]);
      if (dim == 3)
        this->forces_tables[i_boundary].add_value("f_z",
                                                  this->forces_[i_boundary][2]);
      else
        this->forces_tables[i_boundary].add_value("f_z", 0.);

      // Precision
      this->forces_tables[i_boundary].set_precision(
        "f_x", nsparam.forcesParameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "f_y", nsparam.forcesParameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "f_z", nsparam.forcesParameters.output_precision);
      this->forces_tables[i_boundary].set_precision(
        "time", nsparam.forcesParameters.output_precision);
    }
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocessing_torques(
  const VectorType &       evaluation_point,
  const SimulationControl &simulationControl)
{
  TimerOutput::Scope t(this->computing_timer, "calculate_torques");

  this->torques_ = calculate_torques(this->fe,
                                     this->dof_handler,
                                     evaluation_point,
                                     nsparam.physicalProperties,
                                     nsparam.femParameters,
                                     nsparam.boundaryConditions,
                                     mpi_communicator);

  if (nsparam.forcesParameters.verbosity == Parameters::Verbosity::verbose &&
      this->this_mpi_process == 0)
    {
      this->pcout << std::endl;
      TableHandler table;

      for (unsigned int boundary_id = 0;
           boundary_id < nsparam.boundaryConditions.size;
           ++boundary_id)
        {
          table.add_value("Boundary ID", boundary_id);
          table.add_value("T_x", this->torques_[boundary_id][0]);
          table.add_value("T_y", this->torques_[boundary_id][1]);
          table.set_precision("T_x",
                              nsparam.forcesParameters.display_precision);
          table.set_precision("T_y",
                              nsparam.forcesParameters.display_precision);
          table.add_value("T_z", this->torques_[boundary_id][2]);
          table.set_precision("T_z",
                              nsparam.forcesParameters.display_precision);
        }

      std::cout << "+------------------------------------------+" << std::endl;
      std::cout << "|  Torque summary                          |" << std::endl;
      std::cout << "+------------------------------------------+" << std::endl;
      table.write_text(std::cout);
    }

  for (unsigned int boundary_id = 0;
       boundary_id < nsparam.boundaryConditions.size;
       ++boundary_id)
    {
      this->torques_tables[boundary_id].add_value("time",
                                                  simulationControl.getTime());
      this->torques_tables[boundary_id].add_value(
        "T_x", this->torques_[boundary_id][0]);
      this->torques_tables[boundary_id].add_value(
        "T_y", this->torques_[boundary_id][1]);
      this->torques_tables[boundary_id].add_value(
        "T_z", this->torques_[boundary_id][2]);

      // Precision
      this->torques_tables[boundary_id].set_precision(
        "T_x", nsparam.forcesParameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "T_y", nsparam.forcesParameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "T_z", nsparam.forcesParameters.output_precision);
      this->torques_tables[boundary_id].set_precision(
        "time", nsparam.forcesParameters.output_precision);
    }
}

// Find the l2 norm of the error between the finite element sol'n and the exact
// sol'n
template <int dim, typename VectorType, typename DofsType>
std::pair<double, double>
NavierStokesBase<dim, VectorType, DofsType>::calculate_L2_error(
  const VectorType &evaluation_point)
{
  TimerOutput::Scope t(this->computing_timer, "error");

  QGauss<dim>         quadrature_formula(this->degreeQuadrature_ + 1);
  const MappingQ<dim> mapping(this->degreeVelocity_,
                              nsparam.femParameters.qmapping_all);
  FEValues<dim>       fe_values(mapping,
                          this->fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);

  const unsigned int dofs_per_cell =
    this->fe.dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim + 1));

  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  std::vector<double>         local_pressure_values(n_q_points);

  Function<dim> *l_exact_solution = this->exact_solution;

  double pressure_integral       = 0;
  double exact_pressure_integral = 0;

  // loop over elements to calculate average pressure
  {
    typename DoFHandler<dim>::active_cell_iterator cell = this->dof_handler
                                                            .begin_active(),
                                                   endc =
                                                     this->dof_handler.end();
    for (; cell != endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            fe_values[pressure].get_function_values(evaluation_point,
                                                    local_pressure_values);
            // Get the exact solution at all gauss points
            l_exact_solution->vector_value_list(
              fe_values.get_quadrature_points(), q_exactSol);


            // Retrieve the effective "connectivity matrix" for this element
            cell->get_dof_indices(local_dof_indices);

            for (unsigned int q = 0; q < n_q_points; q++)
              {
                pressure_integral +=
                  local_pressure_values[q] * fe_values.JxW(q);
                exact_pressure_integral +=
                  q_exactSol[q][dim] * fe_values.JxW(q);
              }
          }
      }
  }

  pressure_integral =
    Utilities::MPI::sum(pressure_integral, this->mpi_communicator);
  exact_pressure_integral =
    Utilities::MPI::sum(exact_pressure_integral, this->mpi_communicator);
  double average_pressure       = pressure_integral / globalVolume_;
  double average_exact_pressure = exact_pressure_integral / globalVolume_;


  double l2errorU = 0.;
  double l2errorP = 0.;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[velocities].get_function_values(evaluation_point,
                                                    local_velocity_values);
          fe_values[pressure].get_function_values(evaluation_point,
                                                  local_pressure_values);

          // Retrieve the effective "connectivity matrix" for this element
          cell->get_dof_indices(local_dof_indices);

          // Get the exact solution at all gauss points
          l_exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                              q_exactSol);

          for (unsigned int q = 0; q < n_q_points; q++)
            {
              // Find the values of x and u_h (the finite element solution) at
              // the quadrature points
              double ux_sim   = local_velocity_values[q][0];
              double ux_exact = q_exactSol[q][0];

              double uy_sim   = local_velocity_values[q][1];
              double uy_exact = q_exactSol[q][1];

              l2errorU +=
                (ux_sim - ux_exact) * (ux_sim - ux_exact) * fe_values.JxW(q);
              l2errorU +=
                (uy_sim - uy_exact) * (uy_sim - uy_exact) * fe_values.JxW(q);

              if (dim == 3)
                {
                  double uz_sim   = local_velocity_values[q][2];
                  double uz_exact = q_exactSol[q][2];
                  l2errorU += (uz_sim - uz_exact) * (uz_sim - uz_exact) *
                              fe_values.JxW(q);
                }

              double p_sim   = local_pressure_values[q] - average_pressure;
              double p_exact = q_exactSol[q][dim] - average_exact_pressure;
              l2errorP +=
                (p_sim - p_exact) * (p_sim - p_exact) * fe_values.JxW(q);
            }
        }
    }
  l2errorU = Utilities::MPI::sum(l2errorU, this->mpi_communicator);
  l2errorP = Utilities::MPI::sum(l2errorP, this->mpi_communicator);

  return std::make_pair(std::sqrt(l2errorU), std::sqrt(l2errorP));
}

/*
 * Attaches manifold to the boundaries of the mesh
 */
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::create_manifolds()
{
  Parameters::Manifolds manifolds = this->nsparam.manifoldsParameters;

  for (unsigned int i = 0; i < manifolds.types.size(); ++i)
    {
      if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::spherical)
        {
          Point<dim> circleCenter;
          circleCenter = Point<dim>(manifolds.arg1[i], manifolds.arg2[i]);
          static const SphericalManifold<dim> manifold_description(
            circleCenter);
          this->triangulation->set_manifold(manifolds.id[i],
                                            manifold_description);
          this->triangulation->set_all_manifold_ids_on_boundary(
            manifolds.id[i], manifolds.id[i]);
        }
      else if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::iges)
        {
          attach_cad_to_manifold(triangulation,
                                 manifolds.cad_files[i],
                                 manifolds.id[i]);
        }
      else if (manifolds.types[i] == Parameters::Manifolds::ManifoldType::none)
        {}
      else
        throw std::runtime_error("Unsupported manifolds type");
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::finish_simulation()
{
  if (nsparam.forcesParameters.calculate_force)
    {
      this->write_output_forces();
    }

  if (nsparam.forcesParameters.calculate_torque)
    {
      this->write_output_torques();
    }
  if (nsparam.analyticalSolution->calculate_error())
    {
      if (simulationControl.getMethod() ==
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
            nsparam.analyticalSolution->get_filename() + ".dat";
          std::ofstream output(filename.c_str());
          error_table.write_text(output);
          std::vector<std::string> sub_columns;
          if (simulationControl.getMethod() ==
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
  if (this->simulationControl.getMethod() !=
      Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      this->solution_m3 = this->solution_m2;
      this->solution_m2 = this->solution_m1;
      this->solution_m1 = this->present_solution;
      const double CFL  = this->calculate_CFL(this->present_solution);
      this->simulationControl.setCFL(CFL);
    }
  if (this->nsparam.restartParameters.checkpoint &&
      this->simulationControl.getIter() %
          this->nsparam.restartParameters.frequency ==
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

// Do an iteration with the GLS NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial condition
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::iterate(bool firstIteration)
{
  // Carry out the integration normally if this is not a method that needs to be
  // started
  if (!firstIteration || !is_bdf(this->simulationControl.getMethod()))
    {
      if (this->simulationControl.getMethod() ==
          Parameters::SimulationControl::TimeSteppingMethod::sdirk2)
        {
          PhysicsSolver<VectorType>::solve_non_linear_system(
            Parameters::SimulationControl::TimeSteppingMethod::sdirk2_1,
            false,
            false);
          this->solution_m2 = this->present_solution;

          PhysicsSolver<VectorType>::solve_non_linear_system(
            Parameters::SimulationControl::TimeSteppingMethod::sdirk2_2,
            false,
            false);
        }

      else if (this->simulationControl.getMethod() ==
               Parameters::SimulationControl::TimeSteppingMethod::sdirk3)
        {
          PhysicsSolver<VectorType>::solve_non_linear_system(
            Parameters::SimulationControl::TimeSteppingMethod::sdirk3_1,
            false,
            false);

          this->solution_m2 = this->present_solution;

          PhysicsSolver<VectorType>::solve_non_linear_system(
            Parameters::SimulationControl::TimeSteppingMethod::sdirk3_2,
            false,
            false);

          this->solution_m3 = this->present_solution;

          PhysicsSolver<VectorType>::solve_non_linear_system(
            Parameters::SimulationControl::TimeSteppingMethod::sdirk3_3,
            false,
            false);
        }
      else
        {
          PhysicsSolver<VectorType>::solve_non_linear_system(
            this->simulationControl.getMethod(), false, false);
        }
    }
  // This is the first iteration and we are using a BDF scheme
  else
    {
      if (this->simulationControl.getMethod() ==
          Parameters::SimulationControl::TimeSteppingMethod::bdf1)
        PhysicsSolver<VectorType>::solve_non_linear_system(
          this->simulationControl.getMethod(), false, true);

      else if (this->simulationControl.getMethod() ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        {
          Parameters::SimulationControl timeParameters =
            this->simulationControl.getParameters();

          // Start the BDF2 with a single Euler time step with a lower time step
          this->simulationControl.setTimeStep(
            timeParameters.dt * timeParameters.startup_timestep_scaling);
          this->simulationControl.setMethod(
            Parameters::SimulationControl::TimeSteppingMethod::bdf1);
          PhysicsSolver<VectorType>::solve_non_linear_system(
            this->simulationControl.getMethod(), false, true);
          this->solution_m2 = this->solution_m1;
          this->solution_m1 = this->present_solution;

          // Reset the time step and do a bdf 2 newton iteration using the two
          // steps to complete the full step
          this->simulationControl.setMethod(
            Parameters::SimulationControl::TimeSteppingMethod::bdf2);
          this->simulationControl.setTimeStep(
            timeParameters.dt * (1. - timeParameters.startup_timestep_scaling));
          PhysicsSolver<VectorType>::solve_non_linear_system(
            this->simulationControl.getMethod(), false, true);
        }

      else if (this->simulationControl.getMethod() ==
               Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        {
          Parameters::SimulationControl timeParameters =
            this->simulationControl.getParameters();

          // Start the BDF2 with a single Euler time step with a lower time step
          this->simulationControl.setTimeStep(
            timeParameters.dt * timeParameters.startup_timestep_scaling);
          this->simulationControl.setMethod(
            Parameters::SimulationControl::TimeSteppingMethod::bdf1);
          PhysicsSolver<VectorType>::solve_non_linear_system(
            this->simulationControl.getMethod(), false, true);
          this->solution_m2 = this->solution_m1;
          this->solution_m1 = this->present_solution;

          // Reset the time step and do a bdf 2 newton iteration using the two
          // steps
          this->simulationControl.setMethod(
            Parameters::SimulationControl::TimeSteppingMethod::bdf1);
          this->simulationControl.setTimeStep(
            timeParameters.dt * timeParameters.startup_timestep_scaling);
          PhysicsSolver<VectorType>::solve_non_linear_system(
            this->simulationControl.getMethod(), false, true);
          this->solution_m3 = this->solution_m2;
          this->solution_m2 = this->solution_m1;
          this->solution_m1 = this->present_solution;

          // Reset the time step and do a bdf 3 newton iteration using the two
          // steps to complete the full step
          this->simulationControl.setMethod(
            Parameters::SimulationControl::TimeSteppingMethod::bdf3);
          this->simulationControl.setTimeStep(
            timeParameters.dt *
            (1. - 2. * timeParameters.startup_timestep_scaling));
          PhysicsSolver<VectorType>::solve_non_linear_system(
            this->simulationControl.getMethod(), false, true);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::refine_mesh()
{
  if (this->simulationControl.getIter() %
        this->nsparam.meshAdaptation.frequency ==
      0)
    {
      if (this->nsparam.meshAdaptation.type ==
          Parameters::MeshAdaptation::Type::kelly)
        {
          refine_mesh_kelly();
        }
      else if (this->nsparam.meshAdaptation.type ==
               Parameters::MeshAdaptation::Type::uniform)
        {
          refine_mesh_uniform();
        }
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
  const MappingQ<dim> mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);
  const FEValuesExtractors::Vector velocity(0);
  const FEValuesExtractors::Scalar pressure(dim);
  if (this->nsparam.meshAdaptation.variable ==
      Parameters::MeshAdaptation::Variable::pressure)
    {
      KellyErrorEstimator<dim>::estimate(
        mapping,
        this->dof_handler,
        QGauss<dim - 1>(this->degreeQuadrature_ + 1),
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe.component_mask(pressure));
    }
  else if (this->nsparam.meshAdaptation.variable ==
           Parameters::MeshAdaptation::Variable::velocity)
    {
      KellyErrorEstimator<dim>::estimate(
        mapping,
        this->dof_handler,
        QGauss<dim - 1>(this->degreeQuadrature_ + 1),
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        this->present_solution,
        estimated_error_per_cell,
        this->fe.component_mask(velocity));
    }

  if (this->nsparam.meshAdaptation.fractionType ==
      Parameters::MeshAdaptation::FractionType::number)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
      tria,
      estimated_error_per_cell,
      this->nsparam.meshAdaptation.fractionRefinement,
      this->nsparam.meshAdaptation.fractionCoarsening,
      this->nsparam.meshAdaptation.maxNbElements);

  else if (this->nsparam.meshAdaptation.fractionType ==
           Parameters::MeshAdaptation::FractionType::fraction)
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      tria,
      estimated_error_per_cell,
      this->nsparam.meshAdaptation.fractionRefinement,
      this->nsparam.meshAdaptation.fractionCoarsening);

  if (tria.n_levels() > this->nsparam.meshAdaptation.maxRefLevel)
    for (typename Triangulation<dim>::active_cell_iterator cell =
           tria.begin_active(this->nsparam.meshAdaptation.maxRefLevel);
         cell != tria.end();
         ++cell)
      cell->clear_refine_flag();
  for (typename Triangulation<dim>::active_cell_iterator cell =
         tria.begin_active(this->nsparam.meshAdaptation.minRefLevel);
       cell != tria.end_active(this->nsparam.meshAdaptation.minRefLevel);
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
  solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);
  solution_transfer_m1.prepare_for_coarsening_and_refinement(this->solution_m1);
  solution_transfer_m2.prepare_for_coarsening_and_refinement(this->solution_m2);
  solution_transfer_m3.prepare_for_coarsening_and_refinement(this->solution_m3);

  tria.execute_coarsening_and_refinement();
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
  this->nonzero_constraints.distribute(tmp);
  this->nonzero_constraints.distribute(tmp_m1);
  this->nonzero_constraints.distribute(tmp_m2);
  this->nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  this->present_solution = tmp;
  this->solution_m1      = tmp_m1;
  this->solution_m2      = tmp_m2;
  this->solution_m3      = tmp_m3;
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
  solution_transfer.prepare_for_coarsening_and_refinement(
    this->present_solution);
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
  this->nonzero_constraints.distribute(tmp);
  this->nonzero_constraints.distribute(tmp_m1);
  this->nonzero_constraints.distribute(tmp_m2);
  this->nonzero_constraints.distribute(tmp_m3);

  // Fix on the new mesh
  this->present_solution = tmp;
  this->solution_m1      = tmp_m1;
  this->solution_m2      = tmp_m2;
  this->solution_m3      = tmp_m3;
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::postprocess(bool firstIter)
{
  if (this->simulationControl.isOutputIteration())
    this->write_output_results(this->present_solution,
                               this->pvdhandler,
                               this->simulationControl.getOutputFolder(),
                               this->simulationControl.getOuputName(),
                               this->simulationControl.getIter(),
                               this->simulationControl.getTime(),
                               this->simulationControl.getSubdivision(),
                               this->simulationControl.getGroupFiles());

  if (this->nsparam.postProcessingParameters.calculate_enstrophy)
    {
      double enstrophy = calculate_enstrophy(this->fe,
                                             this->dof_handler,
                                             this->present_solution,
                                             nsparam.femParameters,
                                             mpi_communicator);

      this->enstrophy_table.add_value("time",
                                      this->simulationControl.getTime());
      this->enstrophy_table.add_value("enstrophy", enstrophy);

      // Display Enstrophy to screen if verbosity is enabled
      if (this->nsparam.postProcessingParameters.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Enstrophy  : " << enstrophy << std::endl;
        }

      // Output Enstrophy to a text file from processor 0
      if (this->simulationControl.getIter() %
              this->nsparam.postProcessingParameters.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            nsparam.postProcessingParameters.enstrophy_output_name + ".dat";
          std::ofstream output(filename.c_str());
          enstrophy_table.set_precision("time", 12);
          enstrophy_table.set_precision("enstrophy", 12);
          this->enstrophy_table.write_text(output);
        }
    }

  if (this->nsparam.postProcessingParameters.calculate_kinetic_energy)
    {
      TimerOutput::Scope t(this->computing_timer, "kinetic_energy_calculation");
      double             kE = calculate_kinetic_energy(this->fe,
                                           this->dof_handler,
                                           this->present_solution,
                                           nsparam.femParameters,
                                           mpi_communicator);
      this->kinetic_energy_table.add_value("time",
                                           this->simulationControl.getTime());
      this->kinetic_energy_table.add_value("kinetic-energy", kE);
      if (this->nsparam.postProcessingParameters.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Kinetic energy : " << kE << std::endl;
        }

      // Output Kinetic Energy to a text file from processor 0
      if (this->simulationControl.getIter() %
              this->nsparam.postProcessingParameters.output_frequency ==
            0 &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            nsparam.postProcessingParameters.kinetic_energy_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          kinetic_energy_table.set_precision("time", 12);
          kinetic_energy_table.set_precision("kinetic-energy", 12);
          this->kinetic_energy_table.write_text(output);
        }
    }

  if (!firstIter)
    {
      // Calculate forces on the boundary conditions
      if (this->nsparam.forcesParameters.calculate_force)
        {
          if (this->simulationControl.getIter() %
                this->nsparam.forcesParameters.calculation_frequency ==
              0)
            this->postprocessing_forces(this->present_solution,
                                        this->simulationControl);
          if (this->simulationControl.getIter() %
                this->nsparam.forcesParameters.output_frequency ==
              0)
            this->write_output_forces();
        }

      // Calculate torques on the boundary conditions
      if (this->nsparam.forcesParameters.calculate_torque)
        {
          if (this->simulationControl.getIter() %
                this->nsparam.forcesParameters.calculation_frequency ==
              0)
            this->postprocessing_torques(this->present_solution,
                                         this->simulationControl);
          if (this->simulationControl.getIter() %
                this->nsparam.forcesParameters.output_frequency ==
              0)
            this->write_output_torques();
        }

      // Calculate error with respect to analytical solution
      if (this->nsparam.analyticalSolution->calculate_error())
        {
          // Update the time of the exact solution to the actual time
          this->exact_solution->set_time(this->simulationControl.getTime());
          const std::pair<double, double> errors =
            this->calculate_L2_error(this->present_solution);
          const double error_velocity = errors.first;
          const double error_pressure = errors.second;
          if (this->simulationControl.getMethod() ==
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
              this->error_table.add_value("time",
                                          this->simulationControl.getTime());
              this->error_table.add_value("error_velocity", error_velocity);
            }
          if (this->nsparam.analyticalSolution->verbosity ==
              Parameters::Verbosity::verbose)
            {
              this->pcout << "L2 error velocity : " << error_velocity
                          << std::endl;
            }
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "read_checkpoint");
  std::string        prefix = this->nsparam.restartParameters.filename;
  this->simulationControl.read(prefix);
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

  VectorType distributed_system(this->newton_update);
  VectorType distributed_system_m1(this->newton_update);
  VectorType distributed_system_m2(this->newton_update);
  VectorType distributed_system_m3(this->newton_update);
  x_system[0] = &(distributed_system);
  x_system[1] = &(distributed_system_m1);
  x_system[2] = &(distributed_system_m2);
  x_system[3] = &(distributed_system_m3);
  parallel::distributed::SolutionTransfer<dim, VectorType> system_trans_vectors(
    this->dof_handler);
  system_trans_vectors.deserialize(x_system);
  this->present_solution = distributed_system;
  this->solution_m1      = distributed_system_m1;
  this->solution_m2      = distributed_system_m2;
  this->solution_m3      = distributed_system_m3;
}

/*
 * Reads a CFD Mesh from a GMSH file or generates a pre-defined primitive
 */
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::read_mesh()
{
  // GMSH input
  if (this->nsparam.mesh.type == Parameters::Mesh::Type::gmsh)
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(*this->triangulation);
      std::ifstream input_file(this->nsparam.mesh.file_name);
      grid_in.read_msh(input_file);
    }

  // Dealii grids
  else if (this->nsparam.mesh.type == Parameters::Mesh::Type::dealii)
    {
      GridGenerator::generate_from_name_and_arguments(
        *this->triangulation,
        this->nsparam.mesh.grid_type,
        this->nsparam.mesh.grid_arguments);
    }
  else
    throw std::runtime_error(
      "Unsupported mesh type - mesh will not be created");

  const int initialSize = this->nsparam.mesh.initialRefinement;
  this->set_periodicity();
  this->triangulation->refine_global(initialSize);
}

/*
 * Periodicity
 */
template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::set_periodicity()
{
  // Setup parallelism for periodic boundary conditions
  for (unsigned int i_bc = 0; i_bc < nsparam.boundaryConditions.size; ++i_bc)
    {
      if (nsparam.boundaryConditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::periodic)
        {
          std::vector<GridTools::PeriodicFacePair<
            typename Triangulation<dim>::cell_iterator>>
            periodicity_vector;
          GridTools::collect_periodic_faces(
            *dynamic_cast<Triangulation<dim> *>(this->triangulation.get()),
            nsparam.boundaryConditions.id[i_bc],
            nsparam.boundaryConditions.periodic_id[i_bc],
            nsparam.boundaryConditions.periodic_direction[i_bc],
            periodicity_vector);
          this->triangulation->add_periodicity(periodicity_vector);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::set_nodal_values()
{
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  const MappingQ<dim>              mapping(this->degreeVelocity_,
                              this->nsparam.femParameters.qmapping_all);
  VectorTools::interpolate(mapping,
                           this->dof_handler,
                           this->nsparam.initialCondition->uvwp,
                           this->newton_update,
                           this->fe.component_mask(velocities));
  VectorTools::interpolate(mapping,
                           this->dof_handler,
                           this->nsparam.initialCondition->uvwp,
                           this->newton_update,
                           this->fe.component_mask(pressure));
  this->nonzero_constraints.distribute(this->newton_update);
  this->present_solution = this->newton_update;
}


template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_results(
  const VectorType & solution,
  PVDHandler &       pvdhandler,
  const std::string  folder,
  const std::string  solutionName,
  const unsigned int iter,
  const double       time,
  const unsigned int subdivision,
  const unsigned int group_files)
{
  TimerOutput::Scope            t(this->computing_timer, "output");
  const MappingQ<dim>           mapping(this->degreeVelocity_,
                              nsparam.femParameters.qmapping_all);
  vorticity_postprocessor<dim>  vorticity;
  qcriterion_postprocessor<dim> qcriterion;


  SRF_postprocessor<dim>   srf(nsparam.velocitySource.omega_x,
                             nsparam.velocitySource.omega_y,
                             nsparam.velocitySource.omega_z);
  std::vector<std::string> solution_names(dim, "velocity");
  solution_names.push_back("pressure");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);

  DataOutBase::VtkFlags flags;
  if (this->degreeVelocity_ > 1)
    flags.write_higher_order_cells = true;

  DataOut<dim> data_out;
  data_out.set_flags(flags);
  data_out.attach_dof_handler(this->dof_handler);
  data_out.add_data_vector(solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);
  data_out.add_data_vector(solution, vorticity);
  data_out.add_data_vector(solution, qcriterion);
  if (nsparam.velocitySource.type ==
      Parameters::VelocitySource::VelocitySourceType::srf)
    {
      data_out.add_data_vector(solution, srf);
    }
  Vector<float> subdomain(this->triangulation->n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = this->triangulation->locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches(mapping,
                         subdivision,
                         DataOut<dim>::curved_inner_cells);

  const int my_id = Utilities::MPI::this_mpi_process(this->mpi_communicator);

  // Write master files (.pvtu,.pvd,.visit) on the master process
  if (my_id == 0)
    {
      std::vector<std::string> filenames;
      const unsigned int       n_processes =
        Utilities::MPI::n_mpi_processes(this->mpi_communicator);
      const unsigned int n_files =
        (group_files == 0) ? n_processes : std::min(group_files, n_processes);

      for (unsigned int i = 0;
           i <
           n_files; // Utilities::MPI::n_mpi_processes(this->mpi_communicator);
           ++i)
        filenames.push_back(solutionName + "." +
                            Utilities::int_to_string(iter, 4) + "." +
                            Utilities::int_to_string(i, 4) + ".vtu");

      std::string pvtu_filename =
        (solutionName + "." + Utilities::int_to_string(iter, 4) + ".pvtu");
      std::ofstream master_output((folder + pvtu_filename).c_str());

      data_out.write_pvtu_record(master_output, filenames);

      const std::string pvdPrefix = (folder + solutionName);
      pvdhandler.append(time, pvtu_filename);
      std::ofstream pvd_output(pvdPrefix + ".pvd");
      DataOutBase::write_pvd_record(pvd_output, pvdhandler.times_and_names_);
    }

  const unsigned int my_file_id =
    (group_files == 0 ? my_id : my_id % group_files);
  int color = my_id % group_files;

  {
    MPI_Comm comm;
    MPI_Comm_split(this->mpi_communicator, color, my_id, &comm);
    const std::string filename =
      (folder + solutionName + "." + Utilities::int_to_string(iter, 4) + "." +
       Utilities::int_to_string(my_file_id, 4) + ".vtu");
    data_out.write_vtu_in_parallel(filename.c_str(), comm);

    MPI_Comm_free(&comm);
  }

  if (nsparam.postProcessingParameters.output_boundaries)
    {
      DataOutFaces<dim>           data_out_faces;
      boundary_postprocessor<dim> manifold;
      data_out_faces.attach_dof_handler(this->dof_handler);
      data_out_faces.add_data_vector(solution, manifold);
      data_out_faces.build_patches();

      int color = my_id % 1;

      MPI_Comm comm;
      MPI_Comm_split(this->mpi_communicator, color, my_id, &comm);

      const std::string face_filename =
        (folder + "boundaries." + Utilities::int_to_string(iter, 4) + ".vtu");
      data_out_faces.write_vtu_in_parallel(face_filename.c_str(), comm);

      MPI_Comm_free(&comm);
    }
}

template <int dim, typename VectorType, typename DofsType>
void
NavierStokesBase<dim, VectorType, DofsType>::write_output_forces()
{
  TimerOutput::Scope t(this->computing_timer, "output_forces");
  for (unsigned int boundary_id = 0;
       boundary_id < nsparam.boundaryConditions.size;
       ++boundary_id)
    {
      std::string filename = nsparam.forcesParameters.force_output_name + "." +
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
       boundary_id < nsparam.boundaryConditions.size;
       ++boundary_id)
    {
      std::string filename = nsparam.forcesParameters.torque_output_name + "." +
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
  std::string        prefix = this->nsparam.restartParameters.filename;
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    this->simulationControl.save(prefix);
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    this->pvdhandler.save(prefix);

  std::vector<const VectorType *> sol_set_transfer;
  sol_set_transfer.push_back(&this->present_solution);
  sol_set_transfer.push_back(&this->solution_m1);
  sol_set_transfer.push_back(&this->solution_m2);
  sol_set_transfer.push_back(&this->solution_m3);
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
