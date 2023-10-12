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
 * ---------------------------------------------------------------------*/

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/multiphysics.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/mf_navier_stokes.h>

#include <deal.II/base/multithread_info.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

template <int dim>
MFNavierStokesSolver<dim>::MFNavierStokesSolver(
  SimulationParameters<dim> &nsparam)
  : NavierStokesBase<dim, VectorType, IndexSet>(nsparam)
{
  AssertThrow(
    nsparam.fem_parameters.velocity_order ==
      nsparam.fem_parameters.pressure_order,
    dealii::ExcMessage(
      "Matrix free Navier-Stokes does not support different orders for the velocity and the pressure!"));

  this->fe = std::make_shared<FESystem<dim>>(
    FE_Q<dim>(nsparam.fem_parameters.velocity_order), dim + 1);
  if ((nsparam.stabilization.use_default_stabilization == true) ||
      nsparam.stabilization.stabilization ==
        Parameters::Stabilization::NavierStokesStabilization::pspg_supg)
    {
      if (is_bdf(this->simulation_control->get_assembly_method()))
        system_operator = std::make_shared<
          NavierStokesTransientSUPGPSPGOperator<dim, double>>();
      else
        system_operator =
          std::make_shared<NavierStokesSUPGPSPGOperator<dim, double>>();
    }
  else
    throw std::runtime_error(
      "Only SUPG/PSPG stabilization is supported at the moment.");
}

template <int dim>
MFNavierStokesSolver<dim>::~MFNavierStokesSolver()
{
  this->dof_handler.clear();
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve()
{
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);
  this->update_multiphysics_time_average_solution();

  while (this->simulation_control->integrate())
    {
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());

      if ((this->simulation_control->get_step_number() %
               this->simulation_parameters.mesh_adaptation.frequency !=
             0 ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          update_boundary_conditions();
        }

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

      if (!this->simulation_control->is_at_start())
        {
          NavierStokesBase<dim, VectorType, IndexSet>::refine_mesh();
        }

      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          calculate_time_derivative_previous_solutions();
          this->time_derivative_previous_solutions.update_ghost_values();
          this->system_operator->evaluate_time_derivative_previous_solutions(
            this->time_derivative_previous_solutions);
        }

      this->iterate();

      this->postprocess(false);
      this->finish_time_step();
    }

  this->finish_simulation();
}

template <int dim>
void
MFNavierStokesSolver<dim>::setup_dofs_fd()
{
  TimerOutput::Scope t(this->computing_timer, "Setup DoFs");

  // Clear the preconditioners
  ilu_preconditioner.reset();
  gc_multigrid_preconditioner.reset();
  ls_multigrid_preconditioner.reset();

  // Clear matrix free operator
  this->system_operator->clear();

  // Fill the dof handler and initialize vectors
  this->dof_handler.distribute_dofs(*this->fe);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    this->dof_handler.distribute_mg_dofs();

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

  // Non Zero constraints
  define_non_zero_constraints();

  // Zero constraints
  define_zero_constraints();

  // Initialize matrix-free object
  unsigned int mg_level = numbers::invalid_unsigned_int;
  this->system_operator->reinit(
    *this->mapping,
    this->dof_handler,
    this->zero_constraints,
    *this->cell_quadrature,
    &(*this->forcing_function),
    this->simulation_parameters.physical_properties_manager
      .get_kinematic_viscosity_scale(),
    mg_level,
    this->simulation_control);


  // Initialize vectors using operator
  this->system_operator->initialize_dof_vector(this->present_solution);
  this->system_operator->initialize_dof_vector(this->evaluation_point);
  this->system_operator->initialize_dof_vector(this->newton_update);
  this->system_operator->initialize_dof_vector(this->system_rhs);
  this->system_operator->initialize_dof_vector(this->local_evaluation_point);
  this->system_operator->initialize_dof_vector(
    this->time_derivative_previous_solutions);

  // Initialize vectors of previous solutions
  for (auto &solution : this->previous_solutions)
    {
      this->system_operator->initialize_dof_vector(solution);
    }

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->average_velocities->initialize_vectors(
        this->locally_owned_dofs,
        this->locally_relevant_dofs,
        this->fe->n_dofs_per_vertex(),
        this->mpi_communicator);

      if (this->simulation_parameters.restart_parameters.checkpoint)
        {
          this->average_velocities->initialize_checkpoint_vectors(
            this->locally_owned_dofs,
            this->locally_relevant_dofs,
            this->mpi_communicator);
        }
    }

  double global_volume =
    GridTools::volume(*this->triangulation, *this->mapping);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;
}

template <int dim>
void
MFNavierStokesSolver<dim>::update_boundary_conditions()
{
  double time = this->simulation_control->get_current_time();
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
        .u.set_time(time);
      this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
        .v.set_time(time);
      this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
        .w.set_time(time);
      this->simulation_parameters.boundary_conditions.bcPressureFunction[i_bc]
        .p.set_time(time);
    }
  define_non_zero_constraints();
  // Distribute constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim>
void
MFNavierStokesSolver<dim>::set_initial_condition_fd(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->present_solution.update_ghost_values();
      this->finish_time_step();
    }
  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      AssertThrow(
        this->simulation_control->is_steady(),
        dealii::ExcMessage(
          "The lethe-fluid-matrix-free solver does not support viscous initial conditions for transient simulations."));

      // Set the nodal values to have an initial condition that is adequate
      this->set_nodal_values();
      this->present_solution.update_ghost_values();

      // Create a pointer to the current viscosity model
      std::shared_ptr<RheologicalModel> viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      viscosity_model->set_kinematic_viscosity(
        this->simulation_parameters.initial_condition->kinematic_viscosity);

      // Set the kinematic viscosity in the system operator to be the temporary
      // viscosity
      this->system_operator->set_kinematic_viscosity(
        this->simulation_parameters.physical_properties_manager
          .get_kinematic_viscosity_scale());

      // Solve the problem with the temporary viscosity
      PhysicsSolver<LinearAlgebra::distributed::Vector<double>>::
        solve_non_linear_system(false);
      this->finish_time_step();

      // Set the kinematic viscosity in the system operator to be the original
      // viscosity
      this->system_operator->set_kinematic_viscosity(
        this->simulation_parameters.physical_properties_manager
          .get_kinematic_viscosity_scale());
    }
  else if (initial_condition_type == Parameters::InitialConditionType::ramp)
    {
      this->pcout << "*********************************" << std::endl;
      this->pcout << " Initial condition using ramp " << std::endl;
      this->pcout << "*********************************" << std::endl;

      // Set the nodal values to have an initial condition that is adequate
      this->set_nodal_values();
      this->present_solution.update_ghost_values();

      // Create a pointer to the current viscosity model
      std::shared_ptr<RheologicalModel> viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Gather viscosity final parameters
      const double viscosity_end = viscosity_model->get_kinematic_viscosity();

      // Kinematic viscosity ramp parameters
      const int n_iter_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .n_iter;
      double kinematic_viscosity =
        n_iter_viscosity > 0 ?
          this->simulation_parameters.initial_condition->ramp.ramp_viscosity
            .kinematic_viscosity_init :
          viscosity_end;
      const double alpha_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .alpha;

      viscosity_model->set_kinematic_viscosity(kinematic_viscosity);

      // Ramp on kinematic viscosity
      for (int i = 0; i < n_iter_viscosity; ++i)
        {
          this->pcout << std::setprecision(4)
                      << "********* Solution for kinematic viscosity = " +
                           std::to_string(kinematic_viscosity) + " *********"
                      << std::endl;

          viscosity_model->set_kinematic_viscosity(kinematic_viscosity);

          // Set the kinematic viscosity in the system operator to be the
          // temporary viscosity
          this->system_operator->set_kinematic_viscosity(
            this->simulation_parameters.physical_properties_manager
              .get_kinematic_viscosity_scale());

          this->simulation_control->set_assembly_method(
            Parameters::SimulationControl::TimeSteppingMethod::steady);
          // Solve the problem with the temporary viscosity
          PhysicsSolver<LinearAlgebra::distributed::Vector<double>>::
            solve_non_linear_system(false);
          this->finish_time_step();

          kinematic_viscosity +=
            alpha_viscosity * (viscosity_end - kinematic_viscosity);
        }
      // Reset kinematic viscosity to simulation parameters
      viscosity_model->set_kinematic_viscosity(viscosity_end);
      this->system_operator->set_kinematic_viscosity(viscosity_end);
    }
  else
    {
      throw std::runtime_error(
        "Type of initial condition is not supported by MF Navier-Stokes");
    }
}

template <int dim>
void
MFNavierStokesSolver<dim>::assemble_system_matrix()
{
  // Required for compilation but not used for matrix free solvers.
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");
}

template <int dim>
void
MFNavierStokesSolver<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_operator->evaluate_residual(this->system_rhs,
                                           this->evaluation_point);

  this->system_rhs *= -1.0;
}

template <int dim>
void
MFNavierStokesSolver<dim>::update_multiphysics_time_average_solution()
{
  // TODO
}

template <int dim>
void
MFNavierStokesSolver<dim>::calculate_time_derivative_previous_solutions()
{
  this->time_derivative_previous_solutions = 0;

  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);

  for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
    {
      this->time_derivative_previous_solutions.add(bdf_coefs[p + 1],
                                                   this->previous_solutions[p]);
    }
}

template <int dim>
double
MFNavierStokesSolver<dim>::estimate_omega(
  std::shared_ptr<NavierStokesOperatorBase<dim, double>> &mg_operator)
{
  double omega = 0.0;

  using OperatorType               = NavierStokesOperatorBase<dim, double>;
  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using ChebyshevPreconditionerType =
    PreconditionChebyshev<OperatorType, VectorType, SmootherPreconditionerType>;
  typename ChebyshevPreconditionerType::AdditionalData
    chebyshev_additional_data;

  chebyshev_additional_data.preconditioner =
    std::make_shared<SmootherPreconditionerType>();
  mg_operator->compute_inverse_diagonal(
    chebyshev_additional_data.preconditioner->get_vector());
  chebyshev_additional_data.constraints.copy_from(this->zero_constraints);
  chebyshev_additional_data.degree =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .eig_estimation_degree;
  chebyshev_additional_data.smoothing_range =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .eig_estimation_smoothing_range;
  chebyshev_additional_data.eig_cg_n_iterations =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .eig_estimation_cg_n_iterations;
  chebyshev_additional_data.eigenvalue_algorithm = ChebyshevPreconditionerType::
    AdditionalData::EigenvalueAlgorithm::power_iteration;
  chebyshev_additional_data.polynomial_type =
    ChebyshevPreconditionerType::AdditionalData::PolynomialType::fourth_kind;

  auto chebyshev = std::make_shared<ChebyshevPreconditionerType>();
  chebyshev->initialize(*mg_operator, chebyshev_additional_data);

  VectorType vec;
  mg_operator->initialize_dof_vector(vec);

  const auto evs = chebyshev->estimate_eigenvalues(vec);

  const double alpha =
    (chebyshev_additional_data.smoothing_range > 1. ?
       evs.max_eigenvalue_estimate / chebyshev_additional_data.smoothing_range :
       std::min(0.9 * evs.max_eigenvalue_estimate,
                evs.min_eigenvalue_estimate));

  omega = 2.0 / (alpha + evs.max_eigenvalue_estimate);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .eig_estimation_verbose != Parameters::Verbosity::quiet)
    {
      this->pcout << "    - minimum eigenvalue: " << evs.min_eigenvalue_estimate
                  << std::endl;
      this->pcout << "    - maximum eigenvalue: " << evs.max_eigenvalue_estimate
                  << std::endl;
      this->pcout << "    - relaxation parameter: " << omega << std::endl;
    }

  return omega;
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve_with_LSMG(SolverGMRES<VectorType> &solver)
{
  this->computing_timer.enter_subsection("Setup LSMG");

  using OperatorType               = NavierStokesOperatorBase<dim, double>;
  using LSTransferType             = MGTransferMatrixFree<dim, double>;
  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType =
    PreconditionRelaxation<OperatorType, SmootherPreconditionerType>;
  using PreconditionerType = PreconditionMG<dim, VectorType, LSTransferType>;

  // Create level objects
  MGLevelObject<std::shared_ptr<OperatorType>> mg_operators;
  LSTransferType                               mg_transfer;
  MGLevelObject<VectorType>                    mg_solution;
  MGLevelObject<VectorType> mg_time_derivative_previous_solutions;
  MGLevelObject<std::shared_ptr<OperatorType>> mg_interface_in;
  MGLevelObject<std::shared_ptr<OperatorType>> mg_interface_out;
  MGLevelObject<AffineConstraints<double>>     level_constraints;
  MGConstrainedDoFs                            mg_constrained_dofs;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_operators;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_in;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_out;
  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners(
    this->dof_handler.get_triangulation().n_global_levels());

  // Extract min and max levels and resize mg level objects accordingly
  const unsigned int n_h_levels =
    this->dof_handler.get_triangulation().n_global_levels();

  unsigned int minlevel = 0;
  unsigned int maxlevel = n_h_levels - 1;

  mg_operators.resize(0, n_h_levels - 1);
  mg_solution.resize(0, n_h_levels - 1);
  mg_time_derivative_previous_solutions.resize(0, n_h_levels - 1);
  level_constraints.resize(0, n_h_levels - 1);
  ls_mg_interface_in.resize(0, n_h_levels - 1);
  ls_mg_interface_out.resize(0, n_h_levels - 1);
  ls_mg_operators.resize(0, n_h_levels - 1);

  // Fill the constraints
  mg_constrained_dofs.initialize(this->dof_handler);

  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);

  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> no_normal_flux_boundaries;
          no_normal_flux_boundaries.insert(
            this->simulation_parameters.boundary_conditions.id[i_bc]);
          for (auto bid : no_normal_flux_boundaries)
            mg_constrained_dofs.make_no_normal_flux_constraints(
              this->dof_handler, bid, 0);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::periodic)
        {
          // IndexSet locally_relevant_dofs;
          // DoFTools::extract_locally_relevant_dofs(this->dof_handler,
          //                                         locally_relevant_dofs);
          // AffineConstraints<double> temp_constraints;
          // DoFTools::make_hanging_node_constraints(this->dof_handler,
          //                                         temp_constraints);
          // std::vector<GridTools::PeriodicFacePair<
          //   typename DoFHandler<dim>::cell_iterator>>
          //   dof_matched_pairs;
          // GridTools::collect_periodic_faces(
          //   this->dof_handler,
          //   this->simulation_parameters.boundary_conditions.id[i_bc],
          //   this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
          //   this->simulation_parameters.boundary_conditions
          //     .periodic_direction[i_bc],
          //   dof_matched_pairs);
          // DoFTools::make_periodicity_constraints<dim, dim>(dof_matched_pairs,
          //                                                  temp_constraints);

          // DoFTools::make_periodicity_constraints(
          //   this->dof_handler,
          //   this->simulation_parameters.boundary_conditions.id[i_bc],
          //
          // this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
          //   this->simulation_parameters.boundary_conditions
          //     .periodic_direction[i_bc],
          //   temp_constraints);
          // mg_constrained_dofs.add_user_constraints(maxlevel,
          // temp_constraints);
          // for (unsigned int level = minlevel; level <= maxlevel; ++level)
          //   mg_constrained_dofs.add_user_constraints(
          //     level, mg_constrained_dofs.get_level_constraints(level));
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::pressure)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::function_weak)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::partial_slip)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::outlet)
        {
          /*do nothing*/
        }
      else
        {
          std::set<types::boundary_id> dirichlet_boundary_id = {
            this->simulation_parameters.boundary_conditions.id[i_bc]};
          mg_constrained_dofs.make_zero_boundary_constraints(
            this->dof_handler,
            dirichlet_boundary_id,
            this->fe->component_mask(velocities));
        }
    }

  // Create mg operators for each level and additional operators needed only for
  // local smoothing
  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(this->dof_handler, level);

      level_constraints[level].reinit(relevant_dofs);
      DoFTools::make_hanging_node_constraints(this->dof_handler,
                                              level_constraints[level]);
      // level_constraints[level].add_lines(
      //   mg_constrained_dofs.get_boundary_indices(level));
      // mg_constrained_dofs.add_user_constraints(
      //   level, mg_constrained_dofs.get_level_constraints(level));
      // level_constraints[level].merge(
      //   mg_constrained_dofs.get_level_constraints(level));
      // level_constraints[level].merge(
      //   mg_constrained_dofs.get_user_constraint_matrix(level));

      level_constraints[level].close();

      if ((this->simulation_parameters.stabilization
             .use_default_stabilization == true) ||
          this->simulation_parameters.stabilization.stabilization ==
            Parameters::Stabilization::NavierStokesStabilization::pspg_supg)
        {
          if (is_bdf(this->simulation_control->get_assembly_method()))
            mg_operators[level] = std::make_shared<
              NavierStokesTransientSUPGPSPGOperator<dim, double>>();
          else
            mg_operators[level] =
              std::make_shared<NavierStokesSUPGPSPGOperator<dim, double>>();
        }

      mg_operators[level]->reinit(
        *this->mapping,
        this->dof_handler,
        level_constraints[level],
        *this->cell_quadrature,
        &(*this->forcing_function),
        this->simulation_parameters.physical_properties_manager
          .get_kinematic_viscosity_scale(),
        level,
        this->simulation_control);

      mg_operators[level]->initialize_dof_vector(mg_solution[level]);
      mg_operators[level]->initialize_dof_vector(
        mg_time_derivative_previous_solutions[level]);

      ls_mg_operators[level].initialize(*mg_operators[level]);
      ls_mg_interface_in[level].initialize(*mg_operators[level]);
      ls_mg_interface_out[level].initialize(*mg_operators[level]);

      partitioners[level] = mg_operators[level]->get_vector_partitioner();
    }

  // Create transfer operator and transfer solution to mg levels
  mg_transfer.initialize_constraints(mg_constrained_dofs);
  mg_transfer.build(this->dof_handler, partitioners);
  mg_transfer.interpolate_to_mg(this->dof_handler,
                                mg_solution,
                                this->present_solution);

  if (is_bdf(this->simulation_control->get_assembly_method()))
    mg_transfer.interpolate_to_mg(this->dof_handler,
                                  mg_time_derivative_previous_solutions,
                                  this->time_derivative_previous_solutions);

  // Evaluate non linear terms for all mg operators
  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      mg_solution[level].update_ghost_values();
      mg_operators[level]->evaluate_non_linear_term(mg_solution[level]);

      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          mg_time_derivative_previous_solutions[level].update_ghost_values();
          mg_operators[level]->evaluate_time_derivative_previous_solutions(
            mg_time_derivative_previous_solutions[level]);
        }
    }

  mg::Matrix<VectorType> mg_matrix(ls_mg_operators);

  // Create smoother, fill parameters for each level and intialize it
  MGSmootherPrecondition<OperatorType, SmootherType, VectorType> mg_smoother;
  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(minlevel,
                                                                     maxlevel);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>();
      mg_operators[level]->compute_inverse_diagonal(
        smoother_data[level].preconditioner->get_vector());
      smoother_data[level].n_iterations =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_smoother_iterations;
      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_eig_estimation)
        smoother_data[level].relaxation = estimate_omega(mg_operators[level]);
      else
        smoother_data[level].relaxation =
          this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_relaxation;
    }

  mg_smoother.initialize(mg_operators, smoother_data);

  // If multigrid number of levels or minimum number of cells in level are
  // specified, change the min level for the coarse-grid solver and the
  // multigrid object, and print levels with appropriate numbering

  int mg_min_level =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_min_level;

  AssertThrow(
    mg_min_level <= static_cast<int>(MGTools::max_level_for_coarse_mesh(
                      this->dof_handler.get_triangulation())),
    ExcMessage(std::string(
      "The maximum level allowed for the coarse mesh (mg min level) is: " +
      std::to_string(MGTools::max_level_for_coarse_mesh(
        this->dof_handler.get_triangulation())) +
      ".")));

  int mg_level_min_cells =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_level_min_cells;

  AssertThrow(
    mg_level_min_cells <=
      static_cast<int>(this->dof_handler.get_triangulation().n_cells(maxlevel)),
    ExcMessage(
      "The mg level min cells specified are larger than the cells of the finest mg level."));


  if (mg_min_level != -1)
    minlevel = mg_min_level;

  if (mg_level_min_cells != -1)
    {
      for (unsigned int level = minlevel; level <= maxlevel; ++level)
        if (static_cast<int>(this->dof_handler.get_triangulation().n_cells(
              level)) >= mg_level_min_cells)
          {
            minlevel = level;
            break;
          }
    }

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_verbosity != Parameters::Verbosity::quiet)
    for (unsigned int level = minlevel; level <= maxlevel; ++level)
      this->pcout << "   MG Level " << level - minlevel << ": "
                  << this->dof_handler.n_dofs(level) << " DoFs, "
                  << this->dof_handler.get_triangulation().n_cells(level)
                  << " cells" << std::endl;

  // Create coarse-grid GMRES solver and AMG preconditioner
  const int max_iterations =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_max_iterations;
  const double tolerance =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_tolerance;
  const double reduce =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_reduce;
  ReductionControl coarse_grid_solver_control(
    max_iterations, tolerance, reduce, false, false);
  SolverGMRES<VectorType>::AdditionalData solver_parameters;
  solver_parameters.max_n_tmp_vectors =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_max_krylov_vectors;

  SolverGMRES<VectorType> coarse_grid_solver(coarse_grid_solver_control,
                                             solver_parameters);

  std::shared_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  TrilinosWrappers::PreconditionAMG precondition_amg;

  TrilinosWrappers::PreconditionILU precondition_ilu;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_coarse_grid_preconditioner ==
      Parameters::LinearSolver::PreconditionerType::amg)
    {
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.elliptic = false;
      if (this->dof_handler.get_fe().degree > 1)
        amg_data.higher_order_elements = true;
      amg_data.n_cycles =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_n_cycles;
      amg_data.w_cycle =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_w_cycles;
      amg_data.aggregation_threshold =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_aggregation_threshold;
      amg_data.smoother_sweeps =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_smoother_sweeps;
      amg_data.smoother_overlap =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_smoother_overlap;
      amg_data.output_details = false;
      amg_data.smoother_type  = "ILU";
      amg_data.coarse_type    = "ILU";
      // Constant modes for velocity and pressure
      // std::vector<std::vector<bool>> constant_modes;
      // ComponentMask                  components(dim + 1, true);
      // DoFTools::extract_constant_modes(this->dof_handler,
      //                                  components,
      //                                  constant_modes);
      // amg_data.constant_modes = constant_modes;

      Teuchos::ParameterList              parameter_ml;
      std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
      amg_data.set_parameters(parameter_ml,
                              distributed_constant_modes,
                              mg_operators[minlevel]->get_system_matrix());
      const double ilu_fill =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_fill;
      const double ilu_atol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_atol;
      const double ilu_rtol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_rtol;
      parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
      parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
      parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

      parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
      parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
      parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);

      precondition_amg.initialize(mg_operators[minlevel]->get_system_matrix(),
                                  parameter_ml);

      mg_coarse = std::make_shared<
        MGCoarseGridIterativeSolver<VectorType,
                                    SolverGMRES<VectorType>,
                                    OperatorType,
                                    decltype(precondition_amg)>>(
        coarse_grid_solver, *mg_operators[minlevel], precondition_amg);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_preconditioner ==
           Parameters::LinearSolver::PreconditionerType::ilu)
    {
      int current_preconditioner_fill_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_fill;
      const double ilu_atol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_atol;
      const double ilu_rtol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_rtol;
      TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
        current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

      precondition_ilu.initialize(mg_operators[minlevel]->get_system_matrix(),
                                  preconditionerOptions);

      mg_coarse = std::make_shared<
        MGCoarseGridIterativeSolver<VectorType,
                                    SolverGMRES<VectorType>,
                                    OperatorType,
                                    decltype(precondition_ilu)>>(
        coarse_grid_solver, *mg_operators[minlevel], precondition_ilu);
    }

  // Create interface matrices needed for local smoothing in case of local
  // refinement
  mg::Matrix<VectorType> mg_interface_matrix_in(ls_mg_interface_in);
  mg::Matrix<VectorType> mg_interface_matrix_out(ls_mg_interface_out);

  // Create main MG object
  Multigrid<VectorType> mg(
    mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother, minlevel);

  if (this->dof_handler.get_triangulation().has_hanging_nodes())
    mg.set_edge_matrices(mg_interface_matrix_in, mg_interface_matrix_out);

  // Create MG preconditioner
  ls_multigrid_preconditioner =
    std::make_shared<PreconditionMG<dim, VectorType, LSTransferType>>(
      this->dof_handler, mg, mg_transfer);

  this->computing_timer.leave_subsection("Setup LSMG");

  this->computing_timer.enter_subsection("Solve linear system");

  solver.solve(*(system_operator),
               this->newton_update,
               this->system_rhs,
               *ls_multigrid_preconditioner);

  this->computing_timer.leave_subsection("Solve linear system");
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve_with_GCMG(SolverGMRES<VectorType> &solver)
{
  this->computing_timer.enter_subsection("Setup GCMG");

  using OperatorType   = NavierStokesOperatorBase<dim, double>;
  using GCTransferType = MGTransferGlobalCoarsening<dim, VectorType>;
  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType =
    PreconditionRelaxation<OperatorType, SmootherPreconditionerType>;
  using PreconditionerType = PreconditionMG<dim, VectorType, GCTransferType>;

  // Create level objects
  MGLevelObject<DoFHandler<dim>>                     dof_handlers;
  MGLevelObject<std::shared_ptr<OperatorType>>       mg_operators;
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;
  MGLevelObject<VectorType>                          mg_solution;
  MGLevelObject<VectorType> mg_time_derivative_previous_solutions;
  MGLevelObject<AffineConstraints<typename VectorType::value_type>> constraints;

  // Create triangulations for all levels
  std::vector<std::shared_ptr<const Triangulation<dim>>>
    coarse_grid_triangulations;

  coarse_grid_triangulations =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      this->dof_handler.get_triangulation());

  // Modify the triangulations if multigrid number of levels or minimum number
  // of cells in level are specified
  std::vector<std::shared_ptr<const Triangulation<dim>>> temp;

  int mg_min_level =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_min_level;

  AssertThrow(
    (mg_min_level + 1) <= static_cast<int>(coarse_grid_triangulations.size()),
    ExcMessage(
      "The mg min level specified is higher than the finest mg level."));

  int mg_level_min_cells =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_level_min_cells;

  AssertThrow(
    mg_level_min_cells <=
      static_cast<int>(
        coarse_grid_triangulations[coarse_grid_triangulations.size() - 1]
          ->n_global_active_cells()),
    ExcMessage(
      "The mg level min cells specified are larger than the cells of the finest mg level."));

  // find first relevant coarse-grid triangulation
  auto ptr =
    std::find_if(coarse_grid_triangulations.begin(),
                 coarse_grid_triangulations.end() - 1,
                 [&mg_min_level, &mg_level_min_cells](const auto &tria) {
                   if (mg_min_level != -1) // minimum number of levels
                     {
                       if ((mg_min_level + 1) <=
                           static_cast<int>(tria->n_global_levels()))
                         return true;
                     }
                   else if (mg_level_min_cells != -1) // minimum number of cells
                     {
                       if (static_cast<int>(tria->n_global_active_cells()) >=
                           mg_level_min_cells)
                         return true;
                     }
                   else
                     {
                       return true;
                     }
                   return false;
                 });

  // consider all triangulations from that one
  while (ptr != coarse_grid_triangulations.end())
    temp.push_back(*(ptr++));

  coarse_grid_triangulations = temp;

  // Extract min and max levels and resize mg level objects accordingly
  const unsigned int n_h_levels = coarse_grid_triangulations.size();

  unsigned int minlevel = 0;
  unsigned int maxlevel = n_h_levels - 1;

  dof_handlers.resize(minlevel, maxlevel);
  mg_operators.resize(minlevel, maxlevel);
  transfers.resize(minlevel, maxlevel);
  mg_solution.resize(minlevel, maxlevel);
  mg_time_derivative_previous_solutions.resize(minlevel, maxlevel);
  constraints.resize(minlevel, maxlevel);

  // Distribute DoFs for each level
  for (unsigned int l = minlevel; l <= maxlevel; ++l)
    {
      dof_handlers[l].reinit(*coarse_grid_triangulations[l]);
      dof_handlers[l].distribute_dofs(this->dof_handler.get_fe());
    }

  // Apply constraints and create mg operators for each level
  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      const auto &level_dof_handler = dof_handlers[level];
      auto       &level_constraint  = constraints[level];

      level_constraint.clear();
      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(level_dof_handler);
      level_constraint.reinit(locally_relevant_dofs);

      DoFTools::make_hanging_node_constraints(level_dof_handler,
                                              level_constraint);

      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar pressure(dim);

      for (unsigned int i_bc = 0;
           i_bc < this->simulation_parameters.boundary_conditions.size;
           ++i_bc)
        {
          if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
              BoundaryConditions::BoundaryType::slip)
            {
              std::set<types::boundary_id> no_normal_flux_boundaries;
              no_normal_flux_boundaries.insert(
                this->simulation_parameters.boundary_conditions.id[i_bc]);
              VectorTools::compute_no_normal_flux_constraints(
                level_dof_handler,
                0,
                no_normal_flux_boundaries,
                level_constraint,
                *this->mapping);
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::periodic)
            {
              DoFTools::make_periodicity_constraints(
                level_dof_handler,
                this->simulation_parameters.boundary_conditions.id[i_bc],
                this->simulation_parameters.boundary_conditions
                  .periodic_id[i_bc],
                this->simulation_parameters.boundary_conditions
                  .periodic_direction[i_bc],
                level_constraint);
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::pressure)
            {
              /*do nothing*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::function_weak)
            {
              /*do nothing*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::partial_slip)
            {
              /*do nothing*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::outlet)
            {
              /*do nothing*/
            }
          else
            {
              VectorTools::interpolate_boundary_values(
                *this->mapping,
                level_dof_handler,
                this->simulation_parameters.boundary_conditions.id[i_bc],
                dealii::Functions::ZeroFunction<dim>(dim + 1),
                level_constraint,
                this->fe->component_mask(velocities));
            }
        }

      level_constraint.close();

      if ((this->simulation_parameters.stabilization
             .use_default_stabilization == true) ||
          this->simulation_parameters.stabilization.stabilization ==
            Parameters::Stabilization::NavierStokesStabilization::pspg_supg)
        {
          if (is_bdf(this->simulation_control->get_assembly_method()))
            mg_operators[level] = std::make_shared<
              NavierStokesTransientSUPGPSPGOperator<dim, double>>();
          else
            mg_operators[level] =
              std::make_shared<NavierStokesSUPGPSPGOperator<dim, double>>();
        }

      mg_operators[level]->reinit(
        *this->mapping,
        level_dof_handler,
        level_constraint,
        *this->cell_quadrature,
        &(*this->forcing_function),
        this->simulation_parameters.physical_properties_manager
          .get_kinematic_viscosity_scale(),
        numbers::invalid_unsigned_int,
        this->simulation_control);
    }

  // Create transfer operators and transfer solution to mg levels
  for (unsigned int level = minlevel; level < maxlevel; ++level)
    transfers[level + 1].reinit(dof_handlers[level + 1],
                                dof_handlers[level],
                                constraints[level + 1],
                                constraints[level]);

  MGTransferGlobalCoarsening<dim, VectorType> mg_transfer(
    transfers, [&](const auto l, auto &vec) {
      mg_operators[l]->initialize_dof_vector(vec);
    });

  mg_transfer.interpolate_to_mg(this->dof_handler,
                                mg_solution,
                                this->present_solution);

  if (is_bdf(this->simulation_control->get_assembly_method()))
    mg_transfer.interpolate_to_mg(this->dof_handler,
                                  mg_time_derivative_previous_solutions,
                                  this->time_derivative_previous_solutions);

  // Evaluate non linear terms for all mg operators
  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      mg_solution[level].update_ghost_values();
      mg_operators[level]->evaluate_non_linear_term(mg_solution[level]);

      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          mg_time_derivative_previous_solutions[level].update_ghost_values();
          mg_operators[level]->evaluate_time_derivative_previous_solutions(
            mg_time_derivative_previous_solutions[level]);
        }
    }

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_verbosity != Parameters::Verbosity::quiet)
    for (unsigned int level = minlevel; level <= maxlevel; ++level)
      this->pcout << "   MG Level " << level << ": "
                  << dof_handlers[level].n_dofs() << " DoFs, "
                  << coarse_grid_triangulations[level]->n_global_active_cells()
                  << " cells" << std::endl;

  mg::Matrix<VectorType> mg_matrix(mg_operators);

  // Create smoother, fill parameters for each level and intialize it
  MGSmootherPrecondition<OperatorType, SmootherType, VectorType> mg_smoother;
  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(minlevel,
                                                                     maxlevel);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>();
      mg_operators[level]->compute_inverse_diagonal(
        smoother_data[level].preconditioner->get_vector());
      smoother_data[level].n_iterations =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_smoother_iterations;
      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_eig_estimation)
        smoother_data[level].relaxation = estimate_omega(mg_operators[level]);
      else
        smoother_data[level].relaxation =
          this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_relaxation;
    }

  mg_smoother.initialize(mg_operators, smoother_data);

  // Create coarse-grid GMRES solver and AMG preconditioner
  const int max_iterations =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_max_iterations;
  const double tolerance =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_tolerance;
  const double reduce =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_reduce;
  ReductionControl coarse_grid_solver_control(
    max_iterations, tolerance, reduce, false, false);
  SolverGMRES<VectorType>::AdditionalData solver_parameters;
  solver_parameters.max_n_tmp_vectors =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .mg_coarse_grid_max_krylov_vectors;

  SolverGMRES<VectorType> coarse_grid_solver(coarse_grid_solver_control,
                                             solver_parameters);

  std::shared_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  TrilinosWrappers::PreconditionAMG precondition_amg;

  TrilinosWrappers::PreconditionILU precondition_ilu;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_coarse_grid_preconditioner ==
      Parameters::LinearSolver::PreconditionerType::amg)
    {
      TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
      amg_data.elliptic = false;
      if (this->dof_handler.get_fe().degree > 1)
        amg_data.higher_order_elements = true;
      amg_data.n_cycles =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_n_cycles;
      amg_data.w_cycle =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_w_cycles;
      amg_data.aggregation_threshold =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_aggregation_threshold;
      amg_data.smoother_sweeps =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_smoother_sweeps;
      amg_data.smoother_overlap =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_smoother_overlap;
      amg_data.output_details = false;
      amg_data.smoother_type  = "ILU";
      amg_data.coarse_type    = "ILU";
      // Constant modes for velocity
      std::vector<std::vector<bool>> constant_modes;
      ComponentMask                  components(dim + 1, true);
      DoFTools::extract_constant_modes(dof_handlers[minlevel],
                                       components,
                                       constant_modes);
      amg_data.constant_modes = constant_modes;

      Teuchos::ParameterList              parameter_ml;
      std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
      amg_data.set_parameters(parameter_ml,
                              distributed_constant_modes,
                              mg_operators[minlevel]->get_system_matrix());
      const double ilu_fill =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_fill;
      const double ilu_atol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_atol;
      const double ilu_rtol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_rtol;
      parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
      parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
      parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

      parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
      parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
      parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);

      precondition_amg.initialize(mg_operators[minlevel]->get_system_matrix(),
                                  parameter_ml);

      mg_coarse = std::make_shared<
        MGCoarseGridIterativeSolver<VectorType,
                                    SolverGMRES<VectorType>,
                                    OperatorType,
                                    decltype(precondition_amg)>>(
        coarse_grid_solver, *mg_operators[minlevel], precondition_amg);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_preconditioner ==
           Parameters::LinearSolver::PreconditionerType::ilu)
    {
      int current_preconditioner_fill_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_fill;
      const double ilu_atol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_atol;
      const double ilu_rtol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .ilu_precond_rtol;
      TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
        current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

      precondition_ilu.initialize(mg_operators[minlevel]->get_system_matrix(),
                                  preconditionerOptions);

      mg_coarse = std::make_shared<
        MGCoarseGridIterativeSolver<VectorType,
                                    SolverGMRES<VectorType>,
                                    OperatorType,
                                    decltype(precondition_ilu)>>(
        coarse_grid_solver, *mg_operators[minlevel], precondition_ilu);
    }

  // Create main MG object
  Multigrid<VectorType> mg(
    mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother);

  // Create MG preconditioner
  gc_multigrid_preconditioner =
    std::make_shared<PreconditionMG<dim, VectorType, GCTransferType>>(
      this->dof_handler, mg, mg_transfer);

  this->computing_timer.leave_subsection("Setup GCMG");

  this->computing_timer.enter_subsection("Solve linear system");

  solver.solve(*(system_operator),
               this->newton_update,
               this->system_rhs,
               *gc_multigrid_preconditioner);

  this->computing_timer.leave_subsection("Solve linear system");
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve_with_ILU(SolverGMRES<VectorType> &solver)
{
  this->computing_timer.enter_subsection("Setup ILU");

  int current_preconditioner_fill_level =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_operator->get_system_matrix(),
                                 preconditionerOptions);

  this->computing_timer.leave_subsection("Setup ILU");

  this->computing_timer.enter_subsection("Solve linear system");

  solver.solve(*(system_operator),
               this->newton_update,
               this->system_rhs,
               *ilu_preconditioner);

  this->computing_timer.leave_subsection("Solve linear system");
}

template <int dim>
void
MFNavierStokesSolver<dim>::define_non_zero_constraints()
{
  double time = this->simulation_control->get_current_time();
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  // Non-zero constraints
  auto &nonzero_constraints = this->get_nonzero_constraints();
  {
    nonzero_constraints.clear();
    nonzero_constraints.reinit(this->locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);
    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->simulation_parameters.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              nonzero_constraints,
              *this->mapping);
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .u.set_time(time);
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .v.set_time(time);
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .w.set_time(time);
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .u,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .v,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .w),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
              this->simulation_parameters.boundary_conditions
                .periodic_direction[i_bc],
              nonzero_constraints);
          }
      }
  }

  this->establish_solid_domain(true);

  nonzero_constraints.close();
}

template <int dim>
void
MFNavierStokesSolver<dim>::define_zero_constraints()
{
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  this->zero_constraints.clear();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);
  this->zero_constraints.reinit(this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> no_normal_flux_boundaries;
          no_normal_flux_boundaries.insert(
            this->simulation_parameters.boundary_conditions.id[i_bc]);
          VectorTools::compute_no_normal_flux_constraints(
            this->dof_handler,
            0,
            no_normal_flux_boundaries,
            this->zero_constraints,
            *this->mapping);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            this->simulation_parameters.boundary_conditions.id[i_bc],
            this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
            this->simulation_parameters.boundary_conditions
              .periodic_direction[i_bc],
            this->zero_constraints);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::pressure)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::function_weak)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::partial_slip)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::outlet)
        {
          /*do nothing*/
        }
      else
        {
          VectorTools::interpolate_boundary_values(
            *this->mapping,
            this->dof_handler,
            this->simulation_parameters.boundary_conditions.id[i_bc],
            dealii::Functions::ZeroFunction<dim>(dim + 1),
            this->zero_constraints,
            this->fe->component_mask(velocities));
        }
    }

  this->establish_solid_domain(false);

  this->zero_constraints.close();
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
                                               const bool /* renewed_matrix */)
{
  const double absolute_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .minimum_residual;
  const double relative_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .relative_residual;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .solver == Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step, absolute_residual, relative_residual);
  else
    AssertThrow(false, ExcMessage("This solver is not allowed"));
  this->rescale_pressure_dofs_in_newton_update();
}

template <int dim>
void
MFNavierStokesSolver<dim>::assemble_L2_projection()
{
  // TODO
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
                                              const double absolute_residual,
                                              const double relative_residual)
{
  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverGMRES<VectorType>::AdditionalData solver_parameters;

  solver_parameters.max_n_tmp_vectors =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .max_krylov_vectors;
  solver_parameters.right_preconditioning = true;

  SolverGMRES<VectorType> solver(solver_control, solver_parameters);

  this->present_solution.update_ghost_values();

  this->system_operator->evaluate_non_linear_term(this->present_solution);

  this->newton_update = 0.0;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    solve_with_LSMG(solver);
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    solve_with_GCMG(solver);
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::ilu)
    solve_with_ILU(solver);
  else
    AssertThrow(
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::ilu ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::lsmg ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::gcmg,
      ExcMessage(
        "This linear solver does not support this preconditioner. Only <ilu|lsmg|gcmg> preconditioners are supported."));

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  constraints_used.distribute(this->newton_update);
}

template class MFNavierStokesSolver<2>;
template class MFNavierStokesSolver<3>;
