// SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_vans_matrix_free_h
#define lethe_fluid_dynamics_vans_matrix_free_h

#include <solvers/fluid_dynamics_matrix_free.h>

#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/void_fraction.h>


/**
 * @brief A geometric multigrid preconditioner implementation for
 * incompressible VANS equations.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 */
template <int dim>
class MFNavierStokesVANSPreconditionGMG
  : public MFNavierStokesPreconditionGMG<dim>
{
public:
  /// Useful aliases to improve readability.
  using VectorType = typename MFNavierStokesPreconditionGMG<dim>::VectorType;
  using MGVectorType =
    typename MFNavierStokesPreconditionGMG<dim>::MGVectorType;
  using MGNumber = typename MFNavierStokesPreconditionGMG<dim>::MGNumber;

  /**
   * @brief Constructor
   * @param[in] param Complete set of parameters for the CFD-DEM or VANS
   * simulations.
   * @param[in] dof_handler DoFHandler used to solve the VANS equation.
   * @param[in] dof_handler_fe_q_iso_q1 DoFHandler with degree=1 to be used in
   * hp multigrid.
   */
  MFNavierStokesVANSPreconditionGMG(
    const CFDDEMSimulationParameters<dim> &param,
    const DoFHandler<dim>                 &dof_handler,
    const DoFHandler<dim>                 &dof_handler_fe_q_iso_q1);

  /**
   * @brief Creates the operator on a given level
   * @param[in] level Level on which the operator is created..
   */
  void
  create_level_operator(const unsigned int level) override;

  /**
   * @brief Initializes the preconditioner
   * @param[in] simulation_control Simulation control object which is used to
   * obtain information on the time step.
   * @param[in] flow_control Flow controller which is used to provide volume
   * forcing in the flow.
   * @param[in] present_solution Current solution for the VANS equation.
   * @param[in] time_derivative_previous_solutions Pre-calculated time
   * derivative.
   * @param[in] void_fraction_manager Manage for the void fraction. This is used
   * to obtain the void fraction on the levels of the grid.
   */
  void
  initialize(const std::shared_ptr<SimulationControl> &simulation_control,
             FlowControl<dim>                         &flow_control,
             const VectorType                         &present_solution,
             const VectorType            &time_derivative_previous_solutions,
             const VoidFractionBase<dim> &void_fraction_manager);

private:
  /// Reference to the simulation parameters
  const CFDDEMSimulationParameters<dim> &cfd_dem_simulation_parameters;
};

/**
 * @brief A solver for the Volume-Averaged incompressible Navier-Stokes equations
 *  implemented in a matrix-free fashion.
 *
 * A matrix-free stabilized solver for the Volume-Averaged incompressible
 * Navier-Stokes (VANS) equations implemented in a matrix-free fashion with a
 * sum-factorization approach. It uses a continuous Galerkin discretization and
 * solves the fully-coupled discretized problem in a monolithic way using
 * multigrid preconditioners. This class greatly inherits from the matrix-free
 * implementation of the Navier-Stokes equations.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved.
 */
template <int dim>
class FluidDynamicsVANSMatrixFree : public FluidDynamicsMatrixFree<dim>
{
public:
  /**
   * @brief Constructor that sets the finite element degree and system operator
   * according to simulation parameters. It initializes the CFD-DEM parameters,
   * the particle mapping, the particle handler, and the void fraction manager.
   * It also checks whether there are periodic boundaries.
   *
   * @param[in] param Relevant parameters for the solver.
   */
  FluidDynamicsVANSMatrixFree(CFDDEMSimulationParameters<dim> &param);

  /**
   * @brief Solve the problem defined by simulation parameters by iterating
   * through time or through the mesh refinements.
   */
  virtual void
  solve() override;

protected:
  /**
   * @brief Setup the degree of freedom of the flow solver and the degrees of freedom
   * and constraints for the void fraction.
   */
  virtual void
  setup_dofs() override;

  /**
   * @brief Finish the time-step and manage the time-step end for the void
   * fraction.
   */
  virtual void
  finish_time_step_fd();

  /**
   * @brief Add data vectors to the data_out object for post_processing
   * additional results. In this case, the void fraction field is added.
   */
  virtual void
  output_field_hook(DataOut<dim> &data_out) override;

  /**
   * @brief Create geometric multigrid preconditioner.
   */
  void
  create_GMG() override;

  /**
   * Initializes (or re-initializes) geometric multigrid preconditioneré
   */
  void
  initialize_GMG() override;

  /// Simulation parameters for CFD-DEM simulations
  CFDDEMSimulationParameters<dim> cfd_dem_simulation_parameters;

  /// Mapping used for the particles
  MappingQGeneric<dim> particle_mapping;

  /// Particle-handler used to store particles in CFD-DEM simulations
  Particles::ParticleHandler<dim, dim> particle_handler;

  /// Object that manages the void fraction calculation from functions
  /// or from parameters.
  VoidFractionBase<dim> void_fraction_manager;


  /// Member variables which are used to manage boundary conditions
  bool           has_periodic_boundaries;
  Tensor<1, dim> periodic_offset;
  unsigned int   periodic_direction;
};
#endif
