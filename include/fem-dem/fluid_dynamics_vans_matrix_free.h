// SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#include <solvers/fluid_dynamics_matrix_free.h>

/**
 * @brief A solver for the incompressible Navier-Stokes equations implemented
 * in a matrix-free fashion.
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
   * according to simulation parameters.
   *
   * @param[in] nsparam Relevant parameters for the solver.
   */
  FluidDynamicsVANSMatrixFree(SimulationParameters<dim> &nsparam);

  /**
   * @brief Destructor.
   *
   */
  ~FluidDynamicsVANSMatrixFree();

  /**
   * @brief Solve the problem defined by simulation parameters by iterating
   * through time or through the mesh refinements.
   */
  virtual void
  solve();

private:
};