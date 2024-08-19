/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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

#ifndef lethe_physics_subequations_solver_h
#define lethe_physics_subequations_solver_h

#include <core/physics_solver.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/numerics/data_out.h>

/**
 * Sub-equations solved inside a physic (or auxiliary physic) that is not part
 * of the main equations set are solved through the PhysicsSubequationsSolver
 * object. It contains all the common elements of a physic solver. It creates
 * the non-linear solver as specified by the user using the parameters file and
 * provides all the necessary elements needed by the solver to solve a physics
 * problem.
 *
 * @tparam dim Number of dimensions of the problem
 * @tparam VectorType Type of vector container used to store solutions
 */

template <int dim, typename VectorType>
class PhysicsSubequationsSolver : public PhysicsSolver<VectorType>
{
public:
  /**
   * @brief Constructor for physics sub-equations that require the use of a
   * linear and/or non-linear solver within the span of the physics resolution.
   *
   * @param[in] non_linear_solver_parameters Set of parameters that will be used
   * to construct the non-linear solver.
   */
  PhysicsSubequationsSolver(
    const Parameters::NonLinearSolver non_linear_solver_parameters)
    : PhysicsSolver<VectorType>(non_linear_solver_parameters)
  {}

  /**
   * @brief Base destructor
   */
  virtual ~PhysicsSubequationsSolver()
  {}

  /**
   * @brief Sets up the DofHandler and the degrees of freedom associated with
   * the equation to solve.
   */
  virtual void
  setup_dofs() = 0;

  /**
   * @brief Sets up the initial conditions associated with the equation to
   * solve.
   */
  virtual void
  set_initial_conditions() = 0;
};

#endif
