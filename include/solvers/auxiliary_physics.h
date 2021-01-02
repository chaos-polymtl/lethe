/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020
 */

#ifndef lethe_auxiliary_physics_h
#define lethe_auxiliary_physics_h

#include <deal.II/distributed/tria_base.h>

#include <deal.II/numerics/data_out.h>

#include <core/parameters.h>
#include <solvers/navier_stokes_solver_parameters.h>


/**
 * Au auxiliary physics is defined as a physics that is solved on top of
 * a core physics, the latter being the Navier-Stokes
 * equations. Auxiliary physics are simpler
 * physics around which the entire simulation process does not need to be
 * tailored. Examples of auxiliary physics are temperature and concentration.
 * Generally, the auxiliary physics does not affect the core physics. For
 * example, temperature is advected by the fluid flow, but buoyancy effects are
 * not taken into account.
 *
 * The auxiliary physics are managed by the multiphysics interface of Lethe.
 *
 * This base class is used to establish all of the routines that an auxiliary
 * physics must be able to provide to the Multiphysics interface. These elements
 * are then called at specific moments of a simulation.
 *
 * Auxiliary physics are templated by the dimension of the problem and the
 * vector type and DofsType that is  used to manage their data.
 *
 *
 * Current limitations:
 *
 * - Auxiliary physics are currently expected to be used in conjuction with a
 *flow solver which can provide a velocity field.
 * - Support for feedback from the auxiliary physics to the core physics is, at
 *best, very limited
 * - Support for interaction between auxiliary physics is currently not
 *implemented / expected
 **/

template <int dim>
class AuxiliaryPhysics
{
public:
  /**
   * @brief AuxiliaryPhysics - Constructs an auxiliary physics, but does not set-up its degrees
   * of freedom. This constructor only initialializes the triangulation and
   * class-specific members which are not related to the core solver.
   * @param p_triangulation A shared pointer to the triangulation on which the simulation is being carried out.
   * @param p_pcout A conditional OS stream used by the rest of the simulation.

   */
  AuxiliaryPhysics(
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                             p_triangulation,
    ConditionalOStream &                     p_pcout,
    const NavierStokesSolverParameters<dim> &p_simulation_parameters)
    : triangulation(p_triangulation)
    , pcout(pcout)
    , simulation_parameters(p_simulation_parameters){};

  /**
   * @brief Call for the assembly of the matrix and the right-hand side.
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) = 0;

  /**
   * @brief Call for the assembly of right-hand side only.
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) = 0;

  /**
   * @brief Attach the solution vector to the data out provided. This function
   * enable the auxiliary physics to output their solution via the core solver.
   */
  virtual void
  attach_solution_to_output(DataOut<dim> &data_out) = 0;


  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess() = 0;

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() = 0;

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions() = 0;

  /**
   * @brief Call for the solution of the linear system of equation using a strategy appropriate
   * to the auxiliary physics
   *
   * @param initial_step Provides the linear solver with indication if this solution is the first
   * one for the system of equation or not
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix has been recalculated or not
   */
  virtual void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) = 0;

private:
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;
  const NavierStokesSolverParameters<dim> simulation_parameters;
  ConditionalOStream                      pcout;
};



#endif
