/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * Implementation of heat transfer as an auxiliary physics.
 * This heat equation is weakly coupled to the velocity field.
 * Equation solved:
 * rho * Cp * (dT/dt + u.gradT) = k div(gradT) + nu/rho * (gradu : gradu)
 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020-
 */

#ifndef lethe_heat_transfer_h
#define lethe_heat_transfer_h

#include <solvers/auxiliary_physics.h>


template <int dim>
class HeatTransfer : public AuxiliaryPhysics<dim>
{
public:
  HeatTransfer<dim>(
    const NavierStokesSolverParameters<dim> &p_simulation_parameters,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                        p_triangulation,
    ConditionalOStream &p_pcout)
    : simulation_parameters(p_simulation_parameters)
    , triangulation(p_triangulation)
    , pcout(p_pcout)
  {}

  /**
   * @brief Call for the assembly of the matrix and the right-hand side.
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  /**
   * @brief Call for the assembly of the right-hand side
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  /**
   * @brief Attach the solution vector to the data out provided. This function
   * enable the auxiliary physics to output their solution via the core solver.
   */
  virtual void
  attach_solution_to_output(DataOut<dim> &data_out) override;


  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  virtual void
  finish_simulation() override;

  /**
   * @brief Carry out the operations require to finish a time step correctly. This
   * includes setting the previous values
   */
  virtual void
  finish_time_step() override;

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess() override;

  /**
   * @brief Returns the dof_handler of the heat transfer physics
   */
  virtual const DoFHandler<dim> &
  get_dof_handler() override{};

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() override;

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions() override;

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
                      const bool renewed_matrix = true);



private:
  template <bool assemble_matrix>
  void
  assemble_system(const Parameters::SimulationControl::TimeSteppingMethod
                    time_stepping_method);

  // Simulation parameters
  // TODO : Refactor to a more neutral name
  const NavierStokesSolverParameters<dim> &simulation_parameters;


  // Triangulation. Stored as a shared_ptr since it is shared with the core
  // solver
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;

  // Condition Output Stream. Constructed from the one of the core solver.
  ConditionalOStream pcout;
};


#endif
