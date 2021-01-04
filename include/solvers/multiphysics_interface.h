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
 * This class provides an interface for multiphysics simulations by enabling
 * the solution of multiple auxiliary physics on top of a computational
 * fluid dynamics simulation. The auxiliary physics are stored in a map
 * whose key are the Parameters::PhysicsID int enum.
 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020
 */

#ifndef lethe_multiphysics_interface_h
#define lethe_multiphysics_interface_h

#include <deal.II/base/exceptions.h>

#include <core/multiphysics.h>
#include <core/parameters_multiphysics.h>
#include <solvers/auxiliary_physics.h>
#include <solvers/navier_stokes_solver_parameters.h>

#include <map>
#include <memory>

using namespace dealii;

template <int dim>
class MultiphysicsInterface
{
  /** @brief Construct the Multiphysics interface from the simulation parameters.
   * Depending on which multiphysics element is enabled, the appropraite
   * auxiliary physics is instantiated.
   *
   */
  MultiphysicsInterface(
    const NavierStokesSolverParameters<dim> &nsparam,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                        p_triangulation,
    ConditionalOStream &p_pcout);

  /**
   * @brief Call for the assembly of the matrix and the right-hand side for the appropiate
   * physics.
   *
   * @param physics_id Identifier of the physics being currently solver
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  assemble_matrix_and_rhs(
    const PhysicsID physics_id,
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method)
  {
    AssertThrow(
      std::find(enabled_physics.begin(), enabled_physics.end(), physics_id) !=
        enabled_physics.end(),
      "This physics has not been constructed by the multiphysics interface.");

    physics[physics_id]->assemble_matrix_and_rhs(time_stepping_method);
  }

  /**
   * @brief Call for the assembly of right-hand side for the appropriate physics.
   *
   * @param physics_id Identifier of the physics being currently solver
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  assemble_rhs(const PhysicsID physics_id,
               const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method)
  {
    AssertThrow(
      std::find(enabled_physics.begin(), enabled_physics.end(), physics_id) !=
        enabled_physics.end(),
      "This physics has not been constructed by the multiphysics interface.");

    physics[physics_id]->assemble_rhs(time_stepping_method);
  }

  /**
   * @brief Call the attachment of the solution vector to the data out for enabled
   * auxiliary physics.
   */
  void
  attach_solution_to_output(DataOut<dim> &data_out)

  {
    for (auto &iphys : physics)
      {
        iphys.second->attach_solution_to_output(data_out);
      }
  }

  /**
   * @brief Carry out the operations required to finish a simulation correctly for
   * all auxiliary physics.
   */
  void
  finish_simulation()
  {
    for (auto &iphys : physics)
      {
        iphys.second->finish_simulation();
      }
  }

  /**
   * @brief Carry out the operations require to finish a time step correctly for
   * all auxiliary physics.
   */
  void
  finish_time_step()
  {
    for (auto &iphys : physics)
      {
        iphys.second->finish_time_step();
      }
  }

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  void
  postprocess()
  {
    for (auto &iphys : physics)
      {
        iphys.second->postprocess();
      }
  }



  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs()
  {
    for (auto &iphys : physics)
      {
        iphys.second->setup_dofs();
      }
  };

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions()
  {
    for (auto &iphys : physics)
      {
        iphys.second->set_initial_conditions();
      }
  };

  /**
   * @brief Call for the solution of the linear system of an auxiliary physics.
   *
   * @param initial_step Provides the linear solver with indication if this solution is the first
   * one for the system of equation or not.
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix has been recalculated or not.
   */
  virtual void
  solve_linear_system(const PhysicsID physics_id,
                      const bool      initial_step,
                      const bool      renewed_matrix = true)
  {
    AssertThrow(
      std::find(enabled_physics.begin(), enabled_physics.end(), physics_id) !=
        enabled_physics.end(),
      "This physics has not been constructed by the multiphysics interface.")

      physics[physics_id]
        ->solver_linear_system(initial_step, renewed_matrix);
  };



private:
  const Parameters::Multiphysics multiphysics_parameters;

  // Data structure to store all physics which were enabled
  std::vector<PhysicsID> enabled_physics;

  // Auxiliary physics are stored within a map of shared pointer to ensure
  // proper deallocation.
  // std::shared_ptr<AuxiliaryPhysics<dim>>
  std::map<PhysicsID, std::shared_ptr<AuxiliaryPhysics<dim>>> physics;
};


#endif
