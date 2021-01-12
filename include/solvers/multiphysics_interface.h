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

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <core/multiphysics.h>
#include <core/parameters_multiphysics.h>
#include <core/simulation_control.h>
#include <solvers/auxiliary_physics.h>
#include <solvers/simulation_parameters.h>

#include <map>
#include <memory>

using namespace dealii;

template <int dim>
class MultiphysicsInterface
{
public:
  /** @brief Construct the Multiphysics interface from the simulation parameters.
   * Depending on which multiphysics element is enabled, the appropraite
   * auxiliary physics is instantiated.
   *
   */
  MultiphysicsInterface(
    const SimulationParameters<dim> &nsparam,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control);

  std::vector<PhysicsID>
  get_active_physics()
  {
    return active_physics;
  }

  /**
   * @brief Call for the solution of all physics
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  solve(const Parameters::SimulationControl::TimeSteppingMethod
                   time_stepping_method,
        const bool force_matrix_renewal)
  {
    for (auto &iphys : physics)
      {
        solve_physics(iphys.first, time_stepping_method, force_matrix_renewal);
      }
    for (auto &iphys : block_physics)
      {
        solve_block_physics(iphys.first,
                            time_stepping_method,
                            force_matrix_renewal);
      }
  }


  /**
   * @brief Call for the solution of a single physics
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  solve_physics(const PhysicsID physics_id,
                const Parameters::SimulationControl::TimeSteppingMethod
                           time_stepping_method,
                const bool force_matrix_renewal)
  {
    AssertThrow(std::find(active_physics.begin(),
                          active_physics.end(),
                          physics_id) != active_physics.end(),
                ExcInternalError());

    physics[physics_id]->solve_non_linear_system(time_stepping_method,
                                                 false,
                                                 force_matrix_renewal);
  }

  /**
   * @brief Call for the solution of a single block physics
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  solve_block_physics(const PhysicsID physics_id,
                      const Parameters::SimulationControl::TimeSteppingMethod
                                 time_stepping_method,
                      const bool force_matrix_renewal)
  {
    AssertThrow(std::find(active_physics.begin(),
                          active_physics.end(),
                          physics_id) != active_physics.end(),
                ExcInternalError());

    block_physics[physics_id]->solve_non_linear_system(time_stepping_method,
                                                       false,
                                                       force_matrix_renewal);
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
    for (auto &iphys : block_physics)
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
    for (auto &iphys : block_physics)
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
    for (auto &iphys : block_physics)
      {
        iphys.second->finish_time_step();
      }
  }

  /**
   * @brief Carry out the operations require to percolate the time vectors
   * corectly at the end of a simulation
   */
  void
  percolate_time_vectors()
  {
    for (auto &iphys : physics)
      {
        iphys.second->percolate_time_vectors();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->percolate_time_vectors();
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
  postprocess(bool first_iteration)
  {
    for (auto &iphys : physics)
      {
        iphys.second->postprocess(first_iteration);
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->postprocess(first_iteration);
      }
  }


  /**
   * @brief Prepares the auxiliary physics for mesh adaptation
   */
  void
  prepare_for_mesh_adaptation()
  {
    for (auto &iphys : physics)
      {
        iphys.second->pre_mesh_adaptation();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->pre_mesh_adaptation();
      }
  }

  /**
   * @brief Interpolate solution onto new mesh
   */
  void
  post_mesh_adaptation()
  {
    for (auto &iphys : physics)
      {
        iphys.second->post_mesh_adaptation();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->post_mesh_adaptation();
      }
  }



  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  void
  setup_dofs()
  {
    for (auto &iphys : physics)
      {
        iphys.second->setup_dofs();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->setup_dofs();
      }
  };



  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionnaly support
   * the use of L2 projection or steady-state solutions.
   */
  void
  set_initial_conditions()
  {
    for (auto &iphys : physics)
      {
        iphys.second->set_initial_conditions();
      }

    for (auto &iphys : block_physics)
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
  void
  solve_linear_system(const PhysicsID physics_id,
                      const bool      initial_step,
                      const bool      renewed_matrix = true)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());

    physics[physics_id]->solve_linear_system(initial_step, renewed_matrix);
  };

  /**
   * @brief fluid_dynamics_is_block Verifies if the fluid dynamics solution
   * is stored as a block vector or not.
   *
   * @return boolean value indicating the fluid dynamics is stored as a block vector
   */

  bool
  fluid_dynamics_is_block() const
  {
    return block_physics_solutions.find(PhysicsID::fluid_dynamics) !=
           block_physics_solutions.end();
  }



  /**
   * @brief Request a DOF handler for a given physics ID
   *
   * @param physics_id The physics of the DOF handler being requested
   */
  DoFHandler<dim> *
  get_dof_handler(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return physics_dof_handler[physics_id];
  }



  /**
   * @brief Request the present solution of a given physics
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::Vector *
  get_solution(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return physics_solutions[physics_id];
  }


  /**
   * @brief Request the present block solution of a given physics
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::BlockVector *
  get_block_solution(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return block_physics_solutions[physics_id];
  }


  /**
   * @brief Sets the reference to the DOFHandler of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param dof_handler The dof handler for which the reference is stored
   */
  void
  set_dof_handler(PhysicsID physics_id, DoFHandler<dim> *dof_handler)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_dof_handler[physics_id] = dof_handler;
  }

  /**
   * @brief Sets the reference to the solution of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_solution(PhysicsID                      physics_id,
               TrilinosWrappers::MPI::Vector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_solutions[physics_id] = solution_vector;
  }

  /**
   * @brief Sets the reference to the solution of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_block_solution(PhysicsID                           physics_id,
                     TrilinosWrappers::MPI::BlockVector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    block_physics_solutions[physics_id] = solution_vector;
  }


private:
  const Parameters::Multiphysics multiphysics_parameters;

  // Data structure to store all physics which were enabled
  std::vector<PhysicsID> active_physics;

  // Auxiliary physics are stored within a map of shared pointer to ensure
  // proper deallocation.
  std::map<
    PhysicsID,
    std::shared_ptr<AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>>>
    physics;

  std::map<
    PhysicsID,
    std::shared_ptr<AuxiliaryPhysics<dim, TrilinosWrappers::MPI::BlockVector>>>
    block_physics;


  std::map<PhysicsID, DoFHandler<dim> *>               physics_dof_handler;
  std::map<PhysicsID, TrilinosWrappers::MPI::Vector *> physics_solutions;
  std::map<PhysicsID, TrilinosWrappers::MPI::BlockVector *>
    block_physics_solutions;
};


#endif
