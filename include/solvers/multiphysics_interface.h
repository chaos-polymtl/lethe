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
 * whose keys are the Parameters::PhysicsID int enum.
 */

#ifndef lethe_multiphysics_interface_h
#define lethe_multiphysics_interface_h

#include <core/multiphysics.h>
#include <core/parameters_multiphysics.h>
#include <core/simulation_control.h>
#include <core/solid_base.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

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
    std::shared_ptr<SimulationControl> p_simulation_control,
    ConditionalOStream &               p_pcout);

  std::vector<PhysicsID>
  get_active_physics()
  {
    return active_physics;
  }

  /**
   * @brief Write physic solved in the terminal
   *
   * @param physics_id number associated with auxiliary physics in multiphysics.h
   */
  void
  announce_physics(const PhysicsID physics_id)
  {
    if (physics_id == PhysicsID::heat_transfer)
      {
        pcout << "Heat Transfer" << std::endl << "--------------" << std::endl;
      }
    else if (physics_id == PhysicsID::tracer)
      {
        pcout << "Tracer" << std::endl << "-------" << std::endl;
      }
    else if (physics_id == PhysicsID::VOF)
      {
        pcout << "VOF" << std::endl << "-------------" << std::endl;
      }
  }

  /**
   * @brief Call for the solution of the physics that should be solved BEFORE the fluid dynamics (called "pre-solving")
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  pre_solve(const Parameters::SimulationControl::TimeSteppingMethod
              time_stepping_method)
  {
    // Loop through all the elements in the physics map. Consequently, iphys is
    // an std::pair where iphys.first is the PhysicsID and iphys.second is the
    // AuxiliaryPhysics pointer. This is how the map can be traversed
    // sequentially.
    for (auto &iphys : physics)
      {
        // Ensure that iphys.first is present in solve_pre_fluid map
        if (solve_pre_fluid.count(iphys.first) != 0)
          {
            if (solve_pre_fluid[iphys.first] == true)
              {
                // Announce physic solved (verbosity =
                // non_linear_solver.verbosity)
                if (verbosity != Parameters::Verbosity::quiet)
                  announce_physics(iphys.first);

                solve_physics(iphys.first, time_stepping_method);

                // Embellish the console output
                if (verbosity != Parameters::Verbosity::quiet)
                  pcout << "-----------------------" << std::endl;
              }
          }
        else
          {
            pcout
              << "Solution order is not specified for this auxiliary physic "
              << std::endl
              << "will be solved after the fluid dynamics by default."
              << std::endl
              << "For clarity, please specify a solution order for this physic "
              << std::endl
              << "in the map solve_pre_fluid, member of multiphysics_interface.h"
              << std::endl
              << "-------------" << std::endl;
          }
      }

    for (auto &iphys : block_physics)
      {
        // Ensure that iphys.first is present in solve_pre_fluid map
        if (solve_pre_fluid.count(iphys.first) != 0)
          {
            if (solve_pre_fluid[iphys.first] == true)
              {
                // Announce physic solved (verbosity =
                // non_linear_solver.verbosity)
                if (verbosity != Parameters::Verbosity::quiet)
                  announce_physics(iphys.first);

                solve_block_physics(iphys.first, time_stepping_method);

                // Embellish the console output
                if (verbosity != Parameters::Verbosity::quiet)
                  pcout << "-----------------------" << std::endl;
              }
          }
        else
          {
            pcout
              << "Solution order is not specified for this auxiliary physic "
              << std::endl
              << "will be solved after the fluid dynamics by default."
              << std::endl
              << "For clarity, please specify a solution order for this physic "
              << std::endl
              << "in the map solve_pre_fluid, member of multiphysics_interface.h"
              << std::endl
              << "-------------" << std::endl;
          }
      }
  }

  /**
   * @brief Call for the solution of the physics that should be solved AFTER the fluid dynamics (called "post-solving")
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  post_solve(const Parameters::SimulationControl::TimeSteppingMethod
               time_stepping_method)
  {
    // Loop through all the elements in the physics map. Consequently, iphys is
    // an std::pair where iphys.first is the PhysicsID and iphys.second is the
    // AuxiliaryPhysics pointer. This is how the map can be traversed
    // sequentially.
    for (auto &iphys : physics)
      {
        // If iphys.first should be solved after fluid dynamics OR if is not
        // present in solve_pre_fluid map
        if (solve_pre_fluid[iphys.first] == false ||
            solve_pre_fluid.count(iphys.first) == 0)
          {
            // Announce physic solved (verbosity =
            // non_linear_solver.verbosity)
            if (verbosity != Parameters::Verbosity::quiet)
              {
                // Embellish the console output
                pcout << "-----------------------" << std::endl;
                announce_physics(iphys.first);
              }

            solve_physics(iphys.first, time_stepping_method);
          }
      }

    for (auto &iphys : block_physics)
      {
        // If iphys.first should be solved after fluid dynamics OR if is not
        // present in solve_pre_fluid map
        if (solve_pre_fluid[iphys.first] == false ||
            solve_pre_fluid.count(iphys.first) == 0)
          {
            // Announce physic solved (verbosity =
            // non_linear_solver.verbosity)
            if (verbosity != Parameters::Verbosity::quiet)
              {
                // Embellish the console output
                pcout << "-----------------------" << std::endl;
                announce_physics(iphys.first);
              }

            solve_block_physics(iphys.first, time_stepping_method);
          }
      }
  }

  /**
   * @brief Call for the solution of a single physic
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  solve_physics(const PhysicsID physics_id,
                const Parameters::SimulationControl::TimeSteppingMethod
                  time_stepping_method)
  {
    AssertThrow(std::find(active_physics.begin(),
                          active_physics.end(),
                          physics_id) != active_physics.end(),
                ExcInternalError());

    physics[physics_id]->time_stepping_method = time_stepping_method;
    physics[physics_id]->solve_non_linear_system(false);
    physics[physics_id]->modify_solution();
  }

  /**
   * @brief Call for the solution of a single block physic
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  solve_block_physics(const PhysicsID physics_id,
                      const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method)
  {
    AssertThrow(std::find(active_physics.begin(),
                          active_physics.end(),
                          physics_id) != active_physics.end(),
                ExcInternalError());

    block_physics[physics_id]->time_stepping_method = time_stepping_method;
    block_physics[physics_id]->solve_non_linear_system(false);
    block_physics[physics_id]->modify_solution();
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
   * @brief Carry out the operations required to finish a time step correctly for
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
   * @brief Carry out the operations required to percolate the time vectors
   * correctly at the end of a simulation
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
   * the calculation of all derived quantities using the solution vector
   * of the physics. It does not concern the output of the solution using
   * the DataOutObject, which is accomplished through the
   * attach_solution_to_output function
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
   * @brief Prepare the auxiliary physics for mesh adaptation
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
   * @brief Set-up the DofHandler and the degree of freedom associated with the physics.
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
   * only support imposing nodal values, but some physics additionnaly
   * support the use of L2 projection or steady-state solutions.
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
   * @brief fluid_dynamics_is_block Verify if the fluid dynamics solution
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
   * @brief Request the time-average solution of a given physics
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::Vector *
  get_time_average_solution(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return physics_time_average_solutions[physics_id];
  }

  /**
   * @brief Request the present block average solution of a given physics
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::BlockVector *
  get_block_time_average_solution(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return block_physics_time_average_solutions[physics_id];
  }

  /**
   * @brief Request the reynolds_stress solution of a given physics.
   * WIP for an upcoming PR, not yet implemented in the solver.
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::Vector *
  get_reynolds_stress_solution(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return reynolds_stress_solutions;
  }

  /**
   @brief Request the solid objects. Used an auxiliary physics
    * needs to apply a boundary condition on a solid through
    * Nitsche immersed boundary method.
    *
    * NB: this method is called only in
    * HeatTransfer<dim>::assemble_nitsche_heat_restriction,
    * which is itself called only if number_solids > 0
    */
  std::vector<std::shared_ptr<SolidBase<dim, dim>>> *
  get_solids()
  {
    // Called only if number_solids > 0, in the method
    // HeatTransfer<dim>::assemble_matrix_and_rhs()
    return solids;
  }


  /**
   * @brief Request the present solution of the filtered phase fraction gradient (PFG)
   */
  TrilinosWrappers::MPI::Vector *
  get_filtered_phase_fraction_gradient_solution();

  /**
   * @brief Request the present solution of the curvature
   */
  TrilinosWrappers::MPI::Vector *
  get_curvature_solution();

  /**
   * @brief Request the filtered curvature DOF handler
   */
  DoFHandler<dim> *
  get_curvature_dof_handler();

  /**
   * @brief Request the filtered phase fraction gradient (pfg) DOF handler
   */
  DoFHandler<dim> *
  get_filtered_phase_fraction_gradient_dof_handler();

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
   * @brief Request the past (minus 1) solution of a given physics
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::Vector *
  get_solution_m1(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return physics_solutions_m1[physics_id];
  }


  /**
   * @brief Request the past (minus 1) block solution of a given physics
   *
   * @param physics_id The physics of the solution being requested
   */
  TrilinosWrappers::MPI::BlockVector *
  get_block_solution_m1(PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return block_physics_solutions_m1[physics_id];
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
   * @brief Sets the reference to the vector of the SolidBase object. This allows the use of the solid base object in multiple physics at the same time.
   *
   * @param solids_input The reference to the vector of solidBase object
   */
  void
  set_solid(std::vector<std::shared_ptr<SolidBase<dim, dim>>> *solids_input)
  {
    solids = solids_input;
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
   * @brief Sets the reference to the time-average solution of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_time_average_solution(PhysicsID                      physics_id,
                            TrilinosWrappers::MPI::Vector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_time_average_solutions[physics_id] = solution_vector;
  }


  /**
   * @brief Sets the reference to the Reynolds stress of the physics in the multiphysics interface.
   * WIP for an upcoming PR, not yet implemented in the solver.
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_reynolds_stress_solutions(PhysicsID                      physics_id,
                                TrilinosWrappers::MPI::Vector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_time_average_solutions[physics_id] = solution_vector;
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

  /**
   * @brief Sets the reference to the time-average solution of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_block_time_average_solution(
    PhysicsID                           physics_id,
    TrilinosWrappers::MPI::BlockVector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    block_physics_time_average_solutions[physics_id] = solution_vector;
  }

  /**
   * @brief Sets the reference to the solution of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_solution_m1(PhysicsID                      physics_id,
                  TrilinosWrappers::MPI::Vector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_solutions_m1[physics_id] = solution_vector;
  }

  /**
   * @brief Sets the reference to the solution of the physics in the multiphysics interface
   *
   * @param physics_id The physics of the DOF handler being requested
   *
   * @param solution_vector The reference to the solution vector of the physics
   */
  void
  set_block_solution_m1(PhysicsID                           physics_id,
                        TrilinosWrappers::MPI::BlockVector *solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    block_physics_solutions_m1[physics_id] = solution_vector;
  }

  /**
   * @brief Mesh refinement according to an auxiliary physic parameter
   * NB : only implemented for VOF for now
   */
  virtual void
  compute_kelly(dealii::Vector<float> &estimated_error_per_cell)
  {
    for (auto &iphys : physics)
      {
        iphys.second->compute_kelly(estimated_error_per_cell);
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->compute_kelly(estimated_error_per_cell);
      }
  };

  /**
   * @brief Prepares auxiliary physics to write simulation checkpoint
   */
  virtual void
  write_checkpoint()
  {
    for (auto &iphys : physics)
      {
        iphys.second->write_checkpoint();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->write_checkpoint();
      }
  };

  /**
   * @brief Read solution from checkpoint from auxiliary physics
   *
   */
  virtual void
  read_checkpoint()
  {
    for (auto &iphys : physics)
      {
        iphys.second->read_checkpoint();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->read_checkpoint();
      }
  };

private:
  const Parameters::Multiphysics multiphysics_parameters;
  const Parameters::Verbosity    verbosity;
  ConditionalOStream             pcout;

  // Data structure to store all physics which were enabled
  std::vector<PhysicsID> active_physics;


  // Map that states if the physics are solved before the fluid dynamics
  std::map<PhysicsID, bool> solve_pre_fluid{{fluid_dynamics, false},
                                            {VOF, false},
                                            {heat_transfer, false},
                                            {tracer, false}};

  // Auxiliary physics are stored within a map of shared pointer to ensure
  // proper memory management.
  std::map<
    PhysicsID,
    std::shared_ptr<AuxiliaryPhysics<dim, TrilinosWrappers::MPI::Vector>>>
    physics;

  std::map<
    PhysicsID,
    std::shared_ptr<AuxiliaryPhysics<dim, TrilinosWrappers::MPI::BlockVector>>>
    block_physics;


  std::map<PhysicsID, DoFHandler<dim> *> physics_dof_handler;

  std::vector<std::shared_ptr<SolidBase<dim, dim>>> *solids;



  // present solution
  std::map<PhysicsID, TrilinosWrappers::MPI::Vector *> physics_solutions;

  // average solution
  std::map<PhysicsID, TrilinosWrappers::MPI::Vector *>
    physics_time_average_solutions;

  // average solution
  std::map<PhysicsID, TrilinosWrappers::MPI::BlockVector *>
    block_physics_time_average_solutions;

  // reynolds stress solution. This is WIP and is not yet implemented in the
  // solver.
  TrilinosWrappers::MPI::Vector *reynolds_stress_solutions;

  std::map<PhysicsID, TrilinosWrappers::MPI::BlockVector *>
    block_physics_solutions;

  // past (minus 1) solution
  std::map<PhysicsID, TrilinosWrappers::MPI::Vector *> physics_solutions_m1;
  std::map<PhysicsID, TrilinosWrappers::MPI::BlockVector *>
    block_physics_solutions_m1;

  // Checks the required dependencies between multiphase models and handles the
  // corresponding assertions
  void
  inspect_multiphysics_models_dependencies(
    const SimulationParameters<dim> &nsparam);
};


#endif
