// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/*
 * This class provides an interface for multiphysics simulations by enabling
 * the solution of multiple auxiliary physics on top of a computational
 * fluid dynamics simulation. The auxiliary physics are stored in a map
 * whose keys are the Parameters::PhysicsID int enum.
 */

#ifndef lethe_multiphysics_interface_h
#define lethe_multiphysics_interface_h

#include <core/exceptions.h>
#include <core/multiphysics.h>
#include <core/parameters_multiphysics.h>
#include <core/simulation_control.h>
#include <core/solid_base.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/distributed/tria_base.h>

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
   * Depending on which multiphysics element is enabled, the appropriate
   * auxiliary physics is instantiated.
   *
   */
  MultiphysicsInterface(
    const SimulationParameters<dim> &nsparam,
    std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                       p_triangulation,
    std::shared_ptr<SimulationControl> p_simulation_control,
    ConditionalOStream                &p_pcout);

  /**
   * @brief Default destructor.
   */
  virtual ~MultiphysicsInterface() = default;

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
        announce_string(pcout, "Heat Transfer");
      }
    else if (physics_id == PhysicsID::tracer)
      {
        announce_string(pcout, "Tracer");
      }
    else if (physics_id == PhysicsID::VOF)
      {
        announce_string(pcout, "VOF");
      }
    else if (physics_id == PhysicsID::cahn_hilliard)
      {
        announce_string(pcout, "Cahn-Hilliard");
      }
  }

  /**
   * @brief Call for the solution of the physics that should be solved
   *
   * @param fluid_dynamics_has_been_solved Boolean that states if the fluid dynamics has been
   * already solved or not. See the map `solve_pre_fluid` to know which
   * subphysics are solved before the fluid dynamics.
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  void
  solve(const bool fluid_dynamics_has_been_solved,
        const Parameters::SimulationControl::TimeSteppingMethod
          time_stepping_method)
  {
    // Loop through all the elements in the physics map. Consequently, iphys is
    // an std::pair where iphys.first is the PhysicsID and iphys.second is the
    // AuxiliaryPhysics pointer. This is how the map can be traversed
    // sequentially.
    for (auto &iphys : physics)
      {
        // If iphys.first should be solved BEFORE fluid dynamics
        if (!fluid_dynamics_has_been_solved && solve_pre_fluid[iphys.first])
          solve_physics(iphys.first, time_stepping_method);

        // If iphys.first should be solved AFTER fluid dynamics OR if is not
        // present in solve_pre_fluid map
        else if (fluid_dynamics_has_been_solved &&
                 (!solve_pre_fluid[iphys.first] ||
                  solve_pre_fluid.count(iphys.first) == 0))
          solve_physics(iphys.first, time_stepping_method);
      }

    for (auto &iphys : block_physics)
      {
        // If iphys.first should be solved BEFORE fluid dynamics
        if (!fluid_dynamics_has_been_solved && solve_pre_fluid[iphys.first])
          solve_block_physics(iphys.first, time_stepping_method);

        // If iphys.first should be solved AFTER fluid dynamics OR if is not
        // present in solve_pre_fluid map
        else if (fluid_dynamics_has_been_solved &&
                 (!solve_pre_fluid[iphys.first] ||
                  solve_pre_fluid.count(iphys.first) == 0))
          solve_block_physics(iphys.first, time_stepping_method);
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
    // Announce physic solved (verbosity = non_linear_solver.verbosity)
    if (verbosity.at(physics_id) != Parameters::Verbosity::quiet)
      announce_physics(physics_id);

    AssertThrow(std::find(active_physics.begin(),
                          active_physics.end(),
                          physics_id) != active_physics.end(),
                ExcInternalError());

    physics[physics_id]->time_stepping_method = time_stepping_method;
    physics[physics_id]->solve_non_linear_system();
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
    // Announce physic solved (verbosity = non_linear_solver.verbosity)
    if (verbosity.at(physics_id) != Parameters::Verbosity::quiet)
      announce_physics(physics_id);

    AssertThrow(std::find(active_physics.begin(),
                          active_physics.end(),
                          physics_id) != active_physics.end(),
                ExcInternalError());

    block_physics[physics_id]->time_stepping_method = time_stepping_method;
    block_physics[physics_id]->solve_non_linear_system();
    block_physics[physics_id]->modify_solution();
  }

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files. This is a variant for GlobalVectorType.
   */
  std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook_global_vector()
  {
    std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;
    for (auto &iphys : physics)
      {
        std::vector<OutputStruct<dim, GlobalVectorType>> output_structs =
          iphys.second->gather_output_hook();
        for (auto &output_struct : output_structs)
          solution_output_structs.push_back(output_struct);
      }

    return solution_output_structs;
  }

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files. This is a variant for GlobalBlockVectorType.
   */
  std::vector<OutputStruct<dim, GlobalBlockVectorType>>
  gather_output_hook_global_block_vector()
  {
    std::vector<OutputStruct<dim, GlobalBlockVectorType>>
      solution_output_structs;
    for (auto &iphys : block_physics)
      {
        std::vector<OutputStruct<dim, GlobalBlockVectorType>> output_structs =
          iphys.second->gather_output_hook();
        for (auto &output_struct : output_structs)
          solution_output_structs.push_back(output_struct);
      }

    return solution_output_structs;
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
   * @brief Rearrange vector solution correctly for transient simulations for
   * all auxiliary physics.
   *
   * @param fluid_dynamics_has_been_solved Boolean that states if the fluid dynamics has been
   * already solved or not. See the map `solve_pre_fluid` to know which
   * subphysics are solved before the fluid dynamics.
   */
  void
  percolate_time_vectors(const bool fluid_dynamics_has_been_solved)
  {
    // Loop through all the elements in the physics map. Consequently, iphys is
    // an std::pair where iphys.first is the PhysicsID and iphys.second is the
    // AuxiliaryPhysics pointer. This is how the map can be traversed
    // sequentially.
    for (auto &iphys : physics)
      {
        // If iphys.first should be percolated BEFORE fluid dynamics is solved
        if (!fluid_dynamics_has_been_solved && solve_pre_fluid[iphys.first])
          iphys.second->percolate_time_vectors();

        // If iphys.first should be percolated AFTER fluid dynamics is solved OR
        // if is not present in solve_pre_fluid map
        else if (fluid_dynamics_has_been_solved &&
                 (!solve_pre_fluid[iphys.first] ||
                  solve_pre_fluid.count(iphys.first) == 0))
          iphys.second->percolate_time_vectors();
      }
    for (auto &iphys : block_physics)
      {
        // If iphys.first should be percolated BEFORE fluid dynamics is solved
        if (!fluid_dynamics_has_been_solved && solve_pre_fluid[iphys.first])
          iphys.second->percolate_time_vectors();

        // If iphys.first should be percolated AFTER fluid dynamics is solved OR
        // if is not present in solve_pre_fluid map
        else if (fluid_dynamics_has_been_solved &&
                 (!solve_pre_fluid[iphys.first] ||
                  solve_pre_fluid.count(iphys.first) == 0))
          iphys.second->percolate_time_vectors();
      }
  }

  /**
   * @param Update the boundary conditions of the auxiliary physics if they are time-dependent
   */
  void
  update_boundary_conditions()
  {
    for (auto &iphys : physics)
      {
        iphys.second->update_boundary_conditions();
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->update_boundary_conditions();
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
   * only support imposing nodal values, but some physics additionally
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
   */
  void
  solve_linear_system(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());

    physics[physics_id]->solve_linear_system();
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
   *
   * @return Reference to the DOF handler of the requested physics
   */
  const DoFHandler<dim> &
  get_dof_handler(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());

    return *physics_dof_handler[physics_id];
  }

  /**
   * @brief Request the reference to the present solution of a given physics
   *
   * @param physics_id The physics ID of the solution being requested
   *
   * @return Reference to the solution vector of the requested physics
   */
  const GlobalVectorType &
  get_solution(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *physics_solutions[physics_id];
  }

  /**
   * @brief Request the reference to the present filtered solution of a given
   * physics (used in VOF or CahnHilliard physics for STF calculation in the
   * momentum balance)
   *
   * @param[in] physics_id ID of the physics for which the filtered solution is
   * being requested
   *
   * @return Reference to the requested filtered solution vector.
   */
  const GlobalVectorType &
  get_filtered_solution(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *physics_filtered_solutions[physics_id];
  }

  /**
   * @brief Request the reference to the present block solution of a given
   * physics
   *
   * @param[in] physics_id The physics ID of the solution being requested
   */
  const GlobalBlockVectorType &
  get_block_solution(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *block_physics_solutions[physics_id];
  }


  /**
   * @brief Request the reference to the time-average solution of a given physics
   *
   * @param[in] physics_id The physics ID of the solution being requested
   */
  const GlobalVectorType &
  get_time_average_solution(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *physics_time_average_solutions[physics_id];
  }

  /**
   * @brief Request the reference to the present block average solution of a
   * given physics
   *
   * @param[in] physics_id The physics ID of the solution being requested
   */
  const GlobalBlockVectorType &
  get_block_time_average_solution(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *block_physics_time_average_solutions[physics_id];
  }

  /**
  * @brief Request the solid objects. Used an auxiliary physics
  * needs to apply a boundary condition on a solid through
  * Nitsche immersed boundary method.
  *
  * @param[in] number_solids The number of solids declared in the parameter file.
  * The value is used to ensure that at least one solid has been declared.
  *
  * @note The method is called only in
  * HeatTransfer<dim>::assemble_nitsche_heat_restriction,
  * which is itself called only if number_solids > 0
  */
  const std::vector<std::shared_ptr<SolidBase<dim, dim>>> &
  get_solids([[maybe_unused]]const int number_solids)
  {
    Assert(number_solids > 0, NoSolidWarning("the"));
    AssertThrow(solids != nullptr,
                dealii::ExcMessage("solids is not initialized"));
    return *solids;
  }

  /**
   * @brief Request the reference to the present solution vector of the
   * projected phase fraction gradient (PFG)
   */
  const GlobalVectorType &
  get_projected_phase_fraction_gradient_solution();

  /**
   * @brief Request the reference to the present solution of the curvature
   */
  const GlobalVectorType &
  get_curvature_solution();

  /**
   * @brief Request the reference to the projected curvature DOF handler
   */
  const DoFHandler<dim> &
  get_curvature_dof_handler();

  /**
   * @brief Request the reference to the projected phase fraction gradient (PFG)
   * DOF handler
   */
  const DoFHandler<dim> &
  get_projected_phase_fraction_gradient_dof_handler();

  /**
   * @brief Request shared pointer to immersed solid shape
   */
  std::shared_ptr<Shape<dim>>
  get_immersed_solid_shape();

  /**
   * @brief Share immersed solid shape
   *
   * @param[in] shape The reference to the shared pointer pointing to the
   * immersed solid shape
   */
  void
  set_immersed_solid_shape(const std::shared_ptr<Shape<dim>> &shape);

  /**
   * @brief Request the reference to the vector of previous solutions of a given
   * physics
   *
   * @param[in] physics_id The physics ID of the solution being requested
   */
  const std::vector<GlobalVectorType> &
  get_previous_solutions(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *physics_previous_solutions[physics_id];
  }


  /**
   * @brief Request the reference to the vector of previous solutions of a given
   * block physics
   *
   * @param[in] physics_id The physics ID of the solution being requested
   */
  const std::vector<GlobalBlockVectorType> &
  get_block_previous_solutions(const PhysicsID physics_id)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    return *block_physics_previous_solutions[physics_id];
  }

  /**
   * @brief Sets the shared pointer to the DOFHandler of the physics in the
   * multiphysics interface
   *
   * @param[in] physics_id The physics of the DOF handler being requested
   *
   * @param[in] dof_handler Shared pointer to the dof handler for which the
   * reference is stored
   */
  void
  set_dof_handler(const PhysicsID                  physics_id,
                  std::shared_ptr<DoFHandler<dim>> dof_handler)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_dof_handler[physics_id] = dof_handler;
  }

  /**
   * @brief Sets the shared pointer to the vector of the SolidBase object. This
   * allows the use of the solid base object in multiple physics at the same
   * time.
   *
   * @param[in] solids_input Shared pointer to the vector of solidBase object
   */
  void
  set_solid(std::shared_ptr<std::vector<std::shared_ptr<SolidBase<dim, dim>>>>
              solids_input)
  {
    solids = solids_input;
  }

  /**
   * @brief Sets the shared pointer to the solution of the physics in the
   * multiphysics interface
   *
   * @param[in] physics_id The physics ID the present solution being set
   *
   * @param[in] solution_vector Shared pointer to the solution vector of the
   * physics
   */
  void
  set_solution(const PhysicsID                   physics_id,
               std::shared_ptr<GlobalVectorType> solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_solutions[physics_id] = solution_vector;
  }

  /**
   * @brief Sets the shared pointer to the filtered solution of the physics in
   * the multiphysics interface (used in VOF or CahnHilliard physics for STF
   * calculation in the momentum balance)
   *
   * @param[in] physics_id ID of the physics for which the filtered solution is
   * being set
   *
   * @param[in] filtered_solution_vector Shared pointer to the filtered solution
   * vector of the physics; this was implemented for VOF and CahnHilliard
   * physics
   */
  void
  set_filtered_solution(
    const PhysicsID                   physics_id,
    std::shared_ptr<GlobalVectorType> filtered_solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_filtered_solutions[physics_id] = filtered_solution_vector;
  }


  /**
   * @brief Sets the shared pointer to the time-average solution of the physics in the multiphysics interface
   *
   * @param[in] physics_id The physics ID of the time averaged solution being
   * set
   *
   * @param[in] solution_vector The shared pointer to the time averaged solution
   * vector of the physics
   */
  void
  set_time_average_solution(const PhysicsID                   physics_id,
                            std::shared_ptr<GlobalVectorType> solution_vector)
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
   * @param[in] physics_id The physics ID of the DOF handler being requested
   *
   * @param[in] solution_vector The shared pointer to the solution vector of the
   * requested physics
   */
  void
  set_block_solution(const PhysicsID                        physics_id,
                     std::shared_ptr<GlobalBlockVectorType> solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    block_physics_solutions[physics_id] = solution_vector;
  }

  /**
   * @brief Sets the shared pointer to the time-average solution of the block
   * physics in the multiphysics interface
   *
   * @param[in] physics_id The physics ID of the time averaged block vector
   * solution being set
   *
   * @param[in] solution_vector The shared pointer to the block solution vector
   * of the physics
   */
  void
  set_block_time_average_solution(
    const PhysicsID                        physics_id,
    std::shared_ptr<GlobalBlockVectorType> solution_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    block_physics_time_average_solutions[physics_id] = solution_vector;
  }

  /**
   * @brief Sets the pointer to the vector of previous solutions of the physics in the multiphysics interface
   *
   * @param[in] physics_id The physics of the DOF handler
   *
   * @param[in] previous_solutions_vector The shared pointer to the vector of
   * previous solutions
   */
  void
  set_previous_solutions(
    const PhysicsID                                physics_id,
    std::shared_ptr<std::vector<GlobalVectorType>> previous_solutions_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    physics_previous_solutions[physics_id] = previous_solutions_vector;
  }

  /**
   * @brief Sets the pointer to the vector of previous solutions of the block
   * physics in the multiphysics interface
   *
   * @param[in] physics_id The physics of the DOF handler
   *
   * @param[in] previous_solutions_vector The shared pointer to the vector of
   * previous block vector solutions
   */
  void
  set_block_previous_solutions(
    const PhysicsID physics_id,
    std::shared_ptr<std::vector<GlobalBlockVectorType>>
      previous_solutions_vector)
  {
    AssertThrow((std::find(active_physics.begin(),
                           active_physics.end(),
                           physics_id) != active_physics.end()),
                ExcInternalError());
    block_physics_previous_solutions[physics_id] = previous_solutions_vector;
  }

  /**
   * @brief Mesh refinement according to an auxiliary physic parameter
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell
   */
  virtual void
  compute_kelly(const std::pair<const Variable,
                                Parameters::MultipleAdaptationParameters> &ivar,
                dealii::Vector<float> &estimated_error_per_cell)
  {
    for (auto &iphys : physics)
      {
        iphys.second->compute_kelly(ivar, estimated_error_per_cell);
      }
    for (auto &iphys : block_physics)
      {
        iphys.second->compute_kelly(ivar, estimated_error_per_cell);
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
  const Parameters::Multiphysics             multiphysics_parameters;
  std::map<PhysicsID, Parameters::Verbosity> verbosity;
  ConditionalOStream                         pcout;

  // Data structure to store all physics which were enabled
  std::vector<PhysicsID> active_physics;


  // Map that states if the physics are solved before the fluid dynamics
  std::map<PhysicsID, bool> solve_pre_fluid{{fluid_dynamics, false},
                                            {VOF, true},
                                            {heat_transfer, false},
                                            {tracer, false},
                                            {cahn_hilliard, true}};

  // Auxiliary physics are stored within a map of shared pointer to ensure
  // proper memory management.
  std::map<PhysicsID, std::shared_ptr<AuxiliaryPhysics<dim, GlobalVectorType>>>
    physics;

  std::map<PhysicsID,
           std::shared_ptr<AuxiliaryPhysics<dim, GlobalBlockVectorType>>>
    block_physics;

  /// Map of physics and shared pointers to their respective DoFHandler
  std::map<PhysicsID, std::shared_ptr<DoFHandler<dim>>> physics_dof_handler;

  /// Shared pointer to the vector containing shared pointers to solid objects
  std::shared_ptr<std::vector<std::shared_ptr<SolidBase<dim, dim>>>> solids;

  /// Map of physics and shared pointers to their respective solutions.
  std::map<PhysicsID, std::shared_ptr<GlobalVectorType>> physics_solutions;

  /**
   * Map of physics and shared pointers to their respective solutions.
   * Same as MultiphysicsInterface::physics_solutions, but used with
   * BlockVector.
   */
  std::map<PhysicsID, std::shared_ptr<GlobalBlockVectorType>>
    block_physics_solutions;

  /**
   * Map of physics and shared pointers to their respective filtered solutions.
   * These solutions are used with both VOF and Cahn-Hilliard multiphase flow
   * approaches for surface tension force calculation in the momentum equation.
   */
  std::map<PhysicsID, std::shared_ptr<GlobalVectorType>>
    physics_filtered_solutions;

  /**
   * Map of physics and shared pointers to their respective filtered solutions.
   * Same as MultiphysicsInterface::physics_filtered_solutions, but used with
   * BlockVector.
   */
  std::map<PhysicsID, std::shared_ptr<GlobalVectorType>>
    block_physics_filtered_solutions;

  /**
   * Map of physics and shared pointers to their respective vector of previous
   * solutions.
   */
  std::map<PhysicsID, std::shared_ptr<std::vector<GlobalVectorType>>>
    physics_previous_solutions;

  /**
   * Map of physics and shared pointers to their respective vector of previous
   * solutions.
   * Same as MultiphysicsInterface::physics_previous_solutions, but used with
   * BlockVector.
   */
  std::map<PhysicsID, std::shared_ptr<std::vector<GlobalBlockVectorType>>>
    block_physics_previous_solutions;

  /**
   * Map of physics and shared pointers to their respective time-averaged
   * solutions.
   */
  std::map<PhysicsID, std::shared_ptr<GlobalVectorType>>
    physics_time_average_solutions;

  /**
   * Map of physics and shared pointers to their respective time-averaged
   * solutions.
   * Same as MultiphysicsInterface::physics_time_average_solutions, but used
   * with BlockVector.
   */
  std::map<PhysicsID, std::shared_ptr<GlobalBlockVectorType>>
    block_physics_time_average_solutions;

  /// Shared pointer to immersed solid shapes to be used by auxiliary physics
  std::shared_ptr<Shape<dim>> immersed_solid_shape;

  /**
   * Map of physics and shared pointers to their respective previous (n-1)
   * solutions.
   */
  std::map<PhysicsID, std::shared_ptr<GlobalVectorType>> physics_solutions_m1;
  std::map<PhysicsID, std::shared_ptr<GlobalBlockVectorType>>
    block_physics_solutions_m1;

  // Checks the required dependencies between multiphase models and handles the
  // corresponding assertions
  void
  inspect_multiphysics_models_dependencies(
    const SimulationParameters<dim> &nsparam);
};


#endif
