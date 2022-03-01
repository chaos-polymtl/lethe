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
 */

#ifndef lethe_auxiliary_physics_h
#define lethe_auxiliary_physics_h

#include <core/parameters.h>
#include <core/physics_solver.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>


/**
 * Au auxiliary physics is defined as a physics that is solved on top of
 * a core physics, the latter being the Navier-Stokes equations.
 * Auxiliary physics are simpler physics around which the entire simulation
 * process does not need to be tailored. Examples of auxiliary physics are
 * temperature and concentration. Generally, the auxiliary physics does not
 * affect the core physics. For example, temperature is advected by the fluid
 * flow, but buoyancy effects are not taken into account.
 *
 * The auxiliary physics are managed by the multiphysics interface of Lethe.
 *
 * This base class is used to establish all of the routines that an auxiliary
 * physics must be able to provide to the Multiphysics interface. These elements
 * are then called at specific moments of a simulation.
 *
 * Auxiliary physics are templated by the dimension of the problem and the
 * vector type that is used to manage their data.
 *
 *
 * Current limitatitations:
 *
 * - Auxiliary physics are currently expected to be used in conjuction with a
 * flow solver which can provide a velocity field.
 * - Support for feedback from the auxiliary physics to the core physics is
 * there but presently not used anywhere
 * - Support for interaction between auxiliary physics is supported but has not
 * been tested
 **/
template <int dim, typename VectorType>
class AuxiliaryPhysics : public PhysicsSolver<VectorType>
{
public:
  /**
   * @brief AuxiliaryPhysics - Base constructor for Auxiliary physics. At the present
   * moment this is an interface with nothing. The auxiliary physics is a pure
   * virtual class.
   */
  AuxiliaryPhysics(
    const Parameters::NonLinearSolver non_linear_solver_parameters)
    : PhysicsSolver<VectorType>(non_linear_solver_parameters)
  {}

  /**
   * @brief AuxiliaryPhysics - Base destructor. At the present
   * moment this is an interface with nothing.
   */
  virtual ~AuxiliaryPhysics()
  {}

  /**
   * @brief Attach the solution vector to the DataOut provided. This function
   * enable the auxiliary physics to output their solution via the core solver.
   */
  virtual void
  attach_solution_to_output(DataOut<dim> &data_out) = 0;

  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  virtual void
  finish_simulation() = 0;

  /**
   * @brief Carry out the operations required to finish a time step correctly.
   */
  virtual void
  finish_time_step() = 0;

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  virtual void
  percolate_time_vectors() = 0;

  /**
   * @brief Carry out modifications on the auxiliary physic solution.
   * To be defined for some physics only (eg. free surface, see vof.h).
   */
  virtual void
  modify_solution(){};

  /**
   * @brief
   */
  virtual TrilinosWrappers::MPI::Vector *
  get_pfg_solution(){};

  /**
   * @brief
   */
  virtual TrilinosWrappers::MPI::Vector *
  get_curvature_solution(){};

  /**
   * @brief
   */
  virtual DoFHandler<dim> *
  get_pfg_dof_handler(){};

  /**
   * @brief
   */
  virtual DoFHandler<dim> *
  get_curvature_dof_handler(){};

  /**
   * @brief Provide the dof handler associated with an auxiliary physics
   * TODO : delete as the auxiliary physics are supposed to pass their
   * dof_handler to the multiphysics_interface directly
   */
  virtual const DoFHandler<dim> &
  get_dof_handler() = 0;

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess(bool first_iteration) = 0;


  /**
   * @brief pre_mesh_adaption Prepares the auxiliary physics variables for a
   * mesh refinement/coarsening
   */
  virtual void
  pre_mesh_adaptation() = 0;

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  virtual void
  post_mesh_adaptation() = 0;

  /**
   * @brief Prepares auxiliary physics to write checkpoint
   */
  virtual void
  write_checkpoint() = 0;

  /**
   * @brief Set solution vector of Auxiliary Physics using checkpoint
   */
  virtual void
  read_checkpoint() = 0;

  /**
   * @brief Compute the Kelly error estimator used to refine mesh on a auxiliary physic parameter.
   */
  virtual void
  compute_kelly(dealii::Vector<float> &estimated_error_per_cell) = 0;

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
};



#endif
