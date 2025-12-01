// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_h
#define lethe_time_harmonic_maxwell_h

#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/physics_solver.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/numerics/data_out.h>

/**
 * @class TimeHamonicMaxwell Auxiliary physics solver for the time harmonic
 * form of the Maxwell equations.
 * TODO Add more details here.
 *
 **/

/// TODO in this PR:
/// [] Remplir chacune des fonctions
///   [] Constructeur
///   [] Destructeur
///   [] Setup_dofs
///   [] Gather output hook
///   [] Finish simulation
///   [] Percolate_time_vector
///   - modify_solution
///   [] update_boundary_conditions
///   [] get_dof_handler
///   [] postprocess
///   [] pre_mesh_adaptation
///   [] post_mesh_adaptation
///   - write_checkpoint
///   - read_checkpoint
///   - gather_tables()
///   - compute_kelly
///   - set_initial_conditions
///   [] setup_preconditioner
/// [] Physics field
/// [] Multiphysics interface components
/// - Physical Properties
///   - Electrical conductivity
///   - Permetivity
///   - Permeability
/// [] FEM section for DPG
/// - Mesh adaptation

template <int dim, typename VectorType>
class TimeHarmonicMaxwell : public AuxiliaryPhysics<dim, VectorType>
{
public:
  /**
   * @brief TimeHarmonicMaxwell - Base constructor for Auxiliary physics. At the present
   * moment this is an interface with nothing. The auxiliary physics is a pure
   * virtual class.
   */
  TimeHarmonicMaxwell(
    const Parameters::NonLinearSolver non_linear_solver_parameters)
    : AuxiliaryPhysics<dim, VectorType>(non_linear_solver_parameters)
  {}

  /**
   * @brief TimeHarmonicMaxwell - Base destructor. At the present
   * moment this is an interface with nothing.
   */
  virtual ~TimeHarmonicMaxwell()
  {}

  /**
   * @brief Gather and return vector of output structs that are particular to some applications.
   *
   * @return Vector of OutputStructs that will be used to write the output results as VTU files.
   */
  virtual std::vector<OutputStruct<dim, VectorType>>
  gather_output_hook()
  {
    return std::vector<OutputStruct<dim, VectorType>>();
  };

  /**
   * @brief Carry out the operations required to finish a simulation correctly.
   */
  virtual void
  finish_simulation()
  {}

  /**
   * @brief Carry out the operations required to rearrange the values of the
   * previous solution at the end of a time step
   */
  virtual void
  percolate_time_vectors()
  {}

  /**
   * @brief Carry out modifications on the auxiliary physic solution.
   * To be defined for some physics only (eg. free surface, see vof.h).
   */
  virtual void
  modify_solution() {};

  /**
   * @brief Update non zero constraints if the boundary is time-dependent
   */
  virtual void
  update_boundary_conditions() {};

  /**
   * @brief Provide the dof handler associated with an auxiliary physics
   * TODO : delete as the auxiliary physics are supposed to pass their
   * dof_handler to the multiphysics_interface directly
   */
  virtual const DoFHandler<dim> &
  get_dof_handler()
  {}

  /**
   * @brief Postprocess the auxiliary physics results. Post-processing this case implies
   * the calculation of all derived quantities using the solution vector of the
   * physics. It does not concern the output of the solution using the
   * DataOutObject, which is accomplished through the attach_solution_to_output
   * function
   */
  virtual void
  postprocess(bool first_iteration) {};


  /**
   * @brief pre_mesh_adaption Prepares the auxiliary physics variables for a
   * mesh refinement/coarsening
   */
  virtual void
  pre_mesh_adaptation() {};

  /**
   * @brief post_mesh_adaption Interpolates the auxiliary physics variables to the new mesh
   */
  virtual void
  post_mesh_adaptation() {};

  /**
   * @brief Prepares auxiliary physics to write checkpoint
   */
  virtual void
  write_checkpoint() {};

  /**
   * @brief Set solution vector of Auxiliary Physics using checkpoint
   */
  virtual void
  read_checkpoint() {};

  /**
   * @brief Returns a vector of references to TableHandler objects that needs to
   * be serialized/deserialized for a given TimeHarmonicMaxwell solver.
   *
   * @return Structure containing the TableHandler objects and their corresponding file names.
   */
  virtual std::vector<OutputStructTableHandler>
  gather_tables()
  {
    return std::vector<OutputStructTableHandler>();
  };

  /**
   * @brief Compute the Kelly error estimator used to refine mesh on a auxiliary physic parameter.
   *
   * @param ivar The current element of the map simulation_parameters.mesh_adaptation.variables
   *
   * @param estimated_error_per_cell The deal.II vector of estimated_error_per_cell
   */
  virtual void
  compute_kelly(const std::pair<const Variable,
                                Parameters::MultipleAdaptationParameters> &ivar,
                dealii::Vector<float> &estimated_error_per_cell) {};

  /**
   * @brief Sets-up the DofHandler and the degree of freedom associated with the physics.
   */
  virtual void
  setup_dofs() {};

  /**
   * @brief Sets-up the initial conditions associated with the physics. Generally, physics
   * only support imposing nodal values, but some physics additionally support
   * the use of L2 projection or steady-state solutions.
   */
  virtual void
  set_initial_conditions() {};

  /**
   * @brief Set up preconditioner. Not used for the auxiliary physics but
   * needed for the compilation of the non-linear solver.
   */
  void
  setup_preconditioner() override {};

private:
};



#endif
