// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_nitsche_h
#define lethe_fluid_dynamics_nitsche_h

#include <core/solid_base.h>

#include <solvers/fluid_dynamics_matrix_based.h>

using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization
 * and Nitsche immersed boundary method
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim, int spacedim = dim>
class FluidDynamicsNitsche : public FluidDynamicsMatrixBased<spacedim>
{
public:
  FluidDynamicsNitsche(SimulationParameters<spacedim> &nsparam);

  virtual void
  solve() override;

private:
  /**
   * @brief Adds the nitsche restriction to the global matrix and global rhs on
   * the cells surrounding the immersed solid
   */
  template <bool assemble_matrix>
  void
  assemble_nitsche_restriction();

  /**
   * @brief Calculates the force due to the fluid motion on the solid
   * @return Tensor of forces on the solid
   */
  Tensor<1, spacedim>
  calculate_forces_on_solid(const unsigned int i_solid);

  /**
   * @brief Calculates the torque due to the fluid motion on the solid
   * @return Tensor of torque on the solid. This is always a 3D tensor even in
   * 2D
   */
  Tensor<1, 3>
  calculate_torque_on_solid(const unsigned int i_solid);

  /**
   * @brief Calculates the torque due to the fluid motion on the solid
   * @return Tensor of torque on the solid. This is always a 3D tensor even in
   * 2D
   */
  Tensor<1, 3>
  calculate_torque_on_solid();

  /**
   * @brief Post-process for forces on solid after an iteration
   * @param i_solid is the solid index
   * @param first_solid_forces is a boolean set to true for the first call of
   * this method, used for table output formatting
   */
  void
  postprocess_solid_forces(const unsigned int i_solid, bool first_solid_forces);


  /**
   * @brief Post-process for torques on solid after an iteration
   * @param i_solid is the solid index
   * @param first_solid_torques is a boolean set to true for the first call of
   * this method, used for table output formatting
   */
  void
  postprocess_solid_torques(const unsigned int i_solid,
                            bool               first_solid_torques);

  /**
   * @brief Call for the assembly of the matrix
   */
  virtual void
  assemble_system_matrix() override
  {
    assemble_matrix_and_rhs();
  }

  /**
   * @brief Call for the assembly of the right-hand side
   */
  virtual void
  assemble_system_rhs() override
  {
    assemble_rhs();
  }

  /**
   * @brief Same as in fluid_dynamics_matrix_based, but calls
   * assemble_nitsche_restriction() when global matrix and rhs are assembled
   *
   * @deprecated This function is to be deprecated when the new assembly
   * mechanism is integrated to this solver
   */
  void
  assemble_matrix_and_rhs();

  /**
   * @brief Same as in fluid_dynamics_matrix_based, but calls
   * assemble_nitsche_restriction() when rhs is assembled
   *
   * @deprecated This function is to be deprecated when the new assembly
   * mechanism is integrated to this solver
   */
  void
  assemble_rhs();

  /**
   * @brief Outputs a vtu file for each output frequency of the particles
   */
  void
  output_solid_particles(const unsigned int i_solid);
  /**
   * @brief Outputs a vtu file for each output frequency of the solid
   * triangulation
   */
  void
  output_solid_triangulation(const unsigned int i_solid);

  /**
   * @brief Allow the refinement of the mesh according to one of the 2 methods
   * proposed. Overrides the regular mesh adaptation by ensuring that the
   * Nitsche solids are prepared for adaptation.
   */
  virtual void
  refine_mesh() override;

  /**
   * @brief Returns a vector of references to TableHandler objects that needs
   * to be serialized/deserialized for a given solver.
   *
   * @return Structure containing a vector of references to TableHandler objects
   * that needs to be serialized/deserialized for a given solver, and their
   * corresponding file names.
   */
  std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Write a gls_nitsche simulation checkpointing to allow for
   * gls_nitsche simulation restart
   */
  virtual void
  write_checkpoint() override;

  /**
   * @brief Read a gls_nitsche simulation checkpoint and initiate simulation
   * restart
   */
  virtual void
  read_checkpoint() override;

  std::shared_ptr<std::vector<std::shared_ptr<SolidBase<dim, spacedim>>>>
                          solids;
  std::vector<PVDHandler> pvdhandler_solid_triangulation;
  std::vector<PVDHandler> pvdhandler_solid_particles;

  std::vector<TableHandler> solid_forces_table;
  std::vector<TableHandler> solid_torques_table;
};


#endif
