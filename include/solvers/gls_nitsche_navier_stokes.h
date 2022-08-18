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
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique
 Montreal, 2020-
 */

#ifndef lethe_gls_nitsche_navier_stokes_h
#define lethe_gls_nitsche_navier_stokes_h

#include <core/solid_base.h>

#include <solvers/gls_navier_stokes.h>

#include <deal.II/lac/trilinos_vector.h>

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
class GLSNitscheNavierStokesSolver : public GLSNavierStokesSolver<spacedim>
{
public:
  GLSNitscheNavierStokesSolver(SimulationParameters<spacedim> &nsparam);

  virtual void
  solve() override;

private:
  /**
   * @brief Adds the nitsche restriction to the global matrix and global rhs on the cells surrounding the immersed solid
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
   * @return Tensor of torque on the solid. This is always a 3D tensor even in 2D
   */
  Tensor<1, 3>
  calculate_torque_on_solid(const unsigned int i_solid);

  /**
   * @brief Calculates the torque due to the fluid motion on the solid
   * @return Tensor of torque on the solid. This is always a 3D tensor even in 2D
   */
  Tensor<1, 3>
  calculate_torque_on_solid();

  /**
   * @brief Post-process for forces on solid after an iteration
   * @param i_solid is the solid index
   * @param first_solid_forces is a boolean set to true for the first call of this method, used for table output formatting
   */
  void
  postprocess_solid_forces(const unsigned int i_solid, bool first_solid_forces);


  /**
   * @brief Post-process for torques on solid after an iteration
   * @param i_solid is the solid index
   * @param first_solid_torques is a boolean set to true for the first call of this method, used for table output formatting
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
   * @brief Same as in gls_navier_stokes, but calls assemble_nitsche_restriction() when global matrix and rhs are assembled
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver
   */
  void
  assemble_matrix_and_rhs();

  /**
   * @brief Same as in gls_navier_stokes, but calls assemble_nitsche_restriction() when rhs is assembled
   *
   * @deprecated This function is to be deprecated when the new assembly mechanism
   * is integrated to this solver
   */
  void
  assemble_rhs();

  /**
   * @brief Outputs a vtu file for each output frequency of the particles
   */
  void
  output_solid_particles(const unsigned int i_solid);
  /**
   * @brief Outputs a vtu file for each output frequency of the solid triangulation
   */
  void
  output_solid_triangulation(const unsigned int i_solid);

  /**
   * @brief Allow the refinement of the mesh according to one of the 2 methods proposed.
   * Overrides the regular mesh adaptation by ensuring that the Nitsche solids
   * are prepared for adaptation.
   */
  virtual void
  refine_mesh() override;

  /**
   * @brief Write a gls_nitsche simulation checkpointing to allow for gls_nitsche simulation restart
   */
  virtual void
  write_checkpoint() override;

  /**
   * @brief Read a gls_nitsche simulation checkpoint and initiate simulation restart
   */
  virtual void
  read_checkpoint() override;

  std::vector<std::shared_ptr<SolidBase<dim, spacedim>>> solids;
  std::vector<PVDHandler> pvdhandler_solid_triangulation;
  std::vector<PVDHandler> pvdhandler_solid_particles;

  std::vector<TableHandler> solid_forces_table;
  std::vector<TableHandler> solid_torques_table;
};


#endif
