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

#include <deal.II/lac/trilinos_vector.h>

#include <core/solid_base.h>

#include "gls_navier_stokes.h"

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
  calculate_forces_on_solid();

  /**
   * @brief Calculates the torque due to the fluid motion on the solid
   * @return Tensor of torque on the solid
   */
  Tensor<1, 3>
  calculate_torque_on_solid();

  /**
   * @brief Post-process for forces on solid after an iteration
   */
  void
  postprocess_solid_forces();


  /**
   * @brief Post-process for torques on solid after an iteration
   */
  void
  postprocess_solid_torques();

  /**
   * @brief Same has in gls_navier_stokes, but calls assemble_nitsche_restriction() when global matrix and rhs are assembled
   */
  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  /**
   * @brief Same has in gls_navier_stokes, but calls assemble_nitsche_restriction() when rhs is assembled
   */
  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  /**
   * @brief Outputs a vtu file for each output frequency of the particles
   */
  void
  output_solid_particles(
    std::shared_ptr<Particles::ParticleHandler<spacedim>> particle_handler);
  /**
   * @brief Outputs a vtu file for each output frequency of the solid triangulation
   */
  void
  output_solid_triangulation();

  SolidBase<dim, spacedim> solid;
  PVDHandler               pvdhandler_solid_triangulation;
  PVDHandler               pvdhandler_solid_particles;


  TableHandler solid_forces_table;
  TableHandler solid_torques_table;

  /**
   * @brief Temporary - Addition for Heat Transfer solving - test in subfunction before global multiphysics implementation
   */
  IndexSet                       locally_owned_dofs_ht;
  IndexSet                       locally_relevant_dofs_ht;
  FE_Q<spacedim>                 fe_ht;
  DoFHandler<spacedim>           dof_handler_ht;
  TrilinosWrappers::SparseMatrix system_matrix_ht;
  SparsityPattern                sparsity_pattern_ht;
  AffineConstraints<double>      zero_constraints_ht;



  TrilinosWrappers::MPI::Vector solution_ht_m1; // minus 1
  TrilinosWrappers::MPI::Vector solution_ht_m2; // minus 2
  TrilinosWrappers::MPI::Vector solution_ht_m3; // minus 3
};


#endif
