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
  GLSNitscheNavierStokesSolver(NavierStokesSolverParameters<spacedim> &nsparam);

  virtual void
  solve() override;


protected:
  virtual void
  setup_dofs() override;


  virtual void
  set_initial_condition(Parameters::InitialConditionType initial_condition_type,
                        bool restart = false) override;

  virtual void
  solve_linear_system(
    const bool initial_step,
    const bool renewed_matrix = true) override; // Interface function


private:
  /**
   * @brief Temporary - Adds heat transfer solving - test in subfunction before global multiphysics implementation
   */
  virtual void
  setup_dofs_ht();

  virtual void
  assemble_matrix_and_rhs_ht(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method);

  virtual void
  solve_linear_system_ht();

  virtual void
  finish_time_step_ht();

  virtual void
  set_initial_condition_ht();

  void
  postprocess_ht();

  void
  finish_ht();

  double
  calculate_l2_error_ht(const DoFHandler<spacedim> &         dof_handler,
                        const TrilinosWrappers::MPI::Vector &evaluation_point,
                        const Function<spacedim> &           exact_solution,
                        const Parameters::FEM &              fem_parameters,
                        const MPI_Comm &                     mpi_communicator);

  ConvergenceTable error_table_ht;


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
   * @brief Post-process for forces on solid after an iteration
   */
  void
  postprocess_solid_forces();

  /**
   * @brief Same has in gls_navier_stokes, but calls assemble_nitsche_restriction() when global matrix and rhs are assembled
   */
  virtual void
  assemble_matrix_and_rhs_fd(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  /**
   * @brief Same has in gls_navier_stokes, but calls assemble_nitsche_restriction() when rhs is assembled
   */
  virtual void
  assemble_rhs_fd(const Parameters::SimulationControl::TimeSteppingMethod
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

  /**
   * @brief a function for adding data vectors to the data_out object for
   * post_processing additional results
   */
  virtual void
  output_field_hook(DataOut<spacedim> &data_out) override;


  SolidBase<dim, spacedim> solid;
  PVDHandler               pvdhandler_solid_triangulation;
  PVDHandler               pvdhandler_solid_particles;


  TableHandler solid_forces_table;

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
