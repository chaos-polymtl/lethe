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
 * ---------------------------------------------------------------------*/

#ifndef lethe_gls_navier_stokes_h
#define lethe_gls_navier_stokes_h

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_base.h>
#include <solvers/navier_stokes_scratch_data.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>



using namespace dealii;

/**
 * A solver class for the Navier-Stokes equation using GLS stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim>
class GLSNavierStokesSolver
  : public NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>
{
public:
  GLSNavierStokesSolver(SimulationParameters<dim> &nsparam);
  ~GLSNavierStokesSolver();

  virtual void
  solve();

protected:
  virtual void
  setup_dofs_fd();



  virtual void
  set_initial_condition_fd(
    Parameters::InitialConditionType initial_condition_type,
    bool                             restart = false) override;

  void
  set_solution_vector(double value);

protected:
  /*
   *  Assemble the matrix associated with the solver
   */
  void
  assemble_system_matrix();

  /*
   *  Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs();


  /*
   * Assemble the local matrix for a given cell
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data);

  /*
   * Assemble the local rhs for a given cell
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data);

  /*
   * Sets up the vector of assembler functions
   */
  void
  setup_assemblers();


  /*
   * Copy local cell information to global matrix
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data);

  /*
   * Copy local cell rhs information to global rhs
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data);

  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme,
            Parameters::VelocitySource::VelocitySourceType    velocity_source>
  void
  assembleGLS();

  // assembleGLS specialized in the case of two fluids with free surface
  // simulation
  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme,
            Parameters::VelocitySource::VelocitySourceType    velocity_source>
  void
  assembleGLSFreeSurface();

  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true);

private:
  void
  assemble_L2_projection();



  /**
   * GMRES solver with ILU(N) preconditioning
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual,
                     const bool   renewed_matrix);

  /**
   * BiCGStab solver with ILU(N) preconditioning
   */
  void
  solve_system_BiCGStab(const bool   initial_step,
                        const double absolute_residual,
                        const double relative_residual,
                        const bool   renewed_matrix);

  /**
   * AMG preconditioner with ILU smoother and coarsener and GMRES final solver
   */
  void
  solve_system_AMG(const bool   initial_step,
                   const double absolute_residual,
                   const double relative_residual,
                   const bool   renewed_matrix);

  /**
   * Direct solver
   */
  void
  solve_system_direct(const bool   initial_step,
                      const double absolute_residual,
                      const double relative_residual,
                      const bool   renewed_matrix);

  /**
   * TFQMR solver with ILU preconditioner
   */
  void
  solve_system_TFQMR(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual,
                     const bool   renewed_matrix);

  /**
   * Set-up AMG preconditioner
   */
  void
  setup_AMG(const int current_amg_ilu_preconditioner_fill_level);

  /**
   * Set-up ILU preconditioner
   */
  void
  setup_ILU(const int current_ilu_preconditioner_fill_level);


  /**
   * Members
   */
protected:
  TrilinosWrappers::SparseMatrix system_matrix;

private:
  SparsityPattern                                    sparsity_pattern;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;
  std::shared_ptr<TrilinosWrappers::PreconditionAMG> amg_preconditioner;

  const bool   SUPG        = true;
  const double GLS_u_scale = 1;
};


#endif
