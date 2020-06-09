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
 * Author: Bruno Blais, Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020-
 */

#ifndef lethe_gls_nitsche_navier_stokes_h
#define lethe_gls_nitsche_navier_stokes_h

#include "navier_stokes_base.h"

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

template <int dim>
class GLSNitscheNavierStokesSolver
  : public NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>
{
public:
  GLSNitscheNavierStokesSolver(NavierStokesSolverParameters<dim> &nsparam,
                        const unsigned int                 degreeVelocity,
                        const unsigned int                 degreePressure);
  ~GLSNitscheNavierStokesSolver();

  void
  solve();

protected:
  virtual void
  setup_dofs();
  void
  set_initial_condition(Parameters::InitialConditionType initial_condition_type,
                        bool                             restart = false);

  void
  set_solution_vector(double value);

private:
  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme,
            Parameters::VelocitySource::VelocitySourceType    velocity_source>
  void
  assembleGLS();

  void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override;

  void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  void
  assemble_L2_projection();

  /**
   * Interface for the solver for the linear system of equations
   */

  void
  solve_linear_system(
    const bool initial_step,
    const bool renewed_matrix = true) override; // Interface function

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
   * Set-up AMG preconditioner
   */
  void
  setup_AMG();

  /**
   * Set-up ILU preconditioner
   */
  void
  setup_ILU();


  /**
   * Members
   */
private:
  SparsityPattern                                    sparsity_pattern;
  TrilinosWrappers::SparseMatrix                     system_matrix;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;
  std::shared_ptr<TrilinosWrappers::PreconditionAMG> amg_preconditioner;

  const bool   SUPG        = true;
  const double GLS_u_scale = 1;
};


#endif
