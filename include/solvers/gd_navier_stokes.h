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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef LETHE_GDNS_H
#define LETHE_GDNS_H

#include <deal.II/lac/trilinos_block_sparse_matrix.h>

#include "core/bdf.h"
#include "navier_stokes_base.h"

using namespace dealii;

/**
 * A solver class for the Steady-Sate  Navier-Stokes equation using Grad-Div
 * stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim>
class GDNavierStokesSolver
  : public NavierStokesBase<dim,
                            TrilinosWrappers::MPI::BlockVector,
                            std::vector<IndexSet>>
{
public:
  GDNavierStokesSolver(NavierStokesSolverParameters<dim> &nsparam,
                       const unsigned int                 degreeVelocity,
                       const unsigned int                 degreePressure);
  ~GDNavierStokesSolver();

  void
  solve();

private:
  void
  assemble_matrix_and_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method) override;

  void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override;

  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme>
  void
  assembleGD();
  void
  assemble_L2_projection();

  void
  set_nodal_values();

  virtual void
  setup_dofs();

  void
  set_initial_condition(Parameters::InitialConditionType initial_condition_type,
                        bool                             restart = false);
  void
  set_solution_vector(double value);

  /**
   * Solver for the L2 Projection linear system
   */

  void
  solve_L2_system(bool   initial_step,
                  double relative_residual,
                  double minimum_residual);

  /**
   * Solver for the NS linear system of equations
   */

  void
  solve_linear_system(bool   initial_step,
                      double relative_residual,
                      double minimum_residual) override;

  /**
   * GMRES solver with ILU preconditioning
   */
  void
  solve_system_GMRES(bool   initial_step,
                     double relative_residual,
                     double minimum_residual);
  /**
   * GMRES solver with AMG preconditioning
   */
  void
  solve_system_AMG(bool   initial_step,
                   double relative_residual,
                   double minimum_residual);


  /**
   * Members
   */
  TrilinosWrappers::BlockSparsityPattern sparsity_pattern;
  TrilinosWrappers::BlockSparseMatrix    system_matrix;
  TrilinosWrappers::SparseMatrix         pressure_mass_matrix;

  std::vector<types::global_dof_index> dofs_per_block;

  const double gamma = 1;
};

template <class BSPreconditioner>
class BlockSchurPreconditioner : public Subscriptor
{
public:
  BlockSchurPreconditioner(double                                     gamma,
                           double                                     viscosity,
                           const TrilinosWrappers::BlockSparseMatrix &S,
                           const TrilinosWrappers::SparseMatrix &     P,
                           const BSPreconditioner * p_amat_preconditioner,
                           const BSPreconditioner * p_pmass_preconditioner,
                           Parameters::LinearSolver solver_parameters);

  void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const;

private:
  const double                               gamma;
  const double                               viscosity;
  const Parameters::LinearSolver             linear_solver_parameters;
  const TrilinosWrappers::BlockSparseMatrix &stokes_matrix;
  const TrilinosWrappers::SparseMatrix &     pressure_mass_matrix;
  const BSPreconditioner *                   amat_preconditioner;
  const BSPreconditioner *                   pmass_preconditioner;
};

// We can notice that the initialization of the inverse of the matrix at the
// top left corner is completed in the constructor. If so, every application
// of the preconditioner then no longer requires the computation of the
// matrix factors.
template <typename BSPreconditioner>
BlockSchurPreconditioner<BSPreconditioner>::BlockSchurPreconditioner(
  double                                     gamma,
  double                                     viscosity,
  const TrilinosWrappers::BlockSparseMatrix &S,
  const TrilinosWrappers::SparseMatrix &     P,
  const BSPreconditioner *                   p_amat_preconditioner,
  const BSPreconditioner *                   p_pmass_preconditioner,

  Parameters::LinearSolver p_solver_parameters)
  : gamma(gamma)
  , viscosity(viscosity)
  , linear_solver_parameters(p_solver_parameters)
  , stokes_matrix(S)
  , pressure_mass_matrix(P)
  , amat_preconditioner(p_amat_preconditioner)
  , pmass_preconditioner(p_pmass_preconditioner)
{}

template <class BSPreconditioner>
void
BlockSchurPreconditioner<BSPreconditioner>::vmult(
  TrilinosWrappers::MPI::BlockVector &      dst,
  const TrilinosWrappers::MPI::BlockVector &src) const
{
  //    TimerOutput computing_timer(std::cout,
  //                                TimerOutput::summary,
  //                                TimerOutput::wall_times);

  TrilinosWrappers::MPI::Vector utmp(src.block(0));
  {
    //        computing_timer.enter_section("Pressure");
    SolverControl solver_control(
      linear_solver_parameters.max_iterations,
      std::max(
        1e-3 * src.block(1).l2_norm(),
        // linear_solver_parameters.relative_residual*src.block(0).l2_norm(),
        linear_solver_parameters.minimum_residual));
    TrilinosWrappers::SolverCG cg(solver_control);

    dst.block(1) = 0.0;
    cg.solve(pressure_mass_matrix,
             dst.block(1),
             src.block(1),
             *pmass_preconditioner);
    dst.block(1) *= -(viscosity + gamma);
    //        computing_timer.exit_section("Pressure");
  }

  {
    //        computing_timer.enter_section("Operations");
    stokes_matrix.block(0, 1).vmult(utmp, dst.block(1));
    utmp *= -1.0;
    utmp += src.block(0);
    //        computing_timer.exit_section("Operations");
  }
  {
    //        computing_timer.enter_section("A Matrix");
    SolverControl solver_control(
      linear_solver_parameters.max_iterations,
      std::max(1e-1 * src.block(0).l2_norm(),
               linear_solver_parameters.minimum_residual));



    // TrilinosWrappers::SolverBicgstab solver(solver_control);
    TrilinosWrappers::SolverGMRES solver(solver_control);

    // A_inverse.solve(stokes_matrix.block(0, 0),dst.block(0), utmp,
    // mp_preconditioner);
    solver.solve(stokes_matrix.block(0, 0),
                 dst.block(0),
                 utmp,
                 *amat_preconditioner);
    //        computing_timer.exit_section("A Matrix");
  }
}



template <class PreconditionerMp>
class BlockDiagPreconditioner : public Subscriptor
{
public:
  BlockDiagPreconditioner(const TrilinosWrappers::BlockSparseMatrix &S,
                          const PreconditionerMp &Mppreconditioner,
                          SolverControl &         A_parameters);

  void
  vmult(TrilinosWrappers::MPI::BlockVector &      dst,
        const TrilinosWrappers::MPI::BlockVector &src) const;

private:
  const TrilinosWrappers::BlockSparseMatrix & stokes_matrix;
  TrilinosWrappers::PreconditionILU           amat_preconditioner;
  TrilinosWrappers::PreconditionILU           pmat_preconditioner;
  const PreconditionerMp &                    mp_preconditioner;
  SolverFGMRES<TrilinosWrappers::MPI::Vector> A_inverse;
};

template <class PreconditionerMp>
BlockDiagPreconditioner<PreconditionerMp>::BlockDiagPreconditioner(
  const TrilinosWrappers::BlockSparseMatrix &S,
  const PreconditionerMp &                   Mppreconditioner,
  SolverControl &                            A_parameters)
  : stokes_matrix(S)
  , mp_preconditioner(Mppreconditioner)
  , A_inverse(A_parameters)
{
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(0,
                                                                          1e-12,
                                                                          1,
                                                                          0);
  amat_preconditioner.initialize(stokes_matrix.block(0, 0),
                                 preconditionerOptions);

  TrilinosWrappers::PreconditionILU::AdditionalData pmat_preconditionerOptions(
    0, 1e-12, 1, 0);
  TrilinosWrappers::PreconditionILU pmat_preconditioner;
  pmat_preconditioner.initialize(stokes_matrix.block(1, 1),
                                 pmat_preconditionerOptions);
}

template <class PreconditionerMp>
void
BlockDiagPreconditioner<PreconditionerMp>::vmult(
  TrilinosWrappers::MPI::BlockVector &      dst,
  const TrilinosWrappers::MPI::BlockVector &src) const
{
  {
    SolverControl              solver_control(100000, 1e-12);
    TrilinosWrappers::SolverCG cg(solver_control);

    cg.solve(stokes_matrix.block(1, 1),
             dst.block(1),
             src.block(1),
             pmat_preconditioner);
  }
  {
    SolverControl solver_control(10000, 1e-12);

    TrilinosWrappers::SolverGMRES solver(solver_control);
    solver.solve(stokes_matrix.block(0, 0),
                 dst.block(0),
                 src.block(0),
                 amat_preconditioner);
  }
}



#endif
