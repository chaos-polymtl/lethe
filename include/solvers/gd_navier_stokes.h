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

#ifndef lethe_gd_navier_stokes_h
#define lethe_gd_navier_stokes_h

#include <core/exceptions.h>
#include <core/multiphysics.h>
#include <core/parameters.h>

#include <solvers/navier_stokes_base.h>

#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

using namespace dealii;

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
  const Parameters::LinearSolver             &linear_solver_parameters;
  const TrilinosWrappers::BlockSparseMatrix &stokes_matrix;
  const TrilinosWrappers::SparseMatrix &     pressure_mass_matrix;
  const BSPreconditioner *                   amat_preconditioner;
  const BSPreconditioner *                   pmass_preconditioner;
};

/**
 * A solver class for the Navier-Stokes equation using Grad-Div
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
  GDNavierStokesSolver(SimulationParameters<dim> &nsparam);
  ~GDNavierStokesSolver();

  void
  solve();

private:
  /**
   *  @brief Assembles the matrix associated with the solver
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief  Update the average velocity field solution in the multiphysics interface
   */
  virtual void
  update_multiphysics_time_average_solution() override;


  /**
   * @brief Assemble the local matrix for a given cell.
   *
   * This function is used by the WorkStream class to assemble
   * the system matrix. It is a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data);

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data);

  /**
   * @brief sets up the vector of assembler functions
   */
  void
  setup_assemblers();


  /**
   * @brief Copy local cell information to global matrix
   */

  void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data);

  /**
   * @brief Copy local cell rhs information to global rhs
   */

  void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data);

  template <bool                                              assemble_matrix,
            Parameters::SimulationControl::TimeSteppingMethod scheme>
  void
  assembleGD();

  void
  assemble_L2_projection();

  virtual void
  setup_dofs_fd() override;

  virtual void
  set_initial_condition_fd(
    Parameters::InitialConditionType initial_condition_type,
    bool                             restart = false) override;

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
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix) override;

  /**
   * GMRES solver with ILU preconditioning
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double relative_residual,
                     const double minimum_residual,
                     const bool   renewed_matrix);
  /**
   * GMRES solver with AMG preconditioning
   */
  void
  solve_system_AMG(const bool   initial_step,
                   const double relative_residual,
                   const double minimum_residual,
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
  TrilinosWrappers::BlockSparsityPattern sparsity_pattern;
  TrilinosWrappers::BlockSparseMatrix    system_matrix;
  TrilinosWrappers::SparseMatrix         pressure_mass_matrix;

  std::vector<types::global_dof_index> dofs_per_block;

  std::shared_ptr<TrilinosWrappers::PreconditionILU>
    velocity_ilu_preconditioner;
  std::shared_ptr<TrilinosWrappers::PreconditionAMG>
    velocity_amg_preconditioner;
  std::shared_ptr<TrilinosWrappers::PreconditionILU>
    pressure_ilu_preconditioner;
  std::shared_ptr<TrilinosWrappers::PreconditionAMG>
    pressure_amg_preconditioner;

  std::shared_ptr<BlockSchurPreconditioner<TrilinosWrappers::PreconditionILU>>
    system_ilu_preconditioner;

  std::shared_ptr<BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG>>
    system_amg_preconditioner;

  const double gamma = 1;
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
  Parameters::LinearSolver                   p_solver_parameters)
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
      linear_solver_parameters.max_iterations.at(PhysicsID::fluid_dynamics),
      std::max(
        1e-3 * src.block(1).l2_norm(),
        // linear_solver_parameters.relative_residual*src.block(0).l2_norm(),
        linear_solver_parameters.minimum_residual.at(
          PhysicsID::fluid_dynamics)));
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
      linear_solver_parameters.max_iterations.at(PhysicsID::fluid_dynamics),
      std::max(1e-1 * src.block(0).l2_norm(),
               linear_solver_parameters.minimum_residual.at(
                 PhysicsID::fluid_dynamics)));



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
