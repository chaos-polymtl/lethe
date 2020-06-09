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

#ifndef LETHE_GLSSHARPNS_H
#define LETHE_GLSSHARPNS_H

#include "navier_stokes_base.h"

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
class GLSSharpNavierStokesSolver
  : public NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>
{
public:
  GLSSharpNavierStokesSolver(NavierStokesSolverParameters<dim> &nsparam,
                        const unsigned int                 degreeVelocity,
                        const unsigned int                 degreePressure);
  ~GLSSharpNavierStokesSolver();

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
  vertices_cell_mapping();

  void
  define_particules();

  void
  force_on_ib();

  void
  sharp_edge(const bool initial_step);


  double
  calculate_L2_error_particules();
  void
  finish_time_step_particules();



  void
  refine_ib();

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

  void
  solve_system_direct(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual,
                     const bool   renewed_matrix);
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
    std::vector<std::vector<typename DoFHandler<dim>::active_cell_iterator>> vertices_to_cell;
    const bool   SUPG        = false;
    const bool   PSPG        = true;
    const double GLS_u_scale = 1;
    double radius=0.21;
    double radius_2=0.6;
    bool couette= false;
    std::vector<std::vector<double>> particules;
    bool initial_step_bool;
    unsigned int iter_ib=0;
    Vector<double> ib_dof;

    TableHandler table_f;
    TableHandler table_t;

};


#endif
