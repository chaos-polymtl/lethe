// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef tests_utilities_h
#define tests_utilities_h

#include <core/parameters.h>

Parameters::LinearSolver
make_default_linear_solver()
{
  Parameters::LinearSolver linear_solver_parameters;
  linear_solver_parameters.minimum_residual  = 1e-15;
  linear_solver_parameters.relative_residual = 1e-15;
  linear_solver_parameters.solver = Parameters::LinearSolver::SolverType::gmres;
  linear_solver_parameters.max_iterations     = 100;
  linear_solver_parameters.max_krylov_vectors = 100;
  linear_solver_parameters.ilu_precond_atol   = 1e-14;
  linear_solver_parameters.ilu_precond_fill   = 0.;
  linear_solver_parameters.ilu_precond_rtol   = 1.;
  linear_solver_parameters.verbosity          = Parameters::Verbosity::verbose;
  linear_solver_parameters.preconditioner =
    Parameters::LinearSolver::PreconditionerType::ilu;
  return linear_solver_parameters;
}

std::shared_ptr<Parameters::VoidFractionParameters<3>>
make_default_void_fraction_parameters()
{
  std::shared_ptr<Parameters::VoidFractionParameters<3>>
    void_fraction_parameters =
      std::make_shared<Parameters::VoidFractionParameters<3>>();

  void_fraction_parameters->mode = Parameters::VoidFractionMode::qcm;
  void_fraction_parameters->l2_smoothing_length          = 0.;
  void_fraction_parameters->n_quadrature_points          = 3;
  void_fraction_parameters->qcm_sphere_diameter          = 0.5;
  void_fraction_parameters->qcm_sphere_equal_cell_volume = false;
  void_fraction_parameters->quadrature_rule =
    Parameters::VoidFractionQuadratureRule::gauss;
  void_fraction_parameters->project_particle_velocity  = false;
  void_fraction_parameters->read_dem                   = false;
  void_fraction_parameters->particle_refinement_factor = 1;
  return void_fraction_parameters;
}

#endif
