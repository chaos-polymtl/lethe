// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Locks in the coupling between a Parameters::* struct's in-class
 * default member initializers and the default strings its
 * declare_parameters() derives from them.
 *
 * For a representative subset of structs in include/core/parameters.h, this
 * default-constructs an instance (using the in-class initializers), calls
 * declare_parameters() into a fresh ParameterHandler, then parse_parameters()
 * the untouched defaults into a second, independently default-constructed
 * instance. If declare_parameters() and the header defaults ever drift apart,
 * the two instances stop matching and this test fails.
 */

// Lethe
#include <core/parameters.h>

// Tests (with common definitions)
#include <../tests/tests.h>

template <typename T>
void
check(const std::string &label, const T &computed, const T &expected)
{
  if (computed == expected)
    deallog << "  OK  " << label << std::endl;
  else
    deallog << "  FAIL " << label << std::endl;
}

void
check(const std::string &label, const double computed, const double expected)
{
  constexpr double tol = 1e-12;
  if (std::abs(computed - expected) <= tol * std::max(1.0, std::abs(expected)))
    deallog << "  OK  " << label << std::endl;
  else
    deallog << "  FAIL " << label << std::endl;
}

void
test_timer()
{
  deallog << "--- Timer ---" << std::endl;
  const Parameters::Timer expected;
  ParameterHandler        prm;
  Parameters::Timer::declare_parameters(prm);
  Parameters::Timer actual;
  actual.parse_parameters(prm);

  check("type", actual.type, expected.type);
  check("write_time_in_error_table",
        actual.write_time_in_error_table,
        expected.write_time_in_error_table);
}

void
test_non_linear_solver()
{
  deallog << "--- NonLinearSolver ---" << std::endl;
  const Parameters::NonLinearSolver expected;
  ParameterHandler                  prm;
  Parameters::NonLinearSolver::declare_parameters(prm, "fluid dynamics");
  Parameters::NonLinearSolver actual;
  actual.parse_parameters(prm, "fluid dynamics");

  check("verbosity", actual.verbosity, expected.verbosity);
  check("solver", actual.solver, expected.solver);
  check("kinsol_strategy", actual.kinsol_strategy, expected.kinsol_strategy);
  check("tolerance", actual.tolerance, expected.tolerance);
  check("max_iterations", actual.max_iterations, expected.max_iterations);
  check("force_rhs_calculation",
        actual.force_rhs_calculation,
        expected.force_rhs_calculation);
  check("matrix_tolerance", actual.matrix_tolerance, expected.matrix_tolerance);
  check("step_tolerance", actual.step_tolerance, expected.step_tolerance);
  check("reuse_matrix", actual.reuse_matrix, expected.reuse_matrix);
  check("reuse_preconditioner",
        actual.reuse_preconditioner,
        expected.reuse_preconditioner);
  check("abort_at_convergence_failure",
        actual.abort_at_convergence_failure,
        expected.abort_at_convergence_failure);
}

void
test_mesh()
{
  deallog << "--- Mesh ---" << std::endl;
  const Parameters::Mesh expected;
  ParameterHandler       prm;
  Parameters::Mesh::declare_parameters(prm);
  Parameters::Mesh actual;
  actual.parse_parameters(prm);

  check("type", actual.type, expected.type);
  check("file_name", actual.file_name, expected.file_name);
  check("grid_type", actual.grid_type, expected.grid_type);
  check("grid_arguments", actual.grid_arguments, expected.grid_arguments);
  check("initial_refinement",
        actual.initial_refinement,
        expected.initial_refinement);
  check("initial_refinement_at_boundaries",
        actual.initial_refinement_at_boundaries,
        expected.initial_refinement_at_boundaries);
  check("refine_until_target_size",
        actual.refine_until_target_size,
        expected.refine_until_target_size);
  check("simplex", actual.simplex, expected.simplex);
  check("target_size", actual.target_size, expected.target_size);
  check("check_for_diamond_cells",
        actual.check_for_diamond_cells,
        expected.check_for_diamond_cells);
  check("expand_particle_wall_contact_search",
        actual.expand_particle_wall_contact_search,
        expected.expand_particle_wall_contact_search);
  check("translation", actual.translation, expected.translation);
  check("rotation_axis", actual.rotation_axis, expected.rotation_axis);
  check("rotation_angle", actual.rotation_angle, expected.rotation_angle);
  check("scale", actual.scale, expected.scale);
}

void
test_linear_solver()
{
  deallog << "--- LinearSolver ---" << std::endl;
  const Parameters::LinearSolver expected;
  ParameterHandler               prm;
  Parameters::LinearSolver::declare_parameters(prm, "fluid dynamics");
  Parameters::LinearSolver actual;
  actual.parse_parameters(prm, "fluid dynamics");

  check("solver", actual.solver, expected.solver);
  check("verbosity", actual.verbosity, expected.verbosity);
  check("rescale_residual_by_volume",
        actual.rescale_residual_by_volume,
        expected.rescale_residual_by_volume);
  check("relative_residual",
        actual.relative_residual,
        expected.relative_residual);
  check("minimum_residual", actual.minimum_residual, expected.minimum_residual);
  check("max_iterations", actual.max_iterations, expected.max_iterations);
  check("max_krylov_vectors",
        actual.max_krylov_vectors,
        expected.max_krylov_vectors);
  check("enable_hessians_jacobian",
        actual.enable_hessians_jacobian,
        expected.enable_hessians_jacobian);
  check("enable_hessians_residual",
        actual.enable_hessians_residual,
        expected.enable_hessians_residual);
  check("preconditioner", actual.preconditioner, expected.preconditioner);
  check("ilu_precond_fill", actual.ilu_precond_fill, expected.ilu_precond_fill);
  check("ilu_precond_atol", actual.ilu_precond_atol, expected.ilu_precond_atol);
  check("amg_n_cycles", actual.amg_n_cycles, expected.amg_n_cycles);
  check("mg_coarsening_type",
        actual.mg_coarsening_type,
        expected.mg_coarsening_type);
  check("mg_p_coarsening_type",
        actual.mg_p_coarsening_type,
        expected.mg_p_coarsening_type);
  check("mg_smoother_preconditioner_type",
        actual.mg_smoother_preconditioner_type,
        expected.mg_smoother_preconditioner_type);
  check("mg_coarse_grid_solver",
        actual.mg_coarse_grid_solver,
        expected.mg_coarse_grid_solver);
  check("mg_gmres_preconditioner",
        actual.mg_gmres_preconditioner,
        expected.mg_gmres_preconditioner);
  check("mg_verbosity", actual.mg_verbosity, expected.mg_verbosity);
}

void
test_stabilization()
{
  deallog << "--- Stabilization ---" << std::endl;
  const Parameters::Stabilization expected;
  ParameterHandler                prm;
  Parameters::Stabilization::declare_parameters(prm);
  Parameters::Stabilization actual;
  actual.parse_parameters(prm);

  check("use_default_stabilization",
        actual.use_default_stabilization,
        expected.use_default_stabilization);
  check("heat_transfer_dcdd_stabilization",
        actual.heat_transfer_dcdd_stabilization,
        expected.heat_transfer_dcdd_stabilization);
  check("cls_dcdd_stabilization",
        actual.cls_dcdd_stabilization,
        expected.cls_dcdd_stabilization);
  check("dcdd_diffusion_coeff",
        actual.dcdd_diffusion_coeff,
        expected.dcdd_diffusion_coeff);
  check("pressure_scaling_factor",
        actual.pressure_scaling_factor,
        expected.pressure_scaling_factor);
  check("stabilization", actual.stabilization, expected.stabilization);
  check("scalar_limiter", actual.scalar_limiter, expected.scalar_limiter);
}

void
test_mesh_box_refinement()
{
  deallog << "--- MeshBoxRefinement ---" << std::endl;
  const Parameters::MeshBoxRefinement expected;
  ParameterHandler                    prm;
  // declare_parameters() is a non-static member function here (it sizes
  // per-box subsections from instance state), so any fresh instance's
  // in-class defaults can be used to declare the entries.
  Parameters::MeshBoxRefinement declare_source;
  declare_source.declare_parameters(prm);
  Parameters::MeshBoxRefinement actual;
  actual.parse_parameters(prm);

  check("number_of_refinement_boxes",
        actual.number_of_refinement_boxes,
        expected.number_of_refinement_boxes);
}

void
test()
{
  test_timer();
  test_non_linear_solver();
  test_mesh();
  test_linear_solver();
  test_stabilization();
  test_mesh_box_refinement();
}

int
main()
{
  try
    {
      initlog();
      test();
      deallog << "OK" << std::endl;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
