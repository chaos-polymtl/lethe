// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Locks in the coupling between a Parameters::* struct's in-class
 * default member initializers and the default strings its
 * declare_parameters() derives from them.
 *
 * For a representative subset of structs across include/core/parameters.h,
 * parameters_cfd_dem.h, parameters_multiphysics.h and
 * parameters_lagrangian.h, this default-constructs an instance
 * (using the in-class initializers), calls declare_parameters() into a fresh
 * ParameterHandler, then parse_parameters() the untouched defaults into a
 * second, independently default-constructed instance. If declare_parameters()
 * and the header defaults ever drift apart, the two instances stop matching
 * and this test fails.
 */

// Lethe
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/parameters_lagrangian.h>
#include <core/parameters_multiphysics.h>

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
  const Parameters::Mesh<2> expected;
  ParameterHandler          prm;
  Parameters::Mesh<2>::declare_parameters(prm);
  Parameters::Mesh<2> actual;
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
  check("boundaries_to_refine",
        actual.boundaries_to_refine,
        expected.boundaries_to_refine);
  check("boundaries_to_refine",
        actual.boundaries_to_refine,
        expected.boundaries_to_refine);
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
  check("ilu_precond_rtol", actual.ilu_precond_rtol, expected.ilu_precond_rtol);
  check("amg_precond_ilu_fill",
        actual.amg_precond_ilu_fill,
        expected.amg_precond_ilu_fill);
  check("amg_precond_ilu_atol",
        actual.amg_precond_ilu_atol,
        expected.amg_precond_ilu_atol);
  check("amg_precond_ilu_rtol",
        actual.amg_precond_ilu_rtol,
        expected.amg_precond_ilu_rtol);
  check("amg_aggregation_threshold",
        actual.amg_aggregation_threshold,
        expected.amg_aggregation_threshold);
  check("ilu_precond_rtol", actual.ilu_precond_rtol, expected.ilu_precond_rtol);
  check("amg_precond_ilu_fill",
        actual.amg_precond_ilu_fill,
        expected.amg_precond_ilu_fill);
  check("amg_precond_ilu_atol",
        actual.amg_precond_ilu_atol,
        expected.amg_precond_ilu_atol);
  check("amg_precond_ilu_rtol",
        actual.amg_precond_ilu_rtol,
        expected.amg_precond_ilu_rtol);
  check("amg_aggregation_threshold",
        actual.amg_aggregation_threshold,
        expected.amg_aggregation_threshold);
  check("amg_n_cycles", actual.amg_n_cycles, expected.amg_n_cycles);
  check("amg_w_cycles", actual.amg_w_cycles, expected.amg_w_cycles);
  check("amg_smoother_sweeps",
        actual.amg_smoother_sweeps,
        expected.amg_smoother_sweeps);
  check("amg_smoother_overlap",
        actual.amg_smoother_overlap,
        expected.amg_smoother_overlap);
  check("force_linear_solver_continuation",
        actual.force_linear_solver_continuation,
        expected.force_linear_solver_continuation);
  check("mg_min_level", actual.mg_min_level, expected.mg_min_level);
  check("mg_level_min_cells",
        actual.mg_level_min_cells,
        expected.mg_level_min_cells);
  check("mg_int_level", actual.mg_int_level, expected.mg_int_level);
  check("mg_enable_hessians_jacobian",
        actual.mg_enable_hessians_jacobian,
        expected.mg_enable_hessians_jacobian);
  check("mg_smoother_iterations",
        actual.mg_smoother_iterations,
        expected.mg_smoother_iterations);
  check("mg_smoother_relaxation",
        actual.mg_smoother_relaxation,
        expected.mg_smoother_relaxation);
  check("mg_smoother_chebyshev_degree",
        actual.mg_smoother_chebyshev_degree,
        expected.mg_smoother_chebyshev_degree);
  check("mg_smoother_chebyshev_smoothing_range",
        actual.mg_smoother_chebyshev_smoothing_range,
        expected.mg_smoother_chebyshev_smoothing_range);
  check("mg_smoother_eig_estimation",
        actual.mg_smoother_eig_estimation,
        expected.mg_smoother_eig_estimation);
  check("eig_estimation_smoothing_range",
        actual.eig_estimation_smoothing_range,
        expected.eig_estimation_smoothing_range);
  check("eig_estimation_cg_n_iterations",
        actual.eig_estimation_cg_n_iterations,
        expected.eig_estimation_cg_n_iterations);
  check("eig_estimation_verbose",
        actual.eig_estimation_verbose,
        expected.eig_estimation_verbose);
  check("mg_use_fe_q_iso_q1",
        actual.mg_use_fe_q_iso_q1,
        expected.mg_use_fe_q_iso_q1);
  check("mg_p_min_coarsening_degree",
        actual.mg_p_min_coarsening_degree,
        expected.mg_p_min_coarsening_degree);
  check("amg_w_cycles", actual.amg_w_cycles, expected.amg_w_cycles);
  check("amg_smoother_sweeps",
        actual.amg_smoother_sweeps,
        expected.amg_smoother_sweeps);
  check("amg_smoother_overlap",
        actual.amg_smoother_overlap,
        expected.amg_smoother_overlap);
  check("force_linear_solver_continuation",
        actual.force_linear_solver_continuation,
        expected.force_linear_solver_continuation);
  check("mg_min_level", actual.mg_min_level, expected.mg_min_level);
  check("mg_level_min_cells",
        actual.mg_level_min_cells,
        expected.mg_level_min_cells);
  check("mg_int_level", actual.mg_int_level, expected.mg_int_level);
  check("mg_enable_hessians_jacobian",
        actual.mg_enable_hessians_jacobian,
        expected.mg_enable_hessians_jacobian);
  check("mg_smoother_iterations",
        actual.mg_smoother_iterations,
        expected.mg_smoother_iterations);
  check("mg_smoother_relaxation",
        actual.mg_smoother_relaxation,
        expected.mg_smoother_relaxation);
  check("mg_smoother_chebyshev_degree",
        actual.mg_smoother_chebyshev_degree,
        expected.mg_smoother_chebyshev_degree);
  check("mg_smoother_chebyshev_smoothing_range",
        actual.mg_smoother_chebyshev_smoothing_range,
        expected.mg_smoother_chebyshev_smoothing_range);
    check("mg_smoother_chebyshev_eig_cg_n_iterations",  
        actual.mg_smoother_chebyshev_eig_cg_n_iterations,  
        expected.mg_smoother_chebyshev_eig_cg_n_iterations); 
  check("mg_smoother_eig_estimation",
        actual.mg_smoother_eig_estimation,
        expected.mg_smoother_eig_estimation);
  check("eig_estimation_smoothing_range",
        actual.eig_estimation_smoothing_range,
        expected.eig_estimation_smoothing_range);
  check("eig_estimation_cg_n_iterations",
        actual.eig_estimation_cg_n_iterations,
        expected.eig_estimation_cg_n_iterations);
  check("eig_estimation_verbose",
        actual.eig_estimation_verbose,
        expected.eig_estimation_verbose);
  check("mg_use_fe_q_iso_q1",
        actual.mg_use_fe_q_iso_q1,
        expected.mg_use_fe_q_iso_q1);
  check("mg_p_min_coarsening_degree",
        actual.mg_p_min_coarsening_degree,
        expected.mg_p_min_coarsening_degree);
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
  check("mg_gmres_max_iterations",
        actual.mg_gmres_max_iterations,
        expected.mg_gmres_max_iterations);
  check("mg_gmres_tolerance",
        actual.mg_gmres_tolerance,
        expected.mg_gmres_tolerance);
  check("mg_gmres_reduce", actual.mg_gmres_reduce, expected.mg_gmres_reduce);
  check("mg_gmres_max_krylov_vectors",
        actual.mg_gmres_max_krylov_vectors,
        expected.mg_gmres_max_krylov_vectors);
  check("mg_amg_use_default_parameters",
        actual.mg_amg_use_default_parameters,
        expected.mg_amg_use_default_parameters);
  check("mg_gmres_max_iterations",
        actual.mg_gmres_max_iterations,
        expected.mg_gmres_max_iterations);
  check("mg_gmres_tolerance",
        actual.mg_gmres_tolerance,
        expected.mg_gmres_tolerance);
  check("mg_gmres_reduce", actual.mg_gmres_reduce, expected.mg_gmres_reduce);
  check("mg_gmres_max_krylov_vectors",
        actual.mg_gmres_max_krylov_vectors,
        expected.mg_gmres_max_krylov_vectors);
  check("mg_amg_use_default_parameters",
        actual.mg_amg_use_default_parameters,
        expected.mg_amg_use_default_parameters);
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
  const Parameters::MeshBoxRefinement<2> expected;
  ParameterHandler                       prm;
  // declare_parameters() is a non-static member function here (it sizes
  // per-box subsections from instance state), so any fresh instance's
  // in-class defaults can be used to declare the entries.
  Parameters::MeshBoxRefinement<2> declare_source;
  declare_source.declare_parameters(prm);
  Parameters::MeshBoxRefinement<2> actual;
  actual.parse_parameters(prm);

  check("number_of_refinement_boxes",
        actual.number_of_refinement_boxes,
        expected.number_of_refinement_boxes);
    check("max_number_of_refinement_boxes",
        actual.max_number_of_refinement_boxes,
        expected.max_number_of_refinement_boxes);
}

void
test_cfddem()
{
  deallog << "--- CFDDEM ---" << std::endl;
  const Parameters::CFDDEM expected;
  ParameterHandler         prm;
  Parameters::CFDDEM::declare_parameters(prm);
  Parameters::CFDDEM actual;
  actual.parse_parameters(prm);

  check("grad_div", actual.grad_div, expected.grad_div);
  check("void_fraction_time_derivative",
        actual.void_fraction_time_derivative,
        expected.void_fraction_time_derivative);
  check("interpolated_void_fraction",
        actual.interpolated_void_fraction,
        expected.interpolated_void_fraction);
  check("drag_force", actual.drag_force, expected.drag_force);
  check("buoyancy_force", actual.buoyancy_force, expected.buoyancy_force);
  check("shear_force", actual.shear_force, expected.shear_force);
  check("pressure_force", actual.pressure_force, expected.pressure_force);
  check("saffman_lift_force",
        actual.saffman_lift_force,
        expected.saffman_lift_force);
  check("magnus_lift_force",
        actual.magnus_lift_force,
        expected.magnus_lift_force);
  check("rotational_viscous_torque",
        actual.rotational_viscous_torque,
        expected.rotational_viscous_torque);
  check("vortical_viscous_torque",
        actual.vortical_viscous_torque,
        expected.vortical_viscous_torque);
  check("void_fraction_time_derivative",
        actual.void_fraction_time_derivative,
        expected.void_fraction_time_derivative);
  check("interpolated_void_fraction",
        actual.interpolated_void_fraction,
        expected.interpolated_void_fraction);
  check("drag_force", actual.drag_force, expected.drag_force);
  check("buoyancy_force", actual.buoyancy_force, expected.buoyancy_force);
  check("shear_force", actual.shear_force, expected.shear_force);
  check("pressure_force", actual.pressure_force, expected.pressure_force);
  check("saffman_lift_force",
        actual.saffman_lift_force,
        expected.saffman_lift_force);
  check("magnus_lift_force",
        actual.magnus_lift_force,
        expected.magnus_lift_force);
  check("rotational_viscous_torque",
        actual.rotational_viscous_torque,
        expected.rotational_viscous_torque);
  check("vortical_viscous_torque",
        actual.vortical_viscous_torque,
        expected.vortical_viscous_torque);
  check("drag_model", actual.drag_model, expected.drag_model);
  check("drag_coupling", actual.drag_coupling, expected.drag_coupling);
  check("dem_iteration_control",
        actual.dem_iteration_control,
        expected.dem_iteration_control);
  check("dem_iteration_control",
        actual.dem_iteration_control,
        expected.dem_iteration_control);
  check("vans_model", actual.vans_model, expected.vans_model);
  check("coupling_frequency",
        actual.coupling_frequency,
        expected.coupling_frequency);
  check("fraction_of_rayleigh_time",
        actual.fraction_of_rayleigh_time,
        expected.fraction_of_rayleigh_time);
  check("cstar", actual.cstar, expected.cstar);
  check("implicit_stabilization",
        actual.implicit_stabilization,
        expected.implicit_stabilization);
  check("implicit_stabilization",
        actual.implicit_stabilization,
        expected.implicit_stabilization);
  check("particle_statistics",
        actual.particle_statistics,
        expected.particle_statistics);
  check("project_particle_forces",
        actual.project_particle_forces,
        expected.project_particle_forces);
  check("project_particle_forces",
        actual.project_particle_forces,
        expected.project_particle_forces);
}

void
test_multiphysics()
{
  deallog << "--- Multiphysics<2> ---" << std::endl;
  const Parameters::Multiphysics<2> expected;
  ParameterHandler                  prm;
  expected.declare_parameters(prm);
  Parameters::Multiphysics<2>      actual;
  const Parameters::Dimensionality dimensions;
  actual.parse_parameters(prm, dimensions);

  check("fluid_dynamics", actual.fluid_dynamics, expected.fluid_dynamics);
  check("heat_transfer", actual.heat_transfer, expected.heat_transfer);
  check("tracer", actual.tracer, expected.tracer);
  check("CLS", actual.CLS, expected.CLS);
  check("cahn_hilliard", actual.cahn_hilliard, expected.cahn_hilliard);
  check("electromagnetics", actual.electromagnetics, expected.electromagnetics);
  check("viscous_dissipation",
        actual.viscous_dissipation,
        expected.viscous_dissipation);
  check("thermal_buoyancy_force",
        actual.thermal_buoyancy_force,
        expected.thermal_buoyancy_force);
  check("microwave_heating",
        actual.microwave_heating,
        expected.microwave_heating);
  check("cls_parameters.viscous_dissipative_fluid",
        actual.cls_parameters.viscous_dissipative_fluid,
        expected.cls_parameters.viscous_dissipative_fluid);
  check("cls_parameters.diffusivity",
        actual.cls_parameters.diffusivity,
        expected.cls_parameters.diffusivity);
  check("cls_parameters.compressible",
        actual.cls_parameters.compressible,
        expected.cls_parameters.compressible);
  check("cahn_hilliard_parameters.potential_smoothing_coefficient",
        actual.cahn_hilliard_parameters.potential_smoothing_coefficient,
        expected.cahn_hilliard_parameters.potential_smoothing_coefficient);
  check("cahn_hilliard_parameters.epsilon_set_method",
        actual.cahn_hilliard_parameters.epsilon_set_method,
        expected.cahn_hilliard_parameters.epsilon_set_method);
  check("cahn_hilliard_parameters.epsilon",
        actual.cahn_hilliard_parameters.epsilon,
        expected.cahn_hilliard_parameters.epsilon);
  check("cahn_hilliard_parameters.epsilon_verbosity",
        actual.cahn_hilliard_parameters.epsilon_verbosity,
        expected.cahn_hilliard_parameters.epsilon_verbosity);
  check("time_harmonic_maxwell_parameters.time_coupling_strategy",
        actual.time_harmonic_maxwell_parameters.time_coupling_strategy,
        expected.time_harmonic_maxwell_parameters.time_coupling_strategy);
  check("time_harmonic_maxwell_parameters.coupling_iteration",
        actual.time_harmonic_maxwell_parameters.coupling_iteration,
        expected.time_harmonic_maxwell_parameters.coupling_iteration);
  check("time_harmonic_maxwell_parameters.coupling_time",
        actual.time_harmonic_maxwell_parameters.coupling_time,
        expected.time_harmonic_maxwell_parameters.coupling_time);
  check("time_harmonic_maxwell_parameters.coupling_threshold",
        actual.time_harmonic_maxwell_parameters.coupling_threshold,
        expected.time_harmonic_maxwell_parameters.coupling_threshold);
  // When parsing, the electromagnetic_frequency is also scaled
  check("time_harmonic_maxwell_parameters.electromagnetic_frequency",
        actual.time_harmonic_maxwell_parameters.electromagnetic_frequency,
        expected.time_harmonic_maxwell_parameters.electromagnetic_frequency *
          dimensions.electromagnetic_frequency_scaling);
  check("time_harmonic_maxwell_parameters.electromagnetic_scaling_type",
        actual.time_harmonic_maxwell_parameters.electromagnetic_scaling_type,
        expected.time_harmonic_maxwell_parameters.electromagnetic_scaling_type);
  check("time_harmonic_maxwell_parameters.electric_field_amplitude",
        actual.time_harmonic_maxwell_parameters.electric_field_amplitude,
        expected.time_harmonic_maxwell_parameters.electric_field_amplitude);
  check("time_harmonic_maxwell_parameters.magnetic_field_amplitude",
        actual.time_harmonic_maxwell_parameters.magnetic_field_amplitude,
        expected.time_harmonic_maxwell_parameters.magnetic_field_amplitude);
  check("time_harmonic_maxwell_parameters.number_of_waveguide_inlets",
        actual.time_harmonic_maxwell_parameters.number_of_waveguide_inlets,
        expected.time_harmonic_maxwell_parameters.number_of_waveguide_inlets);
  check("time_harmonic_maxwell_parameters.waveguide_mode",
        actual.time_harmonic_maxwell_parameters.waveguide_mode,
        expected.time_harmonic_maxwell_parameters.waveguide_mode);
  check("time_harmonic_maxwell_parameters.mode_order_m",
        actual.time_harmonic_maxwell_parameters.mode_order_m,
        expected.time_harmonic_maxwell_parameters.mode_order_m);
  check("time_harmonic_maxwell_parameters.mode_order_n",
        actual.time_harmonic_maxwell_parameters.mode_order_n,
        expected.time_harmonic_maxwell_parameters.mode_order_n);
  check("time_harmonic_maxwell_parameters.waveguide_corners",
        actual.time_harmonic_maxwell_parameters.waveguide_corners,
        expected.time_harmonic_maxwell_parameters.waveguide_corners);
}

void
test_lagrangian_model_parameters()
{
  deallog << "--- Lagrangian::ModelParameters<2> ---" << std::endl;
  const Parameters::Lagrangian::ModelParameters<2> expected;
  ParameterHandler                                 prm;
  Parameters::Lagrangian::ModelParameters<2>::declare_parameters(prm);
  Parameters::Lagrangian::ModelParameters<2> actual;
  actual.parse_parameters(prm);

  check("load_balance_method",
        actual.load_balance_method,
        expected.load_balance_method);
  check("load_balance_step",
        actual.load_balance_step,
        expected.load_balance_step);
  check("load_balance_frequency",
        actual.load_balance_frequency,
        expected.load_balance_frequency);
  check("load_balance_threshold",
        actual.load_balance_threshold,
        expected.load_balance_threshold);
  check("load_balance_threshold",
        actual.load_balance_threshold,
        expected.load_balance_threshold);
  check("dynamic_load_balance_check_frequency",
        actual.dynamic_load_balance_check_frequency,
        expected.dynamic_load_balance_check_frequency);
  check("load_balance_particle_weight",
        actual.load_balance_particle_weight,
        expected.load_balance_particle_weight);
  check("active_load_balancing_factor",
        actual.active_load_balancing_factor,
        expected.active_load_balancing_factor);
  check("inactive_load_balancing_factor",
        actual.inactive_load_balancing_factor,
        expected.inactive_load_balancing_factor);
  check("contact_detection_method",
        actual.contact_detection_method,
        expected.contact_detection_method);
  check("contact_detection_frequency",
        actual.contact_detection_frequency,
        expected.contact_detection_frequency);
  check("dynamic_contact_search_factor",
        actual.dynamic_contact_search_factor,
        expected.dynamic_contact_search_factor);
  check("neighborhood_threshold",
        actual.neighborhood_threshold,
        expected.neighborhood_threshold);
  check("particle_particle_contact_force_model",
        actual.particle_particle_contact_force_model,
        expected.particle_particle_contact_force_model);
  check("particle_wall_contact_force_method",
        actual.particle_wall_contact_force_method,
        expected.particle_wall_contact_force_method);
  check("dmt_cut_off_threshold",
        actual.dmt_cut_off_threshold,
        expected.dmt_cut_off_threshold);
  check("rolling_resistance_method",
        actual.rolling_resistance_method,
        expected.rolling_resistance_method);
  check("f_coefficient_epsd",
        actual.f_coefficient_epsd,
        expected.f_coefficient_epsd);
  check("integration_method",
        actual.integration_method,
        expected.integration_method);
  check("solver_type", actual.solver_type, expected.solver_type);
  check("sparse_particle_contacts",
        actual.sparse_particle_contacts,
        expected.sparse_particle_contacts);
  check("advect_particles", actual.advect_particles, expected.advect_particles);
  check("granular_temperature_threshold",
        actual.granular_temperature_threshold,
        expected.granular_temperature_threshold);
  check("solid_fraction_threshold",
        actual.solid_fraction_threshold,
        expected.solid_fraction_threshold);
  check("disable_position_integration",
        actual.disable_position_integration,
        expected.disable_position_integration);
}

void
test_insertion_info()
{
  deallog << "--- Lagrangian::InsertionInfo<3> ---" << std::endl;
  const Parameters::Lagrangian::InsertionInfo<3> expected;
  ParameterHandler                               prm;
  Parameters::Lagrangian::InsertionInfo<3>::declare_parameters(prm);
  Parameters::Lagrangian::InsertionInfo<3> actual;
  actual.parse_parameters(prm);

  check("insertion_method", actual.insertion_method, expected.insertion_method);
  check("inserted_this_step",
        actual.inserted_this_step,
        expected.inserted_this_step);
  check("insertion_frequency",
        actual.insertion_frequency,
        expected.insertion_frequency);
  check("removing_particles_in_region",
        actual.removing_particles_in_region,
        expected.removing_particles_in_region);
  check("clear_box_point_1",
        actual.clear_box_point_1,
        expected.clear_box_point_1);
  check("clear_box_point_2",
        actual.clear_box_point_2,
        expected.clear_box_point_2);
  check("list_of_input_files",
        actual.list_of_input_files,
        expected.list_of_input_files);
  check("insertion_plane_point",
        actual.insertion_plane_point,
        expected.insertion_plane_point);
  check("insertion_plane_normal_vector",
        actual.insertion_plane_normal_vector,
        expected.insertion_plane_normal_vector);
  check("list_x", actual.list_x, expected.list_x);
  check("list_y", actual.list_y, expected.list_y);
  check("list_z", actual.list_z, expected.list_z);
  check("list_vx", actual.list_vx, expected.list_vx);
  check("list_vy", actual.list_vy, expected.list_vy);
  check("list_vz", actual.list_vz, expected.list_vz);
  check("list_wx", actual.list_wx, expected.list_wx);
  check("list_wy", actual.list_wy, expected.list_wy);
  check("list_wz", actual.list_wz, expected.list_wz);
  check("list_d", actual.list_d, expected.list_d);
  check("list_T", actual.list_T, expected.list_T);
  // parse_parameters() has to overwrite the in-class default sequence rather
  // than append to it; appending leaves the parsed axes past index dim - 1,
  // where no consumer of direction_sequence ever reads them.
  check("direction_sequence",
        actual.direction_sequence,
        expected.direction_sequence);
  check("insertion_box_point_1",
        actual.insertion_box_point_1,
        expected.insertion_box_point_1);
  check("insertion_box_point_2",
        actual.insertion_box_point_2,
        expected.insertion_box_point_2);
  check("distance_threshold",
        actual.distance_threshold,
        expected.distance_threshold);
  check("insertion_maximum_offset",
        actual.insertion_maximum_offset,
        expected.insertion_maximum_offset);
  check("seed_for_insertion",
        actual.seed_for_insertion,
        expected.seed_for_insertion);
  check("initial_vel", actual.initial_vel, expected.initial_vel);
  check("initial_omega", actual.initial_omega, expected.initial_omega);
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
  test_cfddem();
  test_multiphysics();
  test_lagrangian_model_parameters();
  test_insertion_info();
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
