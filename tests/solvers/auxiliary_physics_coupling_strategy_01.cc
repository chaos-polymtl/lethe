// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks AuxiliaryPhysics::should_solve_auxiliary_physics(). At
 * the moment, this function is only implemented for TimeHarmonicMaxwell, so we
 * test the behavior of this function for TimeHarmonicMaxwell with different
 * time coupling strategies of the TimeHarmonicMaxwellCouplingStrategy enum:
 *   - none: never solve after step 0
 *   - iteration: solve every N iterations
 *   - time: solve when crossing a time multiple
 *   - threshold: solve when physical properties change by more than threshold
 */

// Deal.II includes
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/multiphysics.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>
#include <solvers/time_harmonic_maxwell.h>

// Tests
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  std::shared_ptr<parallel::distributed::Triangulation<dim>> tria =
    std::make_shared<parallel::distributed::Triangulation<dim>>(
      mpi_communicator,
      typename Triangulation<dim>::MeshSmoothing(
        Triangulation<dim>::smoothing_on_refinement |
        Triangulation<dim>::smoothing_on_coarsening));

  SimulationParameters<dim>     solver_parameters;
  ParameterHandler              dummy_handler;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 1;

  solver_parameters.declare(dummy_handler, size_of_subsections);
  solver_parameters.parse(dummy_handler);

  // Enable only fluid dynamics and electromagnetics
  solver_parameters.multiphysics.fluid_dynamics   = true;
  solver_parameters.multiphysics.electromagnetics = true;

  // Set up transient simulation control with dt = 0.1, time_end = 1.0
  solver_parameters.simulation_control.dt = 0.1;
  solver_parameters.simulation_control.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  solver_parameters.simulation_control.time_end                          = 1.0;
  solver_parameters.simulation_control.time_step_independent_of_end_time = true;
  solver_parameters.simulation_control.adapt_with_capillary_time_step_ratio =
    false;

  // -------------------------------------------------------
  // Test 1: TimeHarmonicMaxwellCouplingStrategy::none
  //         Step 1 -> true, all subsequent steps -> false
  // -------------------------------------------------------
  {
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_strategy =
      Parameters::TimeHarmonicMaxwellCouplingStrategy::none;

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    TimeHarmonicMaxwell<dim> electromagnetics_auxiliary_physics(
      nullptr, solver_parameters, tria, simulation_control);

    deallog << "--- Test: TimeHarmonicMaxwellCouplingStrategy::none ---"
            << std::endl;
    // Advance a few steps
    for (int i = 0; i < 5; ++i)
      {
        simulation_control->integrate();
        deallog << "Step " << simulation_control->get_iteration_number() << ": "
                << (electromagnetics_auxiliary_physics
                        .should_solve_auxiliary_physics() ?
                      "true" :
                      "false")
                << std::endl;
      }
  }

  // -------------------------------------------------------
  // Test 2: TimeHarmonicMaxwellCouplingStrategy::iteration with frequency 3
  //         Step 1 -> true, 2 -> false, 3 -> false,
  //         4 -> true, 5 -> false, 6 -> false, 7 -> true
  // -------------------------------------------------------
  {
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_strategy =
      Parameters::TimeHarmonicMaxwellCouplingStrategy::iteration;
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .coupling_iteration = 3;

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    TimeHarmonicMaxwell<dim> electromagnetics_auxiliary_physics(
      nullptr, solver_parameters, tria, simulation_control);

    deallog
      << "--- Test: TimeHarmonicMaxwellCouplingStrategy::iteration (every 3) ---"
      << std::endl;

    for (int i = 0; i < 7; ++i)
      {
        simulation_control->integrate();
        deallog << "Step " << simulation_control->get_iteration_number() << ": "
                << (electromagnetics_auxiliary_physics
                        .should_solve_auxiliary_physics() ?
                      "true" :
                      "false")
                << std::endl;
      }
  }

  // -------------------------------------------------------
  // Test 3: TimeHarmonicMaxwellCouplingStrategy::time with coupling_time 0.24
  //         dt = 0.1, so we cross 0.24 between step 2→3,
  //         0.48 between step 4→5, 0.72 between step 7→8, etc.
  // -------------------------------------------------------
  {
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_strategy =
      Parameters::TimeHarmonicMaxwellCouplingStrategy::time;
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .coupling_time = 0.24;

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    TimeHarmonicMaxwell<dim> electromagnetics_auxiliary_physics(
      nullptr, solver_parameters, tria, simulation_control);

    deallog
      << "--- Test: TimeHarmonicMaxwellCouplingStrategy::time (every 0.24s, dt=0.1) ---"
      << std::endl;

    for (int i = 0; i < 10; ++i)
      {
        simulation_control->integrate();
        deallog << "Step " << simulation_control->get_iteration_number()
                << " (t=" << simulation_control->get_current_time() << "): "
                << (electromagnetics_auxiliary_physics
                        .should_solve_auxiliary_physics() ?
                      "true" :
                      "false")
                << std::endl;
      }
  }

  // -------------------------------------------------------
  // Test 4: TimeHarmonicMaxwellCouplingStrategy::threshold
  //         Unlike the other three strategies, "threshold" reads the
  //         temperature DoFHandler and solution through the multiphysics
  //         interface, so it cannot be exercised with a null interface or
  //         an empty triangulation the way Tests 1-3 are. This test wires up
  //         the minimum needed for the test to actually run:
  //           - a single-cell mesh (the threshold check loops over cells)
  //           - a real MultiphysicsInterface with heat_transfer active,
  //             carrying a minimal FE_Q(1) temperature DoFHandler/solution
  //             that we register and mutate directly
  //         Electric permittivity is temperature-dependent:
  //         eps_r(T) = -0.01*T + 84.0, eps_i(T) = 0.001*T. With a 1%
  //         threshold:
  //         - Step 1 -> true  (always solve at the first step)
  //         - Step 2 -> false (temperature unchanged since step 1)
  //         - Step 3 -> true  (T: 0 C -> 300 C trigger a > 1% relative change
  //         in |eps|)
  //         - Step 4 -> false (temperature unchanged since step 3)
  // -------------------------------------------------------
  {
    // The threshold strategy loops over mesh cells, so the (until now empty)
    // shared triangulation needs at least one cell to actually run.
    GridGenerator::hyper_cube(*tria, -1., 1.);

    // Temperature-dependent electric permittivity:
    //   eps_r(T) = -0.01*T + 84.0, eps_i(T) = 0.001*T
    Parameters::PhysicalProperties physical_properties;
    physical_properties.number_of_fluids                = 1;
    physical_properties.number_of_solids                = 0;
    physical_properties.number_of_material_interactions = 0;
    physical_properties.fluids.resize(1);
    physical_properties.fluids[0].electric_permittivity_model =
      Parameters::Material::ElectricPermittivityModel::polynomial;
    physical_properties.fluids[0]
      .electric_permittivity_real_polynomial_coefficients = {-0.01, 84.0};
    physical_properties.fluids[0]
      .electric_permittivity_imag_polynomial_coefficients = {0.001, 0.0};
    physical_properties.fluids[0].magnetic_permeability_model =
      Parameters::Material::MagneticPermeabilityModel::constant;
    physical_properties.fluids[0].magnetic_permeability_real = 1.0;
    physical_properties.fluids[0].magnetic_permeability_imag = 0.0;

    PhysicalPropertiesManager physical_properties_manager;
    physical_properties_manager.initialize(physical_properties);
    solver_parameters.physical_properties_manager = physical_properties_manager;

    // heat_transfer must be active for the multiphysics interface to accept
    // a DoFHandler/solution under PhysicsID::heat_transfer below.
    solver_parameters.multiphysics.heat_transfer = true;
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_strategy =
      Parameters::TimeHarmonicMaxwellCouplingStrategy::threshold;
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .coupling_threshold = 0.01; // 1% threshold

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    // should_solve_auxiliary_physics() reads the temperature field through
    // this->multiphysics, so (unlike Tests 1-3) it cannot be null here.
    ConditionalOStream pcout(
      std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0);
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);

    // Register a minimal FE_Q(1) temperature field directly with the
    // interface: this is all should_solve_auxiliary_physics() needs to
    // evaluate the threshold, so there is no need to build/solve a real
    // HeatTransfer physics.
    auto dof_handler_temperature = std::make_shared<DoFHandler<dim>>(*tria);
    FE_Q<dim> fe_temperature(1);
    dof_handler_temperature->distribute_dofs(fe_temperature);
    multiphysics.set_dof_handler(PhysicsID::heat_transfer,
                                 dof_handler_temperature);

    // This is the vector should_solve_auxiliary_physics() will actually read
    // as the "current" temperature: we drive it directly between steps.
    auto temperature_solution = std::make_shared<GlobalVectorType>(
      dof_handler_temperature->locally_owned_dofs(), mpi_communicator);
    *temperature_solution = 0.0; // Celsius
    multiphysics.set_solution(PhysicsID::heat_transfer, temperature_solution);

    TimeHarmonicMaxwell<dim> electromagnetics_auxiliary_physics(
      &multiphysics, solver_parameters, tria, simulation_control);
    // temperature_last_solved_solution starts out default-constructed, but
    // that's fine: step 1 always returns true via the early
    // "iteration_number <= 1" return, before the threshold branch (and thus
    // temperature_last_solved_solution) is ever touched. The loop below
    // initializes it for real right after that first step.

    deallog
      << "--- Test: TimeHarmonicMaxwellCouplingStrategy::threshold (1% threshold) ---"
      << std::endl;

    for (int i = 0; i < 4; ++i)
      {
        simulation_control->integrate();

        // Jump the temperature right before step 3: at 0°C, |eps| = 84.0, at
        // 300°C, |eps| = 81.0 + i0.3 = 84.3, a 3.6% change that is > 1%
        // threshold.
        if (i == 2)
          *temperature_solution = 300.0;

        // Update the temperature solution for the current step.
        const bool solved =
          electromagnetics_auxiliary_physics.should_solve_auxiliary_physics();
        deallog << "Step " << simulation_control->get_iteration_number() << ": "
                << (solved ? "true" : "false") << std::endl;

        // Mirror what the real solver does after it solves: pull the
        // current temperature from the multiphysics interface as the new
        // reference point for future threshold checks.
        if (solved)
          electromagnetics_auxiliary_physics
            .update_material_properties_dependencies();
      }
  }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<3>();
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

  return 0;
}
