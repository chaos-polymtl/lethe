// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the should_solve_electromagnetics() method
 * of MultiphysicsInterface with different time coupling strategies:
 *   - none: never solve after step 0
 *   - iteration: solve every N iterations
 *   - time: solve when crossing a time multiple
 */

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/tria.h>

// Lethe
#include <core/multiphysics.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

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

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // -------------------------------------------------------
  // Test 1: TimeCouplingMethod::none
  //         Step 1 -> true, all subsequent steps -> false
  // -------------------------------------------------------
  {
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_method = Parameters::TimeCouplingMethod::none;

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);

    deallog << "--- Test: TimeCouplingMethod::none ---" << std::endl;
    // Advance a few steps
    for (int i = 0; i < 5; ++i)
      {
        simulation_control->integrate();
        deallog << "Step " << simulation_control->get_step_number() << ": "
                << (multiphysics.should_solve_electromagnetics() ? "true" :
                                                                   "false")
                << std::endl;
      }
  }

  // -------------------------------------------------------
  // Test 2: TimeCouplingMethod::iteration with frequency 3
  //         Step 1 -> true, 2 -> false, 3 -> false,
  //         4 -> true, 5 -> false, 6 -> false, 7 -> true
  // -------------------------------------------------------
  {
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_method = Parameters::TimeCouplingMethod::iteration;
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .coupling_iteration = 3;

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);

    deallog << "--- Test: TimeCouplingMethod::iteration (every 3) ---"
            << std::endl;

    for (int i = 0; i < 7; ++i)
      {
        simulation_control->integrate();
        deallog << "Step " << simulation_control->get_step_number() << ": "
                << (multiphysics.should_solve_electromagnetics() ? "true" :
                                                                   "false")
                << std::endl;
      }
  }

  // -------------------------------------------------------
  // Test 3: TimeCouplingMethod::time with coupling_time 0.25
  //         dt = 0.1, so we cross 0.25 between step 2→3,
  //         0.50 between step 4→5, 0.75 between step 7→8, etc.
  // -------------------------------------------------------
  {
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .time_coupling_method = Parameters::TimeCouplingMethod::time;
    solver_parameters.multiphysics.time_harmonic_maxwell_parameters
      .coupling_time = 0.24;

    std::shared_ptr<SimulationControl> simulation_control =
      std::make_shared<SimulationControlTransient>(
        solver_parameters.simulation_control);

    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);

    deallog << "--- Test: TimeCouplingMethod::time (every 0.24s, dt=0.1) ---"
            << std::endl;

    for (int i = 0; i < 11; ++i)
      {
        simulation_control->integrate();
        deallog << "Step " << simulation_control->get_step_number()
                << " (t=" << simulation_control->get_current_time() << "): "
                << (multiphysics.should_solve_electromagnetics() ? "true" :
                                                                   "false")
                << std::endl;
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
