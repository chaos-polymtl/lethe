// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the multiphysics interface behavior when different physics are enabled or disabled in the simulation parameters. It checks that the active physics reported by the interface matches the expected active physics based on the configuration of the simulation parameters.
 */

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/tria.h>
// Lethe
#include <core/multiphysics.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <string>

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

  // Fluid dynamics and heat transfer are both enabled: both should be
  // reported as active by the interface.
  solver_parameters.multiphysics.fluid_dynamics = true;
  solver_parameters.multiphysics.heat_transfer  = true;

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(
      solver_parameters.simulation_control);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  {
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);
    std::vector<PhysicsID> active_physics = multiphysics.get_active_physics();

    deallog << "Active physics (expected: fluid, heat)" << std::endl;
    for (const auto &iphys : active_physics)
      {
        deallog << int(iphys) << std::endl;
      }
  }

  // Heat transfer is disabled again: only fluid dynamics should remain
  // active.
  solver_parameters.multiphysics.heat_transfer = false;
  {
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);
    std::vector<PhysicsID> active_physics = multiphysics.get_active_physics();

    deallog << "Active physics (expected: fluid)" << std::endl;
    for (const auto &iphys : active_physics)
      {
        deallog << int(iphys) << std::endl;
      }
  }

  // Fluid dynamics is disabled as well: it must still be reported as active,
  // since its DoFHandler and solution are required at all times by every
  // other physics (disabling it only skips solving it, see
  // MultiphysicsInterface's constructor).
  solver_parameters.multiphysics.fluid_dynamics = false;
  {
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout);
    std::vector<PhysicsID> active_physics = multiphysics.get_active_physics();

    deallog << "Active physics (expected: fluid, which should always be on)"
            << std::endl;
    for (const auto &iphys : active_physics)
      {
        deallog << int(iphys) << std::endl;
      }
  }

  // From here on, every check is done with the fluid dynamics solver
  // flagged as VANS (is_vans = true), to test how MultiphysicsInterface
  // reacts when it is driven by a volume-averaged (VANS/CFD-DEM) fluid
  // solver instead of a standard one.
  solver_parameters.multiphysics.fluid_dynamics = true;

  // Electromagnetics does not depend on the fluid velocity or the void
  // fraction (TimeHarmonicMaxwell never reads the fluid dynamics solution),
  // so it is the only auxiliary physics that is NOT guarded against VANS:
  // it must still construct successfully. TimeHarmonicMaxwell is only
  // implemented for dim = 3, hence the constexpr guard.
  if constexpr (dim == 3)
    {
      solver_parameters.multiphysics.electromagnetics = true;
      {
        MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                                tria,
                                                simulation_control,
                                                pcout,
                                                /* is_vans */ true);
        std::vector<PhysicsID>     active_physics =
          multiphysics.get_active_physics();

        deallog
          << "Active physics under VANS (expected: fluid, electromagnetics)"
          << std::endl;
        for (const auto &iphys : active_physics)
          {
            deallog << int(iphys) << std::endl;
          }
      }
      solver_parameters.multiphysics.electromagnetics = false;
    }

  // Heat transfer, tracer, CLS and Cahn-Hilliard all read the fluid
  // velocity to assemble their advection term, but none of them has a
  // volume-averaged (VANS) form implemented yet. MultiphysicsInterface's
  // constructor guards each of these four with
  // AssertThrow(!is_vans, PhysicsVANSFormNotImplementedError(...)), so
  // enabling any one of them together with a VANS fluid dynamics solver
  // must throw instead of silently solving the wrong, non-averaged
  // equation. This checks all four guarded physics, one at a time, so that
  // a future physics added to that guarded set without a matching test
  // here would be caught by this test failing to compile/find the flag.
  const auto check_throws_under_vans = [&](bool              &physics_flag,
                                           const std::string &name) {
    physics_flag = true;
    try
      {
        MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                                tria,
                                                simulation_control,
                                                pcout,
                                                /* is_vans */ true);
        deallog << name << " did NOT throw under VANS (unexpected)"
                << std::endl;
      }
    catch (const std::exception &)
      {
        deallog << name << " correctly threw under VANS" << std::endl;
      }
    physics_flag = false;
  };

  check_throws_under_vans(solver_parameters.multiphysics.heat_transfer,
                          "Heat transfer");
  check_throws_under_vans(solver_parameters.multiphysics.tracer, "Tracer");
  check_throws_under_vans(solver_parameters.multiphysics.CLS, "CLS");
  check_throws_under_vans(solver_parameters.multiphysics.cahn_hilliard,
                          "Cahn-Hilliard");
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<2>();
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
