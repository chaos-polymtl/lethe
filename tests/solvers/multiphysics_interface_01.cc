/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
*/

/**
 * @brief This code tests averaging values in time with Trilinos vectors.
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
  size_of_subsections.boundary_conditions = 0;

  solver_parameters.declare(dummy_handler, size_of_subsections);
  solver_parameters.parse(dummy_handler);

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
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
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
