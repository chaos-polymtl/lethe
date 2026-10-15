// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the registration of the solid phase (Euler-Euler
 * model) auxiliary physics in the MultiphysicsInterface. It checks that:
 * - the solid phase is reported as active when driven by a VANS fluid solver;
 * - the solid phase DoFHandler is shared through the interface and holds
 *   dim + 1 components (solid velocity and solid volume fraction);
 * - enabling the solid phase with a standard Navier-Stokes solver throws.
 */

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <core/multiphysics.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <exception>
#include <memory>
#include <vector>

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
  GridGenerator::hyper_cube(*tria, 0., 1.);
  tria->refine_global(2);

  SimulationParameters<dim>     solver_parameters;
  ParameterHandler              dummy_handler;
  Parameters::SizeOfSubsections size_of_subsections;
  size_of_subsections.boundary_conditions = 1;
  size_of_subsections.manifolds           = 1;

  solver_parameters.declare(dummy_handler, size_of_subsections);
  solver_parameters.parse(dummy_handler);

  solver_parameters.multiphysics.fluid_dynamics = true;
  solver_parameters.multiphysics.solid_phase    = true;

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(
      solver_parameters.simulation_control);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // Under VANS, the solid phase must be active and share its fields
  {
    MultiphysicsInterface<dim> multiphysics(solver_parameters,
                                            tria,
                                            simulation_control,
                                            pcout,
                                            FluidDynamicsFormulation::VANS);
    const std::vector<PhysicsID> active_physics =
      multiphysics.get_active_physics();

    deallog << "Active physics under VANS (expected: fluid, solid phase)"
            << std::endl;
    for (const auto &iphys : active_physics)
      deallog << int(iphys) << std::endl;

    multiphysics.setup_dofs();
    multiphysics.set_initial_conditions();

    const DoFHandler<dim> &solid_dof_handler =
      multiphysics.get_dof_handler(PhysicsID::solid_phase);
    deallog << "Solid phase components: "
            << solid_dof_handler.get_fe().n_components() << std::endl;
    deallog << "Solid phase DoFs: " << solid_dof_handler.n_dofs()
            << std::endl;
  }

  // Without VANS, enabling the solid phase must throw
  try
    {
      MultiphysicsInterface<dim> multiphysics(
        solver_parameters,
        tria,
        simulation_control,
        pcout,
        FluidDynamicsFormulation::standard);
      deallog << "Solid phase did NOT throw without VANS (unexpected)"
              << std::endl;
    }
  catch (const std::exception &)
    {
      deallog << "Solid phase correctly threw without VANS" << std::endl;
    }
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