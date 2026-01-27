// SPDX-FileCopyrightText: Copyright (c) 2024-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the VOF subequations interface.
 */

#include <core/parameters.h>
#include <core/simulation_control.h>

#include <solvers/multiphysics_interface.h>
#include <solvers/physical_properties_manager.h>
#include <solvers/simulation_parameters.h>
#include <solvers/vof_subequations_interface.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/tria.h>

#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> tria =
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

  // For VOF, since the number of fluids must be set to 2, but
  // physical_properties is a private member of SimulationParameters
  Parameters::PhysicalProperties physical_properties;
  physical_properties.declare_parameters(dummy_handler);
  physical_properties.parse_parameters(dummy_handler,
                                       solver_parameters.dimensionality);
  physical_properties.number_of_fluids = 2;
  solver_parameters.physical_properties_manager.initialize(physical_properties);

  std::shared_ptr<SimulationControl> simulation_control =
    std::make_shared<SimulationControlTransient>(
      solver_parameters.simulation_control);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // Construct VOF subequations interface object
  VOFSubequationsInterface<dim> subequations_interface(solver_parameters,
                                                       pcout,
                                                       tria,
                                                       simulation_control);

  // Phase fraction gradient and curvature L2 projection enabled
  {
    // VOF is required when enabling surface_tension_force
    solver_parameters.multiphysics.VOF = true;

    // To test with the phase fraction gradient and the curvature L2 projection
    solver_parameters.multiphysics.vof_parameters.surface_tension_force.enable =
      true;

    subequations_interface.initialize_subequations(solver_parameters,
                                                   tria,
                                                   simulation_control);

    std::vector<VOFSubequationsID> active_subequations =
      subequations_interface.get_active_subequations();

    deallog
      << "Active subequations [expected: phase_gradient_projection (0), curvature_projection (1)]"
      << std::endl;
    for (const auto &subequation_id : active_subequations)
      {
        deallog << int(subequation_id) << std::endl;
      }
  }

  // Disable phase fraction gradient and curvature L2 projection
  {
    solver_parameters.multiphysics.vof_parameters.surface_tension_force.enable =
      false;

    subequations_interface.initialize_subequations(solver_parameters,
                                                   tria,
                                                   simulation_control);

    std::vector<VOFSubequationsID> active_subequations =
      subequations_interface.get_active_subequations();

    deallog << "Active subequations [expected: No active subequation]"
            << std::endl;
    if (active_subequations.size() == 0)
      deallog << "No active subequation" << std::endl;
    else
      for (const auto &subequation_id : active_subequations)
        {
          deallog << int(subequation_id) << std::endl;
        }
  }

  // Enable algebraic interface reinitialization
  solver_parameters.multiphysics.vof_parameters.regularization_method
    .algebraic_interface_reinitialization.enable = true;
  {
    subequations_interface.initialize_subequations(solver_parameters,
                                                   tria,
                                                   simulation_control);

    std::vector<VOFSubequationsID> active_subequations =
      subequations_interface.get_active_subequations();

    deallog
      << "Active subequations [expected: phase_gradient_projection (0), curvature_projection (1), algebraic_interface_reinitialization (2)]"
      << std::endl;
    for (const auto &subequation_id : active_subequations)
      {
        deallog << int(subequation_id) << std::endl;
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
