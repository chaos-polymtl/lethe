// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/base/function_signed_distance.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/interface_tools.h>

#include <../tests/tests.h>

void
test()
{
  /* This test checks the computation of the surface area of the level 0.1 of
  a level-set field and the volume enclosed by it using the 
  InterfaceTools::compute_surface_and_volume function. The level-set field of 
  interest is the one describing a sphere. The surface and volume are computed 
  for 3 mesh refinements and the test checks the error on these metriics and the 
  convergence rate of the method (formally 2).
  */
  Triangulation<3> triangulation;
  DoFHandler<3>    dof_handler;
  FE_Q<3>          fe(1);
  MappingQ<3>      mapping(1);

  const Point<3> p_0 = Point<3>({0, 0, 0});
  const Point<3> p_1 = Point<3>({1, 1, 1});

  GridGenerator::hyper_rectangle(triangulation, p_0, p_1);
  triangulation.refine_global(3);

  const Point<3> sphere_center = Point<3>({0.5, 0.5, 0.5});
  const double   sphere_radius = 0.25;
  const double   iso_level     = 0.1;

  Vector<double> error_volume(3);
  Vector<double> error_surface(3);

  // Loop for the mesh convergence study
  for (unsigned int n = 0; n < 3; n++)
    {
      dof_handler.reinit(triangulation);
      dof_handler.distribute_dofs(fe);

      dealii::TrilinosWrappers::MPI::Vector signed_distance;

      signed_distance.reinit(dof_handler.locally_owned_dofs(),
                                             triangulation.get_communicator());

      // Set the level-set field of the sphere
      VectorTools::interpolate(
        mapping,
        dof_handler,
        Functions::SignedDistance::Sphere<3>(sphere_center, sphere_radius),
        signed_distance);

      // Compute the surface and volume of the sphere with the NonMatching FEValues
      double volume, surface;

      std::tie(volume, surface) =
        InterfaceTools::compute_surface_and_volume(dof_handler,
                                       fe,
                                       signed_distance,
                                       iso_level,
                                       triangulation.get_communicator());
                                        
      // Analytical surface and volume
      const double analytical_surface = 4.0 * M_PI * std::pow(sphere_radius + iso_level, 2);
      const double analytical_volume = 4.0 * M_PI * std::pow(sphere_radius + iso_level, 3)/3.0;

      // Compute and store the surface and volume errors
      error_surface[n] =
        abs(analytical_surface - surface);
    error_volume[n] =
        abs(analytical_volume - volume);
        
      deallog << "The surface error for ref. lev. " << n + 3
              << " is: " << error_surface[n] << std::endl;
      deallog << "The volume error for ref. lev. " << n + 3
              << " is: " << error_volume[n] << std::endl;
      triangulation.refine_global(1);
    }

  // Compute the rates of convergence
  const double convergence_order_surface =
    log(error_surface[2] / error_surface[1]) / log(0.5);
  const double convergence_order_volume =
    log(error_volume[2] / error_volume[1]) / log(0.5);

    deallog << "The convergence rate for the surface is: " << convergence_order_surface << std::endl;

  deallog << "The convergence rate for the volume is: " << convergence_order_volume << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test();
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
