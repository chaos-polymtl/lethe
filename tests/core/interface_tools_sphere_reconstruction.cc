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
  /* This test checks the interface reconstruction of the level 0.1 of a
  level-set field using the InterfaceTools::reconstruct_interface function. The
  level-set field of interest is the one describing a sphere. The reconstruction
  is computed for 3 mesh refinements and the test checks the error on the area
  of the reconstructed surface triangulation and the convergence rate of the
  method (formally 2). The area is computed using dealii GridTools::volume()
  method.
  */
  Triangulation<3> triangulation;

  DoFHandler<3> dof_handler;
  FE_Q<3>       fe(1);
  MappingQ<3>   mapping(1);

  const Point<3> p_0 = Point<3>({0, 0, 0});
  const Point<3> p_1 = Point<3>({1, 1, 1});

  GridGenerator::hyper_rectangle(triangulation, p_0, p_1);
  triangulation.refine_global(3);

  const Point<3> sphere_center = Point<3>({0.5, 0.5, 0.5});
  const double   sphere_radius = 0.25;
  const double   iso_level     = 0.1;

  Vector<double> error_area(3);

  // Loop for the mesh convergence study
  for (unsigned int n = 0; n < 3; n++)
    {
      dof_handler.reinit(triangulation);
      dof_handler.distribute_dofs(fe);

      Vector<double> signed_distance;

      signed_distance.reinit(dof_handler.n_dofs());

      // Set the level-set field of the sphere
      VectorTools::interpolate(
        mapping,
        dof_handler,
        Functions::SignedDistance::Sphere<3>(sphere_center, sphere_radius),
        signed_distance);

      // Initialize data structure for the interface reconstruction
      std::map<types::global_cell_index, std::vector<Point<3>>>
        interface_reconstruction_vertices;
      std::map<types::global_cell_index, std::vector<CellData<2>>>
                                        interface_reconstruction_cells;
      std::set<types::global_dof_index> intersected_dofs;

      /* Reconstruct the interface. This method returns maps containing the
      surface vertices and cells of the reconstruction in volume cell-wise
      fashion. The key of the maps is the volume cell and only the intersected
      cells are store in the maps.*/
      InterfaceTools::reconstruct_interface(mapping,
                                            dof_handler,
                                            fe,
                                            signed_distance,
                                            iso_level,
                                            interface_reconstruction_vertices,
                                            interface_reconstruction_cells,
                                            intersected_dofs);

      double area = 0.0;

      // Loop on the intersected volume cells
      for (auto &intersected_cell : interface_reconstruction_cells)
        {
          const unsigned int cell_index = intersected_cell.first;

          /* Create interface recontruction triangulation (surface
          triangulation) in the intersected volume cell */
          const std::vector<Point<3>> &surface_vertices =
            interface_reconstruction_vertices.at(cell_index);
          const std::vector<CellData<2>> &surface_cells =
            intersected_cell.second;

          Triangulation<2, 3> surface_triangulation;
          surface_triangulation.create_triangulation(surface_vertices,
                                                     surface_cells,
                                                     {});

          // Compute the area of the surface triangulation
          area += GridTools::volume(surface_triangulation);
        }

      // Compute and store the area error
      error_area[n] =
        abs(4.0 * M_PI * std::pow(sphere_radius + iso_level, 2) - area);

      deallog << "The area error for ref. lev. " << n + 3
              << " is: " << error_area[n] << std::endl;
      triangulation.refine_global(1);
    }

  const double convergence_order =
    log(error_area[2] / error_area[1]) / log(0.5);

  deallog << "The convergence is: " << convergence_order << std::endl;
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
