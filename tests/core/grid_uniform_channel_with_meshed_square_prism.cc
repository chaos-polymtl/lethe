// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_uniform_channel_with_meshed_square_prism.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <fstream>
#include <map>

template <int dim>
void
run_test(const std::string &case_name,
         const std::string &grid_arguments,
         const unsigned int global_refinement = 0)
{
  deallog << "==================================================" << std::endl;
  deallog << "Case: " << case_name << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  Triangulation<dim, dim>                           triangulation;
  GridUniformChannelWithMeshedSquarePrism<dim, dim> grid(grid_arguments);
  grid.make_grid(triangulation);

  if (global_refinement > 0)
    triangulation.refine_global(global_refinement);

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of vertices     : " << triangulation.n_vertices()
          << std::endl;
  deallog << "Mesh volume            : " << GridTools::volume(triangulation)
          << std::endl;

  // Count the number of faces per boundary id.
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;

  // Write one VTK file per case for visual checks.
  GridOut           go;
  const std::string vtk_filename =
    "grid_uniform_channel_with_meshed_square_prism_" + case_name + ".vtk";
  std::ofstream vtk_out(vtk_filename);
  go.write_vtk(triangulation, vtk_out);
}


int
main()
{
  try
    {
      initlog();
      run_test<2>("2D half-unit square, no rotation, no padding",
                  "0,0:1,1:0.5,0.5:0.25:0.5");

      run_test<2>("2D 2x1 channel, 10 deg rotation, padded",
                  "0,0:2,1:1,0.5:0.125:0.25:10:1:1:1:1");

      run_test<2>("2D 2x1 channel, 45 deg rotation, padded, refined once",
                  "0,0:2,1:1,0.5:0.125:0.25:45:1:1:1:1",
                  1);

      run_test<3>("3D 2x1x1 channel, 80 deg rotation, padded, boundary ids",
                  "0,0:2,1:1,0.5:0.2:0.3:80:1:1:1:1:1.0:2:true");
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
