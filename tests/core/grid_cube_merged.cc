// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the generation of the cube merged mesh. It combines an extruded
 * plate with hole and a cylinder. The grid is only implemented in 3D, so only
 * dim = spacedim = 3 is tested. The number of active cells, the number of
 * vertices and the mesh volume are reported.
 */

// Deal.II
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_cube_merged.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <fstream>
#include <map>

void
test(const std::string &grid_arguments)
{
  deallog << "==================================================" << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  Triangulation<3, 3>  triangulation;
  GridCubeMerged<3, 3> grid(grid_arguments);
  grid.make_grid(triangulation);

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of vertices     : " << triangulation.n_vertices()
          << std::endl;
  deallog << "Mesh volume            : " << GridTools::volume(triangulation)
          << std::endl;

  // Write in VTK format to a separate file so that it can
  // be opened in ParaView as a visual double-check.
  GridOut           go;
  const std::string vtk_filename = "grid_cube_merged.vtk";
  std::ofstream     vtk_out(vtk_filename);
  go.write_vtk(triangulation, vtk_out);
}

int
main()
{
  try
    {
      initlog();

      // Input arguments: name of the two mesh files to be merged and boundary
      // IDs corresponding to the cylindrical boundaries
      test("../mesh_plate3d.msh:../mesh_cylinder.msh: 4 : 7");
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
