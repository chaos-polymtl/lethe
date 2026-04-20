// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the generation of the Birmingham fluidized bed mesh. The
 * grid is only implemented in 3D, so only dim = spacedim = 3 is exercised.
 * For each configuration, the number of active cells, the number of
 * vertices, the mesh volume and the number of faces carrying each
 * boundary id are reported.
 */

// Deal.II
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_birmingham_fluidized_bed.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <fstream>
#include <map>

void
test(const std::string &grid_arguments, const std::string &case_name)
{
  deallog << "==================================================" << std::endl;
  deallog << "Case: " << case_name << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  Triangulation<3, 3>              triangulation;
  GridBirminghamFluidizedBed<3, 3> grid(grid_arguments);
  grid.make_grid(triangulation);

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of vertices     : " << triangulation.n_vertices()
          << std::endl;
  deallog << "Mesh volume            : " << GridTools::volume(triangulation)
          << std::endl;

  // Count the number of faces per boundary id
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;

  // Write in VTK format to a separate file so that it can
  // be opened in ParaView as a visual double-check.
  GridOut           go;
  const std::string vtk_filename =
    "grid_birmingham_fluidized_bed_" + case_name + ".vtk";
  std::ofstream vtk_out(vtk_filename);
  go.write_vtk(triangulation, vtk_out);
}

int
main()
{
  try
    {
      initlog();

      // Default configuration: chimney enabled, zero inlet offset
      test("", "chimney_default");

      // Chimney enabled, explicit inlet offset
      test("true:0.1", "chimney_offset");

      // Chimney disabled, zero inlet offset
      test("false", "no_chimney_default");

      // Chimney disabled, non-zero inlet offset
      test("false:0.05", "no_chimney_offset");
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
