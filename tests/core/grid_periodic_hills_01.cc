// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the generation of the periodic hills grid. The grid is
 * supported in both 2D and 3D and several combinations of spacing_y,
 * alpha and repetitions are exercised. For each configuration, the
 * number of active cells, the number of vertices, the mesh volume and
 * the number of faces carrying each boundary id are reported.
 */

// Deal.II
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_periodic_hills.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <fstream>
#include <map>

template <int dim>
void
test(const std::string &grid_arguments, const std::string &case_name)
{
  deallog << "==================================================" << std::endl;
  deallog << "Dimension     : " << dim << std::endl;
  deallog << "Case          : " << case_name << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  Triangulation<dim>          triangulation;
  GridPeriodicHills<dim, dim> grid(grid_arguments);
  grid.make_grid(triangulation);

  // The stretched case starts from a very coarse mesh (1 x 1 repetitions);
  // a uniform refinement is applied so that the effect of the stretching
  // transformation on the interior vertices becomes visible in the
  // reported statistics.
  if (case_name == "stretched")
    triangulation.refine_global(2);

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of vertices     : " << triangulation.n_vertices()
          << std::endl;
  deallog << "Mesh volume            : " << GridTools::volume(triangulation)
          << std::endl;

  // Count the number of faces per boundary id
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;

  // Write in VTK format to a separate file so that it can
  // be opened in ParaView as a visual double-check.
  GridOut           go;
  const std::string vtk_filename =
    "grid_periodic_hills_" + std::to_string(dim) + "d_" + case_name + ".vtk";
  std::ofstream vtk_out(vtk_filename);
  go.write_vtk(triangulation, vtk_out);
}

int
main()
{
  try
    {
      initlog();

      // 2D cases: spacing_y ; alpha ; repetitions_x ; repetitions_y
      test<2>("0;1;10;2", "default");
      test<2>("0;1.5;1;1", "stretched");

      // 3D cases: spacing_y ; alpha ; repetitions_x ; repetitions_y ;
      // repetitions_z
      test<3>("0;1;10;2;1", "default");
      test<3>("0;1.5;1;1;1", "stretched");
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
