// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the generation of the impinging-jet mixer mesh. The grid is
 * only implemented in 3D, so only dim = spacedim = 3 is tested. For each
 * configuration, the number of active cells, the number of vertices, the
 * mesh volume and the number of faces carrying each boundary id are
 * reported. The discretisation is fixed internally, so only the parsed
 * geometry dimensions differ between cases: the cell/vertex/face counts are
 * therefore identical across cases while the volume changes, which confirms
 * that the geometry arguments are actually parsed. The common inlet axis is
 * also checked to lie on the z = 0 impingement plane.
 */

// Deal.II
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_impinging_jet_mixer.h>

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

  Triangulation<3, 3>         triangulation;
  GridImpingingJetMixer<3, 3> grid(grid_arguments);
  grid.make_grid(triangulation);

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of vertices     : " << triangulation.n_vertices()
          << std::endl;
  deallog << "Mesh volume            : " << GridTools::volume(triangulation)
          << std::endl;

  // Count the number of faces per boundary id.
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;

  // The two inlets share a common axis that must lie on the z = 0 plane. Check
  // this by taking the mean z of the two inlet-opening (boundary ids 0 and 1)
  // face centres, which is the height of the inlet axis.
  double       inlet_z_sum   = 0.0;
  unsigned int inlet_n_faces = 0;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary() && (cell->face(f)->boundary_id() == 0 ||
                                           cell->face(f)->boundary_id() == 1))
        {
          inlet_z_sum += cell->face(f)->center()[2];
          ++inlet_n_faces;
        }
  const double inlet_axis_z =
    (inlet_n_faces > 0) ? inlet_z_sum / inlet_n_faces : 0.0;
  deallog << "Inlet axis z (mean)    : "
          << ((std::abs(inlet_axis_z) < 1e-12) ? 0.0 : inlet_axis_z)
          << std::endl;

  // Write in VTK format to a separate file so that it can
  // be opened in ParaView as a visual double-check.
  GridOut           go;
  const std::string vtk_filename =
    "grid_impinging_jet_mixer_" + case_name + ".vtk";
  std::ofstream vtk_out(vtk_filename);
  go.write_vtk(triangulation, vtk_out);
}

int
main()
{
  try
    {
      initlog();

      // Default configuration (empty argument string -> default dimensions).
      test("", "default");

      // Explicit dimensions equal to the defaults: must reproduce the default
      // case exactly, confirming the argument string is parsed consistently.
      test("0.05:0.02:0.025:0.16:0.06:0.05:0.08:0.08", "explicit_default");

      // A wider/taller mixer: same topology (fixed discretisation) but a
      // different volume, confirming the geometry arguments are parsed.
      test("0.06:0.025:0.03:0.2:0.08:0.06:0.1:0.1", "wide");
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
