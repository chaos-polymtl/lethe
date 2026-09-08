// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the generation of the Rushton mixer mesh. The grid is only
 * implemented in 3D, so only dim = spacedim = 3 is tested. For each
 * configuration, the rotor and the stator halves are generated and their
 * number of active cells, number of vertices, volume and number of faces
 * carrying each boundary id are reported.
 *
 * The two halves are then checked against the requirements of the mortar
 * implementation, which is the reason this generator produces both of them
 * from the same argument string: the two sides of the interface must carry the
 * same number of faces (asserted when the halves are merged) and must
 * discretize the height identically, since the number of circumferential cells
 * is read from the rotor side of the interface while the axial levels are read
 * from the stator side.
 *
 * The degenerate configuration in which the disc is as thick as the blades is
 * also covered, since the disc and blade levels then collapse onto each other.
 */

// Deal.II
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_rushton_mixer.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <vector>

/**
 * @brief Generate one half of the mixer and report its main characteristics.
 *
 * @param[out] triangulation The triangulation to fill.
 * @param[in] grid_type Either rushton_rotor or rushton_stator.
 * @param[in] grid_arguments The grid arguments passed to the generator.
 */
void
generate_and_report(Triangulation<3, 3> &triangulation,
                    const std::string   &grid_type,
                    const std::string   &grid_arguments)
{
  GridRushtonMixer<3, 3> grid(grid_type, grid_arguments);
  grid.make_grid(triangulation);

  deallog << "  " << grid_type << std::endl;
  deallog << "    Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "    Number of vertices     : " << triangulation.n_vertices()
          << std::endl;

  // Use a Q2 mapping so that the volume is integrated consistently with the
  // curved CylindricalManifold attached to the mesh; a linear mapping would
  // underestimate the volume of the cells touching a cylindrical surface.
  const MappingQ<3> mapping(2);
  deallog << "    Mesh volume            : "
          << GridTools::volume(triangulation, mapping) << std::endl;

  // Count the number of faces per boundary id.
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "    Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;
}

/**
 * @brief Collect the faces of one side of the mortar interface.
 *
 * @param[in] triangulation The triangulation of one half of the mixer.
 * @param[in] interface_boundary_id The boundary id of the mortar interface in
 * that half.
 * @param[out] n_faces The number of faces at the interface.
 * @param[out] heights The distinct heights of the vertices at the interface,
 * rounded so that they can be compared between the two halves.
 */
void
collect_interface(const Triangulation<3, 3> &triangulation,
                  const types::boundary_id   interface_boundary_id,
                  unsigned int              &n_faces,
                  std::set<long long int>   &heights)
{
  n_faces = 0;
  heights.clear();

  for (const auto &cell : triangulation.active_cell_iterators())
    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
      {
        const auto face = cell->face(f);
        if (!face->at_boundary() ||
            face->boundary_id() != interface_boundary_id)
          continue;

        ++n_faces;
        for (const auto vertex_no : face->vertex_indices())
          heights.insert(static_cast<long long int>(
            std::round(face->vertex(vertex_no)[2] * 1e9)));
      }
}

void
test(const std::string &grid_arguments, const std::string &case_name)
{
  deallog << "==================================================" << std::endl;
  deallog << "Case: " << case_name << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  Triangulation<3, 3> rotor;
  Triangulation<3, 3> stator;
  generate_and_report(rotor, "rushton_rotor", grid_arguments);
  generate_and_report(stator, "rushton_stator", grid_arguments);

  // The rotor interface is boundary id 5 and the stator interface boundary
  // id 0; the shift applied when the two halves are merged is not relevant
  // here since each half is generated on its own.
  unsigned int            n_rotor_faces  = 0;
  unsigned int            n_stator_faces = 0;
  std::set<long long int> rotor_heights;
  std::set<long long int> stator_heights;
  collect_interface(rotor, 5, n_rotor_faces, rotor_heights);
  collect_interface(stator, 0, n_stator_faces, stator_heights);

  deallog << "  Mortar interface faces, rotor  : " << n_rotor_faces
          << std::endl;
  deallog << "  Mortar interface faces, stator : " << n_stator_faces
          << std::endl;
  deallog << "  Interface face counts match    : "
          << (n_rotor_faces == n_stator_faces ? "true" : "false") << std::endl;
  deallog << "  Interface heights match        : "
          << (rotor_heights == stator_heights ? "true" : "false") << std::endl;

  // Write in VTK format to separate files so that they can be opened in
  // ParaView as a visual double-check.
  GridOut       go;
  std::ofstream rotor_out("grid_rushton_mixer_rotor_" + case_name + ".vtk");
  go.write_vtk(rotor, rotor_out);
  std::ofstream stator_out("grid_rushton_mixer_stator_" + case_name + ".vtk");
  go.write_vtk(stator, stator_out);
}

int
main()
{
  try
    {
      initlog();

      // Default configuration (empty argument string -> default dimensions).
      test("", "default");

      // A coarse configuration, obtained by lowering the number of angular
      // cells and by prescribing the radial and axial cell size explicitly.
      test(
        "1 : 1 : 0.3333333333333333 : 0.16666666666666666 : 0.25 : 0.02 : 0.06666666666666667 : 0.1 : 0.3333333333333333 : 0.1 : 0.5 : 6 : 4 : 24 : 1 : 1 : 0.15",
        "coarse");

      // A disc as thick as the blades: the disc and blade levels collapse onto
      // each other, which the generation of the axial grid must handle.
      test(
        "1 : 1 : 0.3333333333333333 : 0.16666666666666666 : 0.25 : 0.06666666666666667 : 0.06666666666666667 : 0.1 : 0.3333333333333333 : 0.1 : 0.5 : 6 : 4 : 24 : 1 : 1 : 0.15",
        "thick_disc");
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
