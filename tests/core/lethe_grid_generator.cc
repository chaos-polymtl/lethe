// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Lethe
#include <core/uniform_channel_with_meshed_square_prism_grid.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <array>
#include <set>
#include <string>

namespace
{
  template <int dim>
  void
  print_faces(const Triangulation<dim> &tria)
  {
    if constexpr (dim == 2)
      {
        std::set<std::array<unsigned int, 4>> faces;
        for (const auto &cell : tria.active_cell_iterators())
          {
            std::array<unsigned int, 4> face = {{cell->vertex_index(0),
                                                 cell->vertex_index(1),
                                                 cell->vertex_index(2),
                                                 cell->vertex_index(3)}};
            std::ranges::sort(face);
            faces.insert(face);
          }

        deallog << "Triangulation has : " << faces.size() << " faces."
                << std::endl;
        const auto  &vertices = tria.get_vertices();
        unsigned int face_id  = 0;
        for (const auto &face : faces)
          deallog << "  Face " << face_id++ << " has vertices : ("
                  << vertices[face[0]] << "), (" << vertices[face[1]] << "), ("
                  << vertices[face[2]] << "), (" << vertices[face[3]] << ")"
                  << std::endl;
      }
    else if constexpr (dim == 3)
      {
        std::set<std::array<unsigned int, 4>> faces;
        for (const auto &cell : tria.active_cell_iterators())
          {
            for (const auto f : cell->face_indices())
              {
                std::array<unsigned int, 4> face = {
                  {cell->face(f)->vertex_index(0),
                   cell->face(f)->vertex_index(1),
                   cell->face(f)->vertex_index(2),
                   cell->face(f)->vertex_index(3)}};
                std::ranges::sort(face);
                faces.insert(face);
              }
          }

        deallog << "Triangulation has : " << faces.size() << " faces."
                << std::endl;
        const auto  &vertices = tria.get_vertices();
        unsigned int face_id  = 0;
        for (const auto &face : faces)
          deallog << "  Face " << face_id++ << " has vertices : ("
                  << vertices[face[0]] << "), (" << vertices[face[1]] << "), ("
                  << vertices[face[2]] << "), (" << vertices[face[3]] << ")"
                  << std::endl;
      }
  }


  template <int dim>
  void
  run_test(const std::string &name,
           const std::string &grid_arguments,
           const unsigned int global_refinement = 0)
  {
    Triangulation<dim>                                triangulation;
    UniformChannelWithMeshedSquarePrismGrid<dim, dim> grid(grid_arguments);
    grid.make_grid(triangulation);

    if (global_refinement > 0)
      triangulation.refine_global(global_refinement);

    deallog << name << std::endl;
    deallog << "Active cells: " << triangulation.n_active_cells() << std::endl;
    print_faces(triangulation);
  }
} // namespace


void
test()
{
  deallog << std::setprecision(12);
  deallog << "Beginning" << std::endl;

  deallog << "Case 1: Unsupported specialization" << std::endl;
  try
    {
      Triangulation<2, 3>                           triangulation;
      UniformChannelWithMeshedSquarePrismGrid<2, 3> grid(
        "0,0:1,1:0.5,0.5:0.25:0.3");
      grid.make_grid(triangulation);
      deallog << "  Unexpected success" << std::endl;
    }
  catch (const std::exception &exc)
    {
      deallog << "  Caught expected exception: " << exc.what() << std::endl;
    }

  run_test<2>("Case 2: 2D half-unit square, no rotation, no padding",
              "0,0:1,1:0.5,0.5:0.25:0.5");

  run_test<2>("Case 3: 2D 2x1 channel, 10 deg rotation, padded",
              "0,0:2,1:1,0.5:0.125:0.25:10:1:1:1:1");

  run_test<2>("Case 4: 2D 2x1 channel, 45 deg rotation, padded, refined once",
              "0,0:2,1:1,0.5:0.125:0.25:45:1:1:1:1",
              1);

  run_test<3>("Case 5: 3D 2x1x1 channel, 80 deg rotation, padded, material ids",
              "0,0:2,1:1,0.5:0.2:0.3:80:1:1:1:1:1.0:2:true");

  deallog << "OK" << std::endl;
}


int
main()
{
  try
    {
      initlog();
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
