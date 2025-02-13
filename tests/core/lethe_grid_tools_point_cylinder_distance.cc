// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/grid/grid_out.h>

// Lethe
#include <core/lethe_grid_tools.h>

// This include is only required for 9.7 and above.
#if DEAL_II_VERSION_GTE(9, 7, 0)
#  include <deal.II/grid/cell_data.h>
#endif

#include <../tests/tests.h>

void
test()
{
  
  Triangulation<3> triangulation;
  
  const Tensor<1,3> rotation_axis({0,0,1});
  const double rotation_angle = M_PI/4;
  
  GridGenerator::subdivided_cylinder(triangulation, 3, 1, 1);
  triangulation.refine_global(3);
  
  GridTools::rotate(rotation_axis, M_PI/2 + rotation_angle, triangulation);
  
  Tensor<1,3> cylinder_axis({0,0,1});
  
  const Tensor<2, 3> rotation_matrix =
    Physics::Transformations::Rotations::rotation_matrix_3d(
      rotation_axis, rotation_angle);
      
  cylinder_axis = rotation_matrix * cylinder_axis;
  
  // LetheGridTools::find_point_line_distance(const Point<dim>     &line_origin,
  //                                        const Tensor<1, dim> &line_direction,
  //                                        const Point<dim>     &point)
  
  std::ofstream out("grid.vtu");
  GridOut       grid_out;
  grid_out.write_vtu(triangulation, out);
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