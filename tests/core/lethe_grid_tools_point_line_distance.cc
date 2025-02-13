// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Lethe
#include <core/lethe_grid_tools.h>

#include <../tests/tests.h>

void
test()
{
  Triangulation<3> triangulation;

  const Tensor<1, 3> rotation_axis({0, 0, 1});
  const double       rotation_angle = M_PI / 6;

  Tensor<1, 3> cylinder_axis({0, 0, 1});
  Point<3>     cylinder_origin({0, 0, 0});

  const Tensor<2, 3> rotation_matrix =
    Physics::Transformations::Rotations::rotation_matrix_3d(rotation_axis,
                                                            rotation_angle);

  cylinder_axis = rotation_matrix * cylinder_axis;

  const double cylinder_radius = 1.0;
  Point<3>     x_prime;
  Point<3>     x;

  for (int nz = 0; nz < 6; ++nz)
    {
      x_prime[2] = -1.0 + 0.2 * nz;
      for (int nx = 0; nx < 6; ++nx)
        {
          x_prime[0] = -1.0 + 0.2 * nx;
          x_prime[1] = std::sqrt(cylinder_radius * cylinder_radius -
                                 x_prime[0] * x_prime[0]);

          x = rotation_matrix * x_prime;

          const double distance =
            LetheGridTools::find_point_line_distance(cylinder_origin,
                                                     cylinder_axis,
                                                     x);

          deallog << "The distance is " << distance << std::endl;
        }
    }
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
