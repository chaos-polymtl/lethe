// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Lethe
#include <core/lethe_grid_tools.h>

#include <../tests/tests.h>

void
test()
{
  /* This test checks the minimum distance computation between a point and a
  line by the method LetheGridTools::find_point_line_distance(). We consider a
  cylinder in 3D described by the equation:

    x'^2 + y'^2 = R^2 for any z'

  The cylinder is then rotated of an angle theta around the z-axis using a
  rotation matrix M_rot,z:

   (x,y,z) = M_rot,z(theta)*(x',y',z')

  Then, the distance between each point (x,y,z) and the cylinder axis is
  computed using the method LetheGridTools::find_point_line_distance(). It
  should always return the radius R of the cylinder, set to 1 in this test.
   */

  Triangulation<3> triangulation;

  const Tensor<1, 3> rotation_axis({0, 0, 1});
  const double       rotation_angle = M_PI / 6;

  // Define the cylinder
  Tensor<1, 3>   cylinder_axis({0, 0, 1});
  const Point<3> cylinder_origin({0, 0, 0});
  const double   cylinder_radius = 1.0;

  // Get the rotation matrix
  const Tensor<2, 3> rotation_matrix =
    Physics::Transformations::Rotations::rotation_matrix_3d(rotation_axis,
                                                            rotation_angle);
  // Rotate the cylinder axis
  cylinder_axis = rotation_matrix * cylinder_axis;

  // Point in the cylinder frame of reference
  Point<3> x_prime;

  // Point in the x,y,z frame of reference
  Point<3> x;

  // Loop to create planes along the z' axis from -1 to 1
  for (int nz = 0; nz < 6; ++nz)
    {
      x_prime[2] = -1.0 + 0.2 * nz;

      // Loop to create x' and y' coordinates
      for (int nx = 0; nx < 6; ++nx)
        {
          // x' coordinate
          x_prime[0] = -1.0 + 0.2 * nx;

          // Compute the y' coordinate using the equation of a disk
          x_prime[1] = std::sqrt(cylinder_radius * cylinder_radius -
                                 x_prime[0] * x_prime[0]);

          // Rotate the (x',y',z') in the (x,y,z) frame of reference
          x = rotation_matrix * x_prime;

          // Compute the distance between the cylinder_axis and the point
          // (x,y,z). It should always return 1.
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
