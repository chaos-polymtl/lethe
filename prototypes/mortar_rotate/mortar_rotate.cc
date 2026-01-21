// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief MyMortarManagerCircle: check generated points.
 */

#include "./generate_grids.h"
#include "./mortar_operator.h"

void
test_01()
{
  const unsigned int              dim                 = 3;
  const std::vector<unsigned int> n_subdivisions      = {10, 3};
  const unsigned int              n_quadrature_points = 3;
  const std::vector<double>       radius              = {1.2, 0.4};

  // cell angle variation
  const double delta = 2 * numbers::PI / n_subdivisions[0];

  // rotate inner mesh using random scaling factors
  for (const double scale : {0.0})
    {
      const double rotate = delta * scale;
      std::cout << "Rotation angle: " << rotate << std::endl;

      const MortarManagerCircle<dim> manager(n_subdivisions,
                                             radius,
                                             QGauss<dim>(n_quadrature_points),
                                             rotate);

      std::cout << manager.get_n_mortars() << std::endl;
      std::cout << manager.get_n_total_mortars() << std::endl;
      std::cout << manager.get_n_points() << std::endl;
      std::cout << manager.get_n_total_points() << std::endl;
      std::cout << std::endl;

      const auto print = [&](const double shift) {
        for (unsigned int j = 0; j < n_subdivisions[1]; ++j)
          for (unsigned int i = 0; i < n_subdivisions[0]; ++i)
            {
              // center point of each cell (in radians)
              const double rad = delta * (i + shift + 0.5);

              std::cout << "Shift: " << shift << std::endl;
              std::cout << "Cell center at " << rad << ": " << std::endl;

              Point<dim> p;
              p[0] = radius[0] * std::cos(rad);
              p[1] = radius[0] * std::sin(rad);

              p[dim - 1] = radius[1] / n_subdivisions[1] * (j + 0.5);

              const auto indices = manager.get_mortar_indices(p, rotate != 0.0);

              for (const auto ii : indices)
                std::cout << ii << " ";

              std::cout << std::endl;
              std::cout << std::endl;
            }
        std::cout << std::endl;
      };

      // print generated points for aligned mesh (scale = 0.0)
      print(0.0);
    }
}

void
test_02()
{
  const unsigned int dim                 = 2;
  const unsigned int n_subdivisions      = 10;
  const unsigned int n_quadrature_points = 3;
  const double       radius              = 1.2;

  // cell angle variation
  const double delta = 2 * numbers::PI / n_subdivisions;

  // rotate inner mesh using random scaling factors
  for (const double scale :
       {0.0, 0.1, 0.5, 0.9, 1.1, 2.1, n_subdivisions - 0.9})
    {
      const double rotate = delta * scale;
      std::cout << "Rotation angle: " << rotate << std::endl;

      const MortarManagerCircle<dim> manager(n_subdivisions,
                                             radius,
                                             QGauss<dim>(n_quadrature_points),
                                             rotate);

      const auto print = [&](const double shift, const bool is_inner) {
        for (unsigned int i = 0; i < n_subdivisions; ++i)
          {
            // center point of each cell (in radians)
            const double rad = delta * (i + shift + 0.5);

            std::cout << "Shift: " << shift << std::endl;
            std::cout << "Cell center at " << rad << ": " << std::endl;

            Point<dim> p;
            p[0] = radius * std::cos(rad);
            p[1] = radius * std::sin(rad);

            const auto indices = manager.get_mortar_indices(p, is_inner);

            const auto weights = manager.get_weights(p, is_inner);
            const auto points  = manager.get_points(p, is_inner);
            const auto normals = manager.get_normals(p, is_inner);

            for (unsigned int i = 0; i < weights.size(); ++i)
              {
                printf("%2d: %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                       (i < n_quadrature_points ? indices[0] : indices[1]) *
                           n_quadrature_points +
                         (i % n_quadrature_points),
                       weights[i],
                       points[i][0],
                       points[i][1],
                       normals[i][0],
                       normals[i][1]);
              }
            std::cout << std::endl;
          }
        std::cout << std::endl;
      };

      // print generated points for aligned mesh (scale = 0.0)
      print(0.0, false);

      // print generated points for non-aligned mesh
      if (scale != 0.0)
        print(std::floor(scale), false); // outer (fixed) mesh

      if (scale != 0.0)
        print(scale, true); // inner (rotated) mesh

      std::cout << std::endl << std::endl << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      test_01();
      test_02();
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
