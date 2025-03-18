// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief MortarManager: check generated points.
 */

#include <core/mortar_coupling_manager.h>

int
main()
{
  const unsigned int dim                 = 2;
  const unsigned int n_subdivisions      = 10;
  const unsigned int n_quadrature_points = 3;
  const double       radius              = 1.2;

  const double delta = 2 * numbers::PI / n_subdivisions;

  for (const double scale :
       {0.0, 0.1, 0.5, 0.9, 1.1, 2.1, n_subdivisions - 0.9})
    {
      const double rotate = delta * scale;

      const MortarManager<dim> manager(n_subdivisions,
                                       n_quadrature_points,
                                       radius,
                                       rotate);

      const auto print = [&](const double shift) {
        for (unsigned int i = 0; i < n_subdivisions; ++i)
          {
            const double rad = delta * (i + shift + 0.5);

            std::cout << rad << ": " << std::endl;

            const auto indices = manager.get_indices(rad);
            const auto weights = manager.get_weights(rad);
            const auto points  = manager.get_points(rad);
            const auto normals = manager.get_normals(rad);

            for (unsigned int i = 0; i < indices.size(); ++i)
              {
                printf("%2d: %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                       indices[i],
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

      print(0.0);

      if (scale != 0.0)
        print(scale - std::floor(scale));

      if (scale != 0.0)
        print(scale);

      std::cout << std::endl << std::endl << std::endl;
    }
}
