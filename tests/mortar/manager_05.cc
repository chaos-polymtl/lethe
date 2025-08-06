// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief MyMortarManagerCircle: check generated points.
 */

#include <core/mortar_coupling_manager.h>

#include "./tests.h"


int
main()
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

      const MyMortarManagerCircle<dim> manager(n_subdivisions,
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

              auto p     = radius_to_point<dim>(radius[0], rad);
              p[dim - 1] = radius[1] / n_subdivisions[1] * (j + 0.5);

              const auto indices = manager.get_mortar_indices(p);

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
