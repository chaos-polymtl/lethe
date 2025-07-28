// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief MortarManagerCircle: check aligned mesh and output points.
 */

#include <core/mortar_coupling_manager.h>

int
main()
{
  const unsigned int N                   = 100;
  const unsigned int dim                 = 2;
  const unsigned int n_subdivisions      = 10;
  const unsigned int n_quadrature_points = 3;
  const double       radius              = 1.0;

  for (unsigned int i = 0; i <= N; ++i)
    {
      const double rotate = 2 * numbers::PI / N * i;

      const MortarManagerCircle<dim> manager(n_subdivisions,
                                             radius,
                                             QGauss<dim>(n_quadrature_points),
                                             rotate);

      std::cout << rotate << " "
                << static_cast<unsigned int>(manager.is_mesh_aligned()) << " "
                << manager.get_n_total_points() << std::endl;
    }
}
