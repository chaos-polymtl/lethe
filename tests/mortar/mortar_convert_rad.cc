// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Convert rad to point and back.
 */

#include <core/mortar_coupling_manager.h>

template <int dim>
void
test()
{
  const unsigned int N      = 100;
  const double       radius = 1.0;

  for (unsigned int i = 0; i <= N; ++i)
    {
      const double rad_in = 2 * numbers::PI / N * i;

      const auto point = radius_to_point<dim>(radius, rad_in);

      const double rad_out = point_to_angle(point);

      std::cout << rad_in << " " << point << " " << rad_out << std::endl;
    }

  std::cout << std::endl;
}

int
main()
{
  test<2>();
  test<3>();
}
