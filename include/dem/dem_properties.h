/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020-
 */

#include <string>
#include <tuple>
#include <vector>

#ifndef LETHE_dem_properties_H
#define LETHE_dem_properties_H

namespace DEM {
/* This enum class is reponsible for the handling the specific indexes of the
 * particles properties within the property pool and also to generate
 * the associative names for the properties
 * This is the only part where we should use a classical enum because it is
 * used as an integer within the code I think this is a temporary solution for
 * now
 */

enum PropertiesIndex : int {
  id = 0,
  type = 1,
  dp = 2,
  rho = 3,
  x = 4, // Why are we storing this information as a property?
  y = 5, // Why are we storing this information as a property?
  z = 6, // Why are we storing this information as a property?
  v_x = 7,
  v_y = 8,
  v_z = 9,
  acc_x = 10,
  acc_y = 11,
  acc_z = 12,
  force_x = 13,
  force_y = 14,
  force_z = 15,
  omega_x = 16,
  omega_y = 17,
  omega_z = 18,
  mass = 19,
  mom_inertia = 20,
  M_x = 21,
  M_y = 22,
  M_z = 23,
  n_properties = 24
};

unsigned int get_number_properties();

template <int dim> class DEMProperties {
public:
  std::vector<std::tuple<std::string, int>> get_properties_name();
};

} // namespace DEM

#endif
