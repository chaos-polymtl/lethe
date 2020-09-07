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

#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#ifndef Lethe_DEM_properties_h
#  define Lethe_DEM_properties_h

namespace DEM
{
  /* This enum class is reponsible for the handling the specific indexes of the
   * particles properties within the property pool and also to generate
   * the associative names for the properties
   * This is the only part where we should use a classical enum because it is
   * used as an integer within the code I think this is a temporary solution for
   * now
   */

  enum PropertiesIndex : int
  {
    id           = 0,
    type         = 1,
    dp           = 2,
    rho          = 3,
    v_x          = 4,
    v_y          = 5,
    v_z          = 6,
    acc_x        = 7,
    acc_y        = 8,
    acc_z        = 9,
    force_x      = 10,
    force_y      = 11,
    force_z      = 12,
    omega_x      = 13,
    omega_y      = 14,
    omega_z      = 15,
    mass         = 16,
    mom_inertia  = 17,
    M_x          = 18,
    M_y          = 19,
    M_z          = 20,
    n_properties = 21
  };

  unsigned int
  get_number_properties();

  template <int dim>
  class DEMProperties
  {
  public:
    std::vector<std::pair<std::string, int>>
    get_properties_name();
  };

} // namespace DEM

#endif
