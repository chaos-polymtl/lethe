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
    type         = 0,
    dp           = 1,
    v_x          = 2,
    v_y          = 3,
    v_z          = 4,
    acc_x        = 5,
    acc_y        = 6,
    acc_z        = 7,
    omega_x      = 8,
    omega_y      = 9,
    omega_z      = 10,
    mass         = 11,
    n_properties = 12,
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
