// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_properties_h
#define lethe_dem_properties_h

#include <string>
#include <vector>

namespace DEM
{
  /* This enum class is responsible for the handling the specific indexes of the
   * particles properties within the property pool and also to generate
   * the associative names for the properties
   * This is the only part where we should use a classical enum because it is
   * used as an integer within the code I think this is a temporary solution for
   * now
   */

  enum PropertiesIndex : int
  {
    type                    = 0,
    dp                      = 1,
    v_x                     = 2,
    v_y                     = 3,
    v_z                     = 4,
    omega_x                 = 5,
    omega_y                 = 6,
    omega_z                 = 7,
    fem_force_x             = 8,
    fem_force_y             = 9,
    fem_force_z             = 10,
    fem_torque_x            = 11,
    fem_torque_y            = 12,
    fem_torque_z            = 13,
    mass                    = 14,
    volumetric_contribution = 15,
    n_properties            = 16,
  };

  unsigned int
  get_number_properties();

  template <int dim>
  class DEMProperties
  {
  public:
    static std::vector<std::pair<std::string, int>>
    get_properties_name();
  };

} // namespace DEM

#endif
