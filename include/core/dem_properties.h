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
    mass                    = 2,
    moment_of_inertia       = 3,
    v_x                     = 4,
    v_y                     = 5,
    v_z                     = 6,
    omega_x                 = 7,
    omega_y                 = 8,
    omega_z                 = 9,
    force_x                 = 10,
    force_y                 = 11,
    force_z                 = 12,
    torque_x                = 13,
    torque_y                = 14,
    torque_z                = 15,
    fem_force_x             = 16,
    fem_force_y             = 17,
    fem_force_z             = 18,
    fem_torque_x            = 19,
    fem_torque_y            = 20,
    fem_torque_z            = 21,
    volumetric_contribution = 22,
    n_properties            = 23,
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
