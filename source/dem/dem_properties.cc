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

#include <dem/dem_properties.h>

namespace DEM
{
  template <int dim>
  std::vector<std::pair<std::string, int>>
  DEMProperties<dim>::get_properties_name()
  {
    std::vector<std::pair<std::string, int>> properties(
      PropertiesIndex::n_properties);
    properties[PropertiesIndex::type]    = std::make_pair("Type", 1);
    properties[PropertiesIndex::dp]      = std::make_pair("Diameter", 1);
    properties[PropertiesIndex::rho]     = std::make_pair("Density", 1);
    properties[PropertiesIndex::v_x]     = std::make_pair("Velocity", dim);
    properties[PropertiesIndex::v_y]     = std::make_pair("Velocity", 1);
    properties[PropertiesIndex::v_z]     = std::make_pair("Velocity", 1);
    properties[PropertiesIndex::acc_x]   = std::make_pair("Acceleration", dim);
    properties[PropertiesIndex::acc_y]   = std::make_pair("Acceleration", 1);
    properties[PropertiesIndex::acc_z]   = std::make_pair("Acceleration", 1);
    properties[PropertiesIndex::force_x] = std::make_pair("Force", dim);
    properties[PropertiesIndex::force_y] = std::make_pair("Force", 1);
    properties[PropertiesIndex::force_z] = std::make_pair("Force", 1);
    properties[PropertiesIndex::omega_x] = std::make_pair("Omega", dim);
    properties[PropertiesIndex::omega_y] = std::make_pair("Omega", 1);
    properties[PropertiesIndex::omega_z] = std::make_pair("Omega", 1);
    properties[PropertiesIndex::mass]    = std::make_pair("mass", 1);
    properties[PropertiesIndex::mom_inertia] = std::make_pair("MOI", 1);
    properties[PropertiesIndex::M_x]         = std::make_pair("M", dim);
    properties[PropertiesIndex::M_y]         = std::make_pair("M", 1);
    properties[PropertiesIndex::M_z]     = std::make_pair("M", 1);
    properties[PropertiesIndex::displacement] =
      std::make_pair("Displacement", dim);

    return properties;
  }

  unsigned int
  get_number_properties()
  {
    return PropertiesIndex::n_properties;
  }

  template class DEMProperties<2>;
  template class DEMProperties<3>;
} // namespace DEM
