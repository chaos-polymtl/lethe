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

#include <core/dem_properties.h>

namespace DEM
{
  template <int dim>
  std::vector<std::pair<std::string, int>>
  DEMProperties<dim>::get_properties_name()
  {
    std::vector<std::pair<std::string, int>> properties(
      PropertiesIndex::n_properties);
    properties[PropertiesIndex::type]        = std::make_pair("type", 1);
    properties[PropertiesIndex::dp]          = std::make_pair("diameter", 1);
    properties[PropertiesIndex::v_x]         = std::make_pair("velocity", dim);
    properties[PropertiesIndex::v_y]         = std::make_pair("velocity", 1);
    properties[PropertiesIndex::v_z]         = std::make_pair("velocity", 1);
    properties[PropertiesIndex::omega_x]     = std::make_pair("omega", dim);
    properties[PropertiesIndex::omega_y]     = std::make_pair("omega", 1);
    properties[PropertiesIndex::omega_z]     = std::make_pair("omega", 1);
    properties[PropertiesIndex::fem_force_x] = std::make_pair("fem_force", dim);
    properties[PropertiesIndex::fem_force_y] = std::make_pair("fem_force", 1);
    properties[PropertiesIndex::fem_force_z] = std::make_pair("fem_force", 1);
    properties[PropertiesIndex::fem_torque_x] =
      std::make_pair("fem_torque", dim);
    properties[PropertiesIndex::fem_torque_y] = std::make_pair("fem_torque", 1);
    properties[PropertiesIndex::fem_torque_z] = std::make_pair("fem_torque", 1);
    properties[PropertiesIndex::mass]         = std::make_pair("mass", 1);
    properties[PropertiesIndex::volumetric_contribution] =
      std::make_pair("volumetric_contribution", 1);

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
