// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <type_traits>

namespace DEM
{
  template <int dim, typename PropertiesIndex>
  std::vector<std::pair<std::string, int>>
  ParticleProperties<dim, PropertiesIndex>::get_properties_name()
  {
    std::vector<std::pair<std::string, int>> properties(
      PropertiesIndex::n_properties);

    properties[PropertiesIndex::type]    = std::make_pair("type", 1);
    properties[PropertiesIndex::dp]      = std::make_pair("diameter", 1);
    properties[PropertiesIndex::v_x]     = std::make_pair("velocity", 3);
    properties[PropertiesIndex::v_y]     = std::make_pair("velocity", 1);
    properties[PropertiesIndex::v_z]     = std::make_pair("velocity", 1);
    properties[PropertiesIndex::omega_x] = std::make_pair("omega", 3);
    properties[PropertiesIndex::omega_y] = std::make_pair("omega", 1);
    properties[PropertiesIndex::omega_z] = std::make_pair("omega", 1);
    properties[PropertiesIndex::mass]    = std::make_pair("mass", 1);

    if constexpr (std::is_same_v<PropertiesIndex,
                                 DEM::CFDDEMProperties::PropertiesIndex>)
      {
        properties[PropertiesIndex::fem_force_two_way_coupling_x] =
          std::make_pair("fem_force", 3);
        properties[PropertiesIndex::fem_force_two_way_coupling_y] =
          std::make_pair("fem_force", 1);
        properties[PropertiesIndex::fem_force_two_way_coupling_z] =
          std::make_pair("fem_force", 1);
        properties[PropertiesIndex::fem_force_one_way_coupling_x] =
          std::make_pair("fem_force_particle_only", 3);
        properties[PropertiesIndex::fem_force_one_way_coupling_y] =
          std::make_pair("fem_force_particle_only", 1);
        properties[PropertiesIndex::fem_force_one_way_coupling_z] =
          std::make_pair("fem_force_particle_only", 1);
        properties[PropertiesIndex::fem_drag_x] = std::make_pair("fem_drag", 3);
        properties[PropertiesIndex::fem_drag_y] = std::make_pair("fem_drag", 1);
        properties[PropertiesIndex::fem_drag_z] = std::make_pair("fem_drag", 1);
        properties[PropertiesIndex::fem_torque_x] =
          std::make_pair("fem_torque", 3);
        properties[PropertiesIndex::fem_torque_y] =
          std::make_pair("fem_torque", 1);
        properties[PropertiesIndex::fem_torque_z] =
          std::make_pair("fem_torque", 1);
        properties[PropertiesIndex::volumetric_contribution] =
          std::make_pair("volumetric_contribution", 1);
        properties[PropertiesIndex::momentum_transfer_coefficient] =
          std::make_pair("momentum_transfer_coefficient", 1);
      }

    if constexpr (std::is_same_v<PropertiesIndex,
                                 DEM::DEMMPProperties::PropertiesIndex>)
      {
        properties[PropertiesIndex::T] = std::make_pair("temperature", 1);
        properties[PropertiesIndex::specific_heat] =
          std::make_pair("specific_heat", 1);
      }
    return properties;
  }

  template class ParticleProperties<2, DEM::DEMProperties::PropertiesIndex>;
  template class ParticleProperties<2, DEM::CFDDEMProperties::PropertiesIndex>;
  template class ParticleProperties<2, DEM::DEMMPProperties::PropertiesIndex>;
  template class ParticleProperties<3, DEM::DEMProperties::PropertiesIndex>;
  template class ParticleProperties<3, DEM::CFDDEMProperties::PropertiesIndex>;
  template class ParticleProperties<3, DEM::DEMMPProperties::PropertiesIndex>;

} // namespace DEM
