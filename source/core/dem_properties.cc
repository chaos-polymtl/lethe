// SPDX-FileCopyrightText: Copyright (c) 2020-2023 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

namespace DEM
{
  template <int dim, SolverType solver_type>
  std::vector<std::pair<std::string, int>>
  DEMProperties<dim, solver_type>::get_properties_name()
  {
    if constexpr (solver_type == dem)
      {
        std::vector<std::pair<std::string, int>> properties(
          PropertiesIndex<dem>::n_properties);
        properties[PropertiesIndex<dem>::type] = std::make_pair("type", 1);
        properties[PropertiesIndex<dem>::dp] =
          std::make_pair("diameter", 1);
        properties[PropertiesIndex<dem>::v_x] =
          std::make_pair("velocity", dim);
        properties[PropertiesIndex<dem>::v_y] =
          std::make_pair("velocity", 1);
        properties[PropertiesIndex<dem>::v_z] =
          std::make_pair("velocity", 1);
        properties[PropertiesIndex<dem>::omega_x] =
          std::make_pair("omega", dim);
        properties[PropertiesIndex<dem>::omega_y] =
          std::make_pair("omega", 1);
        properties[PropertiesIndex<dem>::omega_z] =
          std::make_pair("omega", 1);
        properties[PropertiesIndex<dem>::mass] = std::make_pair("mass", 1);

        return properties;
      }

    if constexpr (solver_type == cfd_dem)
      {
        std::vector<std::pair<std::string, int>> properties(
          PropertiesIndex<cfd_dem>::n_properties);
        properties[PropertiesIndex<cfd_dem>::type] =
          std::make_pair("type", 1);
        properties[PropertiesIndex<cfd_dem>::dp] =
          std::make_pair("diameter", 1);
        properties[PropertiesIndex<cfd_dem>::v_x] =
          std::make_pair("velocity", dim);
        properties[PropertiesIndex<cfd_dem>::v_y] =
          std::make_pair("velocity", 1);
        properties[PropertiesIndex<cfd_dem>::v_z] =
          std::make_pair("velocity", 1);
        properties[PropertiesIndex<cfd_dem>::omega_x] =
          std::make_pair("omega", dim);
        properties[PropertiesIndex<cfd_dem>::omega_y] =
          std::make_pair("omega", 1);
        properties[PropertiesIndex<cfd_dem>::omega_z] =
          std::make_pair("omega", 1);
        properties[PropertiesIndex<cfd_dem>::fem_force_x]
          = std::make_pair("fem_force", dim);
        properties[PropertiesIndex<cfd_dem>::fem_force_y] = std::make_pair("fem_force", 1);
        properties[PropertiesIndex<cfd_dem>::fem_force_z] = std::make_pair("fem_force", 1);
        properties[PropertiesIndex<cfd_dem>::fem_torque_x] =
          std::make_pair("fem_torque", dim);
        properties[PropertiesIndex<cfd_dem>::fem_torque_y] = std::make_pair("fem_torque", 1);
        properties[PropertiesIndex<cfd_dem>::fem_torque_z] = std::make_pair("fem_torque", 1);
        properties[PropertiesIndex<cfd_dem>::mass] =
          std::make_pair("mass", 1);

        return properties;
      }
  }

  template <SolverType solver_type>
  unsigned int
  get_number_properties()
  {
    return PropertiesIndex<solver_type>::n_properties;
  }

  template class DEMProperties<2, SolverType::dem>;
  template class DEMProperties<2, SolverType::cfd_dem>;
  template class DEMProperties<3, SolverType::dem>;
  template class DEMProperties<3, SolverType::cfd_dem>;

  template unsigned int DEM::get_number_properties<DEM::SolverType::dem>();
  template unsigned int DEM::get_number_properties<DEM::SolverType::cfd_dem>();
} // namespace DEM
