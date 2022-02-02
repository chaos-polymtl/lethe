/*---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * This class manages the physical property models for a simulation
 */

#ifndef lethe_physical_properties_manager_h
#define lethe_physical_properties_manager_h

#include <core/density_model.h>
#include <core/rheological_model.h>
#include <core/specific_heat_model.h>
#include <core/thermal_conductivity_model.h>

using namespace dealii;


class PhysicalPropertiesManager
{
public:
  /**
   * @brief Default constructor that initialized none of the physical properties
   */
  PhysicalPropertiesManager()
  {}

  PhysicalPropertiesManager(Parameters::PhysicalProperties physical_properties);

  void
  initialize(Parameters::PhysicalProperties physical_properties);

  std::vector<std::shared_ptr<DensityModel>>             density;
  std::vector<std::shared_ptr<SpecificHeatModel>>        specific_heat;
  std::vector<std::shared_ptr<ThermalConductivityModel>> thermal_conductivity;
  std::vector<std::shared_ptr<RheologicalModel>>         rheology;

  bool non_newtonian_flow;

private:
  unsigned int number_of_fluids;
};

#endif
