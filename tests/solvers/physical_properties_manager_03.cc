/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*/

/**
 * @brief This code tests the PhysicalPropertiesManager class and its capacity to instantiate the various models for physical properties
 * for the case when there are 2 fluids and 2 solids. This test exceeds the
 * actual capacity of the solver (2 fluids, 1 solid).
 */

// Lethe
#include <core/parameters.h>
#include <core/physical_property_model.h>

#include <solvers/physical_properties_manager.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  // Create a physical property manager
  Parameters::PhysicalProperties physical_properties;
  physical_properties.number_of_fluids                = 2;
  physical_properties.number_of_solids                = 2;
  physical_properties.number_of_material_interactions = 5;

  physical_properties.fluids.resize(2);
  physical_properties.solids.resize(2);
  physical_properties.material_interactions.resize(5);

  for (int i = 0; i < 2; ++i)
    {
      // Generate fluid properties
      physical_properties.fluids[i].density_model =
        Parameters::Material::DensityModel::constant;
      physical_properties.fluids[i].specific_heat_model =
        Parameters::Material::SpecificHeatModel::constant;
      physical_properties.fluids[i].thermal_conductivity_model =
        Parameters::Material::ThermalConductivityModel::constant;

      physical_properties.fluids[i].density              = 1 + 10 * i;
      physical_properties.fluids[i].thermal_conductivity = 3 + 10 * i;
      physical_properties.fluids[i].specific_heat        = 2 + 10 * i;


      // Generate solid properties
      physical_properties.solids[i].density_model =
        Parameters::Material::DensityModel::constant;
      physical_properties.solids[i].specific_heat_model =
        Parameters::Material::SpecificHeatModel::constant;
      physical_properties.solids[i].thermal_conductivity_model =
        Parameters::Material::ThermalConductivityModel::constant;

      physical_properties.solids[i].density              = 10 + 100 * i;
      physical_properties.solids[i].thermal_conductivity = 30 + 100 * i;
      physical_properties.solids[i].specific_heat        = 20 + 100 * i;
    }

  for (int i = 0; i < 5; ++i)
    {
      // Generate fluid-fluid and fluid-solid interaction properties
      physical_properties.material_interactions[i].surface_tension_model =
        Parameters::MaterialInteractions::SurfaceTensionModel::constant;
      physical_properties.material_interactions[i]
        .surface_tension_parameters.surface_tension_coefficient = 10 * (i + 1);
    }

  PhysicalPropertiesManager physical_properties_manager;
  physical_properties_manager.initialize(physical_properties);

  std::map<field, double> dummy_fields;

  for (int i = 0; i < 2; ++i)
    {
      deallog << "Testing PhysicalPropertiesManager - Fluid " << i << std::endl;
      deallog << "Density              : "
              << physical_properties_manager.get_density(i, 0)->value(
                   dummy_fields)
              << std::endl;
      deallog << "Specific heat        : "
              << physical_properties_manager.get_specific_heat(i, 0)->value(
                   dummy_fields)
              << std::endl;
      deallog << "Thermal conductivity : "
              << physical_properties_manager.get_thermal_conductivity(i, 0)
                   ->value(dummy_fields)
              << std::endl;

      deallog << "Testing PhysicalPropertiesManager - Solid " << i << std::endl;
      deallog << "Density              : "
              << physical_properties_manager.get_density(0, i + 1)->value(
                   dummy_fields)
              << std::endl;
      deallog << "Specific heat        : "
              << physical_properties_manager.get_specific_heat(0, i + 1)->value(
                   dummy_fields)
              << std::endl;
      deallog << "Thermal conductivity : "
              << physical_properties_manager.get_thermal_conductivity(0, i + 1)
                   ->value(dummy_fields)
              << std::endl;
    }
  for (int i = 0; i < 5; ++i)
    {
      deallog << "Testing PhysicalPropertiesManager - Material interaction "
              << i << std::endl;
      deallog << "Surface tension      : "
              << physical_properties_manager.get_surface_tension(i)->value(
                   dummy_fields)
              << std::endl;
    }
}

int
main()
{
  try
    {
      initlog();
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
