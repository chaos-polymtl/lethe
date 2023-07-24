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
  physical_properties.number_of_fluids                = 1;
  physical_properties.number_of_solids                = 0;
  physical_properties.number_of_material_interactions = 0;
  physical_properties.fluids.resize(1);
  physical_properties.fluids[0].density_model =
    Parameters::Material::DensityModel::constant;
  physical_properties.fluids[0].specific_heat_model =
    Parameters::Material::SpecificHeatModel::constant;
  physical_properties.fluids[0].thermal_conductivity_model =
    Parameters::Material::ThermalConductivityModel::constant;

  physical_properties.fluids[0].density              = 1;
  physical_properties.fluids[0].thermal_conductivity = 3;
  physical_properties.fluids[0].specific_heat        = 2;

  PhysicalPropertiesManager physical_properties_manager;
  physical_properties_manager.initialize(physical_properties);

  std::map<field, double> dummy_fields;

  deallog << "Testing PhysicalPropertiesManager" << std::endl;
  deallog << "Density              : "
          << physical_properties_manager.get_density()->value(dummy_fields)
          << std::endl;
  deallog << "Specific heat        : "
          << physical_properties_manager.get_specific_heat()->value(
               dummy_fields)
          << std::endl;
  deallog << "Thermal conductivity : "
          << physical_properties_manager.get_thermal_conductivity()->value(
               dummy_fields)
          << std::endl;
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
