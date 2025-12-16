// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the PhysicalPropertiesManager class and its capacity to instantiate the various models for physical properties
 * for the case when there is also a solid region.
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
  physical_properties.number_of_solids                = 1;
  physical_properties.number_of_material_interactions = 1;
  physical_properties.reference_temperature           = 0;


  // Generate fluid properties
  physical_properties.fluids.resize(1);
  physical_properties.fluids[0].density_model =
    Parameters::Material::DensityModel::constant;
  physical_properties.fluids[0].specific_heat_model =
    Parameters::Material::SpecificHeatModel::constant;
  physical_properties.fluids[0].thermal_conductivity_model =
    Parameters::Material::ThermalConductivityModel::constant;
  physical_properties.fluids[0].electric_conductivity_model =
    Parameters::Material::ElectricConductivityModel::constant;
  physical_properties.fluids[0].electric_permittivity_model =
    Parameters::Material::ElectricPermittivityModel::constant;
  physical_properties.fluids[0].magnetic_permeability_model =
    Parameters::Material::MagneticPermeabilityModel::constant;

  // Generate solid properties
  physical_properties.solids.resize(1);
  physical_properties.solids[0].density_model =
    Parameters::Material::DensityModel::constant;
  physical_properties.solids[0].specific_heat_model =
    Parameters::Material::SpecificHeatModel::constant;
  physical_properties.solids[0].thermal_conductivity_model =
    Parameters::Material::ThermalConductivityModel::constant;
  physical_properties.solids[0].electric_conductivity_model =
    Parameters::Material::ElectricConductivityModel::constant;
  physical_properties.solids[0].electric_permittivity_model =
    Parameters::Material::ElectricPermittivityModel::constant;
  physical_properties.solids[0].magnetic_permeability_model =
    Parameters::Material::MagneticPermeabilityModel::constant;

  // Generate fluid-solid interaction properties
  physical_properties.material_interactions.resize(1);
  physical_properties.material_interactions[0].surface_tension_model =
    Parameters::MaterialInteractions::SurfaceTensionModel::constant;
  physical_properties.material_interactions[0].mobility_cahn_hilliard_model =
    Parameters::MaterialInteractions::MobilityCahnHilliardModel::constant;

  // Fix fluid properties
  physical_properties.fluids[0].density                    = 1;
  physical_properties.fluids[0].thermal_conductivity       = 3;
  physical_properties.fluids[0].specific_heat              = 2;
  physical_properties.fluids[0].electric_conductivity      = 0.5;
  physical_properties.fluids[0].electric_permittivity_real = 1.5;
  physical_properties.fluids[0].electric_permittivity_imag = 0.1;
  physical_properties.fluids[0].magnetic_permeability_real = 2.5;
  physical_properties.fluids[0].magnetic_permeability_imag = 0.2;

  // Fix solid properties
  physical_properties.solids[0].density                    = 10;
  physical_properties.solids[0].thermal_conductivity       = 30;
  physical_properties.solids[0].specific_heat              = 20;
  physical_properties.solids[0].electric_conductivity      = 5.;
  physical_properties.solids[0].electric_permittivity_real = 15;
  physical_properties.solids[0].electric_permittivity_imag = 1;
  physical_properties.solids[0].magnetic_permeability_real = 25;
  physical_properties.solids[0].magnetic_permeability_imag = 2;

  // Fix fluid-solid interaction properties
  physical_properties.material_interactions[0]
    .surface_tension_parameters.surface_tension_coefficient = 70;
  physical_properties.material_interactions[0]
    .mobility_cahn_hilliard_parameters.mobility_cahn_hilliard_constant = 10;

  PhysicalPropertiesManager physical_properties_manager;
  physical_properties_manager.initialize(physical_properties);

  std::map<field, double> dummy_fields;

  deallog << "Testing PhysicalPropertiesManager - Fluid 0 " << std::endl;
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
  deallog << "Electric conductivity: "
          << physical_properties_manager.get_electric_conductivity()->value(
               dummy_fields)
          << std::endl;
  deallog << "Electric permittivity (real): "
          << physical_properties_manager.get_electric_permittivity_real()
               ->value(dummy_fields)
          << std::endl;
  deallog << "Electric permittivity (imag): "
          << physical_properties_manager.get_electric_permittivity_imag()
               ->value(dummy_fields)
          << std::endl;
  deallog << "Magnetic permeability (real): "
          << physical_properties_manager.get_magnetic_permeability_real()
               ->value(dummy_fields)
          << std::endl;
  deallog << "Magnetic permeability (imag): "
          << physical_properties_manager.get_magnetic_permeability_imag()
               ->value(dummy_fields)
          << std::endl;

  deallog << "Testing PhysicalPropertiesManager - Solid 0 " << std::endl;
  deallog << "Density              : "
          << physical_properties_manager.get_density(0, 1)->value(dummy_fields)
          << std::endl;
  deallog << "Specific heat        : "
          << physical_properties_manager.get_specific_heat(0, 1)->value(
               dummy_fields)
          << std::endl;
  deallog << "Thermal conductivity : "
          << physical_properties_manager.get_thermal_conductivity(0, 1)->value(
               dummy_fields)
          << std::endl;
  deallog << "Electric conductivity: "
          << physical_properties_manager.get_electric_conductivity(0, 1)->value(
               dummy_fields)
          << std::endl;
  deallog << "Electric permittivity (real): "
          << physical_properties_manager.get_electric_permittivity_real(0, 1)
               ->value(dummy_fields)
          << std::endl;
  deallog << "Electric permittivity (imag): "
          << physical_properties_manager.get_electric_permittivity_imag(0, 1)
               ->value(dummy_fields)
          << std::endl;
  deallog << "Magnetic permeability (real): "
          << physical_properties_manager.get_magnetic_permeability_real(0, 1)
               ->value(dummy_fields)
          << std::endl;
  deallog << "Magnetic permeability (imag): "
          << physical_properties_manager.get_magnetic_permeability_imag(0, 1)
               ->value(dummy_fields)
          << std::endl;

  deallog << "Testing PhysicalPropertiesManager - Material interaction 0 "
          << std::endl;
  deallog << "Surface tension      : "
          << physical_properties_manager.get_surface_tension(0)->value(
               dummy_fields)
          << std::endl;
  deallog << "Mobility     : "
          << physical_properties_manager.get_mobility_cahn_hilliard(0)->value(
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
