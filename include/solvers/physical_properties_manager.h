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
#include <core/thermal_expansion_model.h>
#include <core/tracer_diffusivity_model.h>

using namespace dealii;


DeclException1(
  RequiresConstantDensity,
  std::string,
  "The following assembler or post-processing utility "
    << arg1
    << " that you are trying to use does not support the use of a non-constant density model."
       " Modifications to Lethe are required to take this into account.");


DeclException1(
  RequiresConstantViscosity,
  std::string,
  "The following assembler or post-processing utility "
    << arg1
    << " that you are trying to use does not support the use of a non-constant viscosity model. "
       "Modifications to Lethe are required to take this into account.");



/** @class The PhysicalPropertiesManager class manages the physical properties
 * model which are required to calculate the various physical properties
 * This centralizes the place where the models are created.
 * The class can be constructed empty and initialized from parameters
 * or it can be constructed directly with a Parameters section.
 *
 */

class PhysicalPropertiesManager
{
public:
  /**
   * @brief Default constructor that initialized none of the physical properties
   */
  PhysicalPropertiesManager()
    : is_initialized(false)

  {}

  /**
   * @brief PhysicalPropertiesManager Constructor that initializes the class
   * @param physical_properties Parameters for the physical properties
   */

  PhysicalPropertiesManager(Parameters::PhysicalProperties physical_properties);

  /**
   * @brief initialize Initializes the physical properties using the parameters
   * @param physical_properties Parameters for the physical properties
   */

  void
  initialize(Parameters::PhysicalProperties physical_properties);



  inline unsigned int
  get_number_of_fluids() const
  {
    return number_of_fluids;
  }

  // Getters for the physical property models
  std::shared_ptr<DensityModel>
  get_density(const unsigned int fluid_id = 0) const
  {
    return density[fluid_id];
  }

  std::shared_ptr<SpecificHeatModel>
  get_specific_heat(const unsigned int fluid_id = 0) const
  {
    return specific_heat[fluid_id];
  }

  std::shared_ptr<ThermalConductivityModel>
  get_thermal_conductivity(const unsigned int fluid_id = 0) const
  {
    return thermal_conductivity[fluid_id];
  }

  std::shared_ptr<RheologicalModel>
  get_rheology(const unsigned int fluid_id = 0) const
  {
    return rheology[fluid_id];
  }

  std::shared_ptr<ThermalExpansionModel>
  get_thermal_expansion(const unsigned int fluid_id = 0) const
  {
    return thermal_expansion[fluid_id];
  }

  std::shared_ptr<TracerDiffusivityModel>
  get_tracer_diffusivity(const unsigned int fluid_id = 0) const
  {
    return tracer_diffusivity[fluid_id];
  }

  // Vector Getters for the physical property models
  std::vector<std::shared_ptr<DensityModel>>
  get_density_vector() const
  {
    return density;
  }

  std::vector<std::shared_ptr<SpecificHeatModel>>
  get_specific_heat_vector() const
  {
    return specific_heat;
  }

  std::vector<std::shared_ptr<ThermalConductivityModel>>
  get_thermal_conductivity_vector() const
  {
    return thermal_conductivity;
  }

  std::vector<std::shared_ptr<RheologicalModel>>
  get_rheology_vector() const
  {
    return rheology;
  }

  std::vector<std::shared_ptr<ThermalExpansionModel>>
  get_thermal_expansion_vector() const
  {
    return thermal_expansion;
  }

  std::vector<std::shared_ptr<TracerDiffusivityModel>>
  get_tracer_diffusivity_vector() const
  {
    return tracer_diffusivity;
  }

  void
  set_rheology(std::shared_ptr<RheologicalModel> p_rheology,
               const unsigned int                fluid_id = 0)
  {
    rheology[fluid_id] = p_rheology;
  }

  bool
  is_non_newtonian() const
  {
    return non_newtonian_flow;
  }

  bool
  field_is_required(const field id) const
  {
    return required_fields.at(id);
  }

  bool
  density_is_constant() const
  {
    return constant_density;
  }

private:
  void
  establish_fields_required_by_model(PhysicalPropertyModel &model);


  // Temporary scaling variables. This will be deprecated once the migration is
  // well finished.
public:
  // Viscosity and density scaling used for some post-processing capabilities.
  double viscosity_scale;
  double density_scale;

  bool is_initialized;


private:
  std::vector<std::shared_ptr<DensityModel>>             density;
  std::vector<std::shared_ptr<SpecificHeatModel>>        specific_heat;
  std::vector<std::shared_ptr<ThermalConductivityModel>> thermal_conductivity;
  std::vector<std::shared_ptr<RheologicalModel>>         rheology;
  std::vector<std::shared_ptr<ThermalExpansionModel>>    thermal_expansion;
  std::vector<std::shared_ptr<TracerDiffusivityModel>>   tracer_diffusivity;

  std::map<field, bool> required_fields;

  bool non_newtonian_flow;
  bool constant_density;

  unsigned int number_of_fluids;
};

#endif
