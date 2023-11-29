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
#include <core/mobility_cahn_hilliard_model.h>
#include <core/rheological_model.h>
#include <core/specific_heat_model.h>
#include <core/surface_tension_model.h>
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


enum material_interactions_type
{
  fluid_fluid,
  fluid_solid
};


/** @brief Class that manages the physical properties
 * model which are required to calculate the various physical properties
 * This centralizes the place where the models are created.
 * The class can be constructed empty and initialized from parameters
 * or it can be constructed directly with a Parameters section.
 * The architecture of the PhysicalPropertiesManager is made to ensure that it
 * can store as many physical properties models as required. It first stores
 * the physical properties of the fluids, which are always assumed to lie in
 * material_id=0, and then stores the physical properties of the solid, which
 * may be in cells for which material_id>0.
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

  inline unsigned int
  get_number_of_solids() const
  {
    return number_of_solids;
  }

  inline unsigned int
  get_number_of_material_interactions() const
  {
    return number_of_material_interactions;
  }

  // Getters for the physical property models
  std::shared_ptr<DensityModel>
  get_density(const unsigned int fluid_id    = 0,
              const unsigned int material_id = 0) const
  {
    return density[calculate_global_id(fluid_id, material_id)];
  }

  std::shared_ptr<SpecificHeatModel>
  get_specific_heat(const unsigned int fluid_id    = 0,
                    const unsigned int material_id = 0) const
  {
    return specific_heat[calculate_global_id(fluid_id, material_id)];
  }

  std::shared_ptr<ThermalConductivityModel>
  get_thermal_conductivity(const unsigned int fluid_id    = 0,
                           const unsigned int material_id = 0) const
  {
    return thermal_conductivity[calculate_global_id(fluid_id, material_id)];
  }

  std::shared_ptr<RheologicalModel>
  get_rheology(const unsigned int fluid_id    = 0,
               const unsigned int material_id = 0) const
  {
    return rheology[calculate_global_id(fluid_id, material_id)];
  }

  std::shared_ptr<ThermalExpansionModel>
  get_thermal_expansion(const unsigned int fluid_id    = 0,
                        const unsigned int material_id = 0) const
  {
    return thermal_expansion[calculate_global_id(fluid_id, material_id)];
  }

  std::shared_ptr<TracerDiffusivityModel>
  get_tracer_diffusivity(const unsigned int fluid_id    = 0,
                         const unsigned int material_id = 0) const
  {
    return tracer_diffusivity[calculate_global_id(fluid_id, material_id)];
  }

  std::shared_ptr<SurfaceTensionModel>
  get_surface_tension(const unsigned int material_interaction_id = 0) const
  {
    return surface_tension[material_interaction_id];
  }


  std::shared_ptr<MobilityCahnHilliardModel>
  get_mobility_cahn_hilliard(
    const unsigned int material_interaction_id = 0) const
  {
    return mobility_cahn_hilliard[material_interaction_id];
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

  std::vector<std::shared_ptr<SurfaceTensionModel>>
  get_surface_tension_vector() const
  {
    return surface_tension;
  }


  std::vector<std::shared_ptr<MobilityCahnHilliardModel>>
  get_mobility_cahn_hilliard_vector() const
  {
    return mobility_cahn_hilliard;
  }


  double
  get_kinematic_viscosity_scale() const
  {
    return rheology[0]->get_kinematic_viscosity_scale();
  }

  double
  get_density_scale() const
  {
    return density[0]->get_density_ref();
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

  bool
  surface_tension_is_constant() const
  {
    return constant_surface_tension;
  }

  unsigned int
  get_material_interaction_id(
    const material_interactions_type material_interaction_type,
    const unsigned int               material_0_id,
    const unsigned int               material_1_id) const
  {
    if (material_interaction_type == material_interactions_type::fluid_fluid)
      {
        std::pair<unsigned int, unsigned int> fluid_fluid_material_interaction(
          material_0_id, material_1_id);
        return fluid_fluid_interactions_with_material_interaction_ids
          .find(fluid_fluid_material_interaction)
          ->second;
      }
    else if (material_interaction_type ==
             material_interactions_type::fluid_solid)
      {
        std::pair<unsigned int, unsigned int> fluid_solid_material_interaction(
          material_0_id, material_1_id);
        return fluid_solid_interactions_with_material_interaction_ids
          .find(fluid_solid_material_interaction)
          ->second;
      }
    else
      throw(std::runtime_error(
        "Invalid type of material interaction. The choices are <fluid-fluid|fluid-solid>"));
  }

private:
  void
  establish_fields_required_by_model(PhysicalPropertyModel &model);

  void
  establish_fields_required_by_model(InterfacePropertyModel &model);

  /** @brief Calculates the global id of the physical property. By default. Lethe stores all
   *  the properties of the fluids, then the properties of the solid. Fluids
   * need to have a material id of 0, then solids are consequential.
   *
   * @param fluid_id the id of the fluid (0 for single phase simulations, 0 or 1 for VOF simulations)
   *
   * @param material_id the material id of the cell. Cells with material_id=0 are fluid cells (0 or 1) and cells with material_id>0 are solid cells
   */
  unsigned int
  calculate_global_id(const unsigned int fluid_id,
                      const unsigned int material_id) const
  {
    if (material_id < 1)
      return fluid_id;
    else
      {
        return number_of_fluids + material_id - 1;
      }
  }

public:
  bool is_initialized;

private:
  std::vector<std::shared_ptr<DensityModel>>             density;
  std::vector<std::shared_ptr<SpecificHeatModel>>        specific_heat;
  std::vector<std::shared_ptr<ThermalConductivityModel>> thermal_conductivity;
  std::vector<std::shared_ptr<RheologicalModel>>         rheology;
  std::vector<std::shared_ptr<ThermalExpansionModel>>    thermal_expansion;
  std::vector<std::shared_ptr<TracerDiffusivityModel>>   tracer_diffusivity;
  std::vector<std::shared_ptr<SurfaceTensionModel>>      surface_tension;
  std::vector<std::shared_ptr<MobilityCahnHilliardModel>>
    mobility_cahn_hilliard;

  std::map<field, bool> required_fields;

  bool non_newtonian_flow;
  bool constant_density;
  bool constant_surface_tension;

  unsigned int number_of_fluids;
  unsigned int number_of_solids;
  unsigned int number_of_material_interactions;

  // For material interactions
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
    fluid_fluid_interactions_with_material_interaction_ids;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>
    fluid_solid_interactions_with_material_interaction_ids;
};

#endif
