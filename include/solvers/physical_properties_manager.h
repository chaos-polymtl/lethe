// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physical_properties_manager_h
#define lethe_physical_properties_manager_h

#include <core/density_model.h>
#include <core/electric_conductivity_model.h>
#include <core/electric_permittivity_model.h>
#include <core/magnetic_permeability_model.h>
#include <core/mobility_cahn_hilliard_model.h>
#include <core/rheological_model.h>
#include <core/specific_heat_model.h>
#include <core/surface_tension_model.h>
#include <core/thermal_conductivity_model.h>
#include <core/thermal_expansion_model.h>
#include <core/tracer_diffusivity_model.h>
#include <core/tracer_reaction_model.h>

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

  void
  provide_simulation_control(
    std::shared_ptr<SimulationControl> &simulation_control)
  {
    // At the present moment, only SpecificHeatModel can be time-dependent
    // Consequently we only pass the SimulationControl object to the
    // SpecificHeatModel.
    for (unsigned int f = 0; f < number_of_fluids; ++f)
      specific_heat[f]->provide_simulation_control(simulation_control);
    for (unsigned int s = 0; s < number_of_solids; ++s)
      specific_heat[s + number_of_fluids]->provide_simulation_control(
        simulation_control);
  }

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

  std::shared_ptr<TracerReactionPrefactorModel>
  get_tracer_reaction_prefactor(const unsigned int fluid_id    = 0,
                                const unsigned int material_id = 0) const
  {
    return tracer_reaction_prefactor[calculate_global_id(fluid_id,
                                                         material_id)];
  }

  std::shared_ptr<SurfaceTensionModel>
  get_surface_tension(const unsigned int material_interaction_id = 0) const
  {
    return surface_tension[material_interaction_id];
  }

  /**
   * @brief Returns the mobility model for a given material_interaction_id
   * (integer)
   * @param material_interaction_id
   * @return A shared pointer to a MobilityCahnHilliardModel which inherits from
   * the InterfacePropertyModel class.
   */
  std::shared_ptr<MobilityCahnHilliardModel>
  get_mobility_cahn_hilliard(
    const unsigned int material_interaction_id = 0) const
  {
    return mobility_cahn_hilliard[material_interaction_id];
  }

  /**
   * @brief Returns the electric conductivity model for a given fluid_id and material_id
   * @param fluid_id
   * @param material_id
   * @return A shared pointer to a ElectricConductivityModel which inherits from
   * the PhysicalPropertyModel class.
   */
  std::shared_ptr<ElectricConductivityModel>
  get_electric_conductivity(const unsigned int fluid_id    = 0,
                            const unsigned int material_id = 0) const
  {
    return electric_conductivity[calculate_global_id(fluid_id, material_id)];
  }

  /**
   * @brief Returns the electric permittivity model of the real part for a given fluid_id and material_id
   * @param fluid_id
   * @param material_id
   * @return A shared pointer to a ElectricPermittivityModel which inherits from
   * the PhysicalPropertyModel class.
   */
  std::shared_ptr<ElectricPermittivityModel>
  get_electric_permittivity_real(const unsigned int fluid_id    = 0,
                                 const unsigned int material_id = 0) const
  {
    return electric_permittivity_real[calculate_global_id(fluid_id,
                                                          material_id)];
  }

  /**
   * @brief Returns the electric permittivity model of the imaginary part for a given fluid_id and material_id
   * @param fluid_id
   * @param material_id
   * @return A shared pointer to a ElectricPermittivityModel which inherits from
   * the PhysicalPropertyModel class.
   */
  std::shared_ptr<ElectricPermittivityModel>
  get_electric_permittivity_imag(const unsigned int fluid_id    = 0,
                                 const unsigned int material_id = 0) const
  {
    return electric_permittivity_imag[calculate_global_id(fluid_id,
                                                          material_id)];
  }

  /**
   * @brief Returns the magnetic permeability model of the real part for a given fluid_id and material_id
   * @param fluid_id
   * @param material_id
   * @return A shared pointer to a MagneticPermeabilityModel which inherits from
   * the PhysicalPropertyModel class.
   */
  std::shared_ptr<MagneticPermeabilityModel>
  get_magnetic_permeability_real(const unsigned int fluid_id    = 0,
                                 const unsigned int material_id = 0) const
  {
    return magnetic_permeability_real[calculate_global_id(fluid_id,
                                                          material_id)];
  }

  /**
   * @brief Returns the magnetic permeability model of the imaginary part for a given fluid_id and material_id
   * @param fluid_id
   * @param material_id
   * @return A shared pointer to a MagneticPermeabilityModel which inherits from
   * the PhysicalPropertyModel class.
   */
  std::shared_ptr<MagneticPermeabilityModel>
  get_magnetic_permeability_imag(const unsigned int fluid_id    = 0,
                                 const unsigned int material_id = 0) const
  {
    return magnetic_permeability_imag[calculate_global_id(fluid_id,
                                                          material_id)];
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

  std::vector<Parameters::PhaseChange>
  get_phase_change_parameters_vector() const
  {
    return phase_change_parameters;
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

  /**
   * @brief Returns the vector of all mobility models defined in the problem.
   * @return A vector of shared pointer, each one them pointing to a
   * MobilityCahnHilliardModel.
   */
  std::vector<std::shared_ptr<MobilityCahnHilliardModel>>
  get_mobility_cahn_hilliard_vector() const
  {
    return mobility_cahn_hilliard;
  }

  /**
   * @brief Returns the vector of all electric conductivity models defined in the problem.
   * @return A vector of shared pointer, each one them pointing to a
   * ElectricConductivityModel.
   */
  std::vector<std::shared_ptr<ElectricConductivityModel>>
  get_electric_conductivity_vector() const
  {
    return electric_conductivity;
  }

  /**
   * @brief Returns the vector of all real part electric permittivity models defined in the problem.
   * @return A vector of shared pointer, each one them pointing to a
   * ElectricPermittivityModel.
   */
  std::vector<std::shared_ptr<ElectricPermittivityModel>>
  get_electric_permittivity_real_vector() const
  {
    return electric_permittivity_real;
  }

  /**
   * @brief Returns the vector of all imaginary part electric permittivity models defined in the problem.
   * @return A vector of shared pointer, each one them pointing to a
   * ElectricPermittivityModel.
   */
  std::vector<std::shared_ptr<ElectricPermittivityModel>>
  get_electric_permittivity_imag_vector() const
  {
    return electric_permittivity_imag;
  }

  /**
   * @brief Returns the vector of all real part magnetic permeability models defined in the problem.
   * @return A vector of shared pointer, each one them pointing to a
   * ElectricPermittivityModel.
   */
  std::vector<std::shared_ptr<MagneticPermeabilityModel>>
  get_magnetic_permeability_real_vector() const
  {
    return magnetic_permeability_real;
  }

  /**
   * @brief Returns the vector of all imaginary part magnetic permeability models defined in the problem.
   * @return A vector of shared pointer, each one them pointing to a
   * ElectricPermittivityModel.
   */
  std::vector<std::shared_ptr<MagneticPermeabilityModel>>
  get_magnetic_permeability_imag_vector() const
  {
    return magnetic_permeability_imag;
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

  Parameters::PhysicalProperties
  get_physical_properties_parameters() const
  {
    return physical_properties_parameters;
  }

  bool
  has_phase_change() const
  {
    return phase_change;
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

  double
  get_reference_temperature() const
  {
    return reference_temperature;
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
  std::vector<std::shared_ptr<TracerReactionPrefactorModel>>
                                                    tracer_reaction_prefactor;
  std::vector<double>                               tracer_reaction_order;
  std::vector<std::shared_ptr<SurfaceTensionModel>> surface_tension;
  std::vector<std::shared_ptr<MobilityCahnHilliardModel>>
                                       mobility_cahn_hilliard;
  std::vector<Parameters::PhaseChange> phase_change_parameters;
  std::vector<std::shared_ptr<ElectricConductivityModel>> electric_conductivity;
  std::vector<std::shared_ptr<ElectricPermittivityModel>>
    electric_permittivity_real;
  std::vector<std::shared_ptr<ElectricPermittivityModel>>
    electric_permittivity_imag;
  std::vector<std::shared_ptr<MagneticPermeabilityModel>>
    magnetic_permeability_real;
  std::vector<std::shared_ptr<MagneticPermeabilityModel>>
    magnetic_permeability_imag;

  std::map<field, bool> required_fields;

  bool non_newtonian_flow;
  bool constant_density;
  bool constant_surface_tension;
  bool phase_change;
  /*
   * Reference temperature used for the calculation of all physical properties
   * of all materials. Currently, this is only used for the thermal expansion
   * model.
   */
  double reference_temperature;

  // Internal copy of the parameters used to build the manager
  Parameters::PhysicalProperties physical_properties_parameters;

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
