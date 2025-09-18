// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_heat_transfer_h
#define lethe_particle_heat_transfer_h

#include <dem/contact_type.h>



/**
 * @brief Calculate the corrected contact radius for a particle-particle or
 * particle-wall pair in contact.
 *
 * @param[in] effective_radius Effective radius of the pair in contact.
 * @param[in] effective_youngs_modulus Effective young's modulus of the pair in
 * contact.
 * @param[in] effective_real_youngs_modulus Effective real young's modulus of
 * the pair in contact.
 * @param[in] normal_force_norm Norm of the normal contact force for the pair in
 * contact.
 * @return Corrected contact radius of the pair in contact.
 */
double
calculate_corrected_contact_radius(const double effective_radius,
                                   const double effective_youngs_modulus,
                                   const double effective_real_youngs_modulus,
                                   const double normal_force_norm);

/**
 * @brief Calculate the macrocontact thermal resistance for a particle-particle or
 * particle-wall pair in contact.
 *
 * @param[in] harmonic_conductivity Harmonic mean of the thermal
 * conductivities for a particle-particle contact or twice the particle thermal
 * conductivity for a particle-wall contact.
 * @param[in] contact_radius Contact radius of the pair in contact.
 * @return Macrocontact thermal resistance for the pair in contact.
 */
double
calculate_macrocontact_resistance(const double harmonic_conductivity,
                                  const double contact_radius);

/**
 * @brief Calculate the microcontact thermal resistance for a particle-particle or
 * particle-wall pair in contact.
 *
 * @param[in] equivalent_surface_slope Equivalent surface slope of the pair in
 * contact.
 * @param[in] equivalent_surface_roughness Equivalent surface roughness of the
 * pair in contact.
 * @param[in] effective_microhardness Effective microhardness of the pair in
 * contact.
 * @param[in] contact_radius_squared Squared contact radius of the pair in
 * contact.
 * @param[in] harmonic_conductivity Harmonic mean of the thermal
 * conductivities of the pair in contact.
 * @param[in] maximum_pressure Maximum pressure for hertzian contacts.
 * @return Microcontact thermal resistance for the pair in contact.
 */
double
calculate_microcontact_resistance(const double equivalent_surface_slope,
                                  const double equivalent_surface_roughness,
                                  const double effective_microhardness,
                                  const double contact_radius_squared,
                                  const double harmonic_conductivity,
                                  const double maximum_pressure);

/**
 * @brief Calculate the solid macrogap thermal resistance for one particle.
 *
 * @param[in] radius Radius of the particle.
 * @param[in] thermal_conductivity Thermal conductivity of the particle.
 * @param[in] contact_radius_squared Squared contact radius of the pair in
 * contact.
 * @return Solid macrogap thermal resistance for one particle.
 */
double
calculate_solid_macrogap_resistance(const double radius,
                                    const double thermal_conductivity,
                                    const double contact_radius_squared);

/**
 * @brief Calculate the interstitial gas microgap thermal resistance for a
 * particle-particle or particle-wall pair in contact.
 *
 * @param[in] equivalent_surface_roughness Equivalent surface roughness of the
 * pair in contact.
 * @param[in] contact_radius_squared Squared contact radius of the pair in
 * contact.
 * @param[in] gas_parameter_m Gas parameter.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] maximum_pressure Maximum pressure for hertzian contacts.
 * @param[in] effective_microhardness Effective microhardness of the pair in
 * contact.
 *
 * @return Interstitial gas microgap thermal resistance for the pair in contact.
 */
double
calculate_interstitial_gas_microgap_resistance(
  const double equivalent_surface_roughness,
  const double contact_radius_squared,
  const double gas_parameter_m,
  const double thermal_conductivity_gas,
  const double maximum_pressure,
  const double effective_microhardness);

/**
 * @brief Calculate the interstitial gas macrogap thermal resistance for a
 * particle-particle or particle-wall pair in contact.
 *
 * @param[in] harmonic_radius Harmonic mean of the particles' radii for a
 * particle-particle contact, radius of the particle for a particle-wall
 * contact.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] contact_radius_squared Squared contact radius of the pair in
 * contact.
 * @param[in] gas_parameter_m Gas parameter.
 *
 * @return Interstitial gas microgap thermal resistance for the pair in contact.
 */
double
calculate_interstitial_gas_macrogap_resistance(
  const double harmonic_radius,
  const double thermal_conductivity_gas,
  const double contact_radius_squared,
  const double gas_parameter_m);

/**
 * @brief Calculate the total thermal conductance for particle-particle
 * or particle-wall contacts.
 * The resistance_macrocontact, resistance_solid_macrogap and
 * resistance_gas_macrogap differ from particle-particle contacts to
 * particle-wall contacts.
 *
 * @tparam contact_type Type of contact. Only particle-particle and
 * particle-floating-mesh contacts are accepted for now.
 *
 * @param[in] radius_one Radius of particle one.
 * @param[in] radius_two Radius of particle two for particle-particle contacts.
 * Unused for particle-wall contacts.
 * @param[in] effective_youngs_modulus Effective young's modulus of the contact
 * pair.
 * @param[in] effective_real_youngs_modulus Effective real young's modulus of
 * the contact pair.
 * @param[in] equivalent_surface_roughness Equivalent surface roughness of the
 * contact pair.
 * @param[in] equivalent_surface_slope Equivalent surface slope of the contact
 * pair.
 * @param[in] effective_microhardness Effective microhardness of the contact
 * pair.
 * @param[in] thermal_conductivity_one Thermal conductivity of particle one.
 * @param[in] thermal_conductivity_two Thermal conductivity of particle two or
 * of wall.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] gas_parameter_m Gas parameter.
 * @param[in] normal_overlap Normal overlap of the contact pair.
 * @param[in] normal_force_norm Norm of the normal contact force of the contact
 * pair.
 * @param[out] thermal_conductance Total thermal conductance for the contact
 * pair.
 */
template <ContactType contact_type>
void
calculate_contact_thermal_conductance(
  const double                  radius_one,
  [[maybe_unused]] const double radius_two,
  const double                  effective_youngs_modulus,
  const double                  effective_real_youngs_modulus,
  const double                  equivalent_surface_roughness,
  const double                  equivalent_surface_slope,
  const double                  effective_microhardness,
  const double                  thermal_conductivity_one,
  const double                  thermal_conductivity_two,
  const double                  thermal_conductivity_gas,
  const double                  gas_parameter_m,
  const double                  normal_overlap,
  const double                  normal_force_norm,
  double                       &thermal_conductance);

/**
 * @brief Apply the heat transfer to a local-local particle pair in contact.
 *
 * @param[in] temperature_one Temperature of particle one.
 * @param[in] temperature_two Temperature of particle two.
 * @param[in] thermal_conductance Total thermal conductance between the two
 * particles.
 * @param[in,out] particle_one_heat_transfer_rate Heat transfer rate due to
 * contact applied to particle one.
 * @param[in,out] particle_two_heat_transfer_rate Heat transfer rate due to
 * contact applied to particle two.
 */
void
apply_heat_transfer_on_local_particles(const double temperature_one,
                                       const double temperature_two,
                                       const double thermal_conductance,
                                       double &particle_one_heat_transfer_rate,
                                       double &particle_two_heat_transfer_rate);

/**
 * @brief Apply the heat transfer to a local particle. This function is used for
 * local-ghost particle pairs or particle-wall pairs in contact.
 * The heat transfer is only applied to the local particle, thus particle one
 * should be the local particle here.
 *
 * @param[in] temperature_one Temperature of particle one.
 * @param[in] temperature_two Temperature of particle two or of wall.
 * @param[in] thermal_conductance Total thermal conductance of the pair.
 * @param[in,out] particle_one_heat_transfer_rate Heat transfer rate due to
 * contact applied to particle one.
 */
void
apply_heat_transfer_on_single_local_particle(
  const double temperature_one,
  const double temperature_two,
  const double thermal_conductance,
  double      &particle_one_heat_transfer_rate);

#endif
