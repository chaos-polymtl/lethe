// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef particle_particle_heat_transfer_h
#define particle_particle_heat_transfer_h

#include <core/auxiliary_math_functions.h>

#include <boost/math/special_functions/erf.hpp>

#include <float.h>

#include <cmath>
#include <iostream>


/**
 * @brief Calculate the corrected contact radius between two particles.
 *
 * @param[in] effective_radius Effective radius of the two particles.
 * @param[in] effective_youngs_modulus Effective young's modulus of the two
 * particles.
 * @param[in] effective_real_youngs_modulus Effective real young's modulus of
 * the two particles.
 * @param[in] normal_force_norm Norm of the normal contact force between the two
 * particles.
 * @return Corrected contact radius of the two particles.
 */
double
calculate_corrected_contact_radius(const double effective_radius,
                                   const double effective_youngs_modulus,
                                   const double effective_real_youngs_modulus,
                                   const double normal_force_norm);

/**
 * @brief Calculate the macrocontact thermal resistance between two particles.
 *
 * @param[in] harmonic_particle_conductivity Harmonic mean of the two particles'
 * thermal conductivities.
 * @param[in] contact_radius Contact radius of the two particles.
 * @return Macrocontact thermal resistance between the two particles.
 */
double
calculate_macrocontact_resistance(const double harmonic_particle_conductivity,
                                  const double contact_radius);

/**
 * @brief Calculate the microcontact thermal resistance between two particles.
 *
 * @param[in] equivalent_surface_slope Equivalent value of the two particles'
 * surface slope.
 * @param[in] equivalent_surface_roughness Equivalent value of the two
 * particles' surface roughness.
 * @param[in] effective_microhardness Effective microhardness of the two
 * particles.
 * @param[in] contact_radius_squared Squared contact radius of the two
 * particles.
 * @param[in] harmonic_particle_conductivity Harmonic mean of the two particles'
 * thermal conductivities.
 * @param[in] maximum_pressure Maximum pressure for hertzian contacts.
 * @return Microcontact thermal resistance between the two particles.
 */
double
calculate_microcontact_resistance(const double equivalent_surface_slope,
                                  const double equivalent_surface_roughness,
                                  const double effective_microhardness,
                                  const double contact_radius_squared,
                                  const double harmonic_particle_conductivity,
                                  const double maximum_pressure);

/**
 * @brief Calculate the solid macrogap thermal resistance between two particles.
 *
 * @param[in] radius_one Radius of particle one.
 * @param[in] radius_two Radius of particle two.
 * @param[in] thermal_conductivity_one Thermal conductivity of particle one.
 * @param[in] thermal_conductivity_two Thermal conductivity of particle two.
 * @param[in] contact_radius_squared Squared contact radius of the two
 * particles.
 * @return Solid macrogap thermal resistance between the two particles.
 */
double
calculate_solid_macrogap_resistance(const double radius_one,
                                    const double radius_two,
                                    const double thermal_conductivity_one,
                                    const double thermal_conductivity_two,
                                    const double contact_radius_squared);

/**
 * @brief Calculate the interstitial gas microgap thermal resistance between two particles.
 *
 * @param[in] equivalent_surface_roughness Equivalent value of the two
 * particles' surface roughness.
 * @param[in] contact_radius_squared Squared contact radius of the two
 * particles.
 * @param[in] gas_parameter_m Gas parameter.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] maximum_pressure Maximum pressure for hertzian contacts.
 * @param[in] effective_microhardness Effective microhardness of the two
 * particles.
 *
 * @return Interstitial gas microgap thermal resistance between the two particles.
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
 * @brief Calculate the interstitial gas macrogap thermal resistance between two particles.
 *
 * @param[in] harmonic_radius Harmonic mean of the particles' radii.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] contact_radius_squared Squared contact radius of the two
 * particles.
 * @param[in] gas_parameter_m Gas parameter.
 *
 * @return Interstitial gas microgap thermal resistance between the two particles.
 */
double
calculate_interstitial_gas_macrogap_resistance(
  const double harmonic_radius,
  const double thermal_conductivity_gas,
  const double contact_radius_squared,
  const double gas_parameter_m);

/**
 * @brief Calculate the total thermal conductance between two particles.
 *
 * @param[in] radius_one Radius of particle one.
 * @param[in] radius_two Radius of particle two.
 * @param[in] effective_youngs_modulus Effective young's modulus of the two
 * particles.
 * @param[in] effective_real_youngs_modulus Effective real young's modulus of
 * the two particles.
 * @param[in] equivalent_surface_roughness Equivalent value of the two
 * particles' surface roughness.
 * @param[in] equivalent_surface_slope Equivalent value of the two particles'
 * surface slope.
 * @param[in] effective_microhardness Effective microhardness of the two
 * particles.
 * @param[in] thermal_conductivity_one Thermal conductivity of particle one.
 * @param[in] thermal_conductivity_two Thermal conductivity of particle two.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] gas_parameter_m Gas parameter.
 * @param[in] normal_overlap Normal overlap between the two particles.
 * @param[in] normal_force_norm Norm of the normal contact force between the two
 * particles.
 * @param[out] thermal_conductance Total thermal conductance between the two
 * particles.
 */
void
calculate_contact_thermal_conductance(
  const double radius_one,
  const double radius_two,
  const double effective_youngs_modulus,
  const double effective_real_youngs_modulus,
  const double equivalent_surface_roughness,
  const double equivalent_surface_slope,
  const double effective_microhardness,
  const double thermal_conductivity_one,
  const double thermal_conductivity_two,
  const double thermal_conductivity_gas,
  const double gas_parameter_m,
  const double normal_overlap,
  const double normal_force_norm,
  double      &thermal_conductance);

/**
 * @brief Apply the heat transfer to the local-local particle pair in contact.
 *
 * @param[in] temperature_one Temperature of particle one.
 * @param[in] temperature_two Temperature of particle two.
 * @param[in] thermal_conductance Total thermal conductance between the two
 * particles.
 * @param[in,out] particle_one_heat_transfer_rate Heat transfer rate due to contact
 * applied to particle one.
 * @param[in,out] particle_two_heat_transfer_rate Heat transfer rate due to contact
 * applied to particle two.
 *
 */
void
apply_heat_transfer_on_local_particles(const double temperature_one,
                                       const double temperature_two,
                                       const double thermal_conductance,
                                       double      &particle_one_heat_transfer_rate,
                                       double      &particle_two_heat_transfer_rate);

/**
 * @brief Apply the heat transfer to the local-ghost particle pair in contact.
 * The heat transfer is only applied to the local particle, thus particle one
 * should be the local particle here.
 *
 * @param[in] temperature_one Temperature of particle one.
 * @param[in] temperature_two Temperature of particle two.
 * @param[in] thermal_conductance Total thermal conductance between the two
 * particles.
 * @param[in,out] particle_one_heat_transfer_rate Heat transfer rate due to contact
 * applied to particle one.
 *
 */
void
apply_heat_transfer_on_single_local_particle(
  const double temperature_one,
  const double temperature_two,
  const double thermal_conductance,
  double      &particle_one_heat_transfer_rate);

#endif
