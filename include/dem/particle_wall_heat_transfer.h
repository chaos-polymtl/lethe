// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_particle_heat_transfer.h>

/**
 * @brief Calculate the total thermal conductance for particle-wall contact.
 *
 * @param[in] radius_particle Radius of the particle.
 * @param[in] effective_youngs_modulus Effective young's modulus of the wall and
 * the particle.
 * @param[in] effective_real_youngs_modulus Effective real young's modulus of
 * the wall and the particle.
 * @param[in] equivalent_surface_roughness Equivalent value of the wall and the
 * particle's surface roughness.
 * @param[in] equivalent_surface_slope Equivalent value of the wall and the
 * particle's surface slope.
 * @param[in] effective_microhardness Effective microhardness of the wall and
 * the particle.
 * @param[in] thermal_conductivity_particle Thermal conductivity of the
 * particle.
 * @param[in] thermal_conductivity_wall Thermal conductivity of the wall.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] gas_parameter_m Gas parameter.
 * @param[in] normal_overlap Normal overlap between the wall and the particle.
 * @param[in] normal_force_norm Norm of the normal contact force between the
 * wall and the particle.
 * @param[out] thermal_conductance Total thermal conductance between the wall
 * and the particle.
 */
void
calculate_particle_wall_thermal_conductance(
  const double radius_particle,
  const double effective_youngs_modulus,
  const double effective_real_youngs_modulus,
  const double equivalent_surface_roughness,
  const double equivalent_surface_slope,
  const double effective_microhardness,
  const double thermal_conductivity_particle,
  const double thermal_conductivity_wall,
  const double thermal_conductivity_gas,
  const double gas_parameter_m,
  const double normal_overlap,
  const double normal_force_norm,
  double      &thermal_conductance)
{
  const double harmonic_conductivity =
    harmonic_mean(thermal_conductivity_particle, thermal_conductivity_wall);

  // Calculation of contact radius
  const double contact_radius =
    calculate_corrected_contact_radius(radius_particle,
                                       effective_youngs_modulus,
                                       effective_real_youngs_modulus,
                                       normal_force_norm);

  // Squared contact radius, often used in resistances calculations.
  const double contact_radius_squared = contact_radius * contact_radius;

  // In the same way as the contact radius, if the effective Young's modulus is
  // underestimated, the normal overlap is overestimated so a correctional
  // factor is applied here.
  const double corrected_normal_overlap =
    normal_overlap *
    pow(effective_youngs_modulus / effective_real_youngs_modulus, 2.0 / 3.0);
  const double maximum_pressure =
    (2.0 * effective_real_youngs_modulus * corrected_normal_overlap) /
    (M_PI * contact_radius + DBL_MIN);

  // Calculation of each thermal resistance
  const double resistance_macrocontact =
    calculate_macrocontact_resistance(2 * thermal_conductivity_particle,
                                      contact_radius);

  const double resistance_microcontact =
    calculate_microcontact_resistance(equivalent_surface_slope,
                                      equivalent_surface_roughness,
                                      effective_microhardness,
                                      contact_radius_squared,
                                      harmonic_conductivity,
                                      maximum_pressure);

  const double resistance_solid_macrogap =
    calculate_solid_macrogap_resistance(radius_particle,
                                        thermal_conductivity_particle,
                                        contact_radius_squared);

  const double resistance_gas_microgap =
    calculate_interstitial_gas_microgap_resistance(equivalent_surface_roughness,
                                                   contact_radius_squared,
                                                   gas_parameter_m,
                                                   thermal_conductivity_gas,
                                                   maximum_pressure,
                                                   effective_microhardness);

  const double resistance_gas_macrogap =
    0.5 *
    calculate_interstitial_gas_macrogap_resistance(radius_particle,
                                                   thermal_conductivity_gas,
                                                   contact_radius_squared,
                                                   gas_parameter_m);

  // Calculation of the final thermal conductance (1 / total resistance)
  // conductance = 1 / resistance contact area pathway + 1 / resistance magrogap
  // pathway)
  thermal_conductance =
    1.0 / (resistance_macrocontact + 1.0 / (1.0 / resistance_microcontact +
                                            1.0 / resistance_gas_microgap)) +
    1.0 / (resistance_solid_macrogap + resistance_gas_macrogap);
}
