// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_heat_transfer.h>

double
calculate_corrected_contact_radius(const double effective_radius,
                                   const double effective_youngs_modulus,
                                   const double effective_real_youngs_modulus,
                                   const double normal_force_norm)
{
  const double contact_radius = pow((3 * normal_force_norm * effective_radius) /
                                      (4 * effective_youngs_modulus),
                                    (1.0 / 3.0));

  // In certain cases, the simulation effective Young’s modulus
  // can be chosen a few orders lower than the real effective Young’s modulus
  // to ensure simulation stability (related to Rayleigh's time).
  // As a consequence, the contact radius is overestimated, so a correctional
  // factor is applied here.
  return contact_radius *
         pow(effective_youngs_modulus / effective_real_youngs_modulus,
             1.0 / 5.0);
}

double
calculate_macrocontact_resistance(const double harmonic_conductivity,
                                  const double contact_radius)
{
  return 0.5 / (contact_radius * harmonic_conductivity + DBL_MIN);

  // G. K. Batchelor and R. W. O’Brien, “Thermal or electrical conduction
  // through a granular material,” Proc. R. Soc. Lond. A Math. Phys. Sci., vol.
  // 355, no. 1682, pp. 313–333, Jul. 1977
}

double
calculate_microcontact_resistance(const double equivalent_surface_slope,
                                  const double equivalent_surface_roughness,
                                  const double effective_microhardness,
                                  const double contact_radius_squared,
                                  const double harmonic_conductivity,
                                  const double maximum_pressure)
{
  // 1.184 = (1+0.96/2)/1.25
  return 1.184 / (M_PI * harmonic_conductivity * contact_radius_squared) *
         (equivalent_surface_roughness / equivalent_surface_slope) *
         pow(effective_microhardness / (maximum_pressure + DBL_MIN), 0.96);

  // Van Lew, J. T. (2016). On thermal characterization of breeder pebble beds
  // with microscale numerical modeling of thermofluid and pebble-pebble
  // interactions (Doctoral dissertation, University of California, Los Angeles)
  // Equations 3.25 to 3.27.
}

double
calculate_solid_macrogap_resistance(const double radius,
                                    const double thermal_conductivity,
                                    const double contact_radius_squared)
{
  // resistance = characteristic length parallel to heat flux /
  // (thermal_conductivity
  // * characteristic area perpendicular to the heat flux)

  return 0.25 * M_PI * radius /
         (M_PI * (radius * radius - contact_radius_squared) *
          thermal_conductivity);

  // C. Beaulieu, “Impact de la ségrégation granulaire sur le transfert de
  // chaleur dans un lit rotatif”, Ph.D. thesis, Polytechnique Montréal, 2020.
}

double
calculate_interstitial_gas_microgap_resistance(
  const double equivalent_surface_roughness,
  const double contact_radius_squared,
  const double gas_parameter_m,
  const double thermal_conductivity_gas,
  const double maximum_pressure,
  const double effective_microhardness)
{
  const double x_1 = 2.0 * maximum_pressure / effective_microhardness;
  const double x_2 = 0.03 * maximum_pressure / effective_microhardness;

  // If the values are out of the valid range ]0, 2[ of the erfc^-1 function,
  // this resistance has no physical meaning. Particles are too close or not
  // close enough for a microgap between them to make sense. This resistance is
  // used in parallel with another one, so putting it to infinity ignores it in
  // the total resistance.
  if (x_1 >= 2.0 || x_1 <= 0.0)
    {
      return INFINITY;
    }

  const double a_1 = boost::math::erfc_inv(x_1);
  const double a_2 = boost::math::erfc_inv(x_2) - a_1;

  // 2.82842712475 = 2*sqrt(2)
  return (2.82842712475 * equivalent_surface_roughness * a_2) /
         (M_PI * thermal_conductivity_gas * contact_radius_squared *
          std::log(abs(1 + a_2 / (a_1 + gas_parameter_m /
                                          (2.82842712475 *
                                           equivalent_surface_roughness)))));

  // M. Bahrami, M. M. Yovanovich, and J. R. Culham, “Effective thermal
  // conductivity of rough spherical packed beds,” Int. J. Heat Mass Transfer,
  // vol. 49, no. 19–20, pp. 3691–3701, Sep. 2006
}

double
calculate_interstitial_gas_macrogap_resistance(
  const double harmonic_radius,
  const double thermal_conductivity_gas,
  const double contact_radius_squared,
  const double gas_parameter_m)
{
  const double A =
    2. * sqrt(harmonic_radius * harmonic_radius - contact_radius_squared);
  const double S = 2. * harmonic_radius -
                   contact_radius_squared / harmonic_radius + gas_parameter_m;

  return 2.0 / (M_PI * thermal_conductivity_gas * (S * log(S / (S - A)) - A));

  // M. Bahrami, M. M. Yovanovich, and J. R. Culham, “Effective thermal
  // conductivity of rough spherical packed beds,” Int. J. Heat Mass Transfer,
  // vol. 49, no. 19–20, pp. 3691–3701, Sep. 2006
}

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
  double                       &thermal_conductance)
{
  const double harmonic_conductivity =
    harmonic_mean(thermal_conductivity_one, thermal_conductivity_two);
  // For particle-wall contacts, it is as if the radius of the wall is infinite.
  const double harmonic_radius =
    (contact_type == ContactType::particle_floating_mesh) ?
      2. * radius_one :
      harmonic_mean(radius_one, radius_two);

  // Calculation of contact radius
  // Hertz contact radius
  const double contact_radius =
    calculate_corrected_contact_radius(harmonic_radius * 0.5,
                                       effective_youngs_modulus,
                                       effective_real_youngs_modulus,
                                       normal_force_norm);

  // In the same way as the contact radius, if the effective Young's modulus is
  // underestimated, the normal overlap is overestimated so a correctional
  // factor is applied here.
  const double corrected_normal_overlap =
    normal_overlap *
    pow(effective_youngs_modulus / effective_real_youngs_modulus, 2.0 / 3.0);

  // Geometric contact radius
  // (from analytical contact area between overlapping spheres)
  // const double contact_radius = [&]() {
  //   if constexpr (contact_type == ContactType::particle_floating_mesh)
  //     return sqrt(corrected_normal_overlap *
  //                 (2. * radius_one - corrected_normal_overlap));
  //   else
  //     return sqrt(
  //       corrected_normal_overlap *
  //       (2. * radius_one - corrected_normal_overlap) *
  //       (2. * radius_two - corrected_normal_overlap) *
  //       (2. * radius_one + 2. * radius_two - corrected_normal_overlap) /
  //       ((radius_one + radius_two - corrected_normal_overlap) *
  //        (radius_one + radius_two - corrected_normal_overlap)));
  // }();

  // Squared contact radius, often used in resistances calculations.
  const double contact_radius_squared = contact_radius * contact_radius;

  const double maximum_pressure =
    (2.0 * effective_real_youngs_modulus * corrected_normal_overlap) /
    (M_PI * contact_radius + DBL_MIN);

  // Calculation of each thermal resistance
  // For particle-wall contacts, only the macrocontact resistance of the
  // particle is considered.
  const double resistance_macrocontact = [&]() {
    if constexpr (contact_type == ContactType::particle_floating_mesh)
      return calculate_macrocontact_resistance(2. * thermal_conductivity_one,
                                               contact_radius);
    else
      return calculate_macrocontact_resistance(harmonic_conductivity,
                                               contact_radius);
  }();

  const double resistance_microcontact =
    calculate_microcontact_resistance(equivalent_surface_slope,
                                      equivalent_surface_roughness,
                                      effective_microhardness,
                                      contact_radius_squared,
                                      harmonic_conductivity,
                                      maximum_pressure);

  double resistance_solid_macrogap =
    calculate_solid_macrogap_resistance(radius_one,
                                        thermal_conductivity_one,
                                        contact_radius_squared);
  // For particle-wall contacts, only the solid macrogap resistance of the
  // particle is considered.
  if constexpr (contact_type != ContactType::particle_floating_mesh)
    resistance_solid_macrogap +=
      calculate_solid_macrogap_resistance(radius_two,
                                          thermal_conductivity_two,
                                          contact_radius_squared);

  const double resistance_gas_microgap =
    calculate_interstitial_gas_microgap_resistance(equivalent_surface_roughness,
                                                   contact_radius_squared,
                                                   gas_parameter_m,
                                                   thermal_conductivity_gas,
                                                   maximum_pressure,
                                                   effective_microhardness);

  double resistance_gas_macrogap =
    calculate_interstitial_gas_macrogap_resistance(harmonic_radius,
                                                   thermal_conductivity_gas,
                                                   contact_radius_squared,
                                                   gas_parameter_m);

  // For particle-wall contacts, there is only a macrogap around the particle
  // and not the wall, so only half of the resistance is kept.
  if constexpr (contact_type == ContactType::particle_floating_mesh)
    resistance_gas_macrogap *= 0.5;

  // Calculation of the final thermal conductance (1 / total resistance)
  // conductance = 1 / resistance contact area pathway + 1 / resistance magrogap
  // pathway)
  thermal_conductance =
    1.0 / (resistance_macrocontact + 1.0 / (1.0 / resistance_microcontact +
                                            1.0 / resistance_gas_microgap)) +
    1.0 / (resistance_solid_macrogap + resistance_gas_macrogap);
}

void
apply_heat_transfer_on_local_particles(const double temperature_one,
                                       const double temperature_two,
                                       const double thermal_conductance,
                                       double &particle_one_heat_transfer_rate,
                                       double &particle_two_heat_transfer_rate)
{
  const double heat_transfer_rate =
    thermal_conductance * (temperature_two - temperature_one);
  particle_one_heat_transfer_rate += heat_transfer_rate;
  particle_two_heat_transfer_rate -= heat_transfer_rate;
}

void
apply_heat_transfer_on_single_local_particle(
  const double temperature_one,
  const double temperature_two,
  const double thermal_conductance,
  double      &particle_one_heat_transfer_rate)
{
  particle_one_heat_transfer_rate +=
    thermal_conductance * (temperature_two - temperature_one);
}

// only particle-particle and particle-floating-mesh contacts
template void
calculate_contact_thermal_conductance<ContactType::local_particle_particle>(
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

template void
calculate_contact_thermal_conductance<ContactType::ghost_particle_particle>(
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

template void
calculate_contact_thermal_conductance<
  ContactType::local_periodic_particle_particle>(
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


template void
calculate_contact_thermal_conductance<
  ContactType::ghost_periodic_particle_particle>(
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

template void
calculate_contact_thermal_conductance<
  ContactType::ghost_local_periodic_particle_particle>(
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

template void
calculate_contact_thermal_conductance<ContactType::particle_floating_mesh>(
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
