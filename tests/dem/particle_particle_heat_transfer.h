// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef particle_particle_heat_transfer_h
#define particle_particle_heat_transfer_h

#define _USE_MATH_DEFINES

#include <dem/particle_particle_contact_force.h>

#include <cmath>


double
calculate_corrected_contact_radius(const double effective_radius,
                                   const double effective_youngs_modulus,
                                   const double effective_real_youngs_modulus,
                                   const double normal_force_norm)
{
  // Using the normal force
  double contact_radius = pow((3 * normal_force_norm * effective_radius) /
                                (4 * effective_youngs_modulus),
                              1 / 3);
  contact_radius        = contact_radius *
                   pow(effective_youngs_modulus / effective_real_youngs_modulus,
                       1 / 5); // apply correction
  return contact_radius;

  // Using the distance between particles
  // double distance = radius_one + radius_two -
  // normal_overlap; double contact_area =
  // - 1 / 4 * (distance - radius_one - radius_two) *
  // (distance + radius_one - radius_two) * (distance -
  // radius_one + radius_two) * (distance +
  // radius_one + radius_two) / (distance * distance); double
  // contact_radius = sqrt(contact_area); contact_radius = contact_radius * pow(
  // effective_youngs_modulus / effective_real_youngs_modulus , 1/5); // apply
  // correction return contact_radius;
}

double
calculate_macrocontact_spreading_resistance(
  const double effective_particle_conductivity,
  const double contact_radius)
{
  return 1 / (4 * contact_radius * effective_particle_conductivity);

  // see if difference with (On thermal characterisation...) : (1-a/r_h)^(3/2) /
  // 4ka
  // consider if possible to add collision -> how to get contact duration ?
}

double
calculate_microcontact_spreading_resistance(
  const double effective_surface_slope,
  const double effective_surface_roughness,
  const double effective_microhardness,
  const double contact_radius,
  const double effective_particle_conductivity,
  const double maximum_pressure)
{
  return 1.184 / M_PI *
         (effective_surface_roughness / effective_surface_slope) *
         pow(effective_microhardness / maximum_pressure, 0.96) /
         effective_particle_conductivity / (contact_radius * contact_radius);
}

double
calculate_macrogap_solid_layers_resistance(radius_one,
                                           radius_two,
                                           thermal_conductivity_one,
                                           thermal_conductivity_two,
                                           contact_radius)
{
  // R_c_i = characteristic length parallel to heat flux / (thermal_conductivity
  // * characteristic area perpendicular to heat flux)
  double R_c_one =
    (M_PI * radius_one / 4) * 1 /
    (M_PI * (radius_one * radius_one - contact_radius * contact_radius) *
     thermal_conductivity_one);
  double R_c_two =
    (M_PI * radius_two / 4) * 1 /
    (M_PI * (radius_two * radius_two - contact_radius * contact_radius) *
     thermal_conductivity_two);

  return R_c_one + R_c_two;
}

double
erfc_inverse_approximation(x)
{
  if (x < 1e-9 || x > 1.9)
    {
      std::cerr << "Input out of range (10^-9 <= x <= 1.9)" << std::endl;
      return NAN;
    }

  if (x <= 0.02)
    {
      return 1.0 / (0.218 + 0.735 * pow(x, 0.173));
    }
  else if (x <= 0.5)
    {
      return (1.05 * pow(0.175, x)) / pow(x, 0.12);
    }
  else if (x <= 1.9)
    {
      return (1 - x) / (0.707 + 0.862 * x - 0.431 * x * x);
    }
}


double
calculate_microgap_interstitial_gas_resistance(
  const double effective_surface_roughness,
  const double contact_radius,
  const double gas_parameter_m,
  const double thermal_conductivity_gas,
  const double maximum_pressure,
  const double effective_microhardness)
{
  const double a1 =
    erfc_inverse_approximation(2 * maximum_pressure / effective_microhardness);
  const double a2 = erfc_inverse_approximation(0.03 * maximum_pressure /
                                               effective_microhardness) -
                    a1;

  return (2 * sqrt(2) * effective_surface_roughness * a2) /
         (M_PI * thermal_conductivity_gas * contact_radius * contact_radius *
          log(1 + a2 / (a1 + gas_parameter_m /
                               (2 * sqrt(2) * effective_surface_roughness))))
}

double
calculate_macrogap_interstitial_gas_resistance(
  const double effective_radius,
  const double thermal_conductivity_gas,
  const double contact_radius,
  const double gas_parameter_m)
{
  const double A = 2 * sqrt(effective_radius * effective_radius -
                            contact_radius * contact_radius);
  const double S = 2 * (effective_radius - (contact_radius * contact_radius) /
                                             (2 * effective_radius)) +
                   gas_parameter_m;
  return 2 / (M_PI * thermal_conductivity_gas * (S * log(S / (S - A)) - A));

  // // other simplified model
  // const double l_gas = (effective_radius * effective_radius * (1 - M_PI / 4))
  // /
  //                      (effective_radius - contact_radius);
  // const double A_gas = M_PI * (2 * effective_radius * effective_radius -
  //                              contact_radius * contact_radius);
  // return 1 / thermal_conductivity_gas * (l_gas + gas_parameter_m) / A_gas;
}

// better to calculate heat transfer directly so as not to pass properties as
// parameters of apply function ?
inline void
calculate_contact_thermal_conductance(
  const double radius_one,
  const double radius_two,
  const double effective_radius,
  const double effective_youngs_modulus,
  const double effective_real_youngs_modulus,
  const double effective_surface_roughness,
  const double effective_surface_slope,
  const double effective_microhardness,
  const double thermal_conductivity_one,
  const double thermal_conductivity_two,
  const double thermal_conductivity_gas,
  const double gas_parameter_m,
  const double normal_overlap,
  const double normal_force_norm,
  double      &thermal_conductance)
{
  const double effective_particle_conductivity =
    (thermal_conductivity_one * thermal_conductivity_two) /
    (thermal_conductivity_one + thermal_conductivity_two);

  // Calculation of contact radius
  const double contact_radius =
    calculate_corrected_contact_radius(effective_radius,
                                       effective_youngs_modulus,
                                       effective_real_youngs_modulus,
                                       normal_force_norm);
  const double maximum_pressure =
    (2 * effective_real_youngs_modulus * normal_overlap) /
    (M_PI * contact_radius);

  // Calculation of each thermal resistance
  const double R_L =
    calculate_macrocontact_spreading_resistance(effective_particle_conductivity,
                                                contact_radius);

  const double R_s =
    calculate_microcontact_spreading_resistance(effective_surface_slope,
                                                effective_surface_roughness,
                                                effective_microhardness,
                                                contact_radius,
                                                effective_particle_conductivity,
                                                maximum_pressure);

  const double R_c =
    calculate_macrogap_solid_layers_resistance(radius_one,
                                               radius_two,
                                               thermal_conductivity_one,
                                               thermal_conductivity_two,
                                               contact_radius);

  const double R_g =
    calculate_microgap_interstitial_gas_resistance(effective_surface_roughness,
                                                   contact_radius,
                                                   gas_parameter_m,
                                                   thermal_conductivity_gas,
                                                   maximum_pressure,
                                                   effective_microhardness);

  const double R_G =
    calculate_macrogap_interstitial_gas_resistance(effective_radius,
                                                   thermal_conductivity_gas,
                                                   contact_radius,
                                                   gas_parameter_m);

  // Calculation of the final thermal conductance (1 / total resistance)
  // conductance = 1 / resistance contact area pathway + 1 / resistance magrogap
  // pathway)
  thermal_conductance = 1 / (R_L + 1 / (1 / R_s + 1 / R_g)) + 1 / (R_c + R_G);
}

inline void
apply_heat_transfer_on_local_particles(const double temperature_one,
                                       const double temperature_two,
                                       const double thermal_conductance,
                                       double      &particle_one_heat_transfer,
                                       double      &particle_two_heat_transfer)
{
  const double heat_transfer =
    thermal_conductance * (temperature_two - temperature_one);
  particle_one_heat_transfer += heat_transfer;
  particle_two_heat_transfer -= heat_transfer;
}

// one should be local particle
inline void
apply_heat_transfer_on_single_local_particle(const double temperature_one,
                                             const double temperature_two,
                                             const double thermal_conductance,
                                             double &particle_one_heat_transfer)
{
  particle_one_heat_transfer +=
    thermal_conductance * (temperature_two - temperature_one);
}


#endif
