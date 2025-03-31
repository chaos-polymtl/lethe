// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_particle_heat_transfer.h>

inline double
calculate_corrected_contact_radius(const double effective_radius,
                                   const double effective_youngs_modulus,
                                   const double effective_real_youngs_modulus,
                                   const double normal_force_norm)
{
  double contact_radius = pow((3 * normal_force_norm * effective_radius) /
                                (4 * effective_youngs_modulus),
                              (1.0 / 3.0));

  // apply correction
  return contact_radius *
         pow(effective_youngs_modulus / effective_real_youngs_modulus,
             1.0 / 5.0);
}

inline double
calculate_macrocontact_resistance(const double harmonic_particle_conductivity,
                                  const double contact_radius)
{
  return 1 / (2.0 * contact_radius * harmonic_particle_conductivity + DBL_MIN);
}

inline double
calculate_microcontact_resistance(const double equivalent_surface_slope,
                                  const double equivalent_surface_roughness,
                                  const double effective_microhardness,
                                  const double contact_radius,
                                  const double harmonic_particle_conductivity,
                                  const double maximum_pressure)
{
  return 1.184 /
         (M_PI * harmonic_particle_conductivity * contact_radius *
          contact_radius) *
         (equivalent_surface_roughness / equivalent_surface_slope) *
         pow(effective_microhardness / (maximum_pressure + DBL_MIN), 0.96);
}

inline double
calculate_solid_macrogap_resistance(const double radius_one,
                                    const double radius_two,
                                    const double thermal_conductivity_one,
                                    const double thermal_conductivity_two,
                                    const double contact_radius)
{
  // resistance_i = characteristic length parallel to heat flux /
  // (thermal_conductivity
  // * characteristic area perpendicular to the heat flux)
  const double resistance_one =
    (M_PI * radius_one / 4.0) * 1 /
    (M_PI * (radius_one * radius_one - contact_radius * contact_radius) *
     thermal_conductivity_one);
  const double resistance_two =
    (M_PI * radius_two / 4.0) * 1 /
    (M_PI * (radius_two * radius_two - contact_radius * contact_radius) *
     thermal_conductivity_two);

  return resistance_one + resistance_two;
}

inline double
erfc_inverse_approximation(const double x)
{
  if (x < 1e-9 || x > 1.9)
    {
      std::cerr
        << "Error. Input out of range (10^-9 <= x <= 1.9) for erfc-1 function in thermal resistance calculation."
        << std::endl;
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

  return NAN;
}

inline double
calculate_interstitial_gas_microgap_resistance(
  const double equivalent_surface_roughness,
  const double contact_radius,
  const double gas_parameter_m,
  const double thermal_conductivity_gas,
  const double maximum_pressure,
  const double effective_microhardness)
{
  const double a_1 =
    erfc_inverse_approximation(2 * maximum_pressure / effective_microhardness);
  const double a_2 = erfc_inverse_approximation(0.03 * maximum_pressure /
                                                effective_microhardness) -
                     a_1;
  // If the values are out of the range of the erfc^-1 function
  if (a_1 == NAN || a_2 == NAN)
    {
      return INFINITY;
    }
  else
    {
      return (2 * sqrt(2) * equivalent_surface_roughness * a_2) /
             (M_PI * thermal_conductivity_gas * contact_radius *
              contact_radius *
              log(1 + a_2 / (a_1 +
                             gas_parameter_m /
                               (2 * sqrt(2) * equivalent_surface_roughness))));
    }
}

inline double
calculate_interstitial_gas_macrogap_resistance(
  const double harmonic_radius,
  const double thermal_conductivity_gas,
  const double contact_radius,
  const double gas_parameter_m)
{
  const double A = 2 * sqrt(harmonic_radius * harmonic_radius -
                            contact_radius * contact_radius);
  const double S = 2 * (harmonic_radius - (contact_radius * contact_radius) /
                                            (2 * harmonic_radius)) +
                   gas_parameter_m;
  return 2.0 / (M_PI * thermal_conductivity_gas * (S * log(S / (S - A)) - A));
}

inline void
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
  double      &thermal_conductance)
{
  const double harmonic_particle_conductivity =
    harmonic_mean(thermal_conductivity_one, thermal_conductivity_two);
  const double harmonic_radius = harmonic_mean(radius_one, radius_two);

  // Calculation of contact radius
  const double contact_radius =
    calculate_corrected_contact_radius(harmonic_radius * 0.5,
                                       effective_youngs_modulus,
                                       effective_real_youngs_modulus,
                                       normal_force_norm);
  const double corrected_normal_overlap =
    normal_overlap *
    pow(effective_youngs_modulus / effective_real_youngs_modulus, 2.0 / 3.0);
  const double maximum_pressure =
    (2.0 * effective_real_youngs_modulus * corrected_normal_overlap) /
    (M_PI * contact_radius + DBL_MIN);

  // Calculation of each thermal resistance
  const double resistance_macrocontact =
    calculate_macrocontact_resistance(harmonic_particle_conductivity,
                                      contact_radius);

  const double resistance_microcontact =
    calculate_microcontact_resistance(equivalent_surface_slope,
                                      equivalent_surface_roughness,
                                      effective_microhardness,
                                      contact_radius,
                                      harmonic_particle_conductivity,
                                      maximum_pressure);

  const double resistance_solid_macrogap =
    calculate_solid_macrogap_resistance(radius_one,
                                        radius_two,
                                        thermal_conductivity_one,
                                        thermal_conductivity_two,
                                        contact_radius);

  const double resistance_gas_microgap =
    calculate_interstitial_gas_microgap_resistance(equivalent_surface_roughness,
                                                   contact_radius,
                                                   gas_parameter_m,
                                                   thermal_conductivity_gas,
                                                   maximum_pressure,
                                                   effective_microhardness);

  const double resistance_gas_macrogap =
    calculate_interstitial_gas_macrogap_resistance(harmonic_radius,
                                                   thermal_conductivity_gas,
                                                   contact_radius,
                                                   gas_parameter_m);

  // Calculation of the final thermal conductance (1 / total resistance)
  // conductance = 1 / resistance contact area pathway + 1 / resistance magrogap
  // pathway)
  thermal_conductance =
    1.0 / (resistance_macrocontact + 1.0 / (1.0 / resistance_microcontact +
                                            1.0 / resistance_gas_microgap)) +
    1.0 / (resistance_solid_macrogap + resistance_gas_macrogap);
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

inline void
apply_heat_transfer_on_single_local_particle(const double temperature_one,
                                             const double temperature_two,
                                             const double thermal_conductance,
                                             double &particle_one_heat_transfer)
{
  particle_one_heat_transfer +=
    thermal_conductance * (temperature_two - temperature_one);
}
