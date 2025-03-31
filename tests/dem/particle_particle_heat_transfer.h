// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef particle_particle_heat_transfer_h
#define particle_particle_heat_transfer_h

#include <core/auxiliary_math_functions.h>

#include <cmath>
#include <iostream>


/**
 * @brief Calculate the contact radius between two particles.
 *
 * @param[in] effective_radius Effective radius of the two particles.
 * @param[in] effective_youngs_modulus Effective young's modulus of the two
 * particles.
 * @param[in] effective_real_youngs_modulus Effective real young's modulus of
 * the two particles.
 * @param[in] normal_force_norm Norm of the normal contact force between the two
 * particles.
 * @return Contact radius of the two particles.
 */
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

/**
 * @brief Calculate the macrocontact thermal resistance between two particles.
 *
 * @param[in] harmonic_particle_conductivity Harmonic mean of the two particles'
 * thermal conductivities.
 * @param[in] contact_radius Contact radius of the two particles.
 * @return Macrocontact thermal resistance between the two particles.
 */
inline double
calculate_macrocontact_resistance(const double harmonic_particle_conductivity,
                                  const double contact_radius)
{
  return 1 / (2.0 * contact_radius * harmonic_particle_conductivity);
}

/**
 * @brief Calculate the microcontact thermal resistance between two particles.
 *
 * @param[in] equivalent_surface_slope Equivalent value of the two particles'
 * surface slope.
 * @param[in] equivalent_surface_roughness Equivalent value of the two
 * particles' surface roughness.
 * @param[in] effective_microhardness Effective microhardness of the two
 * particles.
 * @param[in] contact_radius Contact radius of the two particles.
 * @param[in] harmonic_particle_conductivity Harmonic mean of the two particles'
 * thermal conductivities.
 * @param[in] maximum_pressure Maximum pressure for hertzian contacts.
 * @return Microcontact thermal resistance between the two particles.
 */
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
         pow(effective_microhardness / maximum_pressure, 0.96);
}

/**
 * @brief Calculate the solid macrogap thermal resistance between two particles.
 *
 * @param[in] radius_one Radius of particle one.
 * @param[in] radius_two Radius of particle two.
 * @param[in] thermal_conductivity_one Thermal conductivity of particle one.
 * @param[in] thermal_conductivity_two Thermal conductivity of particle two.
 * @param[in] contact_radius Contact radius of the two particles.
 * @return Solid macrogap thermal resistance between the two particles.
 */
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

/**
 * @brief Define an approximation of the erfc^-1 function for values between 1e-9 and 1.9.
 *
 * @param[in] x Value to which apply the function.
 * @return erfc^-1(x)
 */
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

/**
 * @brief Calculate the interstitial gas microgap thermal resistance between two particles.
 *
 * @param[in] equivalent_surface_roughness Equivalent value of the two
 * particles' surface roughness.
 * @param[in] contact_radius Contact radius of the two particles.
 * @param[in] gas_parameter_m Gas parameter.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] maximum_pressure Maximum pressure for hertzian contacts.
 * @param[in] effective_microhardness Effective microhardness of the two
 * particles.
 *
 * @return Interstitial gas microgap thermal resistance between the two particles.
 */
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

/**
 * @brief Calculate the interstitial gas macrogap thermal resistance between two particles.
 *
 * @param[in] harmonic_radius Harmonic mean of the particles' radii.
 * @param[in] thermal_conductivity_gas Thermal conductivity of the interstitial
 * gas.
 * @param[in] contact_radius Contact radius of the two particles.
 * @param[in] gas_parameter_m Gas parameter.
 *
 * @return Interstitial gas microgap thermal resistance between the two particles.
 */
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

  // const double length_gas =
  //   (harmonic_radius * harmonic_radius * (1 - M_PI / 4)) /
  //   (harmonic_radius - contact_radius);
  // const double area_gas = M_PI * (2 * harmonic_radius * harmonic_radius -
  //                                 contact_radius * contact_radius);
  // return 1 / thermal_conductivity_gas * (length_gas + gas_parameter_m) /
  //        area_gas;
}

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
    (M_PI * contact_radius);

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


/**
 * @brief Apply the heat transfer to the local-local particle pair in contact.
 *
 * @param[in] temperature_one Temperature of particle one.
 * @param[in] temperature_two Temperature of particle two.
 * @param[in] thermal_conductance Total thermal conductance between the two
 * particles.
 * @param[in,out] particle_one_heat_transfer Heat transfer due to contact
 * applied to particle one.
 * @param[in,out] particle_two_heat_transfer Heat transfer due to contact
 * applied to particle two.
 *
 */
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

/**
 * @brief Apply the heat transfer to the local-ghost particle pair in contact.
 * The heat transfer is only applied to the local particle, so particle one
 * should be the local particle here.
 *
 * @param[in] temperature_one Temperature of particle one.
 * @param[in] temperature_two Temperature of particle two.
 * @param[in] thermal_conductance Total thermal conductance between the two
 * particles.
 * @param[in,out] particle_one_heat_transfer Heat transfer due to contact
 * applied to particle one.
 *
 */
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
