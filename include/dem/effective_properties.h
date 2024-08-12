/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_effective_properties_h
#define lethe_effective_properties_h

#include <core/auxiliary_math_functions.h>

#include <dem/dem_solver_parameters.h>

#include <vector>

template <int dim>
class DEMEffectiveProperties
{
public:
  /**
   * @brief Copy constructor as a delete function to make sure it can not be
   * copied. It will never be used.
   *
   * @param copy The object to be copied
   */
  DEMEffectiveProperties(const DEMEffectiveProperties &copy) = delete;

  /**
   * @brief Copy constructor as a delete function to make sure it can not be
   * assigned. It will never be used.
   *
   * @param copy The object to be copied
   */
  DEMEffectiveProperties &
  operator=(const DEMEffectiveProperties &copy) = delete;

  /**
   * @brief Getter of the unique instance of the DEMActionManager.
   */
  static DEMEffectiveProperties *
  get_effective_properties();

  inline unsigned int
  vec_particle_type_index(const unsigned int i, const unsigned int j)
  {
    return i * n_particle_types + j;
  }

  /**
   * @brief Set every containers needed to carry the particle-particle force
   * calculation.
   *
   * @param dem_parameters DEM parameters declared in the .prm file.
   */
  void
  set_effective_properties(const DEMSolverParameters<dim> &dem_parameters)
  {
    auto properties = dem_parameters.lagrangian_physical_properties;

    n_particle_types = properties.particle_type_number;
    effective_youngs_modulus.resize(n_particle_types * n_particle_types);
    effective_shear_modulus.resize(n_particle_types * n_particle_types);
    effective_coefficient_of_restitution.resize(n_particle_types *
                                                n_particle_types);
    effective_coefficient_of_friction.resize(n_particle_types *
                                             n_particle_types);
    effective_coefficient_of_rolling_friction.resize(n_particle_types *
                                                     n_particle_types);
    model_parameter_beta.resize(n_particle_types * n_particle_types);
    effective_surface_energy.resize(n_particle_types * n_particle_types);
    effective_hamaker_constant.resize(n_particle_types * n_particle_types);

    for (unsigned int i = 0; i < n_particle_types; ++i)
      {
        const double youngs_modulus_i =
          properties.youngs_modulus_particle.at(i);
        const double poisson_ratio_i = properties.poisson_ratio_particle.at(i);
        const double restitution_coefficient_i =
          properties.restitution_coefficient_particle.at(i);
        const double friction_coefficient_i =
          properties.friction_coefficient_particle.at(i);
        const double rolling_friction_coefficient_i =
          properties.rolling_friction_coefficient_particle.at(i);
        const double surface_energy_i =
          properties.surface_energy_particle.at(i);
        const double hamaker_constant_i =
          properties.hamaker_constant_particle.at(i);

        for (unsigned int j = 0; j < n_particle_types; ++j)
          {
            const unsigned int k = i * n_particle_types + j;

            const double youngs_modulus_j =
              properties.youngs_modulus_particle.at(j);
            const double poisson_ratio_j =
              properties.poisson_ratio_particle.at(j);
            const double restitution_coefficient_j =
              properties.restitution_coefficient_particle.at(j);
            const double friction_coefficient_j =
              properties.friction_coefficient_particle.at(j);
            const double rolling_friction_coefficient_j =
              properties.rolling_friction_coefficient_particle.at(j);
            const double surface_energy_j =
              properties.surface_energy_particle.at(j);
            const double hamaker_constant_j =
              properties.hamaker_constant_particle.at(j);

            this->effective_youngs_modulus[k] =
              (youngs_modulus_i * youngs_modulus_j) /
              ((youngs_modulus_j * (1.0 - poisson_ratio_i * poisson_ratio_i)) +
               (youngs_modulus_i * (1.0 - poisson_ratio_j * poisson_ratio_j)) +
               DBL_MIN);

            this->effective_shear_modulus[k] =
              (youngs_modulus_i * youngs_modulus_j) /
              (2.0 * ((youngs_modulus_j * (2.0 - poisson_ratio_i) *
                       (1.0 + poisson_ratio_i)) +
                      (youngs_modulus_i * (2.0 - poisson_ratio_j) *
                       (1.0 + poisson_ratio_j))) +
               DBL_MIN);

            this->effective_coefficient_of_restitution[k] =
              harmonic_mean(restitution_coefficient_i,
                            restitution_coefficient_j);

            this->effective_coefficient_of_friction[k] =
              harmonic_mean(friction_coefficient_i, friction_coefficient_j);

            this->effective_coefficient_of_rolling_friction[k] =
              harmonic_mean(rolling_friction_coefficient_i,
                            rolling_friction_coefficient_j);

            this->effective_surface_energy[k] =
              surface_energy_i + surface_energy_j -
              std::pow(std::sqrt(surface_energy_i) -
                         std::sqrt(surface_energy_j),
                       2);

            this->effective_hamaker_constant[k] =
              0.5 * (hamaker_constant_i + hamaker_constant_j);

            double restitution_coefficient_particle_log =
              std::log(this->effective_coefficient_of_restitution[k]);

            this->model_parameter_beta[k] =
              restitution_coefficient_particle_log /
              sqrt(restitution_coefficient_particle_log *
                     restitution_coefficient_particle_log +
                   9.8696);
          }
      }
  }

  /**
   * @brief Pointer to the unique of the DEMActionManager.
   */
  static DEMEffectiveProperties *instance;

  unsigned int        n_particle_types;
  std::vector<double> effective_youngs_modulus;
  std::vector<double> effective_shear_modulus;
  std::vector<double> effective_coefficient_of_restitution;
  std::vector<double> effective_coefficient_of_friction;
  std::vector<double> effective_coefficient_of_rolling_friction;
  std::vector<double> effective_surface_energy;
  std::vector<double> effective_hamaker_constant;
  std::vector<double> model_parameter_beta;
};



#endif
