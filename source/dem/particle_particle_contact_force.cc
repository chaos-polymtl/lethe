/* ---------------------------------------------------------------------
  *
  * Copyright (C) 2019 - 2020 by the Lethe authors
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
  * Author: Shahab Golshan, Polytechnique Montreal, 2019
  */

#include <core/parameters_lagrangian.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/particle_particle_contact_force.h>

using namespace DEM;

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters)
{
  auto properties  = dem_parameters.lagrangian_physical_properties;
  n_particle_types = properties.particle_type_number;
  effective_youngs_modulus.resize(n_particle_types * n_particle_types);
  effective_shear_modulus.resize(n_particle_types * n_particle_types);
  effective_coefficient_of_restitution.resize(n_particle_types *
                                              n_particle_types);
  effective_coefficient_of_friction.resize(n_particle_types * n_particle_types);
  effective_coefficient_of_rolling_friction.resize(n_particle_types *
                                                   n_particle_types);

  model_parameter_beta.resize(n_particle_types * n_particle_types);

  for (unsigned int i = 0; i < n_particle_types; ++i)
    {
      const double youngs_modulus_i = properties.youngs_modulus_particle.at(i);
      const double poisson_ratio_i  = properties.poisson_ratio_particle.at(i);
      const double restitution_coefficient_i =
        properties.restitution_coefficient_particle.at(i);
      const double friction_coefficient_i =
        properties.friction_coefficient_particle.at(i);
      const double rolling_friction_coefficient_i =
        properties.rolling_friction_coefficient_particle.at(i);

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
            harmonic_mean(restitution_coefficient_i, restitution_coefficient_j);

          this->effective_coefficient_of_friction[k] =
            harmonic_mean(friction_coefficient_i, friction_coefficient_j);

          this->effective_coefficient_of_rolling_friction[k] =
            harmonic_mean(rolling_friction_coefficient_i,
                          rolling_friction_coefficient_j);

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

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  calculate_particle_particle_contact_force(
    DEMContactManager<dim> &   container_manager,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const Tensor<1, dim>       periodic_offset)
{
  auto &local_adjacent_particles = container_manager.local_adjacent_particles;
  auto &ghost_adjacent_particles = container_manager.ghost_adjacent_particles;
  auto &local_periodic_adjacent_particles =
    container_manager.local_periodic_adjacent_particles;
  auto &ghost_periodic_adjacent_particles =
    container_manager.ghost_periodic_adjacent_particles;
  auto &ghost_local_periodic_adjacent_particles =
    container_manager.ghost_local_periodic_adjacent_particles;

  // Contact forces calculations of local-local and local-ghost particle
  // pairs are performed in separate loops

  // Looping over local_adjacent_particles values with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          // Gather information about particle 1 and set it up.
          auto first_contact_info = adjacent_particles_list.begin();
          auto particle_one       = first_contact_info->second.particle_one;
          auto particle_one_properties = particle_one->get_properties();

          types::particle_index particle_one_id =
            particle_one->get_local_index();
          Tensor<1, 3> &particle_one_torque = torque[particle_one_id];
          Tensor<1, 3> &particle_one_force  = force[particle_one_id];

          // Fix particle one location for 2d and 3d
          Point<3> particle_one_location = [&] {
            if constexpr (dim == 3)
              {
                return particle_one->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle_one->get_location()));
              }
          }();

          for (auto &&contact_info :
               adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle 2 in
              // contact with particle 1
              auto particle_two            = contact_info.particle_two;
              auto particle_two_properties = particle_two->get_properties();

              // Get particle 2 location in dimension independent way
              Point<3> particle_two_location = [&] {
                if constexpr (dim == 3)
                  {
                    return particle_two->get_location();
                  }
                if constexpr (dim == 2)
                  {
                    return (point_nd_to_3d(particle_two->get_location()));
                  }
              }();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    this->normal_relative_velocity_value,
                    this->normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }


                  types::particle_index particle_two_id =
                    particle_two->get_local_index();

                  Tensor<1, 3> &particle_two_torque = torque[particle_two_id];
                  Tensor<1, 3> &particle_two_force  = force[particle_two_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_local_particles(
                    this->normal_force,
                    this->tangential_force,
                    this->particle_one_tangential_torque,
                    this->particle_two_tangential_torque,
                    this->rolling_resistance_torque,
                    particle_one_torque,
                    particle_two_torque,
                    particle_one_force,
                    particle_two_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  contact_info.tangential_overlap.clear();
                }
            }
        }
    }

  // Doing the same calculations for local-ghost particle pairs

  // Looping over ghost_adjacent_particles with iterator
  // adjacent_particles_iterator
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          // Gather information about particle 1 and set it up.
          auto first_contact_info = adjacent_particles_list.begin();
          auto particle_one       = first_contact_info->second.particle_one;
          auto particle_one_properties = particle_one->get_properties();

          types::particle_index particle_one_id =
            particle_one->get_local_index();
          Tensor<1, 3> &particle_one_torque = torque[particle_one_id];
          Tensor<1, 3> &particle_one_force  = force[particle_one_id];

          // Fix particle one location for 2d and 3d
          Point<3> particle_one_location = [&] {
            if constexpr (dim == 3)
              {
                return particle_one->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle_one->get_location()));
              }
          }();

          for (auto &&contact_info :
               adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle one
              // and two in contact
              auto particle_two            = contact_info.particle_two;
              auto particle_two_properties = particle_two->get_properties();
              // Get particle 2 location in dimension independent way
              Point<3> particle_two_location = [&] {
                if constexpr (dim == 3)
                  {
                    return particle_two->get_location();
                  }
                if constexpr (dim == 2)
                  {
                    return (point_nd_to_3d(particle_two->get_location()));
                  }
              }();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    this->normal_relative_velocity_value,
                    this->normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }


                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_single_local_particle(
                    this->normal_force,
                    this->tangential_force,
                    this->particle_one_tangential_torque,
                    this->rolling_resistance_torque,
                    particle_one_torque,
                    particle_one_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  contact_info.tangential_overlap.clear();
                }
            }
        }
    }


  for (auto &&periodic_adjacent_particles_list :
       local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      if (!periodic_adjacent_particles_list.empty())
        {
          // Gather information about particle 1 and set it up.
          auto first_contact_info = periodic_adjacent_particles_list.begin();
          auto particle_one       = first_contact_info->second.particle_one;
          auto particle_one_properties = particle_one->get_properties();

          types::particle_index particle_one_id =
            particle_one->get_local_index();
          Tensor<1, 3> &particle_one_torque = torque[particle_one_id];
          Tensor<1, 3> &particle_one_force  = force[particle_one_id];

          // Fix particle one location for 2d and 3d
          Point<3> particle_one_location = [&] {
            if constexpr (dim == 3)
              {
                return particle_one->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle_one->get_location()));
              }
          }();

          for (auto &&contact_info :
               periodic_adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle 2 in
              // contact with particle 1
              auto particle_two            = contact_info.particle_two;
              auto particle_two_properties = particle_two->get_properties();

              // Get particle 2 location in dimension independent way
              Point<3> particle_two_location = [&] {
                if constexpr (dim == 3)
                  {
                    return particle_two->get_location() - periodic_offset;
                  }
                if constexpr (dim == 2)
                  {
                    return (point_nd_to_3d(particle_two->get_location() -
                                           periodic_offset));
                  }
              }();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    this->normal_relative_velocity_value,
                    this->normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }


                  types::particle_index particle_two_id =
                    particle_two->get_local_index();

                  Tensor<1, 3> &particle_two_torque = torque[particle_two_id];
                  Tensor<1, 3> &particle_two_force  = force[particle_two_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_local_particles(
                    this->normal_force,
                    this->tangential_force,
                    this->particle_one_tangential_torque,
                    this->particle_two_tangential_torque,
                    this->rolling_resistance_torque,
                    particle_one_torque,
                    particle_two_torque,
                    particle_one_force,
                    particle_two_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  contact_info.tangential_overlap.clear();
                }
            }
        }
    }

  // Doing the same calculations for local-ghost particle pairs

  // Looping over ghost_adjacent_particles with iterator
  // adjacent_particles_iterator
  for (auto &&periodic_adjacent_particles_list :
       ghost_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      if (!periodic_adjacent_particles_list.empty())
        {
          // Gather information about particle 1 and set it up.
          auto first_contact_info = periodic_adjacent_particles_list.begin();
          auto particle_one       = first_contact_info->second.particle_one;
          auto particle_one_properties = particle_one->get_properties();

          types::particle_index particle_one_id =
            particle_one->get_local_index();
          Tensor<1, 3> &particle_one_torque = torque[particle_one_id];
          Tensor<1, 3> &particle_one_force  = force[particle_one_id];

          // Fix particle one location for 2d and 3d
          Point<3> particle_one_location = [&] {
            if constexpr (dim == 3)
              {
                return particle_one->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle_one->get_location()));
              }
          }();

          for (auto &&contact_info :
               periodic_adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle one
              // and two in contact
              auto particle_two            = contact_info.particle_two;
              auto particle_two_properties = particle_two->get_properties();
              // Get particle 2 location in dimension independent way
              Point<3> particle_two_location = [&] {
                if constexpr (dim == 3)
                  {
                    return particle_two->get_location() - periodic_offset;
                  }
                if constexpr (dim == 2)
                  {
                    return (point_nd_to_3d(particle_two->get_location() -
                                           periodic_offset));
                  }
              }();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    this->normal_relative_velocity_value,
                    this->normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_one_properties,
                        particle_two_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }


                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_single_local_particle(
                    this->normal_force,
                    this->tangential_force,
                    this->particle_one_tangential_torque,
                    this->rolling_resistance_torque,
                    particle_one_torque,
                    particle_one_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  contact_info.tangential_overlap.clear();
                }
            }
        }
    }

  for (auto &&periodic_adjacent_particles_list :
       ghost_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      if (!periodic_adjacent_particles_list.empty())
        {
          // Gather information about particle 1 and set it up.
          auto first_contact_info = periodic_adjacent_particles_list.begin();
          auto particle_one       = first_contact_info->second.particle_one;
          auto particle_one_properties = particle_one->get_properties();

          // Fix particle one location for 2d and 3d
          Point<3> particle_one_location = [&] {
            if constexpr (dim == 3)
              {
                return particle_one->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle_one->get_location()));
              }
          }();

          for (auto &&contact_info :
               periodic_adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle 2 in
              // contact with particle 1
              auto particle_two            = contact_info.particle_two;
              auto particle_two_properties = particle_two->get_properties();

              types::particle_index particle_two_id =
                particle_two->get_local_index();
              Tensor<1, 3> &particle_two_torque = torque[particle_two_id];
              Tensor<1, 3> &particle_two_force  = force[particle_two_id];

              // Get particle 2 location in dimension independent way
              Point<3> particle_two_location = [&] {
                if constexpr (dim == 3)
                  {
                    return particle_two->get_location() - periodic_offset;
                  }
                if constexpr (dim == 2)
                  {
                    return (point_nd_to_3d(particle_two->get_location() -
                                           periodic_offset));
                  }
              }();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    this->normal_relative_velocity_value,
                    this->normal_unit_vector,
                    particle_two_properties,
                    particle_one_properties,
                    particle_two_location,
                    particle_one_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_one_tangential_torque,
                        this->particle_two_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_two_tangential_torque,
                        this->particle_one_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_two_tangential_torque,
                        this->particle_one_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact(
                        contact_info,
                        this->normal_relative_velocity_value,
                        this->normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        this->normal_force,
                        this->tangential_force,
                        this->particle_two_tangential_torque,
                        this->particle_one_tangential_torque,
                        this->rolling_resistance_torque);
                    }

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_single_local_particle(
                    this->normal_force,
                    this->tangential_force,
                    this->particle_two_tangential_torque,
                    this->rolling_resistance_torque,
                    particle_two_torque,
                    particle_two_force);
                }
              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  contact_info.tangential_overlap.clear();
                }
            }
        }
    }
}

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  calculate_IB_particle_particle_contact_force(
    const double                         normal_overlap,
    particle_particle_contact_info<dim> &contact_info,
    Tensor<1, 3> &                       normal_force,
    Tensor<1, 3> &                       tangential_force,
    Tensor<1, 3> &                       particle_one_tangential_torque,
    Tensor<1, 3> &                       particle_two_tangential_torque,
    Tensor<1, 3> &                       rolling_resistance_torque,
    IBParticle<dim> &                    particle_one,
    IBParticle<dim> &                    particle_two,
    const Point<dim> &                   particle_one_location,
    const Point<dim> &                   particle_two_location,
    const double                         dt,
    const double                         particle_one_radius,
    const double                         particle_two_radius,
    const double                         particle_one_mass,
    const double                         particle_two_mass)
{
  Point<3> particle_one_location_3d = [&] {
    if constexpr (dim == 3)
      {
        return particle_one_location;
      }
    if constexpr (dim == 2)
      {
        return (point_nd_to_3d(particle_one_location));
      }
  }();

  Point<3> particle_two_location_3d = [&] {
    if constexpr (dim == 3)
      {
        return particle_two_location;
      }
    if constexpr (dim == 2)
      {
        return (point_nd_to_3d(particle_two_location));
      }
  }();

  auto particle_one_properties = particle_one.get_properties();
  particle_one_properties[PropertiesIndex::mass] = particle_one_mass;
  particle_one_properties[PropertiesIndex::type] = 0;
  particle_one_properties[PropertiesIndex::dp]   = 2.0 * particle_one_radius;

  auto particle_two_properties = particle_one.get_properties();
  particle_two_properties[PropertiesIndex::mass] = particle_two_mass;
  particle_two_properties[PropertiesIndex::type] = 0;
  particle_two_properties[PropertiesIndex::dp]   = 2.0 * particle_two_radius;

  n_particle_types = 1;
  effective_youngs_modulus.resize(n_particle_types * n_particle_types);
  effective_shear_modulus.resize(n_particle_types * n_particle_types);
  effective_coefficient_of_restitution.resize(n_particle_types *
                                              n_particle_types);
  effective_coefficient_of_friction.resize(n_particle_types * n_particle_types);
  effective_coefficient_of_rolling_friction.resize(n_particle_types *
                                                   n_particle_types);
  model_parameter_beta.resize(n_particle_types * n_particle_types);

  // DEM::PropertiesIndex::type is the first (0) property of particles in the
  // DEM solver. For the IB particles, the first property is ID. For force and
  // torque calculations, we need pairwise properties (such as effective
  // Young's modulus, effective coefficient of restitution, etc.) We rewrite
  // these pairwise properties by using the ID of IB particles (using
  // DEM::PropertiesIndex::type) and use them in force calculations.
  const unsigned int particle_one_type =
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];

  if constexpr (contact_model == Parameters::Lagrangian::
                                   ParticleParticleContactForceModel::linear)
    {
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        ((particle_two.youngs_modulus *
          (1.0 - particle_one.poisson_ratio * particle_one.poisson_ratio)) +
         (particle_one.youngs_modulus *
          (1.0 - particle_two.poisson_ratio * particle_two.poisson_ratio)) +
         DBL_MIN);

      this->effective_shear_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        (2.0 *
           ((particle_two.youngs_modulus * (2.0 - particle_one.poisson_ratio) *
             (1.0 + particle_one.poisson_ratio)) +
            (particle_one.youngs_modulus * (2.0 - particle_two.poisson_ratio) *
             (1.0 + particle_two.poisson_ratio))) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.restitution_coefficient,
                      particle_two.restitution_coefficient);
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.friction_coefficient,
                      particle_two.friction_coefficient);
      this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.rolling_friction_coefficient,
                      particle_two.rolling_friction_coefficient);

      // Since the normal overlap is already calculated we update
      // this element of the container here. The rest of information
      // are updated using the following function
      update_contact_information(contact_info,
                                 this->normal_relative_velocity_value,
                                 this->normal_unit_vector,
                                 particle_one_properties,
                                 particle_two_properties,
                                 particle_one_location_3d,
                                 particle_two_location_3d,
                                 dt);


      calculate_linear_contact(contact_info,
                               this->normal_relative_velocity_value,
                               this->normal_unit_vector,
                               normal_overlap,
                               particle_one_properties,
                               particle_two_properties,
                               this->normal_force,
                               this->tangential_force,
                               this->particle_one_tangential_torque,
                               this->particle_two_tangential_torque,
                               this->rolling_resistance_torque);
    }

  if constexpr (contact_model == Parameters::Lagrangian::
                                   ParticleParticleContactForceModel::hertz)
    {
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        ((particle_two.youngs_modulus *
          (1.0 - particle_one.poisson_ratio * particle_one.poisson_ratio)) +
         (particle_one.youngs_modulus *
          (1.0 - particle_two.poisson_ratio * particle_two.poisson_ratio)) +
         DBL_MIN);

      this->effective_shear_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        (2.0 *
           ((particle_two.youngs_modulus * (2.0 - particle_one.poisson_ratio) *
             (1.0 + particle_one.poisson_ratio)) +
            (particle_one.youngs_modulus * (2.0 - particle_two.poisson_ratio) *
             (1.0 + particle_two.poisson_ratio))) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.restitution_coefficient,
                      particle_two.restitution_coefficient);
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.friction_coefficient,
                      particle_two.friction_coefficient);
      this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.rolling_friction_coefficient,
                      particle_two.rolling_friction_coefficient);

      const double restitution_coefficient_particle_log = std::log(
        this->effective_coefficient_of_restitution[vec_particle_type_index(
          particle_one_type, particle_two_type)]);

      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] =
        restitution_coefficient_particle_log /
        sqrt(restitution_coefficient_particle_log *
               restitution_coefficient_particle_log +
             9.8696);

      // Since the normal overlap is already calculated we update
      // this element of the container here. The rest of information
      // are updated using the following function
      this->update_contact_information(contact_info,
                                       this->normal_relative_velocity_value,
                                       this->normal_unit_vector,
                                       particle_one_properties,
                                       particle_two_properties,
                                       particle_one_location_3d,
                                       particle_two_location_3d,
                                       dt);

      calculate_hertz_contact(contact_info,
                              this->normal_relative_velocity_value,
                              this->normal_unit_vector,
                              normal_overlap,
                              particle_one_properties,
                              particle_two_properties,
                              normal_force,
                              tangential_force,
                              particle_one_tangential_torque,
                              particle_two_tangential_torque,
                              rolling_resistance_torque);
    }

  if constexpr (contact_model ==
                Parameters::Lagrangian::ParticleParticleContactForceModel::
                  hertz_mindlin_limit_force)
    {
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        ((particle_two.youngs_modulus *
          (1.0 - particle_one.poisson_ratio * particle_one.poisson_ratio)) +
         (particle_one.youngs_modulus *
          (1.0 - particle_two.poisson_ratio * particle_two.poisson_ratio)) +
         DBL_MIN);

      this->effective_shear_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        (2.0 *
           ((particle_two.youngs_modulus * (2.0 - particle_one.poisson_ratio) *
             (1.0 + particle_one.poisson_ratio)) +
            (particle_one.youngs_modulus * (2.0 - particle_two.poisson_ratio) *
             (1.0 + particle_two.poisson_ratio))) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.restitution_coefficient,
                      particle_two.restitution_coefficient);
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.friction_coefficient,
                      particle_two.friction_coefficient);
      this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.rolling_friction_coefficient,
                      particle_two.rolling_friction_coefficient);

      const double restitution_coefficient_particle_log = std::log(
        this->effective_coefficient_of_restitution[vec_particle_type_index(
          particle_one_type, particle_two_type)]);

      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] =
        restitution_coefficient_particle_log /
        sqrt(restitution_coefficient_particle_log *
               restitution_coefficient_particle_log +
             9.8696);

      // Since the normal overlap is already calculated we update
      // this element of the container here. The rest of information
      // are updated using the following function
      this->update_contact_information(contact_info,
                                       this->normal_relative_velocity_value,
                                       this->normal_unit_vector,
                                       particle_one_properties,
                                       particle_two_properties,
                                       particle_one_location_3d,
                                       particle_two_location_3d,
                                       dt);

      calculate_hertz_mindlin_limit_force_contact(
        contact_info,
        this->normal_relative_velocity_value,
        this->normal_unit_vector,
        normal_overlap,
        particle_one_properties,
        particle_two_properties,
        normal_force,
        tangential_force,
        particle_one_tangential_torque,
        particle_two_tangential_torque,
        rolling_resistance_torque);
    }

  if constexpr (contact_model ==
                Parameters::Lagrangian::ParticleParticleContactForceModel::
                  hertz_mindlin_limit_overlap)
    {
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        ((particle_two.youngs_modulus *
          (1.0 - particle_one.poisson_ratio * particle_one.poisson_ratio)) +
         (particle_one.youngs_modulus *
          (1.0 - particle_two.poisson_ratio * particle_two.poisson_ratio)) +
         DBL_MIN);

      this->effective_shear_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        (particle_one.youngs_modulus * particle_two.youngs_modulus) /
        (2.0 *
           ((particle_two.youngs_modulus * (2.0 - particle_one.poisson_ratio) *
             (1.0 + particle_one.poisson_ratio)) +
            (particle_one.youngs_modulus * (2.0 - particle_two.poisson_ratio) *
             (1.0 + particle_two.poisson_ratio))) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.restitution_coefficient,
                      particle_two.restitution_coefficient);
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.friction_coefficient,
                      particle_two.friction_coefficient);
      this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] =
        harmonic_mean(particle_one.rolling_friction_coefficient,
                      particle_two.rolling_friction_coefficient);

      const double restitution_coefficient_particle_log = std::log(
        this->effective_coefficient_of_restitution[vec_particle_type_index(
          particle_one_type, particle_two_type)]);

      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] =
        restitution_coefficient_particle_log /
        sqrt(restitution_coefficient_particle_log *
               restitution_coefficient_particle_log +
             9.8696);

      // Since the normal overlap is already calculated we update
      // this element of the container here. The rest of information
      // are updated using the following function
      this->update_contact_information(contact_info,
                                       this->normal_relative_velocity_value,
                                       this->normal_unit_vector,
                                       particle_one_properties,
                                       particle_two_properties,
                                       particle_one_location_3d,
                                       particle_two_location_3d,
                                       dt);

      calculate_hertz_mindlin_limit_overlap_contact(
        contact_info,
        this->normal_relative_velocity_value,
        this->normal_unit_vector,
        normal_overlap,
        particle_one_properties,
        particle_two_properties,
        normal_force,
        tangential_force,
        particle_one_tangential_torque,
        particle_two_tangential_torque,
        rolling_resistance_torque);
    }
}



template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;

template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;


template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_force,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::
    hertz_mindlin_limit_overlap,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
