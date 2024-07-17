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
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
{
  set_effective_properties(dem_parameters);
}

template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  calculate_particle_particle_contact_force(
    DEMContactManager<dim>    &container_manager,
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

  // Define local variables which will be used within the contact calculation
  //  Namely: normal and tangential contact forces, tangential and rolling
  //  torques, normal unit vector of the contact and contact relative velocity
  //  in the normal direction
  Tensor<1, 3> normal_unit_vector;
  Tensor<1, 3> normal_force;
  Tensor<1, 3> tangential_force;
  Tensor<1, 3> particle_one_tangential_torque;
  Tensor<1, 3> particle_two_tangential_torque;
  Tensor<1, 3> rolling_resistance_torque;
  double       normal_relative_velocity_value;
  Tensor<1, 3> tangential_relative_velocity;

  // Contact forces calculations of local-local and local-ghost particle
  // pairs are performed in separate loops

  // Set the force_calculation_threshold_distance. This is useful for non-
  // contact cohesive force models such as the DMT force model. This is used in
  // every loop.
  const double force_calculation_threshold_distance =
    set_force_calculation_threshold_distance();

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

              if (normal_overlap > force_calculation_threshold_distance)
                {
                  // This means that the adjacent particles are in contact
                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    tangential_relative_velocity,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);
                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::DMT)
                    {
                      this->calculate_DMT_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(contact_info,
                                               tangential_relative_velocity,
                                               normal_relative_velocity_value,
                                               normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz_JKR)
                    {
                      this->calculate_hertz_JKR_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact<
                        contact_model>(contact_info,
                                       tangential_relative_velocity,
                                       normal_relative_velocity_value,
                                       normal_unit_vector,
                                       normal_overlap,
                                       particle_one_properties,
                                       particle_two_properties,
                                       normal_force,
                                       tangential_force,
                                       particle_one_tangential_torque,
                                       particle_two_tangential_torque,
                                       rolling_resistance_torque);
                    }


                  types::particle_index particle_two_id =
                    particle_two->get_local_index();

                  Tensor<1, 3> &particle_two_torque = torque[particle_two_id];
                  Tensor<1, 3> &particle_two_force  = force[particle_two_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_local_particles(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    particle_two_tangential_torque,
                    rolling_resistance_torque,
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

              if (normal_overlap > force_calculation_threshold_distance)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    tangential_relative_velocity,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(contact_info,
                                               tangential_relative_velocity,
                                               normal_relative_velocity_value,
                                               normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::DMT)
                    {
                      this->calculate_DMT_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz_JKR)
                    {
                      this->calculate_hertz_JKR_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact<
                        contact_model>(contact_info,
                                       tangential_relative_velocity,
                                       normal_relative_velocity_value,
                                       normal_unit_vector,
                                       normal_overlap,
                                       particle_one_properties,
                                       particle_two_properties,
                                       normal_force,
                                       tangential_force,
                                       particle_one_tangential_torque,
                                       particle_two_tangential_torque,
                                       rolling_resistance_torque);
                    }


                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_single_local_particle(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    rolling_resistance_torque,
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

              if (normal_overlap > force_calculation_threshold_distance)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    tangential_relative_velocity,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(contact_info,
                                               tangential_relative_velocity,
                                               normal_relative_velocity_value,
                                               normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::DMT)
                    {
                      this->calculate_DMT_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz_JKR)
                    {
                      this->calculate_hertz_JKR_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact<
                        contact_model>(contact_info,
                                       tangential_relative_velocity,
                                       normal_relative_velocity_value,
                                       normal_unit_vector,
                                       normal_overlap,
                                       particle_one_properties,
                                       particle_two_properties,
                                       normal_force,
                                       tangential_force,
                                       particle_one_tangential_torque,
                                       particle_two_tangential_torque,
                                       rolling_resistance_torque);
                    }


                  types::particle_index particle_two_id =
                    particle_two->get_local_index();

                  Tensor<1, 3> &particle_two_torque = torque[particle_two_id];
                  Tensor<1, 3> &particle_two_force  = force[particle_two_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_local_particles(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    particle_two_tangential_torque,
                    rolling_resistance_torque,
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

              if (normal_overlap > force_calculation_threshold_distance)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    tangential_relative_velocity,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(contact_info,
                                               tangential_relative_velocity,
                                               normal_relative_velocity_value,
                                               normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::DMT)
                    {
                      this->calculate_DMT_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz_JKR)
                    {
                      this->calculate_hertz_JKR_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
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
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact<
                        contact_model>(contact_info,
                                       tangential_relative_velocity,
                                       normal_relative_velocity_value,
                                       normal_unit_vector,
                                       normal_overlap,
                                       particle_one_properties,
                                       particle_two_properties,
                                       normal_force,
                                       tangential_force,
                                       particle_one_tangential_torque,
                                       particle_two_tangential_torque,
                                       rolling_resistance_torque);
                    }


                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_single_local_particle(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    rolling_resistance_torque,
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

              if (normal_overlap > force_calculation_threshold_distance)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated, we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    tangential_relative_velocity,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_two_properties,
                    particle_one_properties,
                    particle_two_location,
                    particle_one_location,
                    dt);
                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::DMT)
                    {
                      this->calculate_DMT_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        normal_force,
                        tangential_force,
                        particle_two_tangential_torque,
                        particle_one_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz)
                    {
                      this->calculate_hertz_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        normal_force,
                        tangential_force,
                        particle_two_tangential_torque,
                        particle_one_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::hertz_JKR)
                    {
                      this->calculate_hertz_JKR_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        normal_force,
                        tangential_force,
                        particle_two_tangential_torque,
                        particle_one_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_force)
                    {
                      calculate_hertz_mindlin_limit_force_contact(
                        contact_info,
                        tangential_relative_velocity,
                        normal_relative_velocity_value,
                        normal_unit_vector,
                        normal_overlap,
                        particle_two_properties,
                        particle_one_properties,
                        normal_force,
                        tangential_force,
                        particle_two_tangential_torque,
                        particle_one_tangential_torque,
                        rolling_resistance_torque);
                    }

                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::
                                    hertz_mindlin_limit_overlap)
                    {
                      calculate_hertz_mindlin_limit_overlap_contact<
                        contact_model>(contact_info,
                                       tangential_relative_velocity,
                                       normal_relative_velocity_value,
                                       normal_unit_vector,
                                       normal_overlap,
                                       particle_two_properties,
                                       particle_one_properties,
                                       normal_force,
                                       tangential_force,
                                       particle_two_tangential_torque,
                                       particle_one_tangential_torque,
                                       rolling_resistance_torque);
                    }
                  if constexpr (contact_model ==
                                Parameters::Lagrangian::
                                  ParticleParticleContactForceModel::linear)
                    {
                      calculate_linear_contact(contact_info,
                                               tangential_relative_velocity,
                                               normal_relative_velocity_value,
                                               normal_unit_vector,
                                               normal_overlap,
                                               particle_two_properties,
                                               particle_one_properties,
                                               normal_force,
                                               tangential_force,
                                               particle_one_tangential_torque,
                                               particle_two_tangential_torque,
                                               rolling_resistance_torque);
                    }

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_single_local_particle(
                    normal_force,
                    tangential_force,
                    particle_two_tangential_torque,
                    rolling_resistance_torque,
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

// No resistance
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
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
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
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
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
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
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
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
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleParticleContactForce<
  2,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::DMT,
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
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::hertz_JKR,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
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
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
  Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance>;
