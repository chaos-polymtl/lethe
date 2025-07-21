// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/log_collision_data.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
log_collision_data(
  const DEMSolverParameters<dim> &parameters,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<dim> &ongoing_collision_log,
  CollisionEventLog<dim>   &collision_event_log)
{
  // particles_in_contact_now will be used to find the particles that ended
  // their collision during this time step, so that we can log the end of the
  // collision.
  std::set<unsigned int> particles_in_contact_now;

  // Looping over all the active particles in particle-wall pairs
  for (auto &&pairs_in_contact_content :
       particle_wall_pairs_in_contact | boost::adaptors::map_values)
    {
      // Iterating over a map which contains the required information for
      // calculation of the contact force for each particle
      for (auto &&contact_info :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          // Defining local variables which will be used within the contact
          // calculation
          auto     particle            = contact_info.particle;
          auto     particle_properties = particle->get_properties();
          Point<3> point_on_boundary   = contact_info.point_on_boundary;
          // Normal vector from the wall to the particle
          Tensor<1, 3> normal_vector = contact_info.normal_vector;

          // Getting particle 3d location
          Point<3> particle_location_3d = get_location(particle);

          // Defining a tensor which connects the point_on_boundary to the
          // center of particle
          Tensor<1, 3> point_to_particle_vector =
            particle_location_3d - point_on_boundary;

          Tensor<1, 3> projected_vector =
            find_projection(point_to_particle_vector, normal_vector);

          // Calculating the particle-wall distance using the projected vector
          double normal_overlap =
            ((particle_properties[PropertiesIndex::dp]) * 0.5) -
            (projected_vector.norm());

          types::boundary_id boundary_id = contact_info.boundary_id;
          if (boundary_id == parameters.model_parameters.wall_boundary_id ||
              parameters.model_parameters.log_collisions_with_all_walls)
            {
              if (normal_overlap > 0)
                {
                  unsigned int particle_id = particle->get_id();
                  particles_in_contact_now.insert(particle_id);

                  if (!ongoing_collision_log.is_in_collision(particle_id))
                    {
                      collision_log<dim> start_log;
                      start_log.particle_id = particle_id;

                      start_log.velocity[0] =
                        particle_properties[PropertiesIndex::v_x];
                      start_log.velocity[1] =
                        particle_properties[PropertiesIndex::v_y];
                      start_log.velocity[2] =
                        particle_properties[PropertiesIndex::v_z];

                      start_log.omega[0] =
                        particle_properties[PropertiesIndex::omega_x];
                      start_log.omega[1] =
                        particle_properties[PropertiesIndex::omega_y];
                      start_log.omega[2] =
                        particle_properties[PropertiesIndex::omega_z];

                      start_log.time        = current_time;
                      start_log.boundary_id = contact_info.boundary_id;

                      ongoing_collision_log.start_collision(start_log);

                      if (parameters.model_parameters.collision_verbosity ==
                          Parameters::Verbosity::verbose)
                        {
                          std::cout << "Collision with boundary "
                                    << start_log.boundary_id
                                    << " started for particle " << particle_id
                                    << std::endl;
                        }
                    }
                }
            }
        }
    }
  for (auto &&pairs_in_contact_content :
       particle_wall_pairs_in_contact | boost::adaptors::map_values)
    {
      // Iterating over a map which contains the required information for
      // calculation of the contact force for each particle
      for (auto &&contact_info :
           pairs_in_contact_content | boost::adaptors::map_values)

        {
          auto         particle    = contact_info.particle;
          unsigned int particle_id = particle->get_id();

          // If the particle is not in contact now, but was in contact before
          if (particles_in_contact_now.find(particle_id) ==
                particles_in_contact_now.end() &&
              ongoing_collision_log.is_in_collision(particle_id))
            {
              collision_log<dim> end_log;
              end_log.particle_id      = particle_id;
              auto particle_properties = particle->get_properties();

              end_log.velocity[0] = particle_properties[PropertiesIndex::v_x];
              end_log.velocity[1] = particle_properties[PropertiesIndex::v_y];
              end_log.velocity[2] = particle_properties[PropertiesIndex::v_z];

              end_log.omega[0] = particle_properties[PropertiesIndex::omega_x];
              end_log.omega[1] = particle_properties[PropertiesIndex::omega_y];
              end_log.omega[2] = particle_properties[PropertiesIndex::omega_z];

              end_log.time = current_time;
              collision_log<dim> start_log;
              bool               ended =
                ongoing_collision_log.end_collision(particle_id, start_log);
              if (ended)
                {
                  end_log.boundary_id = start_log.boundary_id;
                  collision_event<dim> event;
                  event.particle_id = particle_id;
                  event.start_log   = start_log;
                  event.end_log     = end_log;
                  collision_event_log.add_event(event);
                  if (parameters.model_parameters.collision_verbosity ==
                      Parameters::Verbosity::verbose)
                    {
                      std::cout << "Collision with boundary "
                                << end_log.boundary_id << " ended for particle "
                                << particle_id << std::endl;
                    }
                }
            }
        }
    }
}

template void
log_collision_data<2, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &parameters,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
                         &particle_wall_pairs_in_contact,
  const double            current_time,
  OngoingCollisionLog<2> &ongoing_collision_log,
  CollisionEventLog<2>   &collision_event_log);

template void
log_collision_data<3, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &parameters,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
                         &particle_wall_pairs_in_contact,
  const double            current_time,
  OngoingCollisionLog<3> &ongoing_collision_log,
  CollisionEventLog<3>   &collision_event_log);

template void
log_collision_data<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &parameters,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
                         &particle_wall_pairs_in_contact,
  const double            current_time,
  OngoingCollisionLog<2> &ongoing_collision_log,
  CollisionEventLog<2>   &collision_event_log);

template void
log_collision_data<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &parameters,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
                         &particle_wall_pairs_in_contact,
  const double            current_time,
  OngoingCollisionLog<3> &ongoing_collision_log,
  CollisionEventLog<3>   &collision_event_log);

template void
log_collision_data<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &parameters,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
                         &particle_wall_pairs_in_contact,
  const double            current_time,
  OngoingCollisionLog<2> &ongoing_collision_log,
  CollisionEventLog<2>   &collision_event_log);

template void
log_collision_data<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &parameters,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
                         &particle_wall_pairs_in_contact,
  const double            current_time,
  OngoingCollisionLog<3> &ongoing_collision_log,
  CollisionEventLog<3>   &collision_event_log);
