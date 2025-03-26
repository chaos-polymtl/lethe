// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/lethe_grid_tools.h>

#include <dem/particle_wall_dmt_force.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
ParticleWallDMTForce<dim, PropertiesIndex>::ParticleWallDMTForce(
  const DEMSolverParameters<dim>        &dem_parameters,
  const std::vector<types::boundary_id> &boundary_index)
  : ParticleWallNonLinearForce<dim, PropertiesIndex>(dem_parameters)
  , dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
{
  initialize_particle_wall_properties(dem_parameters);
  if (dem_parameters.model_parameters.rolling_resistance_method ==
      Parameters::Lagrangian::RollingResistanceMethod::no_resistance)
    {
      calculate_rolling_resistance_torque =
        &ParticleWallDMTForce<dim, PropertiesIndex>::no_resistance;
    }
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::RollingResistanceMethod::constant_resistance)
    {
      calculate_rolling_resistance_torque =
        &ParticleWallDMTForce<dim, PropertiesIndex>::constant_resistance;
    }
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::RollingResistanceMethod::viscous_resistance)
    {
      calculate_rolling_resistance_torque =
        &ParticleWallDMTForce<dim, PropertiesIndex>::viscous_resistance;
    }
  this->calculate_force_torque_on_boundary =
    dem_parameters.forces_torques.calculate_force_torque;
  this->center_mass_container = dem_parameters.forces_torques.point_center_mass;
  this->boundary_index        = boundary_index;
  this->force_on_walls        = this->initialize();
  this->torque_on_walls       = this->initialize();
}

template <int dim, typename PropertiesIndex>
void
ParticleWallDMTForce<dim, PropertiesIndex>::
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                              &particle_wall_pairs_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force)
{
  constexpr double M_2PI = 6.283185307179586; // 2. * M_PI

  ParticleWallContactForce<dim, PropertiesIndex>::force_on_walls =
    ParticleWallContactForce<dim, PropertiesIndex>::initialize();
  ParticleWallContactForce<dim, PropertiesIndex>::torque_on_walls =
    ParticleWallContactForce<dim, PropertiesIndex>::initialize();

  // Set the force_calculation_threshold_distance. This is useful for non-
  // contact cohesive force models such as the DMT force model.
  const double force_calculation_threshold_distance =
    set_dmt_cut_off_distance();

  // Looping over particle_wall_pairs_in_contact, which means looping over all
  // the active particles with iterator particle_wall_pairs_in_contact_iterator
  for (auto &&pairs_in_contact_content :
       particle_wall_pairs_in_contact | boost::adaptors::map_values)
    {
      // Now an iterator (particle_wall_contact_information_iterator) on each
      // element of the particle_wall_pairs_in_contact vector is defined. This
      // iterator iterates over a map which contains the required information
      // for calculation of the contact force for each particle
      for (auto &&contact_information :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          // Defining the total force of contact, properties of particle as
          // local parameters
          auto particle            = contact_information.particle;
          auto particle_properties = particle->get_properties();

          auto normal_vector     = contact_information.normal_vector;
          auto point_on_boundary = contact_information.point_on_boundary;

          Point<3> particle_location_3d = [&] {
            if constexpr (dim == 3)
              {
                return particle->get_location();
              }
            if constexpr (dim == 2)
              {
                return (point_nd_to_3d(particle->get_location()));
              }
          }();

          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance
          Tensor<1, 3> point_to_particle_vector =
            particle_location_3d - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, 3> projected_vector =
            this->find_projection(point_to_particle_vector, normal_vector);

          double normal_overlap =
            ((particle_properties[PropertiesIndex::dp]) * 0.5) -
            (projected_vector.norm());

          // Minimal delta_star. We know a force has to be computed.
          if (normal_overlap > force_calculation_threshold_distance)
            {
              // i is the particle, j is the wall.
              // we need to put a minus sign in front of the
              // normal_vector to respect the convention (i ->
              // j)
              const unsigned int particle_type =
                particle_properties[PropertiesIndex::type];
              const double effective_radius =
                0.5 * particle_properties[PropertiesIndex::dp];
              const double effective_surface_energy =
                this->effective_surface_energy[particle_type];
              const double effective_hamaker_constant =
                this->effective_hamaker_constant[particle_type];

              const double F_po =
                M_2PI * effective_radius * effective_surface_energy;

              const double delta_0 = -std::sqrt(effective_hamaker_constant *
                                                effective_radius / (6. * F_po));
              // Cohesive force. This will need to be added to the
              // first vector inside the tuple.
              Tensor<1, 3> cohesive_force;

              // This tuple (forces and torques) contains four
              // elements which are: 1, normal force, 2,
              // tangential force, 3, tangential torque and 4,
              // rolling resistance torque, respectively
              std::tuple<Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>>
                forces_and_torques;

              forces_and_torques = std::make_tuple(Tensor<1, 3>(),
                                                   Tensor<1, 3>(),
                                                   Tensor<1, 3>(),
                                                   Tensor<1, 3>());
              // Contact + constant cohesive force.
              if (normal_overlap > 0.)
                {
                  // i is the particle, j is the wall.
                  // we need to put a minus sign infront of the
                  // normal_vector to respect the convention (i ->
                  // j)
                  cohesive_force = -F_po * -normal_vector;

                  contact_information.normal_overlap = normal_overlap;
                  this->update_contact_information(contact_information,
                                                   particle_location_3d,
                                                   particle_properties,
                                                   dt);

                  // This tuple (forces and torques) contains
                  // four elements which are: 1, normal force,
                  // 2, tangential force, 3, tangential torque
                  // and 4, rolling resistance torque,
                  // respectively
                  forces_and_torques =
                    this->calculate_nonlinear_contact_force_and_torque(
                      dt, contact_information, particle_properties);
                }
              // No contact, but still in the constant zone for the cohesive
              // force.
              else if (normal_overlap > delta_0)
                {
                  // i is the particle, j is the wall.
                  // we need to put a minus sign infront of the
                  // normal_vector to respect the convention (i ->
                  // j)
                  cohesive_force = -F_po * -normal_vector;

                  // No contact means there is no tangential
                  // overlap.
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_information.tangential_overlap[d] = 0.;
                    }
                }
              else // between delta_* and delta_0.
                {
                  // i is the particle, j is the wall.
                  // we need to put a minus sign infront of the
                  // normal_vector to respect the convention (i ->
                  // j)
                  const double F_cohesion =
                    effective_hamaker_constant * effective_radius /
                    (6. * Utilities::fixed_power<2>(normal_overlap));

                  // No contact means there is not tangential
                  // overlap.
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_information.tangential_overlap[d] = 0.;
                    }
                  cohesive_force = -F_cohesion * -normal_vector;
                }

              // Get particle's torque and force
              types::particle_index particle_id = particle->get_local_index();

              Tensor<1, 3> &particle_torque = torque[particle_id];
              Tensor<1, 3> &particle_force  = force[particle_id];

              // Added the cohesive term
              std::get<0>(forces_and_torques) += cohesive_force;

              // Apply the calculated forces and torques on the particle
              this->apply_force_and_torque(forces_and_torques,
                                           particle_torque,
                                           particle_force,
                                           point_on_boundary,
                                           contact_information.boundary_id);
            }
          else // No cohesive force
            {
              // No contact means there is not tangential
              // overlap.
              for (int d = 0; d < dim; ++d)
                {
                  contact_information.tangential_overlap[d] = 0.;
                }
            }
        }
    }
}


template <int dim, typename PropertiesIndex>
void
ParticleWallDMTForce<dim, PropertiesIndex>::
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                              &particle_floating_mesh_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids)
{
  constexpr double M_2PI = 6.283185307179586; // 2. * M_PI

  std::vector<Particles::ParticleIterator<dim>> particle_locations;
  std::vector<Point<dim>> triangle(this->vertices_per_triangle);

  // Set the force_calculation_threshold_distance. This is useful for non-
  // contact cohesive force models such as the DMT force model.
  const double force_calculation_threshold_distance =
    set_dmt_cut_off_distance();

  for (unsigned int solid_counter = 0; solid_counter < solids.size();
       ++solid_counter)
    {
      // Get translational and rotational velocities and
      // center of rotation
      Tensor<1, 3> translational_velocity =
        solids[solid_counter]->get_translational_velocity();
      Tensor<1, 3> angular_velocity =
        solids[solid_counter]->get_angular_velocity();
      Point<3> center_of_rotation =
        solids[solid_counter]->get_center_of_rotation();

      auto &particle_floating_mesh_contact_pair =
        particle_floating_mesh_in_contact[solid_counter];

      for (auto &[cut_cell, map_info] : particle_floating_mesh_contact_pair)
        {
          if (!map_info.empty())
            {
              // Clear the particle locations vector for the new cut cell
              particle_locations.clear();
              const unsigned int n_particles = map_info.size();

              // Gather all the particles locations in a vector
              for (auto &&contact_info : map_info | boost::adaptors::map_values)
                {
                  particle_locations.push_back(contact_info.particle);
                }

              // Build triangle vector
              for (unsigned int vertex = 0;
                   vertex < this->vertices_per_triangle;
                   ++vertex)
                {
                  // Find vertex-floating wall distance
                  triangle[vertex] = cut_cell->vertex(vertex);
                }

              // Call find_particle_triangle_projection to get the
              // distance and projection of particles on the triangle
              // (floating mesh cell)
              auto particle_triangle_information = LetheGridTools::
                find_particle_triangle_projection<dim, PropertiesIndex>(
                  triangle, particle_locations, n_particles);

              const std::vector<bool> pass_distance_check =
                std::get<0>(particle_triangle_information);
              const std::vector<Point<3>> projection_points =
                std::get<1>(particle_triangle_information);
              const std::vector<Tensor<1, 3>> normal_vectors =
                std::get<2>(particle_triangle_information);

              unsigned int particle_counter = 0;

              for (auto &&contact_info : map_info | boost::adaptors::map_values)
                {
                  // If particle passes the distance check
                  if (pass_distance_check[particle_counter])
                    {
                      // Define the total force of contact, properties of
                      // particle as local parameters
                      auto &particle            = contact_info.particle;
                      auto  particle_properties = particle->get_properties();

                      const Point<3> &projection_point =
                        projection_points[particle_counter];

                      Point<3> particle_location_3d;

                      if constexpr (dim == 3)
                        particle_location_3d = particle->get_location();

                      if constexpr (dim == 2)
                        particle_location_3d =
                          point_nd_to_3d(particle->get_location());

                      const double particle_triangle_distance =
                        particle_location_3d.distance(projection_point);

                      // Find normal overlap
                      double normal_overlap =
                        ((particle_properties[PropertiesIndex::dp]) * 0.5) -
                        particle_triangle_distance;

                      // We check if a force need to be computed.
                      if (normal_overlap > force_calculation_threshold_distance)
                        {
                          // i is the particle, j is the wall.
                          // we need to put a minus sign in front of the
                          // normal_vector to respect the convention (i ->
                          // j)
                          const unsigned int particle_type =
                            particle_properties[PropertiesIndex::type];
                          const double effective_radius =
                            0.5 * particle_properties[PropertiesIndex::dp];
                          const double effective_surface_energy =
                            this->effective_surface_energy[particle_type];
                          const double effective_hamaker_constant =
                            this->effective_hamaker_constant[particle_type];

                          const double F_po =
                            M_2PI * effective_radius * effective_surface_energy;

                          const double delta_0 =
                            -std::sqrt(effective_hamaker_constant *
                                       effective_radius / (6. * F_po));

                          contact_info.normal_overlap = normal_overlap;

                          contact_info.normal_vector =
                            normal_vectors[particle_counter];

                          contact_info.point_on_boundary = projection_point;

                          contact_info.boundary_id = solid_counter;

                          this
                            ->update_particle_floating_wall_contact_information(
                              contact_info,
                              particle_properties,
                              dt,
                              translational_velocity,
                              angular_velocity,
                              center_of_rotation.distance(
                                particle_location_3d));

                          // Cohesive force. This will need to be added to the
                          // first vector inside the tuple.
                          Tensor<1, 3> cohesive_force;

                          // This tuple (forces and torques) contains four
                          // elements which are: 1, normal force, 2,
                          // tangential force, 3, tangential torque and 4,
                          // rolling resistance torque, respectively
                          std::tuple<Tensor<1, 3>,
                                     Tensor<1, 3>,
                                     Tensor<1, 3>,
                                     Tensor<1, 3>>
                            forces_and_torques;

                          forces_and_torques = std::make_tuple(Tensor<1, 3>(),
                                                               Tensor<1, 3>(),
                                                               Tensor<1, 3>(),
                                                               Tensor<1, 3>());

                          // Contact + constant cohesive force.
                          if (normal_overlap > 0.)
                            {
                              // i is the particle, j is the wall.
                              // we need to put a minus sign infront of the
                              // normal_vector to respect the convention (i ->
                              // j)
                              cohesive_force =
                                -F_po * -normal_vectors[particle_counter];

                              // This tuple (forces and torques) contains
                              // four elements which are: 1, normal force,
                              // 2, tangential force, 3, tangential torque
                              // and 4, rolling resistance torque,
                              // respectively
                              forces_and_torques =
                                this
                                  ->calculate_nonlinear_contact_force_and_torque(
                                    dt, contact_info, particle_properties);
                            }
                          // No contact, but still in the constant zone for the
                          // cohesive force.
                          else if (normal_overlap > delta_0)
                            {
                              // i is the particle, j is the wall.
                              // we need to put a minus sign infront of the
                              // normal_vector to respect the convention (i ->
                              // j)
                              cohesive_force =
                                -F_po * -normal_vectors[particle_counter];

                              // No contact means there is no tangential
                              // overlap.
                              for (int d = 0; d < dim; ++d)
                                {
                                  contact_info.tangential_overlap[d] = 0.;
                                  contact_info
                                    .rolling_resistance_spring_torque[d] = 0.;
                                }
                            }
                          else // between delta_* and delta_0.
                            {
                              // i is the particle, j is the wall.
                              // we need to put a minus sign infront of the
                              // normal_vector to respect the convention (i ->
                              // j)
                              const double F_cohesion =
                                effective_hamaker_constant * effective_radius /
                                (6. *
                                 Utilities::fixed_power<2>(normal_overlap));

                              // No contact means there is not tangential
                              // overlap.
                              for (int d = 0; d < dim; ++d)
                                {
                                  contact_info.tangential_overlap[d] = 0.;
                                  contact_info
                                    .rolling_resistance_spring_torque[d] = 0.;
                                }
                              cohesive_force =
                                -F_cohesion * -normal_vectors[particle_counter];
                            }
                          // Get particle's torque and force
                          types::particle_index particle_id =
                            particle->get_local_index();

                          Tensor<1, 3> &particle_torque = torque[particle_id];
                          Tensor<1, 3> &particle_force  = force[particle_id];

                          // Added the cohesive term
                          std::get<0>(forces_and_torques) += cohesive_force;

                          // Apply the calculated forces and torques on the
                          // particle
                          this->apply_force_and_torque(
                            forces_and_torques,
                            particle_torque,
                            particle_force,
                            projection_point,
                            contact_info.boundary_id);
                        }
                      else // No cohesive force
                        {
                          // No contact means there is not tangential
                          // overlap.
                          for (int d = 0; d < dim; ++d)
                            {
                              contact_info.tangential_overlap[d] = 0.;
                              contact_info.rolling_resistance_spring_torque[d] =
                                0.;
                            }
                        }
                      particle_counter++;
                    }
                }
            }
        }
    }
}
template class ParticleWallDMTForce<2, DEM::DEMProperties::PropertiesIndex>;
template class ParticleWallDMTForce<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class ParticleWallDMTForce<2, DEM::DEMMPProperties::PropertiesIndex>;
template class ParticleWallDMTForce<3, DEM::DEMProperties::PropertiesIndex>;
template class ParticleWallDMTForce<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class ParticleWallDMTForce<3, DEM::DEMMPProperties::PropertiesIndex>;
