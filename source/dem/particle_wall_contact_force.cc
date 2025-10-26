// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/auxiliary_math_functions.h>

#include <dem/contact_type.h>
#include <dem/particle_heat_transfer.h>
#include <dem/particle_wall_contact_force.h>

#include <algorithm> // for std::ranges::find
#include <ranges>
#include <vector>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  ParticleWallContactForce(const DEMSolverParameters<dim> &dem_parameters)
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
  , f_coefficient_epsd(dem_parameters.model_parameters.f_coefficient_epsd)
{
  set_effective_properties(dem_parameters);
  if constexpr (std::is_same_v<PropertiesIndex,
                               DEM::DEMMPProperties::PropertiesIndex>)
    {
      set_multiphysic_properties(dem_parameters);
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
void
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  calculate_particle_wall_contact(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                &particle_wall_pairs_in_contact,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{
  // Getting the threshold distance for contact force, this is useful
  // for non-contact cohesive force models such as the DMT.
  const double force_calculation_threshold_distance =
    get_force_calculation_threshold_distance();

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
          Tensor<1, 3> normal_force;
          Tensor<1, 3> tangential_force;
          Tensor<1, 3> tangential_torque;
          Tensor<1, 3> rolling_resistance_torque;
          double       normal_relative_velocity_value;
          Tensor<1, 3> tangential_relative_velocity;

          // Getting particle 3d location
          Point<3> particle_location_3d = get_location(particle);

          // Defining a tensor which connects the point_on_boundary to the
          // center of particle
          Tensor<1, 3> point_to_particle_vector =
            particle_location_3d - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary
          Tensor<1, 3> projected_vector =
            this->find_projection(point_to_particle_vector, normal_vector);

          // Calculating the particle-wall distance using the projected vector
          double normal_overlap =
            ((particle_properties[PropertiesIndex::dp]) * 0.5) -
            (projected_vector.norm());

          if (normal_overlap > force_calculation_threshold_distance)
            {
              // Updating contact information
              this->update_contact_information(contact_info,
                                               tangential_relative_velocity,
                                               normal_relative_velocity_value,
                                               particle_location_3d,
                                               particle_properties,
                                               dt);

              // Calculating contact force and torque
              this->calculate_contact(contact_info,
                                      tangential_relative_velocity,
                                      normal_relative_velocity_value,
                                      normal_overlap,
                                      dt,
                                      particle_properties,
                                      normal_force,
                                      tangential_force,
                                      tangential_torque,
                                      rolling_resistance_torque);

              // Applying the calculated forces and torques on the particle
              types::particle_index particle_index =
                particle->get_local_index();
              Tensor<1, 3> &particle_torque =
                contact_outcome.torque[particle_index];
              Tensor<1, 3> &particle_force =
                contact_outcome.force[particle_index];

              this->apply_force_and_torque(normal_force,
                                           tangential_force,
                                           tangential_torque,
                                           rolling_resistance_torque,
                                           particle_torque,
                                           particle_force);
            }
          else
            {
              contact_info.tangential_displacement.clear();
              contact_info.rolling_resistance_spring_torque.clear();
            }
        }
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
void
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  calculate_particle_solid_object_contact(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                &particle_floating_mesh_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{
  // Get the threshold distance for contact force, this is
  // useful for non-contact cohesive force models such as
  // the DMT.
  const double force_calculation_threshold_distance =
    get_force_calculation_threshold_distance();

  std::vector<Particles::ParticleIterator<dim>> particle_locations;
  std::vector<Point<dim>> triangle(this->vertices_per_triangle);

  struct particle_iterator_hash
  {
    std::size_t
    operator()(const typename Particles::ParticleIterator<dim> &P) const
    {
      // Use the particle's unique ID as the hash
      return std::hash<std::size_t>()(P->get_id());
    }
  };
  // For each solid surface, we create a map used to store every contact.
  // The key of that map is the particle_iterator.
  // The value is a vector of tuple storing the required information to
  // compute the contact force later on.
  // The information includes the type of contact : face_contact,
  // edge_contact, vetex_contact
  using particle_triangle_contact_description = std::vector<
    std::tuple<typename Triangulation<dim - 1, dim>::active_cell_iterator,
               Point<3>,
               Tensor<1, 3>,
               double,
               LetheGridTools::ContactIndicator>>;
  using particle_triangle_contact_record =
    ankerl::unordered_dense::map<Particles::ParticleIterator<dim>,
                                 particle_triangle_contact_description,
                                 particle_iterator_hash>;

  // Iterating over the solid objects
  for (unsigned int solid_counter = 0; solid_counter < solids.size();
       ++solid_counter)
    {
      const auto [this_solid_cp_es_neighbors,
                  this_solid_cp_vs_neighbors,
                  this_solid_np_es_neighbors,
                  this_solid_np_vs_neighbors] =
        solids[solid_counter]->get_neighbors_maps();

      particle_triangle_contact_record contact_record;

      auto &particle_floating_mesh_contact_pair =
        particle_floating_mesh_in_contact[solid_counter];

      for (auto &[triangle_cell_iterator, map_info] :
           particle_floating_mesh_contact_pair)
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
                  triangle[vertex] = triangle_cell_iterator->vertex(vertex);
                }

              // Getting the projection of particles on the triangle
              // (floating mesh cell)
              auto particle_triangle_information = LetheGridTools::
                find_particle_triangle_projection<dim, PropertiesIndex>(
                  triangle, particle_locations, n_particles);

              const auto &[pass_distance_check,
                           projection_points,
                           normal_vectors,
                           contact_indicators] = particle_triangle_information;

              unsigned int particle_counter = 0;

              for (auto &&contact_info : map_info | boost::adaptors::map_values)
                {
                  // If the particle is close enough to the triangle in
                  // the direction normal to the triangle.
                  if (pass_distance_check[particle_counter])
                    {
                      // Defining local variables which will be used to store
                      // the contact in the contact record
                      auto &particle            = contact_info.particle;
                      auto  particle_properties = particle->get_properties();

                      Point<3> particle_location_3d = get_location(particle);

                      const double particle_triangle_distance =
                        particle_location_3d.distance(
                          projection_points[particle_counter]);

                      // Find normal overlap
                      const double normal_overlap =
                        0.5 * particle_properties[PropertiesIndex::dp] -
                        particle_triangle_distance;

                      // The contact is potentially valid, thus we need to store
                      // it
                      if (normal_overlap > force_calculation_threshold_distance)
                        {
                          contact_record[particle].emplace_back(std::make_tuple(
                            triangle_cell_iterator,
                            projection_points[particle_counter],
                            normal_vectors[particle_counter],
                            normal_overlap,
                            contact_indicators[particle_counter]));
                        }
                    }
                  particle_counter++;
                } // Loop on each solid object
            }
        }
      // Every contact between this solid object and the particles are
      // recorded. Now, we need to remove the invalid contacts from
      // individual record and compute the force associated with the
      // valid contacts.

      // Loop on every particle / contact record
      for (auto &[particle, this_contact_record] : contact_record)
        {
          // Loop on every contact in the contact record
          for (auto contact_1 = this_contact_record.begin();
               contact_1 != this_contact_record.end();)
            {
              // Extract the information of C1
              auto &[T1_cell,
                     projection_point_C1,
                     normal_vector_C1,
                     normal_overlap_C1,
                     contact_indicator_C1] = *contact_1;

              // Assigning the triangle neighboring list of T1;
              std::vector<
                typename Triangulation<dim - 1, dim>::active_cell_iterator>
                T1_cp_es_neighbors, T1_cp_vs_neighbors, T1_np_es_neighbors,
                T1_np_vs_neighbors;
              T1_cp_es_neighbors = this_solid_cp_es_neighbors.at(T1_cell);
              T1_cp_vs_neighbors = this_solid_cp_vs_neighbors.at(T1_cell);
              T1_np_es_neighbors = this_solid_np_es_neighbors.at(T1_cell);
              T1_np_vs_neighbors = this_solid_np_vs_neighbors.at(T1_cell);


              // Initialize variables
              bool erase_contact_1 = false;
              auto contact_2       = std::next(contact_1);
              while (contact_2 != this_contact_record.end())
                {
                  // Extract the information of C2
                  auto &[T2_cell,
                         projection_point_C2,
                         normal_vector_C2,
                         normal_overlap_C2,
                         contact_indicator_C2] = *contact_2;

                  // First, we check if both triangle are neighbors. If they are
                  // not neighbors, C1 and C2 are automaticly valid.
                  // (Disconnected triangles)
                  if (std::ranges::find(T1_cp_es_neighbors, T2_cell) ==
                        T1_cp_es_neighbors.end() &&
                      std::ranges::find(T1_np_es_neighbors, T2_cell) ==
                        T1_np_es_neighbors.end() &&
                      std::ranges::find(T1_cp_vs_neighbors, T2_cell) ==
                        T1_cp_vs_neighbors.end() &&
                      std::ranges::find(T1_np_vs_neighbors, T2_cell) ==
                        T1_np_vs_neighbors.end())
                    continue;

                  // If both triangles are neighbord (connected), we need to
                  // check for double contacts.

                  // If C1 is a face contact
                  if (contact_indicator_C1 ==
                      LetheGridTools::ContactIndicator::face_contact)
                    {
                      // If C1 and C2 are both face contacts, they
                      // are both valid. This case is either can happen with
                      // either edge sharing or vertex sharing triangles that
                      // are non-coplanar.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::face_contact)
                        continue;

                      // If C1 is a face contact and C2 is an edge contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::edge_contact)
                        {
                          // If T1 and T2 are non-coplanar and vertex
                          // sharing, C1 and C2 are valid. This could happen
                          // with two concave triangles.
                          if (std::ranges::find(T1_np_vs_neighbors, T2_cell) !=
                              T1_np_vs_neighbors.end())
                            continue;

                          // Otherwise, C2 is invalid and is erased from the
                          // current contact record.
                          else
                            {
                              contact_2 = this_contact_record.erase(contact_2);
                              continue;
                            }
                        }
                      else // If C1 is a face contact and C2 is a vertex
                        // contact.
                        {
                          // C2 is always invalid
                          contact_2 = this_contact_record.erase(contact_2);
                          continue;
                        }
                    }
                  // If C1 is an edge contact.
                  if (contact_indicator_C1 ==
                      LetheGridTools::ContactIndicator::edge_contact)
                    {
                      // If C1 is an edge contact and C2 is a face contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::face_contact)
                        {
                          // C2 is valid.
                          // C1 is invalid if T1 and T2 are connected. We
                          // already checked at the beginning of the while loop
                          // that this is the case, thus we know if we enter
                          // this if that both triangle are connected.
                          erase_contact_1 = true;
                          break; // must exit the while loop
                        }
                      // If C1 and C2 are edge contacts,
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::edge_contact)
                        {
                          // C1 is always valid.
                          // C2 is invalid if T1 and T2 are edge sharing. This
                          // case means that both conta // if T1 and T2 share an
                          // edge, C2 is invalid. This means that both contacts
                          // occur on the same edge which is connecting both
                          // triangles.
                          if (std::ranges::find(T1_cp_es_neighbors, T2_cell) !=
                                T1_cp_es_neighbors.end() ||
                              std::ranges::find(T1_np_es_neighbors, T2_cell) !=
                                T1_np_es_neighbors.end())
                            {
                              contact_2 = this_contact_record.erase(contact_2);
                              continue;
                            }

                          // If they are vertex sharing, this means that each
                          // contact is occurring on two different edges. This
                          // is not a double contact, C2 is valid
                          continue;
                        }
                    }
                  // If C1 is an vertex contact
                  if (contact_indicator_C1 ==
                      LetheGridTools::ContactIndicator::vertex_contact)
                    {
                      // If C1 is a vertex contact and C2 is a face contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::face_contact)
                        {
                          // C2 is valid.
                          // C1 is invalid if T1 and T2 are connected. We
                          // already checked at the beginning of the while loop
                          // that this is the case, thus we know if we enter
                          // this if that both triangle are connected.
                          erase_contact_1 = true;
                          break; // must exit the while loop
                        }
                      // If C1 is a vertex contact and C2 is an edge contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::edge_contact)
                        {
                          // C2 is valid
                          // C1 is invalid if T1 and T2 are vertex sharing. This
                          // case means that both contacts occur on the same
                          // vertex which is connecting both triangles.
                          if (std::ranges::find(T1_cp_vs_neighbors, T2_cell) !=
                                T1_cp_vs_neighbors.end() ||
                              std::ranges::find(T1_np_vs_neighbors, T2_cell) !=
                                T1_np_vs_neighbors.end())
                            {
                              erase_contact_1 = true;
                              break;
                            }
                        }
                      // If C1 and C2 are vertex contacts,
                      if (contact_indicator_C2 ==
                          LetheGridTools::ContactIndicator::vertex_contact)
                        {
                          // C2 is invalid
                          contact_2 = this_contact_record.erase(contact_2);
                          continue;
                        }
                    }
                  ++contact_2;
                } // While loop

              // Erase C1 from the contact record is marked as true.
              if (erase_contact_1)
                {
                  contact_1 = this_contact_record.erase(contact_1);
                  continue;
                }
              ++contact_1;
            } // For loop


          // At this point, every contact in the contact record are valid,
          // thus we can compute and applied the force and torques.
          // Get translational and rotational velocities and center of
          // rotation of solid object
          Tensor<1, 3> translational_velocity =
            solids[solid_counter]->get_translational_velocity();
          Tensor<1, 3> angular_velocity =
            solids[solid_counter]->get_angular_velocity();
          Point<3> center_of_rotation =
            solids[solid_counter]->get_center_of_rotation();
          Parameters::ThermalBoundaryType thermal_boundary_type =
            solids[solid_counter]->get_thermal_boundary_type();

          // We need the contact info of the current solid surface
          // auto &particle_floating_mesh_contact_pair =
          // particle_floating_mesh_in_contact[solid_counter];



          for (auto contact = this_contact_record.begin();
               contact != this_contact_record.end();
               contact++)
            {
              // Extract the information of contact
              auto &[triangle_cell,
                     projection_point,
                     normal_vector,
                     normal_overlap,
                     contact_indicator] = *contact;

              auto map_p_index_to_c_info =
                particle_floating_mesh_contact_pair[triangle_cell];

              auto &contact_info =
                map_p_index_to_c_info.find(particle->get_local_index())->second;

              // Defining local variables which will be used within the
              // contact calculation
              auto         particle_properties = particle->get_properties();
              Tensor<1, 3> normal_force;
              Tensor<1, 3> tangential_force;
              Tensor<1, 3> tangential_torque;
              Tensor<1, 3> rolling_resistance_torque;
              double       normal_relative_velocity_value;
              Tensor<1, 3> tangential_relative_velocity;

              Point<3> particle_location_3d = get_location(particle);

              contact_info.normal_vector     = normal_vector;
              contact_info.point_on_boundary = projection_point;
              contact_info.boundary_id       = solid_counter;

              // Updating contact information
              this->update_particle_solid_object_contact_information(
                contact_info,
                tangential_relative_velocity,
                normal_relative_velocity_value,
                particle_properties,
                dt,
                translational_velocity,
                angular_velocity,
                center_of_rotation.distance(particle_location_3d));

              // Calculating contact force and torque
              this->calculate_contact(contact_info,
                                      tangential_relative_velocity,
                                      normal_relative_velocity_value,
                                      normal_overlap,
                                      dt,
                                      particle_properties,
                                      normal_force,
                                      tangential_force,
                                      tangential_torque,
                                      rolling_resistance_torque);

              // Applying the calculated forces and torques on
              // the particle
              types::particle_index particle_id = particle->get_local_index();
              Tensor<1, 3>         &particle_torque =
                contact_outcome.torque[particle_id];
              Tensor<1, 3> &particle_force = contact_outcome.force[particle_id];

              this->apply_force_and_torque(normal_force,
                                           tangential_force,
                                           tangential_torque,
                                           rolling_resistance_torque,
                                           particle_torque,
                                           particle_force);

              if constexpr (std::is_same_v<
                              PropertiesIndex,
                              DEM::DEMMPProperties::PropertiesIndex>)
                {
                  double temperature_wall;
                  if (thermal_boundary_type !=
                      Parameters::ThermalBoundaryType::adiabatic)
                    temperature_wall = solids[solid_counter]->get_temperature();
                  if ((thermal_boundary_type !=
                       Parameters::ThermalBoundaryType::adiabatic) &&
                      (normal_overlap > 0))
                    {

                      const unsigned int particle_type =
                        static_cast<unsigned int>(
                          particle_properties[PropertiesIndex::type]);
                      const double temperature_particle =
                        particle_properties[PropertiesIndex::T];
                      double &particle_heat_transfer_rate =
                        contact_outcome
                          .heat_transfer_rate[particle->get_local_index()];
                      double thermal_conductance;

                      calculate_contact_thermal_conductance<
                        ContactType::particle_floating_mesh>(
                        0.5 * particle_properties[PropertiesIndex::dp],
                        0,
                        this->effective_youngs_modulus[particle_type],
                        this->effective_real_youngs_modulus[particle_type],
                        this->equivalent_surface_roughness[particle_type],
                        this->equivalent_surface_slope[particle_type],
                        this->effective_microhardness[particle_type],
                        this->particle_thermal_conductivity[particle_type],
                        this->wall_thermal_conductivity,
                        this->gas_thermal_conductivity,
                        this->gas_parameter_m[particle_type],
                        normal_overlap,
                        normal_force.norm(),
                        thermal_conductance);

                      // Apply the heat transfer to the particle
                      apply_heat_transfer_on_single_local_particle(
                        temperature_particle,
                        temperature_wall,
                        thermal_conductance,
                        particle_heat_transfer_rate);
                    }
                }
            }
          //     else
          //       {
          //         contact_info.tangential_displacement.clear();
          //         contact_info.rolling_resistance_spring_torque.clear();
          //       }
        }
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
void
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  set_effective_properties(const DEMSolverParameters<dim> &dem_parameters)
{
  auto properties = dem_parameters.lagrangian_physical_properties;

  n_particle_types = properties.particle_type_number;
  effective_youngs_modulus.resize(n_particle_types);
  effective_shear_modulus.resize(n_particle_types);
  effective_coefficient_of_restitution.resize(n_particle_types);
  effective_coefficient_of_friction.resize(n_particle_types);
  effective_coefficient_of_rolling_viscous_damping.resize(n_particle_types);
  effective_coefficient_of_rolling_friction.resize(n_particle_types);
  model_parameter_beta.resize(n_particle_types);
  effective_surface_energy.resize(n_particle_types);
  effective_hamaker_constant.resize(n_particle_types);

  // Initialize wall variables and boundary conditions
  this->center_mass_container = dem_parameters.forces_torques.point_center_mass;
  this->boundary_translational_velocity_map =
    dem_parameters.boundary_conditions.boundary_translational_velocity;
  this->boundary_rotational_speed_map =
    dem_parameters.boundary_conditions.boundary_rotational_speed;
  this->boundary_rotational_vector =
    dem_parameters.boundary_conditions.boundary_rotational_vector;
  this->point_on_rotation_vector =
    dem_parameters.boundary_conditions.point_on_rotation_axis;

  // Wall properties
  const double wall_youngs_modulus = properties.youngs_modulus_wall;
  const double wall_poisson_ratio  = properties.poisson_ratio_wall;
  const double wall_restitution_coefficient =
    properties.restitution_coefficient_wall;
  const double wall_friction_coefficient = properties.friction_coefficient_wall;
  const double wall_rolling_friction_coefficient =
    properties.rolling_friction_wall;
  const double wall_rolling_viscous_damping =
    properties.rolling_viscous_damping_wall;
  const double wall_surface_energy   = properties.surface_energy_wall;
  const double wall_hamaker_constant = properties.hamaker_constant_wall;

  for (unsigned int i = 0; i < n_particle_types; ++i)
    {
      // Particle properties
      const double particle_youngs_modulus =
        properties.youngs_modulus_particle.at(i);
      const double particle_poisson_ratio =
        properties.poisson_ratio_particle.at(i);
      const double particle_restitution_coefficient =
        properties.restitution_coefficient_particle.at(i);
      const double particle_friction_coefficient =
        properties.friction_coefficient_particle.at(i);
      const double particle_rolling_friction_coefficient =
        properties.rolling_friction_coefficient_particle.at(i);
      const double particle_rolling_viscous_damping_coefficient =
        properties.rolling_viscous_damping_coefficient_particle.at(i);
      const double particle_surface_energy =
        properties.surface_energy_particle.at(i);
      const double particle_hamaker_constant =
        properties.hamaker_constant_particle.at(i);

      // Effective particle-wall properties.
      this->effective_youngs_modulus[i] =
        (particle_youngs_modulus * wall_youngs_modulus) /
        (wall_youngs_modulus *
           (1. - particle_poisson_ratio * particle_poisson_ratio) +
         particle_youngs_modulus *
           (1. - wall_poisson_ratio * wall_poisson_ratio) +
         DBL_MIN);

      this->effective_shear_modulus[i] =
        (particle_youngs_modulus * wall_youngs_modulus) /
        ((2. * wall_youngs_modulus * (2. - particle_poisson_ratio) *
          (1. + particle_poisson_ratio)) +
         (2. * particle_youngs_modulus * (2. - wall_poisson_ratio) *
          (1. + wall_poisson_ratio)) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[i] =
        harmonic_mean(particle_restitution_coefficient,
                      wall_restitution_coefficient);

      this->effective_coefficient_of_friction[i] =
        harmonic_mean(particle_friction_coefficient, wall_friction_coefficient);

      this->effective_coefficient_of_rolling_friction[i] =
        harmonic_mean(particle_rolling_friction_coefficient,
                      wall_rolling_friction_coefficient);

      this->effective_coefficient_of_rolling_viscous_damping[i] =
        harmonic_mean(particle_rolling_viscous_damping_coefficient,
                      wall_rolling_viscous_damping);

      this->effective_surface_energy[i] =
        particle_surface_energy + wall_surface_energy -
        std::pow(std::sqrt(particle_surface_energy) -
                   std::sqrt(wall_surface_energy),
                 2);

      this->effective_hamaker_constant[i] =
        0.5 * (particle_hamaker_constant + wall_hamaker_constant);

      const double log_coeff_restitution =
        std::log(this->effective_coefficient_of_restitution[i]);
      this->model_parameter_beta[i] =
        log_coeff_restitution /
        sqrt((log_coeff_restitution * log_coeff_restitution) + 9.8696);
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleWallContactForceModel contact_model,
          RollingResistanceMethod       rolling_friction_model>
void
ParticleWallContactForce<dim,
                         PropertiesIndex,
                         contact_model,
                         rolling_friction_model>::
  set_multiphysic_properties(const DEMSolverParameters<dim> &dem_parameters)
{
  auto properties = dem_parameters.lagrangian_physical_properties;

  n_particle_types = properties.particle_type_number;
  effective_real_youngs_modulus.resize(n_particle_types);
  equivalent_surface_roughness.resize(n_particle_types);
  equivalent_surface_slope.resize(n_particle_types);
  effective_microhardness.resize(n_particle_types);
  particle_thermal_conductivity.resize(n_particle_types);
  gas_parameter_m.resize(n_particle_types);
  this->gas_thermal_conductivity = properties.thermal_conductivity_gas;

  // Wall properties
  const double wall_real_youngs_modulus = properties.real_youngs_modulus_wall;
  const double wall_poisson_ratio       = properties.poisson_ratio_wall;
  const double wall_surface_roughness   = properties.surface_roughness_wall;
  const double wall_surface_slope       = properties.surface_slope_wall;
  const double wall_microhardness       = properties.microhardness_wall;
  const double wall_thermal_accommodation =
    properties.thermal_accommodation_wall;
  this->wall_thermal_conductivity = properties.thermal_conductivity_wall;

  for (unsigned int i = 0; i < n_particle_types; ++i)
    {
      // Particle properties
      const double particle_real_youngs_modulus =
        properties.real_youngs_modulus_particle.at(i);
      const double particle_poisson_ratio =
        properties.poisson_ratio_particle.at(i);
      const double particle_surface_roughness =
        properties.surface_roughness_particle.at(i);
      const double particle_surface_slope =
        properties.surface_slope_particle.at(i);
      const double particle_microhardness =
        properties.microhardness_particle.at(i);
      const double particle_thermal_accommodation =
        properties.thermal_accommodation_particle.at(i);
      this->particle_thermal_conductivity[i] =
        properties.thermal_conductivity_particle.at(i);

      // Effective particle-wall properties
      this->effective_real_youngs_modulus[i] =
        (particle_real_youngs_modulus * wall_real_youngs_modulus) /
        (wall_real_youngs_modulus *
           (1. - particle_poisson_ratio * particle_poisson_ratio) +
         particle_real_youngs_modulus *
           (1. - wall_poisson_ratio * wall_poisson_ratio) +
         DBL_MIN);
      this->equivalent_surface_roughness[i] =
        sqrt(particle_surface_roughness * particle_surface_roughness +
             wall_surface_roughness * wall_surface_roughness);
      this->equivalent_surface_slope[i] =
        sqrt(particle_surface_slope * particle_surface_slope +
             wall_surface_slope * wall_surface_slope);
      this->effective_microhardness[i] =
        harmonic_mean(particle_microhardness, wall_microhardness);
      this->gas_parameter_m[i] =
        ((2. - particle_thermal_accommodation) /
           particle_thermal_accommodation +
         (2. - wall_thermal_accommodation) / wall_thermal_accommodation) *
        (2. * properties.specific_heats_ratio_gas) /
        (1. + properties.specific_heats_ratio_gas) *
        properties.molecular_mean_free_path_gas /
        (properties.dynamic_viscosity_gas * properties.specific_heat_gas /
         properties.thermal_conductivity_gas);
    }
}

// dem
// No resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::none>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd>;

// cfd_dem
// No resistance
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::none>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd>;

// dem_mp
//  No resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::none>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::none>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::none>;

// Constant resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::constant>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant>;

// Viscous resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous>;

// EPSD resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd>;
