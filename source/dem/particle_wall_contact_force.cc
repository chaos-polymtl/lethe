// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/auxiliary_math_functions.h>

#include <dem/contact_type.h>
#include <dem/particle_heat_transfer.h>
#include <dem/particle_wall_contact_force.h>

#include <algorithm>
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
                               DEMMPProperties::PropertiesIndex>)
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
    typename dem_data_structures<dim>::particle_wall_in_contact
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
              Tensor<1, 3> tangential_relative_velocity;
              double       normal_relative_velocity_value;
              Tensor<1, 3> normal_force;
              Tensor<1, 3> tangential_force;
              Tensor<1, 3> tangential_torque;
              Tensor<1, 3> rolling_resistance_torque;

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
          else // If there is no contact (or interaction), we need to clear the
               // tangential displacement and rolling resistance torque
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
    typename dem_data_structures<
      dim>::particle_floating_mesh_potentially_in_contact
                &particle_floating_mesh_potentially_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{
  // Get the threshold distance for contact force, this is
  // useful for non-contact cohesive force models such as
  // the DMT.
  const double force_calculation_threshold_distance =
    get_force_calculation_threshold_distance();

  // Initiate containers
  std::vector<Point<dim>> triangle(this->vertices_per_triangle);

  // Iterating over the solid objects
  for (unsigned int solid_counter = 0; solid_counter < solids.size();
       ++solid_counter)
    {
      // Extract the neighboring cell maps of this solid objects
      const auto &[this_solid_es_neighbors, this_solid_vs_neighbors] =
        solids[solid_counter]->get_neighbors_maps();

      // For each solid surface, we create a map, called the contact_record,
      // used to store every contact. The key of that map is the particle local
      // ID. The value is a vector of tuple storing the required information to
      // compute the contact force later on.
      // The information includes:
      // 1. The triangle cell with which the contact is occurring,
      // 2. The normal overlap,
      // 3. The type of contact (face, edge or vertex)
      // 4. The contact info associated.
      particle_triangle_contact_record contact_record;

      typename dem_data_structures<
        dim>::particle_triangle_cell_from_mesh_potentially_in_contact
        &particle_floating_mesh_potential_contact_pair =
          particle_floating_mesh_potentially_in_contact[solid_counter];

      // Loop over every triangle of the solid object
      for (auto &[triangle_cell_iterator, map_info] :
           particle_floating_mesh_potential_contact_pair)
        {
          if (map_info.empty())
            continue;

          // Build triangle vector from the triangle cell vertices
          for (unsigned int vertex = 0; vertex < this->vertices_per_triangle;
               ++vertex)
            triangle[vertex] = triangle_cell_iterator->vertex(vertex);

          // We reserve the contact record of this triangle an arbitrary number
          // of contacts
          contact_record.reserve(map_info.size());

          // Loop on every particle using their contact info
          for (auto &&contact_info : map_info | boost::adaptors::map_values)
            {
              // We check the contact between the triangle and the particle.
              auto particle_triangle_information = LetheGridTools::
                find_particle_triangle_projection<dim, PropertiesIndex>(
                  triangle, contact_info.particle);

              const auto &[pass_distance_check,
                           projection_point,
                           normal_vector,
                           contact_indicator] = particle_triangle_information;

              // If the particle is close enough to the triangle in
              // the direction normal to the triangle.
              if (pass_distance_check)
                {
                  // Defining local variables which will be used to store
                  // the contact in the contact record
                  auto        &particle            = contact_info.particle;
                  auto         particle_properties = particle->get_properties();
                  Point<3>     particle_location_3d = get_location(particle);
                  const double particle_triangle_distance =
                    particle_location_3d.distance(projection_point);

                  // Find normal overlap
                  const double normal_overlap =
                    0.5 * particle_properties[PropertiesIndex::dp] -
                    particle_triangle_distance;

                  // The contact is potentially valid, thus we need to store
                  // it in the contact record.
                  if (normal_overlap > force_calculation_threshold_distance)
                    {
                      // Update information in the contact_info
                      contact_info.normal_vector     = normal_vector;
                      contact_info.point_on_boundary = projection_point;
                      contact_info.boundary_id       = solid_counter;

                      contact_record[particle->get_local_index()].emplace_back(
                        triangle_cell_iterator,
                        normal_overlap,
                        contact_indicator,
                        contact_info);
                    }
                }
            }
        }
      // Every contact between this solid object and the particles are
      // recorded. Now, we need to remove the invalid contacts from
      // the records and compute the forces associated with the valid
      // contacts.

      // Loop on every particle / contact record
      for (auto &[particle_index, this_contact_record] : contact_record)
        {
          // Loop on every contact in the contact record
          for (auto C1 = this_contact_record.begin();
               C1 != this_contact_record.end();)
            {
              // Extract the information of C1
              auto &[T1_cell,
                     normal_overlap_C1,
                     contact_indicator_C1,
                     contact_info_C1] = *C1;

              // Assigning the triangle neighboring list of T1;
              const auto &T1_es_neighbors = this_solid_es_neighbors.at(T1_cell);
              const auto &T1_vs_neighbors = this_solid_vs_neighbors.at(T1_cell);

              // Initialize variables. We need to compare C1 to every other
              // contact (C2) that was not compared yet.
              bool erase_contact_1 = false;
              auto C2              = std::next(C1);
              while (C2 != this_contact_record.end())
                {
                  // Extract the information of C2
                  auto &[T2_cell,
                         normal_overlap_C2,
                         contact_indicator_C2,
                         contact_info_C2] = *C2;

                  // First, we check if both triangle are neighbors. If they
                  // are not neighbors, C1 and C2 are automatically valid.
                  // (Disconnected triangles)
                  if (std::ranges::find(T1_es_neighbors, T2_cell) ==
                        T1_es_neighbors.end() &&
                      std::ranges::find(T1_vs_neighbors, T2_cell) ==
                        T1_vs_neighbors.end())
                    {
                      ++C2;
                      continue;
                    }
                  // If both triangles are neighbors (connected), we need to
                  // check for double contacts.
                  // If C1 is a face contact
                  if (contact_indicator_C1 ==
                      LetheGridTools::ParticleTriangleContactIndicator::
                        face_contact)
                    {
                      // If C1 and C2 are both face contacts, they are both
                      // valid. This case can happen with either edge-sharing or
                      // vertex-sharing triangles that are non-coplanar.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            face_contact)
                        {
                          ++C2;
                          continue;
                        }
                      // If C1 is a face contact and C2 is an edge contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            edge_contact)
                        {
                          // If T1 and T2 are vertex-sharing, C1 and C2 are
                          // valid. This could happen with two concave
                          // triangles.
                          if (std::ranges::find(T1_vs_neighbors, T2_cell) !=
                              T1_vs_neighbors.end())
                            {
                              ++C2;
                              continue;
                            }
                          // Otherwise, C2 is invalid and is erased from the
                          // current contact record.
                          clear_contact_info(contact_info_C2);
                          C2 = this_contact_record.erase(C2);
                          continue;
                        }

                      // If C1 is a face contact and C2 is a vertex
                      // contact.
                      // C2 is always invalid.
                      clear_contact_info(contact_info_C2);
                      C2 = this_contact_record.erase(C2);
                      continue;
                    }

                  // If C1 is an edge contact.
                  if (contact_indicator_C1 ==
                      LetheGridTools::ParticleTriangleContactIndicator::
                        edge_contact)
                    {
                      // If C1 is an edge contact and C2 is a face contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            face_contact)
                        {
                          // C2 is valid.
                          // C1 is invalid if T1 and T2 are connected. We
                          // already checked at the beginning of the while
                          // loop that this is the case, thus we know that
                          // both triangles are connected at this location
                          // in the code.
                          erase_contact_1 = true;
                          break; // must exit the while loop
                        }
                      // If C1 and C2 are both edge contacts.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            edge_contact)
                        {
                          // C1 is always valid.
                          // C2 is invalid if T1 and T2 are edge-sharing.
                          // This case means that both contacts occur on the
                          // edge that is connecting both triangles.
                          if (std::ranges::find(T1_es_neighbors, T2_cell) !=
                              T1_es_neighbors.end())
                            {
                              clear_contact_info(contact_info_C2);
                              C2 = this_contact_record.erase(C2);
                              continue;
                            }
                          else
                            {
                              // If they are vertex-sharing, this means that
                              // each contact is occurring on two different
                              // edges. This is not a double contact, C2 is
                              // also valid.
                              ++C2;
                              continue;
                            }
                        }
                    }
                  // If C1 is a vertex contact
                  if (contact_indicator_C1 ==
                      LetheGridTools::ParticleTriangleContactIndicator::
                        vertex_contact)
                    {
                      // If C1 is a vertex contact and C2 is a face contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            face_contact)
                        {
                          // C2 is valid.
                          // C1 is invalid if T1 and T2 are connected. We
                          // already checked at the beginning of the while
                          // loop that this is the case, thus we know if we
                          // enter this "if" that both triangle are
                          // connected.
                          erase_contact_1 = true;
                          break; // must exit the while loop
                        }
                      // If C1 is a vertex contact and C2 is an edge
                      // contact.
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            edge_contact)
                        {
                          // C2 is valid
                          // C1 is invalid if T1 and T2 are vertex-sharing.
                          // This case means that both contacts occur on the
                          // vertex that is connecting both triangles,
                          // thus C1 and C2 are double contacts.
                          if (std::ranges::find(T1_vs_neighbors, T2_cell) !=
                              T1_vs_neighbors.end())
                            {
                              erase_contact_1 = true;
                              break;
                            }
                          // Otherwise, C1 is valid.
                          ++C2;
                          continue;
                        }
                      // If C1 and C2 are vertex contacts,
                      if (contact_indicator_C2 ==
                          LetheGridTools::ParticleTriangleContactIndicator::
                            vertex_contact)
                        {
                          // C2 is invalid
                          clear_contact_info(contact_info_C2);
                          C2 = this_contact_record.erase(C2);
                          continue;
                        }
                    }
                  ++C2;
                } // While loop

              // Erase C1 from the contact record is marked as true.
              if (erase_contact_1)
                {
                  clear_contact_info(contact_info_C1);
                  C1 = this_contact_record.erase(C1);
                  continue; // skip increment; erase() already advanced
                }
              ++C1;
            } // For loop

          //  At this point, every contact in the contact record are valid,
          //  thus we can compute and applied the force and torques.

          // Get translational and rotational velocities and center of
          // rotation of solid object
          const Tensor<1, 3> translational_velocity =
            solids[solid_counter]->get_translational_velocity();
          const Tensor<1, 3> angular_velocity =
            solids[solid_counter]->get_angular_velocity();
          const Point<3> center_of_rotation =
            solids[solid_counter]->get_center_of_rotation();

          // Multiphysics properties
          const Parameters::ThermalBoundaryType thermal_boundary_type =
            solids[solid_counter]->get_thermal_boundary_type();

          for (auto contact = this_contact_record.begin();
               contact != this_contact_record.end();
               ++contact)
            {
              //  Extract the information of the contact
              auto &[T_cell, normal_overlap, contact_indicator, contact_info] =
                *contact;

              // Defining local variables which will be used within the
              // contact calculation
              auto        &particle             = contact_info.particle;
              auto         particle_properties  = particle->get_properties();
              Point<3>     particle_location_3d = get_location(particle);
              Tensor<1, 3> normal_force;
              Tensor<1, 3> tangential_force;
              Tensor<1, 3> tangential_torque;
              Tensor<1, 3> rolling_resistance_torque;
              double       normal_relative_velocity_value;
              Tensor<1, 3> tangential_relative_velocity;

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

              // Applying the calculated forces and torques on the particle
              this->apply_force_and_torque(
                normal_force,
                tangential_force,
                tangential_torque,
                rolling_resistance_torque,
                contact_outcome.torque[particle_index],
                contact_outcome.force[particle_index]);
              if constexpr (std::is_same_v<
                              PropertiesIndex,
                              DEM::DEMMPProperties::PropertiesIndex>)
                {
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

                      const double temperature_wall =
                        solids[solid_counter]->get_temperature();

                      // Apply the heat transfer to the particle
                      apply_heat_transfer_on_single_local_particle(
                        temperature_particle,
                        temperature_wall,
                        thermal_conductance,
                        particle_heat_transfer_rate);
                    }
                }
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
