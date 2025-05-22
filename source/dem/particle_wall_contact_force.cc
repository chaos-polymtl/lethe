// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_wall_contact_force.h>

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

  // Iterating over the solid objects
  for (unsigned int solid_counter = 0; solid_counter < solids.size();
       ++solid_counter)
    {
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
      double temperature_wall;
      if (thermal_boundary_type != Parameters::ThermalBoundaryType::adiabatic)
        temperature_wall = solids[solid_counter]->get_temperature();

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

              // Getting the projection of particles on the triangle
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
                      // Defining local variables which will be used within the
                      // contact calculation
                      auto &particle            = contact_info.particle;
                      auto  particle_properties = particle->get_properties();
                      Tensor<1, 3> normal_force;
                      Tensor<1, 3> tangential_force;
                      Tensor<1, 3> tangential_torque;
                      Tensor<1, 3> rolling_resistance_torque;
                      double       normal_relative_velocity_value;
                      Tensor<1, 3> tangential_relative_velocity;

                      const Point<3> &projection_point =
                        projection_points[particle_counter];

                      Point<3> particle_location_3d = get_location(particle);

                      const double particle_triangle_distance =
                        particle_location_3d.distance(projection_point);

                      // Find normal overlap
                      double normal_overlap =
                        ((particle_properties[PropertiesIndex::dp]) * 0.5) -
                        particle_triangle_distance;

                      if (normal_overlap > force_calculation_threshold_distance)
                        {
                          contact_info.normal_vector =
                            normal_vectors[particle_counter];

                          contact_info.point_on_boundary = projection_point;

                          contact_info.boundary_id = solid_counter;

                          // Updating contact information
                          this
                            ->update_particle_solid_object_contact_information(
                              contact_info,
                              tangential_relative_velocity,
                              normal_relative_velocity_value,
                              particle_properties,
                              dt,
                              translational_velocity,
                              angular_velocity,
                              center_of_rotation.distance(
                                particle_location_3d));

                          // Calculating contact force and torque
                          this->calculate_contact(
                            contact_info,
                            tangential_relative_velocity,
                            normal_relative_velocity_value,
                            normal_overlap,
                            dt,
                            particle_properties,
                            normal_force,
                            tangential_force,
                            tangential_torque,
                            rolling_resistance_torque);

                          // Applying the calculated forces and torques on the
                          // particle
                          types::particle_index particle_id =
                            particle->get_local_index();
                          Tensor<1, 3> &particle_torque =
                            contact_outcome.torque[particle_id];
                          Tensor<1, 3> &particle_force =
                            contact_outcome.force[particle_id];

                          this->apply_force_and_torque(
                            normal_force,
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

                      if constexpr (std::is_same_v<
                                      PropertiesIndex,
                                      DEM::DEMMPProperties::PropertiesIndex>)
                        {
                          if ((thermal_boundary_type !=
                               Parameters::ThermalBoundaryType::adiabatic) &&
                              (normal_overlap > 0))
                            {
                              const unsigned int particle_type =
                                particle_properties[PropertiesIndex::type];
                              const double temperature_particle =
                                particle_properties[PropertiesIndex::T];
                              double &particle_heat_transfer_rate =
                                contact_outcome.heat_transfer_rate
                                  [particle->get_local_index()];
                              double thermal_conductance;

                              calculate_contact_thermal_conductance<
                                ContactType::particle_floating_mesh>(
                                0.5 * particle_properties[PropertiesIndex::dp],
                                0,
                                this->effective_youngs_modulus[particle_type],
                                this->effective_real_youngs_modulus
                                  [particle_type],
                                this
                                  ->equivalent_surface_roughness[particle_type],
                                this->equivalent_surface_slope[particle_type],
                                this->effective_microhardness[particle_type],
                                this->particle_thermal_conductivity
                                  [particle_type],
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
                  particle_counter++;
                }
            }
        }
    }
}

// dem
// No resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;

// cfd_dem
// No resistance
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;

// dem_mp
//  No resistance
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::DMT,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::JKR,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex,
                                        ParticleWallContactForceModel::linear,
                                        RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleWallContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleWallContactForceModel::nonlinear,
  RollingResistanceMethod::epsd_resistance>;
