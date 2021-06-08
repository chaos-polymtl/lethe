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

#include <dem/dem_solver_parameters.h>
#include <dem/pw_contact_force.h>

// Updates the contact information (contact_info) based on the new information
// of particles pair in the current time step
template <int dim>
void
PWContactForce<dim>::update_contact_information(
  pw_contact_info_struct<dim> &  contact_info,
  const ArrayView<const double> &particle_properties,
  const double &                 dt)
{
  auto               normal_vector = contact_info.normal_vector;
  const unsigned int boundary_id   = contact_info.boundary_id;

  // Using velocity and angular velocity of particle as
  // local vectors
  Tensor<1, dim> particle_velocity;
  particle_velocity[0] = particle_properties[DEM::PropertiesIndex::v_x];
  particle_velocity[1] = particle_properties[DEM::PropertiesIndex::v_y];
  if (dim == 3)
    {
      particle_velocity[2] = particle_properties[DEM::PropertiesIndex::v_z];
    }

  Tensor<1, dim> particle_omega;
  particle_omega[0] = particle_properties[DEM::PropertiesIndex::omega_x];
  particle_omega[1] = particle_properties[DEM::PropertiesIndex::omega_y];
  if (dim == 3)
    {
      particle_omega[2] = particle_properties[DEM::PropertiesIndex::omega_z];
    }

  // Defining relative contact velocity
  Tensor<1, dim> contact_relative_velocity;
  if (dim == 3)
    {
      contact_relative_velocity =
        particle_velocity -
        this->boundary_translational_velocity_map[boundary_id] +
        cross_product_3d((0.5 * particle_properties[DEM::PropertiesIndex::dp] *
                            particle_omega +
                          this->triangulation_radius *
                            this->boundary_rotational_speed_map[boundary_id] *
                            this->boundary_rotational_vector[boundary_id]),
                         normal_vector);
    }
  if (dim == 2)
    {
      contact_relative_velocity =
        particle_velocity - this->triangulation_radius *
                              this->boundary_rotational_speed_map[boundary_id] *
                              cross_product_2d(normal_vector);
    }

  // Calculation of normal relative velocity
  double normal_relative_velocity_value =
    contact_relative_velocity * normal_vector;
  Tensor<1, dim> normal_relative_velocity =
    normal_relative_velocity_value * normal_vector;

  // Calculation of tangential relative velocity
  Tensor<1, dim> tangential_relative_velocity =
    contact_relative_velocity - normal_relative_velocity;

  // Calculation of new tangential_overlap, since this value is
  // history-dependent it needs the value at previous time-step
  // This variable is the main reason that we have iteration over
  // two different vectors (pairs_in_contact and
  // contact_pair_candidates): tangential_overlap of the particles
  // which were already in contact (pairs_in_contact) needs to be
  // modified using its history, while the tangential_overlaps of
  // new particles are equal to zero
  Tensor<1, dim> modified_tangential_overlap =
    contact_info.tangential_overlap + tangential_relative_velocity * dt;

  // Updating the contact_info container based on the new calculated values
  contact_info.normal_relative_velocity     = normal_relative_velocity_value;
  contact_info.tangential_overlap           = modified_tangential_overlap;
  contact_info.tangential_relative_velocity = tangential_relative_velocity;
}

template <int dim>
Tensor<1, dim>
  PWContactForce<dim>::calculation_total_torque(Tensor<1, dim> total_force,
                                                Point<dim>     center_mass,
                                                Point<dim> point_on_boundary)
{
  if (dim == 3)
    {
      return cross_product_3d(total_force, point_on_boundary - center_mass);
    }
}



template <int dim>
void
PWContactForce<dim>::calculate_pw_force_torque(
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
    &                      pw_pairs_in_contact,
  const double &           dt,
  DEMSolverParameters<dim> parameters)
{
  std::map<unsigned int, Tensor<1, dim>> total_force;
  std::map<unsigned int, Tensor<1, dim>> total_torque;
  // Tensor<1, dim> total_force; total_force.clear();
  // Looping over pw_pairs_in_contact, which means looping over all the active
  // particles with iterator pw_pairs_in_contact_iterator
  for (auto &&pairs_in_contact_content :
       pw_pairs_in_contact | boost::adaptors::map_values)
    {
      // Now an iterator (pw_contact_information_iterator) on each element of
      // the pw_pairs_in_contact vector is defined. This iterator iterates over
      // a map which contains the required information for calculation of the
      // contact force for each particle
      for (auto &&contact_information :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          // Defining total force of the contact, properties of particle as
          // local parameters
          auto particle            = contact_information.particle;
          auto particle_properties = particle->get_properties();

          auto normal_vector     = contact_information.normal_vector;
          auto point_on_boundary = contact_information.point_on_boundary;

          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance
          Tensor<1, dim> point_to_particle_vector =
            particle->get_location() - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, dim> projected_vector =
            this->find_projection(point_to_particle_vector, normal_vector);
          double normal_overlap =
            ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
            (projected_vector.norm());

          if (normal_overlap > 0)
            {
              contact_information.normal_overlap = normal_overlap;

              this->update_contact_information(contact_information,
                                               particle_properties,
                                               dt);
              std::tuple<Tensor<1, dim>, Tensor<1, dim>> forces_and_torques;

              // This tuple (forces and torques) contains four elements which
              // are: 1, normal force, 2, tangential force, 3, tangential torque
              // and 4, rolling resistance torque, respectively
              if (parameters.model_parameters.pw_contact_force_method ==
                  Parameters::Lagrangian::ModelParameters::PWContactForceModel::
                    pw_linear)
                {
                  const unsigned int particle_type =
                    particle_properties[DEM::PropertiesIndex::type];

                  // Calculation of normal and tangential spring and dashpot
                  // constants using particle properties
                  double rp_sqrt =
                    sqrt(particle_properties[DEM::PropertiesIndex::dp] * 0.5);

                  double normal_spring_constant =
                    1.0667 * rp_sqrt *
                    this->effective_youngs_modulus[particle_type] *
                    pow((0.9375 *
                         particle_properties[DEM::PropertiesIndex::mass] *
                         contact_information.normal_relative_velocity *
                         contact_information.normal_relative_velocity /
                         (rp_sqrt *
                          this->effective_youngs_modulus[particle_type])),
                        0.2);
                  double tangential_spring_constant =
                    -1.0667 * rp_sqrt *
                      this->effective_youngs_modulus[particle_type] *
                      pow((0.9375 *
                           particle_properties[DEM::PropertiesIndex::mass] *
                           contact_information.tangential_relative_velocity *
                           contact_information.tangential_relative_velocity /
                           (rp_sqrt *
                            this->effective_youngs_modulus[particle_type])),
                          0.2) +
                    DBL_MIN;
                  double normal_damping_constant = sqrt(
                    (4 * particle_properties[DEM::PropertiesIndex::mass] *
                     normal_spring_constant) /
                    (1 +
                     pow((M_PI / (log(this->effective_coefficient_of_restitution
                                        [particle_type]) +
                                  DBL_MIN)),
                         2)));

                  // Calculation of normal force using spring and dashpot normal
                  // forces
                  Tensor<1, dim> normal_force =
                    (normal_spring_constant *
                       contact_information.normal_overlap -
                     normal_damping_constant *
                       contact_information.normal_relative_velocity) *
                    contact_information.normal_vector;
                  ;

                  // Calculation of tangential force
                  Tensor<1, dim> tangential_force =
                    tangential_spring_constant *
                    contact_information.tangential_overlap;

                  double coulomb_threshold =
                    this->effective_coefficient_of_friction[particle_type] *
                    normal_force.norm();

                  // Check for gross sliding
                  if (tangential_force.norm() > coulomb_threshold)
                    {
                      // Gross sliding occurs and the tangential overlap and
                      // tangential force are limited to Coulumb's criterion
                      tangential_force =
                        coulomb_threshold *
                        (tangential_force / tangential_force.norm());

                      contact_information.tangential_overlap =
                        tangential_force /
                        (tangential_spring_constant + DBL_MIN);
                    }


                  forces_and_torques =
                    std::make_tuple(normal_force, tangential_force);
                }
              else
                {
                  const unsigned int particle_type =
                    particle_properties[DEM::PropertiesIndex::type];

                  // Calculation of model parameters (beta, sn and st). These
                  // values are used to consider non-linear relation of the
                  // contact force to the normal overlap
                  double radius_times_overlap_sqrt =
                    sqrt(particle_properties[DEM::PropertiesIndex::dp] * 0.5 *
                         contact_information.normal_overlap);
                  double log_coeff_restitution = log(
                    this->effective_coefficient_of_restitution[particle_type]);
                  double model_parameter_beta =
                    log_coeff_restitution /
                    sqrt((log_coeff_restitution * log_coeff_restitution) +
                         9.8696);
                  double model_parameter_sn =
                    2 * this->effective_youngs_modulus[particle_type] *
                    radius_times_overlap_sqrt;

                  // Calculation of normal and tangential spring and dashpot
                  // constants using particle and wall properties
                  double normal_spring_constant =
                    1.3333 * this->effective_youngs_modulus[particle_type] *
                    radius_times_overlap_sqrt;
                  double normal_damping_constant =
                    1.8257 * model_parameter_beta *
                    sqrt(model_parameter_sn *
                         particle_properties[DEM::PropertiesIndex::mass]);
                  double tangential_spring_constant =
                    -8 * this->effective_shear_modulus[particle_type] *
                      radius_times_overlap_sqrt +
                    DBL_MIN;

                  // Calculation of normal force using spring and dashpot normal
                  // forces
                  Tensor<1, dim> normal_force =
                    (normal_spring_constant *
                       contact_information.normal_overlap +
                     normal_damping_constant *
                       contact_information.normal_relative_velocity) *
                    contact_information.normal_vector;

                  // Calculation of tangential force
                  Tensor<1, dim> tangential_force =
                    tangential_spring_constant *
                    contact_information.tangential_overlap;

                  double coulomb_threshold =
                    this->effective_coefficient_of_friction[particle_type] *
                    normal_force.norm();

                  // Check for gross sliding
                  if (tangential_force.norm() > coulomb_threshold)
                    {
                      // Gross sliding occurs and the tangential overlap and
                      // tangential force are limited to Coulumb's criterion
                      tangential_force =
                        coulomb_threshold *
                        (tangential_force / tangential_force.norm());

                      contact_information.tangential_overlap =
                        tangential_force /
                        (tangential_spring_constant + DBL_MIN);
                    }


                  forces_and_torques =
                    std::make_tuple(normal_force, tangential_force);
                }

              // Apply the calculated forces and torques on the particle pair

              Tensor<1, dim> normal_force     = std::get<0>(forces_and_torques);
              Tensor<1, dim> tangential_force = std::get<1>(forces_and_torques);

              // Calculation of total force
              total_force[contact_information.boundary_id] =
                total_force[contact_information.boundary_id] + normal_force +
                tangential_force;
              total_torque[contact_information.boundary_id] =
                total_torque[contact_information.boundary_id] +
                calculation_total_torque(
                  total_force[contact_information.boundary_id],
                  parameters.forces_torques.point_center_mass,
                  point_on_boundary);
              std::cout << "\n" << contact_information.boundary_id;
            }
        }
    }
  std::cout << "\n"
            << -total_force[2][0] << " " << -total_force[2][1] << " "
            << -total_force[2][2];
  std::cout << "\n"
            << -total_force[3][0] << " " << -total_force[3][1] << " "
            << -total_force[3][2];
  std::cout << "\n"
            << -total_force[4][0] << " " << -total_force[4][1] << " "
            << -total_force[4][2];
  std::cout << "\n"
            << -total_torque[2][0] << " " << -total_torque[2][1] << " "
            << -total_torque[2][2];
  std::cout << "\n"
            << -total_torque[3][0] << " " << -total_torque[3][1] << " "
            << -total_torque[3][2];
  std::cout << "\n"
            << -total_torque[4][0] << " " << -total_torque[4][1] << " "
            << -total_torque[4][2];
}


template class PWContactForce<2>;
template class PWContactForce<3>;
