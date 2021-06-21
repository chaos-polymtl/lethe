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
void
PWContactForce<dim>::get_force_torque()
{
  for (const auto &it : ForceOnWall)
    {
      std::cout << "Boundary " << it.first << " :\n"
                << "F = " << it.second << "\nM = " << TorqueOnWall[it.first]
                << "\n\n";
    }
}

template <int dim>
void
PWContactForce<dim>::update_boundary_velocity(
  const DEMSolverParameters<dim> &dem_parameters)
{
  if (dem_parameters.boundary_motion.motion_method ==
      Parameters::Lagrangian::BoundaryMotion<dim>::MotionMethod::free)
    {
      // Definition of all the parameters needed
      double      dt  = dem_parameters.simulation_control.dt;
      double      m   = dem_parameters.forces_torques.mass_solid;
      std::string rad = dem_parameters.mesh.grid_arguments;
      rad.erase(remove(rad.begin(), rad.end(), ' '), rad.end());
      rad           = rad.substr(rad.find(":") + 1,
                       rad.find(":", rad.find(":") + 1) - rad.find(":") - 1);
      double radius = std::stod(rad);
      double mu_s   = dem_parameters.forces_torques.friction_coefficient_mu_s;
      double mu_k   = dem_parameters.forces_torques.friction_coefficient_mu_k;
      double mu_r   = dem_parameters.forces_torques.friction_coefficient_mu_r;
      double theta  = dem_parameters.forces_torques.angle_theta;
      Tensor<1, dim> g(dem_parameters.physical_properties.g), P(m * g);
      Tensor<1, dim> I;
      I[0] = dem_parameters.forces_torques.moment_inertia_x;
      I[1] = dem_parameters.forces_torques.moment_inertia_y;
      if (dim > 2)
        {
          I[2] = dem_parameters.forces_torques.moment_inertia_z;
        }
      Tensor<1, dim> Sum_F, Sum_M, F_friction, M_friction, N;
      for (const auto &it : ForceOnWall)
        {
          Sum_F = Sum_F + it.second;
          Sum_M = Sum_M + TorqueOnWall[it.first];
        }
      // base_change is used to change original coordinate vectors into a
      // simplified model
      Sum_F = base_change(Sum_F, theta);
      // std::cout << "\nSum_F " << Sum_F << "\n";
      // std::cout << "\nSum_M " << Sum_M << "\n";
      // Second's law of Newton says that the normal reaction is equal to all
      // forces in the other direction
      N[2] = -(Sum_F[2] + P[2]);
      // std::cout << "\nN " << N << "\n";
      bool is_moving = (abs(Sum_F[1] + P[1]) > abs(mu_s * N[2]));
      // std::cout << "\nis_moving " << is_moving << "\n";
      std::string rotation_direction;
      if (is_moving)
        {
          rotation_direction =
            (Sum_F[1] + P[1] > abs(mu_s * N[2])) ? "left" : "right";
        }
      if (is_moving && rotation_direction == "right")
        {
          F_friction[1] = -mu_k * N[2];
        }
      else if (is_moving && rotation_direction == "left")
        {
          F_friction[1] = mu_k * N[2];
        }
      else if (!is_moving)
        {
          F_friction[1] = -(P[1] + Sum_F[1]);
        }
      // std::cout << "\nF_friction " << F_friction << "\n";
      Tensor<1, dim> vector_contact;
      vector_contact[2] = -radius;
      M_friction        = cross_product_3d(vector_contact, F_friction);
      // std::cout << "\nM_friction " << M_friction << "\n";
      Tensor<1, dim> add_F, add_M, transition;
      add_F = (dt / m) * (Sum_F + F_friction + P + N);
      // std::cout << "\nadd_F " << add_F << "\n";
      Tensor<1, dim> vector_rolling_contact;
      vector_rolling_contact[1] = mu_r * radius;
      vector_rolling_contact[2] = -radius;
      Tensor<1, dim> M_rolling(cross_product_3d(vector_rolling_contact, N));
      for (int i = 0; i < dim; i++)
        {
          add_M[i] = (dt / I[i]) * (Sum_M[i] + M_friction[i] + M_rolling[i]);
        }
      // std::cout << "\nadd_M " << add_M << "\n";
      for (const auto &it : boundary_translational_velocity_map)
        {
          boundary_translational_velocity_map[it.first] = it.second + add_F;
          // std::cout << "\nboundary_translational_velocity_map[it.first] " <<
          // boundary_translational_velocity_map[it.first] << "\n";
          transition = boundary_rotational_speed_map[it.first] *
                         boundary_rotational_vector[it.first] +
                       add_M;
          boundary_rotational_speed_map[it.first] = transition.norm();
          // std::cout << "\nboundary_rotational_speed_map[it.first] " <<
          // boundary_rotational_speed_map[it.first] << "\n";
          boundary_rotational_vector[it.first] = transition / transition.norm();
          // std::cout << "\nboundary_rotational_vector[it.first] " <<
          // boundary_rotational_vector[it.first] << "\n";
        }
    }
}

template <int dim>
Tensor<1, dim>
PWContactForce<dim>::base_change(const Tensor<1, dim> tensor,
                                 const double         theta)
{
  Tensor<1, dim> out(tensor);
  out[1] =
    tensor[1] * cos(theta * M_PI / 180) - tensor[2] * sin(theta * M_PI / 180);
  out[2] =
    tensor[2] * cos(theta * M_PI / 180) + tensor[1] * sin(theta * M_PI / 180);
  return out;
}


template class PWContactForce<2>;
template class PWContactForce<3>;
