/* ---------------------------------------------------------------------
*
* Copyright (C) 2019 - 2019 by the Lethe authors
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
* Author: Lucka Barbeau, Shahab Golshan, Bruno Blais Polytechnique Montreal,
2021
*/

#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_linear_force.h>
#include <dem/particle_particle_nonlinear_force.h>
#include <dem/particle_wall_linear_force.h>
#include <dem/particle_wall_nonlinear_force.h>
#include <fem-dem/ib_particles_dem.h>

template <int dim>
void
IBParticlesDEM<dim>::initialize(
  const std::shared_ptr<Parameters::IBParticles<dim>> &p_nsparam,
  const MPI_Comm &                                     mpi_communicator_input,
  const std::vector<IBParticle<dim>> &                 particles)
{
  parameters       = p_nsparam;
  mpi_communicator = mpi_communicator_input;
  dem_particles    = particles;
  boundary_cells.resize(dem_particles.size());

  dem_parameters.model_parameters.particle_particle_contact_force_method =
    Parameters::Lagrangian::ModelParameters::ParticleParticleContactForceModel::
      hertz;
  dem_parameters.model_parameters.particle_wall_contact_force_method =
    Parameters::Lagrangian::ModelParameters::ParticleWallContactForceModel::
      nonlinear;
  particle_particle_contact_force_object =
    std::make_shared<ParticleParticleHertz<dim>>(dem_parameters);
  std::vector<types::boundary_id> boundary_index(0);
  particle_wall_contact_force_object =
    std::make_shared<ParticleWallNonLinearForce<dim>>(
      dem_parameters.boundary_conditions.boundary_translational_velocity,
      dem_parameters.boundary_conditions.boundary_rotational_speed,
      dem_parameters.boundary_conditions.boundary_rotational_vector,
      triangulation_cell_diameter,
      dem_parameters,
      boundary_index);
}
template <int dim>
void
IBParticlesDEM<dim>::update_particles(
  const std::vector<IBParticle<dim>> &particles,
  double                              time)
{
  dem_particles = particles;
  cfd_time      = time;
}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pp_contact_force(
  const double               dt_dem,
  std::vector<Tensor<1, 3>> &contact_force,
  std::vector<Tensor<1, 3>> &contact_torque)
{
  for (auto &particle_one : dem_particles)
    {
      for (auto &particle_two : dem_particles)
        {
          if (particle_one.particle_id != particle_two.particle_id and
              particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;
              particle_particle_contact_info_struct<dim> contact_info;

              // Check if there is already information on the contact of these
              // to particles. If not initialize it in the contact map with 0
              // values.
              try
                {
                  contact_info = pp_contact_map[particle_one.particle_id]
                                               [particle_two.particle_id];
                }
              catch (...)
                {
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d]           = 0;
                      contact_info.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id]
                                [particle_two.particle_id] = contact_info;
                }

              // Calculation of normal overlap
              double normal_overlap =
                (particle_one.radius + particle_two.radius) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0)
                // This means that the adjacent particles are in contact
                {
                  Tensor<1, 3> normal_force;
                  Tensor<1, 3> tangential_force;
                  Tensor<1, 3> particle_one_tangential_torque;
                  Tensor<1, 3> particle_two_tangential_torque;
                  Tensor<1, 3> rolling_resistance_torque;

                  particle_particle_contact_force_object
                    ->calculate_IB_particle_particle_contact_force(
                      normal_overlap,
                      pp_contact_map[particle_one.particle_id]
                                    [particle_two.particle_id],
                      normal_force,
                      tangential_force,
                      particle_one_tangential_torque,
                      particle_two_tangential_torque,
                      rolling_resistance_torque,
                      particle_one,
                      particle_two,
                      particle_one_location,
                      particle_two_location,
                      dt_dem,
                      particle_one.radius,
                      particle_two.radius,
                      particle_one.mass,
                      particle_two.mass);



                  contact_force[particle_one.particle_id] -=
                    (normal_force + tangential_force);
                  contact_force[particle_two.particle_id] +=
                    (normal_force + tangential_force);

                  contact_torque[particle_one.particle_id] =
                    contact_torque[particle_one.particle_id] -
                    particle_one_tangential_torque + rolling_resistance_torque;
                  contact_torque[particle_two.particle_id] =
                    contact_torque[particle_two.particle_id] -
                    particle_two_tangential_torque - rolling_resistance_torque;
                }

              else
                {
                  // if the adjacent pair is not in contact anymore
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d]           = 0;
                      contact_info.tangential_relative_velocity[d] = 0;
                    }
                  pp_contact_map[particle_one.particle_id].erase(
                    particle_two.particle_id);
                }
            }
        }
    }
}

template <int dim>
void
IBParticlesDEM<dim>::calculate_pp_lubrification_force(
  const double               dt_dem,
  const double h,
  std::vector<Tensor<1, 3>> &lubrification_force,
  std::vector<Tensor<1, 3>> &lubrification_torque)
{
  using numbers::PI;

  double mu=0.001;
  for (auto &particle_one : dem_particles)
    {
      for (auto &particle_two : dem_particles)
        {
          if (particle_one.particle_id != particle_two.particle_id and
              particle_one.particle_id < particle_two.particle_id)
            {
              const Point<dim> particle_one_location = particle_one.position;
              const Point<dim> particle_two_location = particle_two.position;
              Tensor<1, 3> raidal_vector;
              double radial_volecity;
              Tensor<1,3> f_lub;

              // Calculation of normal overlap
              double gap =particle_one_location.distance(particle_two_location)-(particle_one.radius + particle_two.radius);
              raidal_vector=tensor_nd_to_3d(particle_one.position-particle_two.position);
              if (gap < h)
                // This means that the adjacent particles are very close
                {
                  if(gap<0.0001*particle_one.radius)
                    {
                      gap = 0.0001 * particle_one.radius;
                    }

                  radial_volecity= scalar_product(-raidal_vector,particle_one.velocity)+scalar_product(raidal_vector,particle_two.velocity);
                  f_lub=3/2*PI*mu*(particle_one.radius*2*particle_two.radius*2/(particle_one.radius*2+particle_two.radius*2))*(particle_one.radius*2*particle_two.radius*2/(particle_one.radius*2+particle_two.radius*2))/gap*radial_volecity*raidal_vector/raidal_vector.norm();

                }
              Tensor<1,3> f_fluid_p1_parallel=(scalar_product(particle_one.fluid_forces,raidal_vector/raidal_vector.norm()))*raidal_vector/raidal_vector.norm();
              Tensor<1,3> f_fluid_p2_parallel=(scalar_product(particle_two.fluid_forces,raidal_vector/raidal_vector.norm()))*raidal_vector/raidal_vector.norm();
              Tensor<1,3> f_fluid_p1_orto= particle_one.fluid_forces- f_fluid_p1_parallel;
              Tensor<1,3> f_fluid_p2_orto= particle_two.fluid_forces- f_fluid_p2_parallel;

              Tensor<1,3> f_lub_max_p1;
              Tensor<1,3> f_lub_max_p2;

              if(f_fluid_p1_parallel.norm()>f_lub.norm()){
                  f_lub_max_p1=0;
              }
              else{
                  f_lub_max_p1=f_lub-f_lub_max_p1;
                }
              if(f_fluid_p2_parallel.norm()>f_lub.norm()){
                  f_lub_max_p2=0;
                }
              else{
                  f_lub_max_p2=-f_lub-f_lub_max_p1;
                }

              lubrification_force[particle_one.particle_id]=f_lub_max_p1;
              lubrification_force[particle_two.particle_id]=f_lub_max_p2;
            }
        }
    }
}

template <int dim>
void
IBParticlesDEM<dim>::update_particles_boundary_contact(
  const std::vector<IBParticle<dim>> &particles,
  const DoFHandler<dim> &             dof_handler,
  const Quadrature<dim - 1> &         face_quadrature_formula,
  const Mapping<dim> &                mapping)
{
  const FESystem<dim, dim> fe = dof_handler.get_fe();
  for (unsigned int p_i = 0; p_i < particles.size(); ++p_i)
    {
      // Clear the last boundary cell candidates.
      boundary_cells[p_i].clear();

      // Find the new cells that are at a boundary and in proximity of the
      // particle.
      auto cells_at_boundary = LetheGridTools::find_boundary_cells_in_sphere(
        dof_handler, particles[p_i].position, particles[p_i].radius * 1.5);

      // Loop over the cells at the boundary.
      for (unsigned int i = 0; i < cells_at_boundary.size(); ++i)
        {
          unsigned int      n_face_q_points = face_quadrature_formula.size();
          FEFaceValues<dim> fe_face_values(mapping,
                                           fe,
                                           face_quadrature_formula,
                                           update_values |
                                             update_quadrature_points |
                                             update_normal_vectors);
          // Loop over the faces of the cell at the boundary.
          for (int face_id = 0;
               face_id < int(GeometryInfo<dim>::faces_per_cell);
               ++face_id)
            {
              // Find the face at the boundary
              if (cells_at_boundary[i]->face(face_id)->at_boundary())
                {
                  fe_face_values.reinit(cells_at_boundary[i], face_id);
                  // Loop over the quadrature point of the face at the boundary
                  // to store information about the location and normals of the
                  // boundary.
                  for (unsigned int f_q_point = 0; f_q_point < n_face_q_points;
                       ++f_q_point)
                    {
                      Tensor<1, dim> normal_vector =
                        -fe_face_values.normal_vector(f_q_point);
                      BoundaryCellsInfo boundary_information;
                      boundary_information.normal_vector = normal_vector;
                      boundary_information.point_on_boundary =
                        fe_face_values.quadrature_point(f_q_point);
                      boundary_cells[p_i].push_back(boundary_information);
                    }
                }
            }
        }

      // Regroup the information of all the processor
      auto global_boundary_cell =
        Utilities::MPI::all_gather(this->mpi_communicator, boundary_cells[p_i]);
      boundary_cells[p_i].clear();
      for (unsigned int i = 0; i < global_boundary_cell.size(); ++i)
        {
          boundary_cells[p_i].insert(boundary_cells[p_i].end(),
                                     global_boundary_cell[i].begin(),
                                     global_boundary_cell[i].end());
        }
    }
}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pw_contact_force(
  const double               dt_dem,
  std::vector<Tensor<1, 3>> &contact_force,
  std::vector<Tensor<1, 3>> &contact_torque)
{
  double wall_youngs_modulus = parameters->wall_youngs_modulus;
  double wall_poisson_ratio  = parameters->wall_poisson_ratio;
  double wall_rolling_friction_coefficient =
    parameters->wall_rolling_friction_coefficient;
  double wall_friction_coefficient = parameters->wall_friction_coefficient;
  double wall_restitution_coefficient =
    parameters->wall_restitution_coefficient;

  // Loop over the particles
  for (auto &particle : dem_particles)
    {
      unsigned int boundary_index = 0;
      double       best_dist      = DBL_MAX;
      // BB FIX
      unsigned int best_index = UINT_MAX;
      // For each particle loop over the point and normal identified as
      // potential contact candidate.
      for (auto &boundary_cell_iter : boundary_cells[particle.particle_id])
        {
          // find the best candidate (the closest point).
          double dist =
            (boundary_cell_iter.point_on_boundary - particle.position).norm();
          if (dist < best_dist)
            {
              best_dist  = dist;
              best_index = boundary_index;
            }
          boundary_index += 1;
        }
      // Do the particle wall contact calculation with the best candidate.
      if (boundary_cells[particle.particle_id].size() > 0)
        {
          auto &boundary_cell =
            boundary_cells[particle.particle_id][best_index];

          auto boundary_cell_information = boundary_cell;
          particle_wall_contact_info_struct<dim> contact_info;

          // Check if there is already information on the contact between this
          // particle and this boundary contact point. If not initialize the
          // contact history with 0 values.
          try
            {
              contact_info =
                pw_contact_map[particle.particle_id][boundary_index];
            }
          catch (...)
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_info.tangential_overlap[d]           = 0;
                  contact_info.tangential_relative_velocity[d] = 0;
                }

              // BB temporary fix for unused variables
              contact_info.global_face_id           = 0;
              contact_info.boundary_id              = 0;
              contact_info.normal_overlap           = 0;
              contact_info.normal_relative_velocity = 0;
              // End BB

              pw_contact_map[particle.particle_id][boundary_index] =
                contact_info;
            }

          Tensor<1, 3> normal =
            tensor_nd_to_3d(boundary_cell_information.normal_vector);
          auto point_on_boundary = boundary_cell_information.point_on_boundary;


          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance

          Point<3> particle_position_3d = point_nd_to_3d(particle.position);
          Point<3> point_on_boundary_3d = point_nd_to_3d(point_on_boundary);

          Tensor<1, 3> point_to_particle_vector =
            particle_position_3d - point_on_boundary_3d;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, 3> projected_vector =
            particle_wall_contact_force_object->find_projection(
              point_to_particle_vector, normal);

          // Find the normal overlap
          double normal_overlap = particle.radius - (projected_vector.norm());

          if (normal_overlap > 0)
            {
              // Do the calculation to evaluate the particle wall contact force.
              contact_info.normal_overlap = normal_overlap;

              Tensor<1, 3> normal_force;
              Tensor<1, 3> tangential_force;
              Tensor<1, 3> tangential_torque;
              Tensor<1, 3> rolling_resistance_torque;

              pw_contact_map[particle.particle_id][boundary_index]
                .normal_overlap = normal_overlap;
              pw_contact_map[particle.particle_id][boundary_index]
                .normal_vector = normal;
              pw_contact_map[particle.particle_id][boundary_index].boundary_id =
                0;

              particle_wall_contact_force_object
                ->calculate_IB_particle_wall_contact_force(
                  pw_contact_map[particle.particle_id][boundary_index],
                  normal_force,
                  tangential_force,
                  tangential_torque,
                  rolling_resistance_torque,
                  particle,
                  wall_youngs_modulus,
                  wall_poisson_ratio,
                  wall_restitution_coefficient,
                  wall_friction_coefficient,
                  wall_rolling_friction_coefficient,
                  dt_dem,
                  particle.mass,
                  particle.radius);
              // Updating the force of particles in the particle handler
              contact_force[particle.particle_id] +=
                normal_force + tangential_force;
              // Updating the torque acting on particles
              contact_torque[particle.particle_id] +=
                tangential_torque + rolling_resistance_torque;
            }
          else
            {
              // Set to 0 the tangential overlap if the particle is not in
              // contact with the wall anymore.
              for (int d = 0; d < dim; ++d)
                {
                  contact_info.tangential_overlap[d] = 0;
                }
            }
        }
    }
}


template <int dim>
void
IBParticlesDEM<dim>::calculate_pw_lubrification_force(
  const double               dt_dem,
  const double h,
  std::vector<Tensor<1, 3>> &lubrification_force,
  std::vector<Tensor<1, 3>> &lubrification_torque)
{
  using numbers::PI;

  double mu=0.001;
  // Loop over the particles
  for (auto &particle : dem_particles)
    {
      unsigned int boundary_index = 0;
      double       best_dist      = DBL_MAX;
      unsigned int best_index;
      // For each particle loop over the point and normal identified as
      // potential contact candidate.
      for (auto &boundary_cell_iter : boundary_cells[particle.particle_id])
        {
          // find the best candidate (the closest point).
          double dist =
            (boundary_cell_iter.point_on_boundary - particle.position).norm();
          if (dist < best_dist)
            {
              best_dist  = dist;
              best_index = boundary_index;
            }
          boundary_index += 1;
        }
      // Do the particle wall contact calculation with the best candidate.
      if (boundary_cells[particle.particle_id].size() > 0)
        {
          auto &boundary_cell =
            boundary_cells[particle.particle_id][best_index];

          auto boundary_cell_information = boundary_cell;
          Tensor<1, 3> normal =
            tensor_nd_to_3d(boundary_cell_information.normal_vector);
          auto point_on_boundary = boundary_cell_information.point_on_boundary;
          const Point<dim> particle_one_location = particle.position;

          // Calculation of normal overlap
          Point<3> particle_position_3d = point_nd_to_3d(particle.position);
          Point<3> point_on_boundary_3d = point_nd_to_3d(point_on_boundary);

          Tensor<1, 3> point_to_particle_vector =
            particle_position_3d - point_on_boundary_3d;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, 3> projected_vector =
            particle_wall_contact_force_object->find_projection(
              point_to_particle_vector, normal);

          // Find the normal overlap
          double gap = (projected_vector.norm())-particle.radius;

          Tensor<1, 3>     raidal_vector;
          double           radial_volecity;
          Tensor<1, 3>     f_lub;


          raidal_vector =projected_vector;
          if (gap < h)
            // This means that the adjacent particles are very close
            {
              if (gap < 0.0001 * particle.radius)
                {
                  gap = 0.0001 * particle.radius;
                }
              radial_volecity =
                scalar_product(-raidal_vector, particle.velocity) ;
              f_lub = 3 / 2 * PI * mu *
                      (particle.radius) *
                      (particle.radius) /
                      gap * radial_volecity * raidal_vector /
                      raidal_vector.norm();
              //std::cout<<"f_lub "<<f_lub <<std::endl;
            }
          Tensor<1, 3> f_fluid_p1_parallel =
            (scalar_product(particle.fluid_forces,
                            raidal_vector / raidal_vector.norm())) *
            raidal_vector / raidal_vector.norm();

          Tensor<1, 3> f_fluid_p1_orto =
            particle.fluid_forces - f_fluid_p1_parallel;

          Tensor<1, 3> f_lub_max_p1;

          if (f_fluid_p1_parallel.norm() > f_lub.norm())
            {
              f_lub_max_p1 = 0;
            }
          else
            {
              f_lub_max_p1 = f_lub - f_lub_max_p1;
            }

          lubrification_force[particle.particle_id] = f_lub_max_p1;

        }
    }

}

template <int dim>
void
IBParticlesDEM<dim>::integrate_particles_motion(const double dt,const double h)
{
  // Initialize local containers and physical variables
  using numbers::PI;
  double rho    = parameters->density;
  double dt_dem = dt / parameters->coupling_frequency;

  std::vector<Tensor<1, 3>> contact_force(dem_particles.size());
  std::vector<Tensor<1, 3>> contact_wall_force(dem_particles.size());
  std::vector<Tensor<1, 3>> contact_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> contact_wall_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> current_fluid_force(dem_particles.size());
  std::vector<Tensor<1, 3>> current_fluid_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrification_force(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrification_torque(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrification_wall_force(dem_particles.size());
  std::vector<Tensor<1, 3>> lubrification_wall_torque(dem_particles.size());

  std::vector<Tensor<1, 3>> velocity(dem_particles.size());
  std::vector<Point<dim>>   position(dem_particles.size());

  // Local time for the dem step.
  double t = 0;
  // The gravitational acceleration.
  Tensor<1, 3> g;
  this->parameters->f_gravity->set_time(cfd_time);
  // The gravitational force on the particle.
  Tensor<1, 3> gravity;
  lubrification_force.clear();
  lubrification_force.resize(dem_particles.size());
  lubrification_torque.clear();
  lubrification_torque.resize(dem_particles.size());
  lubrification_wall_force.clear();
  lubrification_wall_force.resize(dem_particles.size());
  lubrification_wall_torque.clear();
  lubrification_wall_torque.resize(dem_particles.size());
  // Initialize the particles
  for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
    {
      dem_particles[p_i].position  = dem_particles[p_i].previous_positions[0];
      dem_particles[p_i].velocity  = dem_particles[p_i].previous_velocity[0];
      dem_particles[p_i].omega     = dem_particles[p_i].previous_omega[0];
      dem_particles[p_i].impulsion = 0;
      dem_particles[p_i].omega_impulsion         = 0;
      dem_particles[p_i].contact_impulsion       = 0;
      dem_particles[p_i].omega_contact_impulsion = 0;
      // Initialized the gravity at the particle position.
      g[0] = this->parameters->f_gravity->value(dem_particles[p_i].position, 0);
      g[1] = this->parameters->f_gravity->value(dem_particles[p_i].position, 1);
      if (dim == 3)
        g[2] =
          this->parameters->f_gravity->value(dem_particles[p_i].position, 2);
    }

  // Integrate with the sub_time_step
  while (t + 0.5 * dt_dem < dt)
    {
      current_fluid_force.clear();
      current_fluid_force.resize(dem_particles.size());
      current_fluid_torque.clear();
      current_fluid_torque.resize(dem_particles.size());
      contact_torque.clear();
      contact_torque.resize(dem_particles.size());
      contact_force.clear();
      contact_force.resize(dem_particles.size());
      contact_wall_force.clear();
      contact_wall_force.resize(dem_particles.size());
      contact_wall_torque.clear();
      contact_wall_torque.resize(dem_particles.size());


      // Calculate particle-particle and particle-wall contact force
      calculate_pp_contact_force(dt_dem, contact_force, contact_torque);
      calculate_pw_contact_force(dt_dem,
                                 contact_wall_force,
                                 contact_wall_torque);
      calculate_pp_lubrification_force(dt_dem,h,lubrification_force,lubrification_torque);
      calculate_pw_lubrification_force(dt_dem,h,lubrification_wall_force,lubrification_wall_torque);


      for (unsigned int p_i = 0; p_i < dem_particles.size(); ++p_i)
        {
          auto inv_inertia = invert(dem_particles[p_i].inertia);
          if (dim == 2)
            gravity = g * (dem_particles[p_i].mass -
                           dem_particles[p_i].radius *
                             dem_particles[p_i].radius * PI * rho);
          if (dim == 3)
            {
              gravity = g * (dem_particles[p_i].mass -
                             4.0 / 3.0 * dem_particles[p_i].radius *
                               dem_particles[p_i].radius *
                               dem_particles[p_i].radius * PI * rho);
            }

          // We consider only the force at t+dt so the scheme is consistent to a
          // BDFn scheme on the fluid side. If there is no contact.
          current_fluid_force[p_i]  = dem_particles[p_i].fluid_forces;
          current_fluid_torque[p_i] = dem_particles[p_i].fluid_torque;

          // Explicit Euler for the sub_time_stepping. This is exact for the
          // gravity and fluid impulsion integration. In case of contact the
          // scheme becomes first order. This could be improved in the future.

          dem_particles[p_i].velocity =
            dem_particles[p_i].velocity +
            (current_fluid_force[p_i] + contact_force[p_i] +
             contact_wall_force[p_i] + gravity+lubrification_force[p_i]+lubrification_wall_force[p_i]) *
              dt_dem / dem_particles[p_i].mass;

          for (unsigned int d = 0; d < dim; ++d)
            {
              dem_particles[p_i].position[d] =
                dem_particles[p_i].position[d] +
                dem_particles[p_i].velocity[d] * dt_dem;
            }

          dem_particles[p_i].omega =
            dem_particles[p_i].omega +
            inv_inertia *
              (current_fluid_torque[p_i] + contact_torque[p_i] +
               contact_wall_torque[p_i]) *
              dt_dem;

          // Integration of the impulsion applied to the particle.
          // This is what will be transferred to the CFD to integrate the
          // particle.
          dem_particles[p_i].impulsion +=
            (current_fluid_force[p_i] + gravity + contact_wall_force[p_i] +
             contact_force[p_i]+lubrification_force[p_i]+lubrification_wall_force[p_i]) *
            dt_dem;
          dem_particles[p_i].contact_impulsion +=
            (contact_wall_force[p_i] + contact_force[p_i]) * dt_dem;
          dem_particles[p_i].omega_impulsion +=
            (current_fluid_torque[p_i] + contact_torque[p_i] +
             contact_wall_torque[p_i]) *
            dt_dem;
          dem_particles[p_i].omega_contact_impulsion +=
            (contact_torque[p_i] + contact_wall_torque[p_i]) * dt_dem;
        }
      t += dt_dem;
    }
}


template class IBParticlesDEM<2>;
template class IBParticlesDEM<3>;
