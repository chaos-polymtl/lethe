// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#ifndef full_contact_functions_h
#define full_contact_functions_h

#include <../tests/dem/test_particles_functions.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/velocity_verlet_integrator.h>

#include <fstream>
#include <iomanip>
#include <string>

/**
 * @brief Store the initial properties of both particles.
 */
template <int dim>
struct initial_particle_properties
{
  Point<dim>     position[2];
  int            id[2];
  int            type[2];
  Tensor<1, dim> v[2];
  Tensor<1, dim> omega[2];
  double         mass[2];
  double         diameter[2];
};

/**
 * @brief Store the output vectors of the time, the force and torque acting on particle 0,
 * the normal overlap, the tangential overlap and the velocities and angular
 * velocities for both particles.
 */
struct full_contact_output
{
  std::vector<double>       time;
  std::vector<Tensor<1, 3>> force;
  std::vector<Tensor<1, 3>> torque;
  std::vector<double>       normal_overlap;
  std::vector<Tensor<1, 3>> tangential_overlap;
  std::vector<Tensor<1, 3>> velocity0;
  std::vector<Tensor<1, 3>> velocity1;
  std::vector<Tensor<1, 3>> omega0;
  std::vector<Tensor<1, 3>> omega1;
};


/**
 * @brief Set the force and torque to 0 for each particle.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param torque Vector of torques acting on each particle.
 * @param force Vector of forces acting on each particle.
 */
template <int dim>
void
reinitialize_force(Particles::ParticleHandler<dim> &particle_handler,
                   std::vector<Tensor<1, 3>>       &torque,
                   std::vector<Tensor<1, 3>>       &force)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Getting id of particle as local variable
      unsigned int particle_id = particle->get_id();

      // Reinitializing forces and torques of particles in the system
      force[particle_id][0] = 0;
      force[particle_id][1] = 0;
      force[particle_id][2] = 0;

      torque[particle_id][0] = 0;
      torque[particle_id][1] = 0;
      torque[particle_id][2] = 0;
    }
}

/**
 * @brief Update the velocity tensors with the values from the PropertiesIndex.
 *
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param particle0_properties Properties of particle 0.
 * @param particle1_properties Properties of particle 1.
 * @param velocity0, velocity1, omega0, omega1 Velocities of particles 0 and 1.
 *
 */
template <typename PropertiesIndex>
void
update_velocities(const ArrayView<const double> &particle0_properties,
                  const ArrayView<const double> &particle1_properties,
                  Tensor<1, 3>                  &velocity0,
                  Tensor<1, 3>                  &velocity1,
                  Tensor<1, 3>                  &omega0,
                  Tensor<1, 3>                  &omega1)
{
  velocity0[0] = particle0_properties[PropertiesIndex::v_x];
  velocity0[1] = particle0_properties[PropertiesIndex::v_y];
  velocity0[2] = particle0_properties[PropertiesIndex::v_z];
  velocity1[0] = particle1_properties[PropertiesIndex::v_x];
  velocity1[1] = particle1_properties[PropertiesIndex::v_y];
  velocity1[2] = particle1_properties[PropertiesIndex::v_z];
  omega0[0]    = particle0_properties[PropertiesIndex::omega_x];
  omega0[1]    = particle0_properties[PropertiesIndex::omega_y];
  omega0[2]    = particle0_properties[PropertiesIndex::omega_z];
  omega1[0]    = particle1_properties[PropertiesIndex::omega_x];
  omega1[1]    = particle1_properties[PropertiesIndex::omega_y];
  omega1[2]    = particle1_properties[PropertiesIndex::omega_z];
}

using namespace Parameters::Lagrangian;

/**
 * @brief Simulate the full contact between two particles with chosen force model.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @tparam contact_model Model for the contact force.
 * @tparam rolling_friction_model Rolling resistance model.
 * @param triangulation Triangulation to access the information of the cells.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param contact_manager Manages the contact between particles.
 * @param dem_parameters DEM parameters declared.
 * @param p Initial properties of particles.
 * @param g A constant volumetric body force applied to all particles
 * @param dt Time Step.
 * @param output_interval Interval of iterations at which data is written on file and stored for output.
 * @param neighborhood_threshold Threshold value of contact detection.
 * @param cut_off_factor Factor to allow for contact forces even when particles are not in contact.
 * @param filename Name of file where force and overlap are written.
 *
 * @return Struct of vectors of time, force and torque acting on particle 0, normal overlap,
 * tangential overlap, velocities, angular velocities.
 */
template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
full_contact_output
simulate_full_contact(parallel::distributed::Triangulation<dim> &triangulation,
                      Particles::ParticleHandler<dim>         &particle_handler,
                      DEMContactManager<dim, PropertiesIndex> &contact_manager,
                      DEMSolverParameters<dim>                &dem_parameters,
                      initial_particle_properties<dim>        &p,
                      const Tensor<1, 3>                      &g,
                      const double                             dt,
                      const unsigned int                       output_interval,
                      const double neighborhood_threshold,
                      const double cut_off_factor,
                      std::string  filename)

{
  // Clearing particle handler
  particle_handler.clear_particles();

  // Struct for output
  full_contact_output output;

  // Constructing particle iterators from particle positions (inserting
  // particles)
  Particles::ParticleIterator<dim> pit1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, p.position[0], p.id[0]);
  Particles::ParticleIterator<dim> pit2 = construct_particle_iterator<dim>(
    particle_handler, triangulation, p.position[1], p.id[1]);

  // Setting particle properties
  set_particle_properties<dim, PropertiesIndex>(
    pit1, p.type[0], p.diameter[0], p.mass[0], p.v[0], p.omega[0]);
  set_particle_properties<dim, PropertiesIndex>(
    pit2, p.type[1], p.diameter[1], p.mass[1], p.v[1], p.omega[1]);

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;
  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  for (auto &moi_val : MOI)
    moi_val = 1;

  // Force and integrator objects
  VelocityVerletIntegrator<dim, PropertiesIndex> integrator_object;
  ParticleParticleContactForce<dim,
                               PropertiesIndex,
                               contact_model,
                               rolling_friction_model>
    force_object(dem_parameters);

  // Defining local variables
  auto          particle0          = particle_handler.begin();
  auto          particle1          = std::next(particle0);
  bool          contact_is_ongoing = true;
  double        time               = 0;
  const int     max_iteration      = 10000;
  double        normal_overlap     = -1;
  double        distance;
  Point<3>     &position0 = particle0->get_location();
  Point<3>     &position1 = particle1->get_location();
  Tensor<1, 3> &force0    = force[particle0->get_id()];
  Tensor<1, 3> &torque0   = torque[particle0->get_id()];
  Tensor<1, 3>  tangential_overlap{{0, 0, 0}};
  Tensor<1, 3>  velocity0;
  Tensor<1, 3>  velocity1;
  Tensor<1, 3>  omega0;
  Tensor<1, 3>  omega1;
  Tensor<1, 3>  contact_vector;
  Tensor<1, 3>  normal_unit_vector;
  Tensor<1, 3>  relative_velocity;
  Tensor<1, 3>  tangential_relative_velocity;
  const double  force_calculation_threshold_distance =
    cut_off_factor * 0.5 * (p.diameter[0] + p.diameter[1]);

  // Open file and write names of columns
  std::ofstream file(filename + ".dat", std::ios::binary);
  file
    << "time force0_x force0_y force0_z torque0_x torque0_y torque0_z normal_overlap tangential_overlap_x tangential_overlap_y tangential_overlap_z velocity0_x velocity0_y velocity0_z velocity1_x velocity1_y velocity1_z omega0_x omega0_y omega0_z omega1_x omega1_y omega1_z"
    << "\n";

  for (int iteration = 0; iteration < max_iteration; ++iteration)
    {
      // Reinitializing forces
      reinitialize_force(particle_handler, torque, force);

      contact_manager.update_local_particles_in_cells(particle_handler);

      // Dummy Adaptive sparse contacts object and particle-particle broad
      // search
      AdaptiveSparseContacts<dim, PropertiesIndex>
        dummy_adaptive_sparse_contacts;
      contact_manager.execute_particle_particle_broad_search(
        particle_handler, dummy_adaptive_sparse_contacts);

      // Calling fine search
      contact_manager.execute_particle_particle_fine_search(
        neighborhood_threshold);

      // Calling force calculation
      force_object.calculate_particle_particle_contact_force(
        contact_manager.get_local_adjacent_particles(),
        contact_manager.get_ghost_adjacent_particles(),
        contact_manager.get_local_local_periodic_adjacent_particles(),
        contact_manager.get_local_ghost_periodic_adjacent_particles(),
        contact_manager.get_ghost_local_periodic_adjacent_particles(),
        dt,
        torque,
        force);

      distance = (position0).distance(position1);

      // Checking if contact is still ongoing (overlap superior to force
      // calculation threshold)
      if (normal_overlap > force_calculation_threshold_distance &&
          (0.5 * (p.diameter[0] + p.diameter[1]) - distance) <=
            force_calculation_threshold_distance)
        {
          contact_is_ongoing = false;
        }

      // Calculating overlap and tangential overlap
      update_velocities<PropertiesIndex>(particle0->get_properties(),
                                         particle1->get_properties(),
                                         velocity0,
                                         velocity1,
                                         omega0,
                                         omega1);
      normal_overlap     = 0.5 * (p.diameter[0] + p.diameter[1]) - distance;
      contact_vector     = position1 - position0;
      normal_unit_vector = contact_vector / contact_vector.norm();
      relative_velocity  = velocity0 - velocity1 +
                          (cross_product_3d(0.5 * (p.diameter[0] * omega0 +
                                                   p.diameter[1] * omega1),
                                            normal_unit_vector));
      tangential_relative_velocity =
        relative_velocity -
        ((relative_velocity * normal_unit_vector) * normal_unit_vector);
      tangential_overlap += tangential_relative_velocity * dt;

      // Printing on file and storing data at each output interval
      //  and at the first iteration after end of contact
      if (iteration % output_interval == 0 || !contact_is_ongoing)
        {
          file << time << " " << force0[0] << " " << force0[1] << " "
               << force0[2] << " " << torque0[0] << " " << torque0[1] << " "
               << torque0[2] << " " << normal_overlap << " "
               << tangential_overlap[0] << " " << tangential_overlap[1] << " "
               << tangential_overlap[2] << " " << velocity0[0] << " "
               << velocity0[1] << " " << velocity0[2] << " " << velocity1[0]
               << " " << velocity1[1] << " " << velocity1[2] << " " << omega0[0]
               << " " << omega0[1] << " " << omega0[2] << " " << omega1[0]
               << " " << omega1[1] << " " << omega1[2] << "\n";

          output.time.push_back(time);
          output.force.push_back(force0);
          output.torque.push_back(torque0);
          output.normal_overlap.push_back(normal_overlap);
          output.tangential_overlap.push_back(tangential_overlap);
          output.velocity0.push_back(velocity0);
          output.velocity1.push_back(velocity1);
          output.omega0.push_back(omega0);
          output.omega1.push_back(omega1);
        }

      // Return output before integration if contact is done
      if (not contact_is_ongoing)
        {
          file.close();
          return output;
        }

      // Integration
      integrator_object.integrate(particle_handler, g, dt, torque, force, MOI);

      // Update contacts
      contact_manager.update_contacts();

      time += dt;
    }

  // Return empty output if contact does not start or end after
  // max_iteration iterations.
  if (normal_overlap > force_calculation_threshold_distance)
    {
      std::cout << "Error. Contact not ending after ";
    }
  else
    {
      std::cout << "Error. No contact found after ";
    }

  std::cout << max_iteration << " iterations." << std::endl;
  file.close();
  return full_contact_output{};
}



#endif // full_contact_functions_h
