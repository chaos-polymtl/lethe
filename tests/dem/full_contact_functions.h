// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#ifndef full_contact_functions_h
#define full_contact_functions_h

#include <../tests/dem/test_particles_functions.h>
#include <dem/dem_contact_manager.h>
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
  Tensor<1, dim> velocity[2];
  Tensor<1, dim> omega[2];
  double         mass[2];
  double         diameter[2];
};

/**
 * @brief Store the output vectors of the time, the force and torque acting on particle 0,
 * the normal overlap, the tangential displacement and the velocities and
 * angular velocities for both particles.
 */
struct full_contact_output
{
  std::vector<double>       time;
  std::vector<Tensor<1, 3>> force;
  std::vector<Tensor<1, 3>> torque;
  std::vector<double>       normal_overlap;
  std::vector<Tensor<1, 3>> tangential_displacement;
  std::vector<Tensor<1, 3>> velocity_0;
  std::vector<Tensor<1, 3>> velocity_1;
  std::vector<Tensor<1, 3>> omega_0;
  std::vector<Tensor<1, 3>> omega_1;
};

/**
 * @brief Update the velocity tensors with the values from the PropertiesIndex.
 *
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param particle_0_properties Properties of particle 0.
 * @param particle_1_properties Properties of particle 1.
 * @param velocity_0, velocity_1, omega_0, omega_1 Velocities of particles 0 and 1.
 *
 */
template <typename PropertiesIndex>
void
update_velocities(const ArrayView<const double> &particle_0_properties,
                  const ArrayView<const double> &particle_1_properties,
                  Tensor<1, 3>                  &velocity_0,
                  Tensor<1, 3>                  &velocity_1,
                  Tensor<1, 3>                  &omega_0,
                  Tensor<1, 3>                  &omega_1)
{
  velocity_0[0] = particle_0_properties[PropertiesIndex::v_x];
  velocity_0[1] = particle_0_properties[PropertiesIndex::v_y];
  velocity_0[2] = particle_0_properties[PropertiesIndex::v_z];
  velocity_1[0] = particle_1_properties[PropertiesIndex::v_x];
  velocity_1[1] = particle_1_properties[PropertiesIndex::v_y];
  velocity_1[2] = particle_1_properties[PropertiesIndex::v_z];
  omega_0[0]    = particle_0_properties[PropertiesIndex::omega_x];
  omega_0[1]    = particle_0_properties[PropertiesIndex::omega_y];
  omega_0[2]    = particle_0_properties[PropertiesIndex::omega_z];
  omega_1[0]    = particle_1_properties[PropertiesIndex::omega_x];
  omega_1[1]    = particle_1_properties[PropertiesIndex::omega_y];
  omega_1[2]    = particle_1_properties[PropertiesIndex::omega_z];
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
 * @param output_frequency Interval of iterations at which data is written on file and stored for output.
 * @param neighborhood_threshold Threshold value of contact detection.
 * @param cut_off_factor Factor to allow for contact forces even when particles are not in contact.
 * @param filename Name of file where force and overlap are written.
 *
 * @return Struct of vectors of time, force and torque acting on particle 0, normal overlap,
 * tangential displacement, velocities, angular velocities.
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
                      const unsigned int                       output_frequency,
                      const double      neighborhood_threshold,
                      const double      cut_off_factor,
                      const std::string filename)

{
  // Clearing particle handler
  particle_handler.clear_particles();

  // Struct for output
  full_contact_output output;

  // Constructing particle iterators from particle positions (inserting
  // particles)
  Particles::ParticleIterator<dim> pit_0 = construct_particle_iterator<dim>(
    particle_handler, triangulation, p.position[0], p.id[0]);
  Particles::ParticleIterator<dim> pit_1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, p.position[1], p.id[1]);

  // Setting particle properties
  set_particle_properties<dim, PropertiesIndex>(
    pit_0, p.type[0], p.diameter[0], p.mass[0], p.velocity[0], p.omega[0]);
  set_particle_properties<dim, PropertiesIndex>(
    pit_1, p.type[1], p.diameter[1], p.mass[1], p.velocity[1], p.omega[1]);

  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  std::vector<double>                          MOI;
  particle_handler.sort_particles_into_subdomains_and_cells();
  const unsigned int number_of_particles =
    particle_handler.get_max_local_particle_index();
  contact_outcome.resize_interaction_containers(number_of_particles);
  MOI.resize(number_of_particles);
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
  auto               particle_0         = particle_handler.begin();
  auto               particle_1         = std::next(particle_0);
  bool               contact_is_ongoing = true;
  double             time               = 0;
  const unsigned int max_iteration      = 10000;
  double             normal_overlap     = -1;
  double             distance;
  Point<3>          &position_0 = particle_0->get_location();
  Point<3>          &position_1 = particle_1->get_location();
  Tensor<1, 3>      &force_0    = contact_outcome.force[particle_0->get_id()];
  Tensor<1, 3>      &torque_0   = contact_outcome.torque[particle_0->get_id()];
  Tensor<1, 3>       tangential_displacement{{0, 0, 0}};
  Tensor<1, 3>       velocity_0;
  Tensor<1, 3>       velocity_1;
  Tensor<1, 3>       omega_0;
  Tensor<1, 3>       omega_1;
  Tensor<1, 3>       normal_unit_vector;
  Tensor<1, 3>       relative_velocity;
  Tensor<1, 3>       tangential_relative_velocity;
  const double       force_calculation_threshold_distance =
    cut_off_factor * 0.5 * (p.diameter[0] + p.diameter[1]);

  // Open file and write names of columns
  std::ofstream file(filename + ".dat", std::ios::binary);
  file
    << "time force_0_x force_0_y force_0_z torque_0_x torque_0_y torque_0_z normal_overlap tangential_displacement_x tangential_displacement_y tangential_displacement_z velocity_0_x velocity_0_y velocity_0_z velocity_1_x velocity_1_y velocity_1_z omega_0_x omega_0_y omega_0_z omega_1_x omega_1_y omega_1_z"
    << "\n";

  for (unsigned int iteration = 0; iteration < max_iteration; ++iteration)
    {
      // Reinitializing contact outcomes
      reinitialize_contact_outcomes<dim, PropertiesIndex>(particle_handler,
                                                          contact_outcome);

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
      force_object.calculate_particle_particle_contact(
        contact_manager.get_local_adjacent_particles(),
        contact_manager.get_ghost_adjacent_particles(),
        contact_manager.get_local_local_periodic_adjacent_particles(),
        contact_manager.get_local_ghost_periodic_adjacent_particles(),
        contact_manager.get_ghost_local_periodic_adjacent_particles(),
        dt,
        contact_outcome);

      distance = position_0.distance(position_1);

      // Checking if contact is still ongoing (overlap superior to force
      // calculation threshold)
      if (normal_overlap > force_calculation_threshold_distance &&
          (0.5 * (p.diameter[0] + p.diameter[1]) - distance) <=
            force_calculation_threshold_distance)
        {
          contact_is_ongoing = false;
        }

      // Calculating overlap and tangential displacement
      update_velocities<PropertiesIndex>(particle_0->get_properties(),
                                         particle_1->get_properties(),
                                         velocity_0,
                                         velocity_1,
                                         omega_0,
                                         omega_1);
      normal_overlap     = 0.5 * (p.diameter[0] + p.diameter[1]) - distance;
      normal_unit_vector = (position_1 - position_0) / distance;
      relative_velocity  = velocity_0 - velocity_1 +
                          (cross_product_3d(0.5 * (p.diameter[0] * omega_0 +
                                                   p.diameter[1] * omega_1),
                                            normal_unit_vector));
      tangential_relative_velocity =
        relative_velocity -
        ((relative_velocity * normal_unit_vector) * normal_unit_vector);
      tangential_displacement += tangential_relative_velocity * dt;

      // Printing on file and storing data at each output interval
      //  and at the first iteration after end of contact
      if (iteration % output_frequency == 0 || !contact_is_ongoing)
        {
          file << time << " " << force_0[0] << " " << force_0[1] << " "
               << force_0[2] << " " << torque_0[0] << " " << torque_0[1] << " "
               << torque_0[2] << " " << normal_overlap << " "
               << tangential_displacement[0] << " "
               << tangential_displacement[1] << " "
               << tangential_displacement[2] << " " << velocity_0[0] << " "
               << velocity_0[1] << " " << velocity_0[2] << " " << velocity_1[0]
               << " " << velocity_1[1] << " " << velocity_1[2] << " "
               << omega_0[0] << " " << omega_0[1] << " " << omega_0[2] << " "
               << omega_1[0] << " " << omega_1[1] << " " << omega_1[2] << "\n";

          output.time.push_back(time);
          output.force.push_back(force_0);
          output.torque.push_back(torque_0);
          output.normal_overlap.push_back(normal_overlap);
          output.tangential_displacement.push_back(tangential_displacement);
          output.velocity_0.push_back(velocity_0);
          output.velocity_1.push_back(velocity_1);
          output.omega_0.push_back(omega_0);
          output.omega_1.push_back(omega_1);
        }

      // Return output before integration if contact is done
      if (not contact_is_ongoing)
        {
          file.close();
          return output;
        }

      // Integration
      integrator_object.integrate(particle_handler,
                                  g,
                                  dt,
                                  contact_outcome.torque,
                                  contact_outcome.force,
                                  MOI);

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
