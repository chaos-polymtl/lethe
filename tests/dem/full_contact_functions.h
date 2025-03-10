// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#ifndef full_contact_functions_h
#define full_contact_functions_h

#include <../tests/dem/test_particles_functions.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/velocity_verlet_integrator.h>

#include <fstream>
#include <string>

template <int dim>
struct initial_particles_properties
{
  Point<dim>     position[2];
  int            id[2];
  int            type;
  Tensor<1, dim> v[2];
  Tensor<1, dim> omega[2];
  double         mass;
  double         diameter;
};

struct contact_output
{
  std::vector<double> time;
  std::vector<double> force;
  std::vector<double> overlap;
};


/**
 * @brief Set the force and torque at 0 for each particle.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param torque Torque acting on particles
 * @param force Force acting on particles
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


using namespace Parameters::Lagrangian;

/**
 * @brief Simulate the full contact between two particles with chosen force model.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @tparam Model for the contact force.
 * @tparam Rolling resistance model.
 * @param triangulation Triangulation to access the information of the cells.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param contact_manager Manages the contact between particles.
 * @param dem_parameters DEM parameters declared.
 * @param p Initial properties of particles.
 * @param g A constant volumetric body force applied to all particles
 * @param dt Time Step.
 * @param output_step Step at which force and overlap are written on file.
 * @param neighborhood_threshold Threshold value of contact detection.
 * @param filename Name of file where force and overlap are written.
 *
 * @return contact_output Struct with the vectors of time, force and overlap.
 */
template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
struct contact_output
simul_full_contact(parallel::distributed::Triangulation<dim> &triangulation,
                   Particles::ParticleHandler<dim>           &particle_handler,
                   DEMContactManager<dim, PropertiesIndex>   &contact_manager,
                   DEMSolverParameters<dim>                  &dem_parameters,
                   initial_particles_properties<dim>         &p,
                   Tensor<1, 3>                              &g,
                   double                                     dt,
                   unsigned int                               output_step,
                   double      neighborhood_threshold,
                   std::string filename)

{
  contact_output output;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  // Constructing particle iterators from particle positions (inserting
  // particles)
  Particles::ParticleIterator<dim> pit1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, p.position[0], p.id[0]);
  Particles::ParticleIterator<dim> pit2 = construct_particle_iterator<dim>(
    particle_handler, triangulation, p.position[1], p.id[1]);

  // Setting particle properties
  set_particle_properties<dim, PropertiesIndex>(
    pit1, p.type, p.diameter, p.mass, p.v[0], p.omega[0]);
  set_particle_properties<dim, PropertiesIndex>(
    pit2, p.type, p.diameter, p.mass, p.v[1], p.omega[1]);

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


  auto   particle0       = particle_handler.begin();
  auto   particle1       = std::next(particle0);
  bool   CONTACT_ONGOING = true;
  double time            = 0;
  int    iteration       = 0;
  int    max_iteration   = 10000;
  double norm_force      = 0;
  double overlap         = -1;
  double distance;

  // Open file and write names of columns
  std::ofstream file(filename + ".dat", std::ios::binary);
  file << "time,force,overlap"
       << "\n";

  while (CONTACT_ONGOING and iteration < max_iteration)
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


      distance = (particle0->get_location() - particle1->get_location()).norm();

      // Checking if contact is still ongoing (positive overlap)
      if (overlap >= 0 and (p.diameter - distance) < 0)
        {
          CONTACT_ONGOING = false;
        }

      overlap = p.diameter - distance;

      // Printing on file and storing contact force and overlap at each output
      // step and at the end of contact
      if (iteration % output_step == 0 or not CONTACT_ONGOING)
        {
          norm_force = (force[particle0->get_id()]).norm();
          file << time << "," << norm_force << "," << overlap << "\n";
          output.time.push_back(time);
          output.force.push_back(norm_force);
          output.overlap.push_back(overlap);
        }

      if (not CONTACT_ONGOING)
        break;

      // Integration
      integrator_object.integrate(particle_handler, g, dt, torque, force, MOI);

      // Update contacts
      contact_manager.update_contacts();

      time += dt;
      iteration += 1;

      if (iteration >= max_iteration)
        {
          std::cout << "Error. Contact not ending after " << max_iteration
                    << " iterations" << std::endl;
          break;
        }
    }

  file.close();
  return output;
}



#endif // full_contact_functions_h
