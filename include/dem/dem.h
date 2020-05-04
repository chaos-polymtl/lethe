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
 * Author: Bruno Blais, Shahab Golshan, Polytechnique Montreal, 2019-
 */

#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/integrator.h>
#include <dem/non_uniform_insertion.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_contact_force.h>
#include <dem/particle_point_line_fine_search.h>
#include <dem/pp_broad_search.h>
#include <dem/pp_contact_force.h>
#include <dem/pp_contact_info_struct.h>
#include <dem/pp_fine_search.h>
#include <dem/pp_linear_force.h>
#include <dem/pp_nonlinear_force.h>
#include <dem/pw_broad_search.h>
#include <dem/pw_contact_force.h>
#include <dem/pw_contact_info_struct.h>
#include <dem/pw_fine_search.h>
#include <dem/pw_linear_force.h>
#include <dem/pw_nonlinear_force.h>
#include <dem/uniform_insertion.h>
#include <dem/velocity_verlet_integrator.h>
#include <dem/visualization.h>
#include <dem/write_vtu.h>

#include <fstream>
#include <iostream>

#ifndef LETHE_DEM_H
#define LETHE_DEM_H

/**
 * The DEM class which initializes all the required parameters and iterates over
 * the DEM iterator
 */

template <int dim> class DEMSolver {
public:
  DEMSolver(DEMSolverParameters<dim> dem_parameters);

  /**
   * Initialiazes all the required parameters and iterates over the DEM iterator
   * (DEM engine).
   */
  void solve();

private:
  MPI_Comm mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
  ConditionalOStream pcout;
  DEMSolverParameters<dim> parameters;
  parallel::distributed::Triangulation<dim> triangulation;
  Particles::PropertyPool property_pool;
  MappingQGeneric<dim> mapping;
  TimerOutput computing_timer;
  Particles::ParticleHandler<dim, dim> particle_handler;
  int DEM_step = 0;
  double DEM_time = 0.0;
  const int number_of_steps = parameters.simulationControl.final_time_step;
  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
      cell_neighbor_list;
  std::vector<typename Triangulation<dim>::active_cell_iterator>
      boundary_cells_with_faces;
  std::vector<std::tuple<typename Triangulation<dim>::active_cell_iterator,
                         Point<dim>, Point<dim>>>
      boundary_cells_with_lines;
  std::vector<
      std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
      boundary_cells_with_points;
  std::vector<boundary_cells_info_struct<dim>> boundary_cells_information;
  std::map<int, std::pair<Particles::ParticleIterator<dim>,
                          Particles::ParticleIterator<dim>>>
      contact_pair_candidates;
  std::map<int, std::tuple<std::pair<Particles::ParticleIterator<dim>, int>,
                           Tensor<1, dim>, Point<dim>>>
      pw_contact_candidates;
  std::map<int, std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
      particle_point_contact_candidates;
  std::map<int,
           std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
      particle_line_contact_candidates;
  std::map<int, particle_point_line_contact_info_struct<dim>>
      particle_points_in_contact, particle_lines_in_contact;

  std::map<int, Particles::ParticleIterator<dim>> particle_container;
  DEM::DEMProperties<dim> properties_class;

  // Initilization of classes and building objects
  PPBroadSearch<dim> pp_broad_search_object;
  PPFineSearch<dim> pp_fine_search_object;
  PWBroadSearch<dim> pw_broad_search_object;
  ParticlePointLineBroadSearch<dim> particle_point_line_broad_search_object;
  PWFineSearch<dim> pw_fine_search_object;
  ParticlePointLineFineSearch<dim> particle_point_line_fine_search_object;
  ParticlePointLineForce<dim> particle_point_line_contact_force_object;
  std::shared_ptr<Integrator<dim>> integrator_object;
  std::shared_ptr<PPContactForce<dim>> pp_contact_force_object;
  std::shared_ptr<PWContactForce<dim>> pw_contact_force_object;

  /**
   * Defines or reads the mesh based on the information provided by the user.
   * Gmsh files can also be read in this function.
   */
  void read_mesh();

  /**
   * Reinitializes exerted forces and momentums on particles
   *
   * @param particle_handler Particle handler to access all the particles in the
   * system
   */
  void reinitialize_force(Particles::ParticleHandler<dim> &particle_handler);

  /**
   * Updates the iterators to particles in a map of particles
   * (particle_container) after calling sorting particles in cells function
   *
   * @param particle_handler Particle handler to access all the particles in the
   * system
   * @return particle_container A map of particles which is used to update the
   * iterators to particles in pp and pw fine search outputs after calling sort
   * particles into cells function
   */
  std::map<int, Particles::ParticleIterator<dim>> update_particle_container(
      const Particles::ParticleHandler<dim> &particle_handler);

  /**
   * Updates the iterators to particles in pp_contact_container (output of pp
   * fine search)
   *
   * @param pairs_in_contact_info Output of particle-particle fine search
   * @param particle_container Output of update_particle_container function
   */
  void update_pp_contact_container_iterators(
      std::map<int, std::map<int, pp_contact_info_struct<dim>>>
          &pairs_in_contact_info,
      const std::map<int, Particles::ParticleIterator<dim>>
          &particle_container);

  /**
   * Updates the iterators to particles in pw_contact_container (output of pw
   * fine search)
   *
   * @param pw_pairs_in_contact Output of particle-wall fine search
   * @param particle_container Output of update_particle_container function
   */
  void update_pw_contact_container_iterators(
      std::map<int, std::map<int, pw_contact_info_struct<dim>>>
          &pw_pairs_in_contact,
      const std::map<int, Particles::ParticleIterator<dim>>
          &particle_container);

  /**
   * Updates the iterators to particles in particle_points_in_contact and
   * particle_lines_in_contact (output of particle point line fine search)
   *
   * @param particle_points_in_contact Output of particle-point fine search
   * @param particle_lines_in_contact Output of particle-line fine search
   * @param particle_container Output of update_particle_container function
   */
  void update_particle_point_line_contact_container_iterators(
      std::map<int, particle_point_line_contact_info_struct<dim>>
          &particle_points_in_contact,
      std::map<int, particle_point_line_contact_info_struct<dim>>
          &particle_lines_in_contact,
      const std::map<int, Particles::ParticleIterator<dim>>
          &particle_container);

  /**
   * Sets the chosen integration method in the parameter handler file
   *
   * @param dem_parameters DEM parameters
   * @return A pointer to the integration object
   */
  std::shared_ptr<Integrator<dim>>
  set_integrator_type(const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Sets the chosen particle-particle contact force model in the parameter
   * handler file
   *
   * @param dem_parameters DEM parameters
   * @return A pointer to the particle-particle contact force object
   */
  std::shared_ptr<PPContactForce<dim>>
  set_pp_contact_force(const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Sets the chosen particle-wall contact force model in the parameter handler
   * file
   *
   * @param dem_parameters DEM parameters
   * @return A pointer to the particle-wall contact force object
   */
  std::shared_ptr<PWContactForce<dim>>
  set_pw_contact_force(const DEMSolverParameters<dim> &dem_parameters);
};

#endif
