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
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/nonuniform_insertion.h>
#include <dem/pp_broad_search.h>
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
      std::vector<std::map<int, pp_contact_info_struct<dim>>>
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
      std::vector<std::map<int, pw_contact_info_struct<dim>>>
          &pw_pairs_in_contact,
      const std::map<int, Particles::ParticleIterator<dim>>
          &particle_container);
};

#endif
