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

#include <core/pvd_handler.h>
#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/integrator.h>
#include <dem/localize_contacts.h>
#include <dem/locate_ghost_particles.h>
#include <dem/locate_local_particles.h>
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

#include <fstream>
#include <iostream>

#ifndef LETHE_DEM_H
#  define LETHE_DEM_H

/**
 * The DEM class which initializes all the required parameters and iterates over
 * the DEM iterator
 */
template <int dim>
class DEMSolver
{
public:
  DEMSolver(DEMSolverParameters<dim> dem_parameters);

  /**
   * Initialiazes all the required parameters and iterates over the DEM iterator
   * (DEM engine).
   */
  void
  solve();

private:
  /**
   * Defines or reads the mesh based on the information provided by the user.
   * Gmsh files can also be read in this function.
   */
  void
  read_mesh();

  /**
   * Reinitializes exerted forces and momentums on particles
   *
   * @param particle_handler Particle handler to access all the particles in the
   * system
   */
  void
  reinitialize_force(Particles::ParticleHandler<dim> &particle_handler);

  /**
   * @brief Manages the call to the particle insertion. Returns true if
   * particles were inserted
   *
   */
  bool
  insert_particles();

  /**
   * @brief Carries out the broad contact detection search using the
   * background triangulation for particle-walls contact.
   *
   */
  void
  particle_wall_broad_search();

  /**
   * @brief Carries out the fine particled-wall contact detection
   *
   */
  void
  particle_wall_fine_search();

  /**
   * @brief Calculates particles-wall contact forces
   *
   */
  void
  particle_wall_contact_force();

  /**
   * @brief finish_simulation
   * Finishes the simulation by calling all
   * the post-processing elements that are required
   */
  void
  finish_simulation();

  /**
   * Sets the chosen insertion method in the parameter handler file
   *
   * @param dem_parameters DEM parameters
   * @return A pointer to the insertion object
   */
  std::shared_ptr<Insertion<dim>>
  set_insertion_type(const DEMSolverParameters<dim> &dem_parameters);

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

  /**
   * Sets the background degree of freedom used for paralle grid output
   *
   */
  void
  setup_background_dofs();

  /**
   * @brief write_output_results
   * Post-processing as parallel VTU files
   */
  void
  write_output_results();

  MPI_Comm                                  mpi_communicator;
  const unsigned int                        n_mpi_processes;
  const unsigned int                        this_mpi_process;
  ConditionalOStream                        pcout;
  DEMSolverParameters<dim>                  parameters;
  parallel::distributed::Triangulation<dim> triangulation;
  Particles::PropertyPool                   property_pool;
  MappingQGeneric<dim>                      mapping;
  TimerOutput                               computing_timer;
  Particles::ParticleHandler<dim, dim>      particle_handler;

  // Simulation control for time stepping and I/Os
  std::shared_ptr<SimulationControl> simulation_control;

  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cells_local_neighbor_list;
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cells_ghost_neighbor_list;
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    boundary_cells_with_faces;
  std::vector<std::tuple<typename Triangulation<dim>::active_cell_iterator,
                         Point<dim>,
                         Point<dim>>>
    boundary_cells_with_lines;
  std::vector<
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
                                                 boundary_cells_with_points;
  std::map<int, boundary_cells_info_struct<dim>> boundary_cells_information;
  std::map<int, std::vector<int>>                local_contact_pair_candidates;
  std::map<int, std::vector<int>>                ghost_contact_pair_candidates;
  std::map<int, std::map<int, pp_contact_info_struct<dim>>>
    local_adjacent_particles;
  std::map<int, std::map<int, pp_contact_info_struct<dim>>>
                                                            ghost_adjacent_particles;
  std::map<int, std::map<int, pw_contact_info_struct<dim>>> pw_pairs_in_contact;
  std::map<
    std::pair<int, int>,
    std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>
    pw_contact_candidates;
  std::map<int, std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
    particle_point_contact_candidates;
  std::map<int,
           std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
    particle_line_contact_candidates;
  std::map<int, particle_point_line_contact_info_struct<dim>>
    particle_points_in_contact, particle_lines_in_contact;

  std::map<int, Particles::ParticleIterator<dim>> particle_container;
  DEM::DEMProperties<dim>                         properties_class;
  std::vector<std::pair<std::string, int>>        properties =
    properties_class.get_properties_name();
  const double neighborhood_threshold =
    parameters.model_parameters.neighborhood_threshold *
    parameters.physical_properties.diameter;

  // Initilization of classes and building objects
  PPBroadSearch<dim>                   pp_broad_search_object;
  PPFineSearch<dim>                    pp_fine_search_object;
  PWBroadSearch<dim>                   pw_broad_search_object;
  ParticlePointLineBroadSearch<dim>    particle_point_line_broad_search_object;
  PWFineSearch<dim>                    pw_fine_search_object;
  ParticlePointLineFineSearch<dim>     particle_point_line_fine_search_object;
  ParticlePointLineForce<dim>          particle_point_line_contact_force_object;
  std::shared_ptr<Integrator<dim>>     integrator_object;
  std::shared_ptr<Insertion<dim>>      insertion_object;
  std::shared_ptr<PPContactForce<dim>> pp_contact_force_object;
  std::shared_ptr<PWContactForce<dim>> pw_contact_force_object;
  Visualization<dim>                   visualization_object;
  PVDHandler                           particles_pvdhandler;

  // Information for parallel grid processing
  DoFHandler<dim> background_dh;
  PVDHandler      grid_pvdhandler;
};

#endif
