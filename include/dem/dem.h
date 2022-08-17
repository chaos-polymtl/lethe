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
 */

#include <core/dem_properties.h>
#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/grid_motion.h>
#include <dem/insertion.h>
#include <dem/integrator.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/localize_contacts.h>
#include <dem/locate_local_particles.h>
#include <dem/non_uniform_insertion.h>
#include <dem/output_force_torque_calculation.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/particle_particle_fine_search.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_contact_force.h>
#include <dem/particle_point_line_fine_search.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>
#include <dem/periodic_boundaries_manipulator.h>
#include <dem/visualization.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/particles/particle_handler.h>

#include <fstream>
#include <iostream>
#include <unordered_set>

#ifndef Lethe_DEM_h
#  define Lethe_DEM_h

using namespace dem_data_containers;

/**
 * The DEM class which initializes all the required parameters and iterates over
 * the DEM iterator
 */
template <int dim>
class DEMSolver
{
  using FuncPtrType = bool (DEMSolver<dim>::*)();
  FuncPtrType check_contact_search_step;
  FuncPtrType check_load_balance_step;

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
   * The cell_weight() function indicates to the triangulation how much
   * computational work is expected to happen on this cell, and consequently
   * how the domain needs to be partitioned so that every MPI rank receives a
   * roughly equal amount of work (potentially not an equal number of cells).
   * While the function is called from the outside, it is connected to the
   * corresponding signal from inside this class, therefore it can be private.
   * This function is the key component that allow us to dynamically balance the
   * computational load. The function attributes a weight to
   * every cell that represents the computational work on this cell. Here the
   * majority of work is expected to happen on the particles, therefore the
   * return value of this function (representing "work for this cell") is
   * calculated based on the number of particles in the current cell.
   * The function is connected to the cell_weight() signal inside the
   * triangulation, and will be called once per cell, whenever the triangulation
   * repartitions the domain between ranks (the connection is created inside the
   * particles_generation() function of this class).
   */
  unsigned int
  cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
      &                                                                  cell,
    const typename parallel::distributed::Triangulation<dim>::CellStatus status)
    const;

  /**
   * Finds contact search steps for constant contact search method
   */
  inline bool
  check_contact_search_step_constant();

  /**
   * Finds contact search steps for dynamic contact search method
   */
  inline bool
  check_contact_search_step_dynamic();

  /**
   * Finds load-balance step for single-step load-balance
   */
  inline bool
  check_load_balance_once();

  /**
   * For cases where load balance method is equal to none
   */
  inline bool
  no_load_balance();

  /**
   * Finds load-balance step for frequent load-balance
   */
  inline bool
  check_load_balance_frequent();

  /**
   * Finds load-balance step for dynamic load-balance
   */
  inline bool
  check_load_balance_dynamic();

  /**
   * @brief Manages the call to the load balancing. Returns true if
   * load balancing is performed
   *
   */
  void
  load_balance();

  /**
   * @brief Manages the call to the particle insertion. Returns true if
   * particles were inserted
   *
   */
  bool
  insert_particles();

  /**
   * @brief Updates moment of inertia container after sorting particles
   * into subdomains
   *
   */
  void
  update_moment_of_inertia(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    std::vector<double> &                    MOI);

  /**
   * @brief Carries out the broad contact detection search using the
   * background triangulation for particle-walls contact
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
  std::shared_ptr<ParticleParticleContactForce<dim>>
  set_particle_particle_contact_force(
    const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Sets the chosen particle-wall contact force model in the parameter handler
   * file
   *
   * @param dem_parameters DEM parameters
   * @return A pointer to the particle-wall contact force object
   */
  std::shared_ptr<ParticleWallContactForce<dim>>
  set_particle_wall_contact_force(
    const DEMSolverParameters<dim> &dem_parameters);

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

  /**
   * @brief post_process_results
   */
  void
  post_process_results();

  MPI_Comm                                  mpi_communicator;
  const unsigned int                        n_mpi_processes;
  const unsigned int                        this_mpi_process;
  ConditionalOStream                        pcout;
  DEMSolverParameters<dim>                  parameters;
  parallel::distributed::Triangulation<dim> triangulation;

  MappingQGeneric<dim>                 mapping;
  bool                                 particles_insertion_step;
  unsigned int                         contact_build_number;
  TimerOutput                          computing_timer;
  double                               smallest_contact_search_criterion;
  Particles::ParticleHandler<dim, dim> particle_handler;
  bool                                 contact_detection_step;
  bool                                 load_balance_step;
  bool                                 checkpoint_step;
  Tensor<1, 3>                         g;
  double                               triangulation_cell_diameter;

  typename dem_data_structures<dim>::cells_total_neighbor_list
    total_neighbor_list;
  typename dem_data_structures<dim>::floating_mesh_information
    floating_mesh_info;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;
  typename dem_data_structures<dim>::particle_floating_mesh_candidates
    particle_floating_mesh_candidates;
  typename dem_data_structures<dim>::particle_floating_mesh_in_contact
    particle_floating_mesh_in_contact;
  typename dem_data_structures<dim>::particle_floating_wall_candidates
    particle_floating_wall_candidates;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_floating_wall_in_contact;
  typename dem_data_structures<dim>::particle_wall_candidates
    particle_wall_candidates;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_in_contact;
  typename dem_data_structures<dim>::particle_point_candidates
    particle_point_candidates;
  typename dem_data_structures<dim>::particle_line_candidates
    particle_line_candidates;
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_points_in_contact;
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_lines_in_contact;
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_candidates;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_adjacent_particles;
  typename dem_data_structures<dim>::particle_index_iterator_map
    particle_container;
  typename dem_data_structures<dim>::boundary_points_and_normal_vectors
    updated_boundary_points_and_normal_vectors;
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;

  DEM::DEMProperties<dim>                  properties_class;
  std::vector<std::pair<std::string, int>> properties =
    properties_class.get_properties_name();
  double             neighborhood_threshold_squared;
  double             maximum_particle_diameter;
  const unsigned int contact_detection_frequency;
  const unsigned int insertion_frequency;

  // Initilization of classes and building objects
  std::shared_ptr<SimulationControl> simulation_control;
  BoundaryCellsInformation<dim>      boundary_cell_object;
  std::shared_ptr<GridMotion<dim>>   grid_motion_object;
  ParticleParticleBroadSearch<dim>   particle_particle_broad_search_object;
  ParticleParticleFineSearch<dim>    particle_particle_fine_search_object;
  ParticleWallBroadSearch<dim>       particle_wall_broad_search_object;
  ParticlePointLineBroadSearch<dim>  particle_point_line_broad_search_object;
  ParticleWallFineSearch<dim>        particle_wall_fine_search_object;
  ParticlePointLineFineSearch<dim>   particle_point_line_fine_search_object;
  ParticlePointLineForce<dim>        particle_point_line_contact_force_object;
  std::shared_ptr<Integrator<dim>>   integrator_object;
  std::shared_ptr<Insertion<dim>>    insertion_object;
  PeriodicBoundariesManipulator<dim> periodic_boundaries_object;
  std::shared_ptr<ParticleParticleContactForce<dim>>
    particle_particle_contact_force_object;
  std::shared_ptr<ParticleWallContactForce<dim>>
                                particle_wall_contact_force_object;
  Visualization<dim>            visualization_object;
  LagrangianPostProcessing<dim> post_processing_object;
  FindCellNeighbors<dim>        cell_neighbors_object;
  PVDHandler                    particles_pvdhandler;
  const double                  standard_deviation_multiplier;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       displacement;
  std::vector<double>       MOI;

  // Information for parallel grid processing
  DoFHandler<dim> background_dh;
  PVDHandler      grid_pvdhandler;
  bool            floating_mesh;

  // Solid DEM objects
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids;
};

#endif
