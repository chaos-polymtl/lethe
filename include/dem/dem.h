/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_dem_h
#define lethe_dem_h

#include <core/dem_properties.h>
#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/force_chains_visualization.h>
#include <dem/grid_motion.h>
#include <dem/insertion.h>
#include <dem/integrator.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/load_balancing.h>
#include <dem/output_force_torque_calculation.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/particle_point_line_contact_force.h>
#include <dem/particle_wall_contact_force.h>
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

using namespace DEM;

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
   * Initializes all the required parameters and iterates over the DEM iterator
   * (DEM engine).
   */
  void
  solve();

private:
  /**
   * @brief Manages the call to the load balance by first identifying if
   * load balancing is required and then performing the load balance.
   */
  void
  load_balance();

  /**
   * @brief Sets the right iteration check function according to the chosen contact detection method.
   *
   * @return Return a function. This function returns a bool indicating if the contact search should be carried out in the current iteration.
   */
  inline std::function<bool()>
  set_contact_search_iteration_function();

  /**
   * @brief Establish if this is a contact detection iteration using the constant contact detection frequency.
   * If the iteration number is a multiple of the frequency, this iteration is
   * considered to be a contact detection iteration.
   *
   * @return bool indicating if the contact search should be carried out in the current iteration.
   */
  inline bool
  check_contact_search_iteration_constant();

  /**
   * @brief Establish if this is a contact detection iteration using the maximal displacement of the particles.
   * If this particle displacement surpasses a threshold, this iteration is a
   * contact detection iteration.
   *
   * @return bool indicating if the contact search should be carried out in the current iteration.
   */
  inline bool
  check_contact_search_iteration_dynamic();

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
    std::vector<double>                     &MOI);


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
  std::shared_ptr<ParticleParticleContactForceBase<dim>>
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
   * Sets the background degree of freedom used for parallel grid output
   *
   */
  void
  setup_background_dofs();

  /**
   * @brief write_output_results
   * Generates VTU file with particles information for visualization and
   * post-processing
   * Post-processing as parallel VTU files
   */
  void
  write_output_results();

  /**
   * @brief post_process_results Calculates average velocity and other post-processed quantities that need to be outputted to files
   */
  void
  post_process_results();


  /**
   * @brief reports_statistics Calculates statistics on the particles and report them to the terminal. This function is notably used to calculate min/max/avg/total values of the linear kinetic energy, angular kinetic energy, linear velocity, angular velocity
   */
  void
  report_statistics();

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
  double                               smallest_floating_mesh_mapping_criterion;
  Particles::ParticleHandler<dim, dim> particle_handler;
  bool                                 contact_detection_step;
  bool                                 load_balance_iteration;
  bool                                 checkpoint_step;
  Tensor<1, 3>                         g;

  DEM::DEMProperties<dim>                  properties_class;
  std::vector<std::pair<std::string, int>> properties =
    properties_class.get_properties_name();
  double             neighborhood_threshold_squared;
  double             maximum_particle_diameter;
  const unsigned int contact_detection_frequency;
  const unsigned int insertion_frequency;

  // Initialization of classes and building objects
  DEMContactManager<dim>             contact_manager;
  LoadBalancing<dim>                 load_balancing;
  std::shared_ptr<SimulationControl> simulation_control;
  BoundaryCellsInformation<dim>      boundary_cell_object;
  std::shared_ptr<GridMotion<dim>>   grid_motion_object;
  ParticlePointLineForce<dim>        particle_point_line_contact_force_object;
  std::shared_ptr<Integrator<dim>>   integrator_object;
  std::shared_ptr<Insertion<dim>>    insertion_object;
  std::shared_ptr<ParticleParticleContactForceBase<dim>>
    particle_particle_contact_force_object;
  std::shared_ptr<ParticlesForceChainsBase<dim>> particles_force_chains_object;
  std::shared_ptr<ParticleWallContactForce<dim>>
                                particle_wall_contact_force_object;
  Visualization<dim>            visualization_object;
  LagrangianPostProcessing<dim> post_processing_object;
  PVDHandler                    particles_pvdhandler;
  PVDHandler                    particles_pvdhandler_force_chains;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       displacement;
  std::vector<double>       MOI;

  // Mesh and boundary information
  typename dem_data_structures<dim>::solid_surfaces_mesh_information
    solid_surfaces_mesh_info;
  typename dem_data_structures<dim>::solid_volumes_mesh_info
    solid_volumes_mesh_info;
  typename dem_data_structures<dim>::boundary_points_and_normal_vectors
    updated_boundary_points_and_normal_vectors;
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    periodic_boundaries_cells_information;

  // Information for periodic boundaries
  PeriodicBoundariesManipulator<dim> periodic_boundaries_object;
  Tensor<1, dim>                     periodic_offset;
  bool                               has_periodic_boundaries;

  // Information for parallel grid processing
  DoFHandler<dim> background_dh;
  PVDHandler      grid_pvdhandler;
  bool            has_solid_objects;

  // Storage of statistics about time and contact lists
  statistics contact_list;
  statistics simulation_time;

  // Solid DEM objects
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solid_surfaces;
  std::vector<std::shared_ptr<SerialSolid<dim, dim>>>     solid_volumes;

  // Distribution objects
  std::vector<std::shared_ptr<Distribution>> size_distribution_object_container;

  // Adaptive sparce contacts (ASC) in cells object
  AdaptiveSparseContacts<dim> sparse_contacts_object;

  // Flag to indicate if sparse contacts are enabled
  bool has_sparse_contacts;

  // Contraints for the background grid needed for ASC with PBC
  AffineConstraints<double> background_constraints;

  // Load balancing iteration check function
  std::function<bool()> load_balance_iteration_check_function;

  // Contact detection iteration check function
  std::function<bool()> contact_detection_iteration_check_function;
};

#endif
