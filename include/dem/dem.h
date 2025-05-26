// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_h
#define lethe_dem_h

#include <core/checkpoint_control.h>
#include <core/dem_properties.h>
#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem_action_manager.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_contact_detection_step.h>
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
 * @brief Solver using the soft-sphere model of the discrete element method
 * (DEM) to simulate granular systems.
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class DEMSolver
{
public:
  DEMSolver(DEMSolverParameters<dim> dem_parameters);

  /**
   * @brief Engine of the DEM solver. Calls all the necessary functions to set
   * parameters, solve the simulation, and finish the simulation.
   */
  void
  solve();

private:
  /**
   * @brief Set the parameters for the DEM simulation
   */
  void
  setup_parameters();

  /**
   * @brief Initialize the distribution type for the particles, and sets the
   * maximum particle diameter and the neighborhood threshold squared in the
   * process.
   */
  void
  setup_distribution_type();

  /**
   * @brief Set up the solid object containers with serial objects.
   * TODO: Modify this set up once the solid volumes are fully implemented.
   */
  void
  setup_solid_objects();

  /**
   * @brief Set up the pointers to the functions and classes according to the
   * parameters.
   */
  void
  setup_functions_and_pointers();

  /**
   * @brief Set the iteration check function according to the chosen
   * contact detection method.
   *
   * @return Contact detection method function.
   */
  inline std::function<void()>
  set_contact_search_iteration_function();

  /**
   * @brief Set the insertion method.
   *
   * @return The pointer to the insertion object
   */
  std::shared_ptr<Insertion<dim, PropertiesIndex>>
  set_insertion_type();

  /**
   * @brief Set the integration method.
   *
   * @return The pointer to the integration object
   */
  std::shared_ptr<Integrator<dim, PropertiesIndex>>
  set_integrator_type();

  /**
   * @brief Set up the triangulation dependent parameters after the reading of
   * the mesh and/or the simulation restart.
   */
  void
  setup_triangulation_dependent_parameters();

  /**
   * @brief Set up the background degree of freedom used for parallel grid
   * output or mapping of periodic dofs when using periodic boundaries.
   */
  void
  setup_background_dofs();

  /**
   * @brief Check if a load balancing is required according to the load
   * balancing method and perform it if necessary.
   */
  void
  load_balance();

  /**
   * @brief Establish if this is a contact detection iteration using the
   * constant contact detection frequency.
   * If the iteration number is a multiple of the frequency, this iteration is
   * considered to be a contact detection iteration.
   */
  inline void
  check_contact_search_iteration_constant();

  /**
   * @brief Establish if this is a contact detection iteration using the
   * maximal displacement of the particles. If this particle displacement
   * surpasses a threshold, this iteration is a contact detection iteration.
   */
  inline void
  check_contact_search_iteration_dynamic();

  /**
   * @brief Check if particles have to be inserted at this iteration and
   * perform it if necessary.
   */
  void
  insert_particles();

  /**
   * @brief Calculates particles-wall contact forces.
   */
  void
  particle_wall_contact_force();

  /**
   * @brief Move all solid objects.
   */
  void
  move_solid_objects();

  /**
   * @brief Execute the last post-processing at the end of the simulation and
   * output test results if necessary.
   */
  void
  finish_simulation();

  /**
   * @brief Generate VTU file with particles information for visualization.
   */
  void
  write_output_results();

  /**
   * @brief Calculate average velocity and other post-processed quantities that
   * need to be outputted to files.
   */
  void
  post_process_results();

  /**
   * @brief Calculate statistics on the particles and report them to the
   * terminal. This function is notably used to monitor the time min, max and
   * total performed contact searches, and the instant min, max, avg and total
   * values of the linear velocity, angular velocity, linear kinetic energy and
   * angular kinetic.
   */
  void
  report_statistics();

  /**
   * @brief Execute the sorting of particle into subdomains and cells, and
   * reinitialize the containers dependent on the local particle ids.
   */
  void
  sort_particles_into_subdomains_and_cells();

  /**
   * @brief The MPI communicator used for the parallel simulation.
   */
  MPI_Comm mpi_communicator;

  /**
   * @brief The number of MPI processes used for the parallel simulation.
   */
  const unsigned int n_mpi_processes;

  /**
   * @brief The rank of the current MPI process.
   */
  const unsigned int this_mpi_process;

  /**
   * @brief The output stream used for the parallel simulation, only used by
   * the process 0.
   */
  ConditionalOStream pcout;

  /**
   * @brief The parameters of the DEM simulation.
   */
  DEMSolverParameters<dim> parameters;

  /**
   * @brief The checkpoint controller of the DEM simulation.
   */
  CheckpointControl checkpoint_controller;

  /**
   * @brief The distributed triangulation used for the DEM simulation.
   */
  parallel::distributed::Triangulation<dim> triangulation;

  /**
   * @brief The polynomial mapping. The mapping is first order in DEM
   * simulations.
   */
  MappingQGeneric<dim> mapping;

  /**
   * @brief The particle handler that manages the particles.
   */
  Particles::ParticleHandler<dim, dim> particle_handler;

  /**
   * @brief The action manager that manages the actions triggered by events.
   */
  DEMActionManager *action_manager;

  /**
   * @brief The timer that keeps track of the time spent in some functions.
   * Currently theses functions are: load balancing and VTU output.
   */
  TimerOutput computing_timer;

  /**
   * @brief The properties of the DEM simulation.
   * @tparam dim An integer that denotes the number of spatial dimensions.
   * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
   */
  DEM::ParticleProperties<dim, PropertiesIndex> properties_class;

  /**
   * @brief The acceleration acting on the particles.
   */
  Tensor<1, 3> g;

  /**
   * @brief The contact search frequency criterion for the dynamic contact
   * detection method.
   * The value is the smallest value between those two quantities:
   * - The diameter of the smallest active cell of a triangulation minus the
   *   largest particle radius.
   *   \f$D_{c,min} - D_{p,max} / 2\f$
   * - A security factor on the neighboring threshold times the largest particle
   *   radius.
   *   \f$factor * (neighborhood_threshold - 1) * D_{p,max} / 2\f$
   */
  double smallest_contact_search_criterion;

  /**
   * @brief The smallest solid object mapping criterion.
   * The value is \f$2^-0.5 * D_{c,min}\f$
   */
  double smallest_solid_object_mapping_criterion;

  /**
   * @brief The neighborhood threshold squared.
   *        \f$(\text{neighborhood_threshold} * D_{p,max})^2\f$
   */
  double neighborhood_threshold_squared;

  /**
   * @brief The maximum particle diameter.
   */
  double maximum_particle_diameter;

  /**
   * @brief The counter to keep track of the number of contact search.
   */
  unsigned int contact_build_number;

  /**
   * @brief The manager of all the contact search operations.
   */
  DEMContactManager<dim, PropertiesIndex> contact_manager;

  /**
   * @brief The load balancing handler.
   */
  LagrangianLoadBalancing<dim, PropertiesIndex> load_balancing;

  /**
   * @brief The simulation control (DEM Transient).
   */
  std::shared_ptr<SimulationControl> simulation_control;

  /**
   * @brief The boundary cells object.
   */
  BoundaryCellsInformation<dim> boundary_cell_object;

  /**
   * @brief The grid motion object.
   */
  std::shared_ptr<GridMotion<dim>> grid_motion_object;

  /**
   * @brief The particle-particle contact force object.
   */
  std::shared_ptr<ParticleParticleContactForceBase<dim, PropertiesIndex>>
    particle_particle_contact_force_object;

  /**
   * @brief The particle-wall contact force object.
   */
  std::shared_ptr<ParticleWallContactForceBase<dim, PropertiesIndex>>
    particle_wall_contact_force_object;

  /**
   * @brief The particle-point-line contact force object.
   */
  ParticlePointLineForce<dim, PropertiesIndex>
    particle_point_line_contact_force_object;

  /**
   * @brief The particle force chains post-processing object.
   */
  std::shared_ptr<ParticlesForceChainsBase<dim, PropertiesIndex>>
    particles_force_chains_object;

  /**
   * @brief The integrator object.
   */
  std::shared_ptr<Integrator<dim, PropertiesIndex>> integrator_object;

  /**
   * @brief The insertion object.
   */
  std::shared_ptr<Insertion<dim, PropertiesIndex>> insertion_object;

  /**
   * @brief The visualization object.
   */
  Visualization<dim, PropertiesIndex> visualization_object;

  /**
   * @brief The particle PVD handler.
   */
  PVDHandler particles_pvdhandler;

  /**
   * @brief The force chains PVD handler.
   */
  PVDHandler particles_pvdhandler_force_chains;

  /**
   * @brief Class object to store the vectors of force, torque and heat transfer rate applied to particles.
   */
  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;

  /**
   * @brief The vector of torque of particles.
   */
  std::vector<Tensor<1, 3>> &torque = contact_outcome.torque;

  /**
   * @brief The vector of force of particles.
   */
  std::vector<Tensor<1, 3>> &force = contact_outcome.force;

  /**
   * @brief The displacement tracking of particles for the dynamic contact
   * detection.
   */
  std::vector<double> displacement;

  /**
   * @brief The moment of inertia of particles.
   */
  std::vector<double> MOI;

  /**
   * @brief The vector of vector of pairs of the mapping of the background mesh
   * to the solid surface.
   */
  typename dem_data_structures<dim>::solid_surfaces_mesh_information
    solid_surfaces_mesh_info;

  /**
   * @brief The vector of vector of pairs of the mapping of the background mesh
   * to the solid volume.
   */
  typename dem_data_structures<dim>::solid_volumes_mesh_info
    solid_volumes_mesh_info;

  /**
   * @brief The map of the point and normal vector of grid.
   * Only used with grid motion.
   */
  typename dem_data_structures<dim>::boundary_points_and_normal_vectors
    updated_boundary_points_and_normal_vectors;

  /**
   * @brief The boundary information of the forces.
   */
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;

  /**
   * @brief The torques boundary information.
   */
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;

  /**
   * @brief The periodic boundaries cells information map.
   * Only used with a simulation with periodic boundaries.
   */
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    periodic_boundaries_cells_information;

  /**
   * @brief The periodic boundaries manipulator.
   * Only used with a simulation with periodic boundaries.
   */
  PeriodicBoundariesManipulator<dim> periodic_boundaries_object;

  /**
   * @brief The background degree of freedom handler uses for parallel grid
   * processing.
   */
  DoFHandler<dim> background_dh;

  /**
   * @brief The grid PVD handler.
   */
  PVDHandler grid_pvdhandler;

  /**
   * @brief The storage of statistics about time and contact lists
   */
  statistics contact_list;

  /**
   * @brief The container of the solid surfaces.
   */
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solid_surfaces;

  /**
   * @brief The container of the solid volumes.
   */
  std::vector<std::shared_ptr<SerialSolid<dim, dim>>> solid_volumes;

  /**
   * @brief The container of the distribution objects.
   */
  std::vector<std::shared_ptr<Distribution>> size_distribution_object_container;

  /**
   * @brief The object handling the adaptive sparse contacts (ASC).
   * Only used with the ASC method.
   */
  AdaptiveSparseContacts<dim, PropertiesIndex> sparse_contacts_object;

  /**
   * @brief The constraints for the background grid needed for the adaptive sparse.
   */
  AffineConstraints<double> background_constraints;

  /**
   * @brief The pointer to the function that checks if a load balancing is
   * required.
   */
  std::function<bool()> load_balance_iteration_check_function;

  /**
   * @brief The pointer to the contact detection function that checks if a
   * contact search is required.
   */
  std::function<void()> contact_detection_iteration_check_function;

  /**
   * @brief Disable position integration, useful for multiphysic simulations
   * with a packed bed, loaded with another prm.
   */
  bool disable_position_integration;
};

#endif
