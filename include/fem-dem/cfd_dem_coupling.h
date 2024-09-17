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


#ifndef lethe_dem_cfd_coupling_h
#define lethe_dem_cfd_coupling_h

#include <solvers/navier_stokes_scratch_data.h>

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/dem_action_manager.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_contact_detection_step.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/periodic_boundaries_manipulator.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/fluid_dynamics_vans.h>
#include <fem-dem/postprocessing_cfd_dem.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

/**
 * @brief Solver using the GLS Volume-averaged Navier-Stokes (VANS) solver for
 * fluid and the soft-sphere model of the discrete element method (DEM) to
 * simulate solid-fluid flow with two-way coupling.
 */
template <int dim>
class CFDDEMSolver : public FluidDynamicsVANS<dim>
{
public:
  CFDDEMSolver(CFDDEMSimulationParameters<dim> &nsparam);

  ~CFDDEMSolver();

  /**
   * @brief Engine of the CFD-DEM solver. Calls all the necessary functions to
   * set parameters, solve the simulation, and finish the simulation.
   */
  virtual void
  solve() override;

private:
  /**
   * @brief Set up the various parameters related to the DEM.
   */
  void
  dem_setup_parameters();

  /**
   * @brief Initialize the distribution type for the particles, and sets the
   * maximum particle diameter and the neighborhood threshold squared in the
   * process.
   */
  void
  setup_distribution_type();

  /**
   * @brief Set the integration method.
   *
   * @return The pointer to the integration object
   */
  std::shared_ptr<Integrator<dim>>
  set_integrator_type();

  /**
   * @brief Initialize some DEM parameters.
   */
  void
  initialize_dem_parameters();

  /**
   * @brief Read the DEM restart files from a DEM simulation.
   */
  void
  read_dem();

  /**
   * @brief Write the CFD-DEM restart files.
   */
  void
  write_checkpoint() override;

  /**
   * @brief Read the CFD-DEM restart files.
   */
  void
  read_checkpoint() override;

  /**
   * @brief Execute the contact detection method.
   *
   * @param[in] counter The DEM iteration.
   */
  void
  check_contact_detection_method(unsigned int counter);

  /**
   * @brief Check if a load balancing is required according to the load
   * balancing method and perform it if necessary.
   */
  void
  load_balance();

  /**
   * @brief Add fluid-particle interaction force to the "force" container.
   */
  void
  add_fluid_particle_interaction_force();

  /**
   * @brief Add fluid-particle interaction torque to the "torque" container.
   */
  void
  add_fluid_particle_interaction_torque();

  /**
   * @brief Calculate particles-wall contact forces.
   */
  void
  particle_wall_contact_force();

  /**
   * @brief Post-processing as parallel VTU files.
   */
  void
  write_dem_output_results();

  /**
   * @brief Calculate statistics on the particles and report them to the
   * terminal. This function is notably used to monitor the time min, max and
   * total performed contact searches, and the instant min, max, avg and total
   * values of the linear velocity, angular velocity, linear kinetic energy and
   * angular kinetic.
   */
  void
  report_particle_statistics();

  /**
   * @brief Print the final summary of the particles.
   */
  void
  print_particles_summary();

  /**
   * @brief Post-process fluid dynamics after an iteration.
   */
  void
  postprocess_fd(bool first_iteration) override;

  /**
   * @brief Post-process cfd-dem after an iteration.
   */
  void
  postprocess_cfd_dem();

  /**
   * @brief Dynamic flow control calculation that take into account the void
   * fraction for the average velocity calculation.
   */
  void
  dynamic_flow_control() override;

  /**
   * @brief Execute the sorting of particle into subdomains and cells, and
   * reinitialize the containers dependent on the local particle ids.
   */
  void
  sort_particles_into_subdomains_and_cells();

  /**
   * @brief Execute a complete DEM iteration, including the particle-particle
   * and particle-wall contact force calculations, integration.
   */
  void
  dem_iterator(unsigned int counter);

  /**
   * @brief Check if a contact search has to be performed and execute it if so.
   */
  void
  dem_contact_build(unsigned int counter);


  unsigned int                               coupling_frequency;
  Tensor<1, 3>                               g;
  std::vector<Tensor<1, 3>>                  torque;
  std::vector<Tensor<1, 3>>                  force;
  std::vector<double>                        displacement;
  std::vector<double>                        MOI;
  double                                     neighborhood_threshold_squared;
  std::vector<std::shared_ptr<Distribution>> size_distribution_object_container;
  double                                     maximum_particle_diameter;
  double                                     smallest_contact_search_criterion;

  DEMContactManager<dim>           contact_manager;
  LagrangianLoadBalancing<dim>     load_balancing;
  ParticlePointLineForce<dim>      particle_point_line_contact_force_object;
  std::shared_ptr<Integrator<dim>> integrator_object;
  std::shared_ptr<Insertion<dim>>  insertion_object;
  std::shared_ptr<ParticleParticleContactForceBase<dim>>
    particle_particle_contact_force_object;
  std::shared_ptr<ParticleWallContactForce<dim>>
                                particle_wall_contact_force_object;
  Visualization<dim>            visualization_object;
  BoundaryCellsInformation<dim> boundary_cell_object;

  // Mesh and boundary information
  typename dem_data_structures<dim>::solid_surfaces_mesh_information
    solid_surfaces_mesh_info;
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

  // Object handling the sparse contacts
  AdaptiveSparseContacts<dim> sparse_contacts_object;

  // Counter of contact searches in a CFD iteration
  unsigned int contact_search_counter;

  // Total number of contact searches since the beginning of the simulation
  unsigned int contact_search_total_number;

  // Storage of statistics about time and contact lists
  statistics contact_list;

  DEM::DEMProperties<dim> properties_class;

  // Information for parallel grid processing
  PVDHandler grid_pvdhandler;
  PVDHandler particles_pvdhandler;

  DEMSolverParameters<dim> dem_parameters;
  double                   dem_time_step;
  const unsigned int       this_mpi_process;
  const unsigned int       n_mpi_processes;

  /// Post-processing variables to output total fluid volume and total particles
  /// volume
  TableHandler table_phase_volumes;

  DEMActionManager *dem_action_manager;
};
#endif
