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
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_contact_detection_step.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/periodic_boundaries_manipulator.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/gls_vans.h>
#include <fem-dem/postprocessing_cfd_dem.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

template <int dim>
class CFDDEMSolver : public GLSVANSSolver<dim>
{
  using FuncPtrType = bool (CFDDEMSolver<dim>::*)(const unsigned int &counter);
  FuncPtrType check_contact_search_step;

public:
  CFDDEMSolver(CFDDEMSimulationParameters<dim> &nsparam);

  ~CFDDEMSolver();

  virtual void
  solve() override;

  /**
   * @brief Manages the call to the load balancing. Returns true if
   * load balancing is performed
   *
   */
  void
  load_balance();

protected:
private:
  /**
   * @brief Carries out the DEM calculations in the DEM_CFD solver. Particle-particle and particle-wall contact force calculations, integration and update_ghost
   */
  void
  dem_iterator(unsigned int counter);

  /**
   * @brief Carries out the particle-particle and particle-wall contact searches, sort_particles_into_subdomains_and_cells and exchange_ghost
   */
  void
  dem_contact_build(unsigned int counter);

  /**
   * @brief Sets up the various parameters related to the DEM contacts
   */
  void
  dem_setup_contact_parameters();

  /**
   * @brief Carries out the initialization of DEM parameters
   */
  void
  initialize_dem_parameters();

  /**
   * @brief write DEM_output_results
   * Post-processing as parallel VTU files
   */
  void
  write_DEM_output_results();

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
   * @brief Updates moment of inertia container after sorting particles
   * into subdomains
   *
   */
  void
  update_moment_of_inertia(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    std::vector<double>                     &MOI);

  /**
   * Sets the chosen integration method in the parameter handler file
   *
   * @return A pointer to the integration object
   */
  std::shared_ptr<Integrator<dim>>
  set_integrator_type();

  /**
   * Adds fluid-particle interaction force to the "force" container
   *
   */
  void
  add_fluid_particle_interaction_force();

  /**
   * Adds fluid-particle interaction torque to the "torque" container
   *
   */
  void
  add_fluid_particle_interaction_torque();

  /**
   * Sets the chosen particle-particle contact force model in the parameter
   * handler file
   *
   * @return A pointer to the particle-particle contact force object
   */
  std::shared_ptr<ParticleParticleContactForceBase<dim>>
  set_particle_particle_contact_force();

  /**
   * Sets the chosen particle-wall contact force model in the parameter handler
   * file
   *
   * @return A pointer to the particle-wall contact force object
   */
  std::shared_ptr<ParticleWallContactForce<dim>>
  set_particle_wall_contact_force();

  void
  read_dem();

  void
  write_checkpoint() override;

  void
  read_checkpoint() override;

  void
  print_particles_summary();

  /**
   * @brief dem_post_process_results
   */
  void
  dem_post_process_results();

  /**
   * @brief postprocess
   * Post-process fluid dynamics after an iteration
   */
  void
  postprocess_fd(bool first_iteration) override;

  /**
   * @brief postprocess for cfd-dem
   * Post-process cfd-dem after an iteration
   */
  void
  postprocess_cfd_dem();

  /**
   * @brief dynamic_flow_control
   * Dynamic flow control calculation that take into account the void fraction
   * for the average velocity calculation
   */
  void
  dynamic_flow_control() override;

  /**
   * @brief Checks all the conditions that require a contact search step. The
   * check of conditions is done in order of suspected frequency occurrence.
   *
   * @param counter The counter of DEM iterations in a CFD iteration.
   */
  inline bool
  contact_search_step(const unsigned int counter)
  {
    if (contact_detection_step)
      {
        // Contact search step according to the contact detection method
        return true;
      }
    else if (counter == (coupling_frequency - 1))
      {
        // Execute the contact search at the last DEM iteration of the CFD time
        return true;
      }
    else if (has_sparse_contacts && (counter == 1))
      {
        // First mobility status identification of the CFD time step (from the
        // velocity computed at the first DEM time step (counter = 0) of the CFD
        // time step) The contact search is executed to make sure the mobility
        // status of cell match the particles that are in.
        return true;
      }
    else if (load_balance_step)
      {
        // Needs to update contacts since particles/cells may have been
        // distributed to a different subdomain
        return true;
      }
    else if ((this->simulation_control->is_at_start() && (counter == 0)) ||
             checkpoint_step)
      {
        // First contact search of the simulation
        return true;
      }
    else if (has_periodic_boundaries && particle_displaced_in_pbc)
      {
        // Particles have been displaced in periodic boundaries
        return true;
      }
    else
      {
        return false;
      }
  }

  unsigned int                               coupling_frequency;
  bool                                       contact_detection_step;
  bool                                       checkpoint_step;
  bool                                       load_balance_step;
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
  Tensor<1, dim>                     periodic_offset;
  bool                               has_periodic_boundaries;
  bool                               particle_displaced_in_pbc;

  // Object handling the sparse contacts
  AdaptiveSparseContacts<dim> sparse_contacts_object;

  // Flag to indicate if sparse contacts are used
  bool has_sparse_contacts;

  // Counter of contact searches in a CFD iteration
  unsigned int contact_search_counter;

  // Total number of contact searches since the beginning of the simulation
  unsigned int contact_search_total_number;

  // Storage of statistics about time and contact lists
  statistics contact_list;
  statistics simulation_time;

  DEM::DEMProperties<dim> properties_class;

  // Information for parallel grid processing
  PVDHandler grid_pvdhandler;
  PVDHandler particles_pvdhandler;

  DEMSolverParameters<dim>      dem_parameters;
  double                        dem_time_step;
  const unsigned int            this_mpi_process;
  const unsigned int            n_mpi_processes;
  LagrangianPostProcessing<dim> dem_post_processing_object;

  /// Post-processing variables to output total fluid volume and total particles
  /// volume
  TableHandler table_phase_volumes;
};
#endif
