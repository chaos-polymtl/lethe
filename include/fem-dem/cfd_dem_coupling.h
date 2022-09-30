/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_contact_detection_step.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/periodic_boundaries_manipulator.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/gls_vans.h>

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
  FuncPtrType check_load_balance_step;

public:
  CFDDEMSolver(CFDDEMSimulationParameters<dim> &nsparam);

  ~CFDDEMSolver();

  virtual void
  solve() override;

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
   * @brief Updates moment of inertia container after sorting particles
   * into subdomains
   *
   */
  void
  update_moment_of_inertia(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    std::vector<double> &                    MOI);

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
   * Sets the chosen particle-particle contact force model in the parameter
   * handler file
   *
   * @return A pointer to the particle-particle contact force object
   */
  std::shared_ptr<ParticleParticleContactForce<dim>>
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

  void
  post_processing() override;

  unsigned int              coupling_frequency;
  bool                      contact_detection_step;
  bool                      checkpoint_step;
  bool                      load_balance_step;
  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       displacement;
  std::vector<double>       MOI;
  double                    neighborhood_threshold_squared;
  double                    maximum_particle_diameter;
  double                    standard_deviation_multiplier;
  double                    smallest_contact_search_criterion;
  double                    triangulation_cell_diameter;

  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_candidates;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_adjacent_particles;
  typename dem_data_structures<dim>::particle_wall_candidates
    particle_wall_candidates;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_in_contact;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_floating_wall_in_contact;
  typename dem_data_structures<dim>::particle_point_candidates
    particle_point_candidates;
  typename dem_data_structures<dim>::particle_line_candidates
    particle_line_candidates;
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_points_in_contact;
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_lines_in_contact;
  typename dem_data_structures<dim>::particle_floating_wall_candidates
    particle_floating_wall_candidates;
  typename dem_data_structures<dim>::particle_floating_mesh_candidates
    particle_floating_mesh_candidates;
  typename dem_data_structures<dim>::particle_floating_mesh_in_contact
    particle_floating_mesh_in_contact;
  typename dem_data_structures<dim>::particle_index_iterator_map
    particle_container;
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;

  ParticleParticleBroadSearch<dim>   particle_particle_broad_search_object;
  ParticleParticleFineSearch<dim>    particle_particle_fine_search_object;
  ParticleWallBroadSearch<dim>       particle_wall_broad_search_object;
  ParticlePointLineBroadSearch<dim>  particle_point_line_broad_search_object;
  ParticleWallFineSearch<dim>        particle_wall_fine_search_object;
  ParticlePointLineFineSearch<dim>   particle_point_line_fine_search_object;
  ParticlePointLineForce<dim>        particle_point_line_contact_force_object;
  PeriodicBoundariesManipulator<dim> periodic_boundaries_object;
  std::shared_ptr<Integrator<dim>>   integrator_object;
  std::shared_ptr<Insertion<dim>>    insertion_object;
  std::shared_ptr<ParticleParticleContactForce<dim>>
    particle_particle_contact_force_object;
  std::shared_ptr<ParticleWallContactForce<dim>>
                                particle_wall_contact_force_object;
  Visualization<dim>            visualization_object;
  BoundaryCellsInformation<dim> boundary_cell_object;
  FindCellNeighbors<dim>        cell_neighbors_object;

  DEM::DEMProperties<dim> properties_class;

  // Information for parallel grid processing
  DoFHandler<dim> background_dh;
  PVDHandler      grid_pvdhandler;
  PVDHandler      particles_pvdhandler;

  DEMSolverParameters<dim>      dem_parameters;
  double                        dem_time_step;
  const unsigned int            this_mpi_process;
  const unsigned int            n_mpi_processes;
  LagrangianPostProcessing<dim> dem_post_processing_object;

  Triangulation<dim> tria;
};
#endif
