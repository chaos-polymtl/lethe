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

#include <dem/dem.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_contact_detection_step.h>
#include <fem-dem/cfd_dem_simulation_parameters.h>
#include <fem-dem/gls_vans.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

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
  std::shared_ptr<PPContactForce<dim>>
  set_pp_contact_force();

  /**
   * Sets the chosen particle-wall contact force model in the parameter handler
   * file
   *
   * @return A pointer to the particle-wall contact force object
   */
  std::shared_ptr<PWContactForce<dim>>
  set_pw_contact_force();

  void
  read_dem();

  void
  write_checkpoint() override;

  void
  read_checkpoint() override;

  void
  print_particles_summary();

  /**
   * @brief postprocess
   * Post-process fluid dynamics after an iteration
   */
  void
  postprocess_fd(bool first_iteration) override;


  unsigned int                coupling_frequency;
  bool                        contact_detection_step;
  bool                        checkpoint_step;
  bool                        load_balance_step;
  std::vector<Tensor<1, dim>> momentum;
  std::vector<Tensor<1, dim>> force;
  std::vector<double>         displacement;
  std::vector<double>         MOI;
  double                      neighborhood_threshold_squared;
  double                      maximum_particle_diameter;
  double                      standard_deviation_multiplier;
  double                      smallest_contact_search_criterion;
  double                      triangulation_cell_diameter;

  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cells_local_neighbor_list;
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cells_ghost_neighbor_list;
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    local_contact_pair_candidates;
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    ghost_contact_pair_candidates;
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
    local_adjacent_particles;
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
    ghost_adjacent_particles;
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
    pw_pairs_in_contact;
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
    pfw_pairs_in_contact;
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  unsigned int,
                                  unsigned int>>>
    pw_contact_candidates;
  std::unordered_map<types::particle_index,
                     std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
    particle_point_contact_candidates;
  std::unordered_map<
    types::particle_index,
    std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
    particle_line_contact_candidates;
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<dim>>
    particle_points_in_contact, particle_lines_in_contact;
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>>
    pfw_contact_candidates;
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    particle_container;
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    ghost_particle_container;
  std::map<unsigned int, std::map<unsigned int, Tensor<1, dim>>>
    forces_boundary_information;
  std::map<unsigned int, std::map<unsigned int, Tensor<1, dim>>>
    torques_boundary_information;


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
  BoundaryCellsInformation<dim>        boundary_cell_object;
  FindCellNeighbors<dim>               cell_neighbors_object;

  DEM::DEMProperties<dim> properties_class;

  // Information for parallel grid processing
  DoFHandler<dim> background_dh;
  PVDHandler      grid_pvdhandler;
  PVDHandler      particles_pvdhandler;

  DEMSolverParameters<dim> dem_parameters;
  double                   dem_time_step;
  const unsigned int       this_mpi_process =
    Utilities::MPI::this_mpi_process(this->mpi_communicator);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(this->mpi_communicator);

  Triangulation<dim> tria;
};
#endif
