// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ray_tracing_h
#define lethe_ray_tracing_h

#include <core/dem_properties.h>
#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/insertion.h>
#include <dem/load_balancing.h>
#include <dem/ray_tracing_solver_parameters.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>

template <int dim, typename PropertiesIndex>
class RayTracingSolver
{
public:
  RayTracingSolver(RayTracingSolverParameters<dim> parameters);


private:
  /**
   * @brief Set the parameters for the DEM simulation
   */
  void
  setup_parameters();

  /**
   * @brief Set the insertion method.
   *
   * @return The pointer to the particle insertion object.
   */
  std::shared_ptr<Insertion<dim, PropertiesIndex>>
  set_particle_insertion_type();

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
   * @brief Insert particles and photons at the beginning of a ray tracing
   * simulation.
   */
  void
  insert_particles_and_photons();

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
   * @brief Execute the sorting of particle into subdomains and cells, and
   * reinitialize the containers dependent on the local particle ids.
   */
  void
  sort_particles_into_subdomains_and_cells();


  /**
   * @brief Calls all the necessary functions to set parameters, solve the intersection
   * points between the photons and the particles, and finish the simulation.
   */
  void
  solve();


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
   * @brief The parameters of the Ray Tracing simulation.
   */
  RayTracingSolverParameters<dim> parameters;

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
   * @brief The timer that keeps track of the time spent in some functions.
   * Currently theses functions are: load balancing and VTU output.
   */
  TimerOutput computing_timer;

  /**
   * @brief The photon handler that manages the photons.
   */
  Particles::ParticleHandler<dim, dim> photon_handler;

  /**
   * @brief The load balancing handler.
   */
  LagrangianLoadBalancing<dim, PropertiesIndex> load_balancing;

  /**
   * @brief The simulation control (DEM Transient).
   */
  std::shared_ptr<SimulationControl> simulation_control;


  /**
   * @brief The particle insertion object.
   */
  std::shared_ptr<Insertion<dim, PropertiesIndex>> particle_insertion_object;
};
#endif // lethe_ray_tracing_h
