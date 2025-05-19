// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_ray_tracing_h
#define lethe_particle_ray_tracing_h

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>
#include <dem/load_balancing.h>

#include <deal.II/particles/particle_handler.h>

using namespace DEM;

/**
 * @brief
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class ParticleRayTracing
{
public:
  ParticleRayTracing(DEMSolverParameters<dim> dem_parameters);

  /**
   * @brief Set the parameters for the ray tracing simulation
   */
  void
  setup_parameters();

  /*
   * @brief Find locals cells containing particles or with neighboring cells, local
   * and ghost, containing particles.
   *
   * @param[in,out] local_cells_near_particles Contains local cells with
   * particles or with particles in their neighbors.
   */
  void
  find_locally_own_cells_near_particles(
    typename dem_data_structures<dim>::cell_set local_cells_near_particles);


  /*
   * @brief Insert the photon at the appropriate location.
   */
  void
  insert_photon();

  /**
   * @brief Calls all the necessary functions to set parameters, solve the
   * simulation, and finish the simulation.
   */
  void
  solve();


  /**
   * @brief The parameters of the DEM simulation.
   */
  DEMSolverParameters<dim> dem_parameters;

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
   * @brief The parameters of the DEM simulation.
   */
  DEMSolverParameters<dim> parameters;

  /**
   * @brief The output stream used for the parallel simulation, only used by
   * the process 0.
   */
  ConditionalOStream pcout;

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
   * @brief The background degree of freedom handler uses for parallel grid
   * processing.
   */
  DoFHandler<dim> background_dh;

  /**
   * @brief The particle handler that manages the photons.
   */
  Particles::ParticleHandler<dim, dim> photon_handler;

  /**
   * @brief The particle handler that manages the particles.
   */
  Particles::ParticleHandler<dim, dim> particle_handler;

  /**
   * @brief The displacement tensor of a photons during a pseudo time-step.
   */
  Tensor<1, dim> photon_displacement_tensor;

  /*
   * @brief Insertion object used to insert particle.
   */
  std::shared_ptr<Insertion<dim, PropertiesIndex>> insertion_object;

  /**
   * @brief The load balancing handler.
   */
  LagrangianLoadBalancing<dim, PropertiesIndex> load_balancing;
};


#endif // lethe_particle_ray_tracing_h
