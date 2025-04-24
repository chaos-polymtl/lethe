// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_ray_tracing_h
#define lethe_particle_ray_tracing_h

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/particle_ray_tracing_parameters.h>

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
  ParticleRayTracing(ParticleRayTracingParameters<dim> parameters);

  void
  find_locally_own_cells_with_particles(
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_with_particle_and_neighbours);


  /**
   * @brief The parameters of the DEM simulation.
   */
  ParticleRayTracingParameters<dim> parameters;

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
   * @brief The distributed triangulation used for the DEM simulation.
   */
  parallel::distributed::Triangulation<dim> triangulation;

  /**
   * @brief The polynomial mapping. The mapping is first order in DEM
   * simulations.
   */
  MappingQGeneric<dim> mapping;

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
  const Tensor<1, dim> photon_displacement_tensor;

  // Container that shows the local/ghost neighbor cells of all local cells in
  // the triangulation
  // typename dem_data_structures<dim>::cells_total_neighbor_list
  //  total_neighbor_list;

  typename dem_data_structures<dim>::cell_set
    local_and_ghost_cells_with_particles_and_neighbors;
};


#endif // lethe_particle_ray_tracing_h
