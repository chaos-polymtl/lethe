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


#ifndef lethe_load_balancing_h
#define lethe_load_balancing_h

#include <core/parameters_lagrangian.h>

#include <dem/data_containers.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 *
 */
template <int dim>
class LoadBalancing
{
public:
  LoadBalancing();

public:
  void
  set_parameters(
    const Parameters::Lagrangian::ModelParameters &model_parameters)
  {
    particle_weight = model_parameters.load_balance_particle_weight;
    inactive_load_balancing_factor =
      model_parameters.inactive_load_balancing_factor;
    active_load_balancing_factor =
      model_parameters.active_load_balancing_factor;
  }


  /**
   * @brief In order to consider the particles when repartitioning the triangulation
   * the algorithm needs to know three things:
   * 1. How much weight to assign to each cell (how many particles are in
   * there)
   * 2. How to pack the particles before shipping data around
   * 3. How to unpack the particles after repartitioning
   * Attach the correct functions to the signals inside
   * parallel::distributed::Triangulation, which will be called every time the
   * repartition() or refinement functions are called.
   * These connections only need to be created once, so we might as well
   * have set them up in the constructor of this class, but for the purpose
   * of this example we want to group the particle related instructions.
   *
   * @param triangulation The triangulation that will be repartitioned
   * @param particle_handler The particle handler that contains the particles
   */
  void
  connect_weight_signals(
    parallel::distributed::Triangulation<dim> &triangulation,
    const Particles::ParticleHandler<dim>     &particle_handler)
  {
    triangulation.signals.weight.connect(
      [](const typename Triangulation<dim>::cell_iterator &,
         const CellStatus) -> unsigned int { return 1000; });

    triangulation.signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->calculate_total_cell_weight(cell,
                                                 status,
                                                 particle_handler);
      });
  }

  void
  connect_mobility_status_weight_signals(
    parallel::distributed::Triangulation<dim> &triangulation,
    const Particles::ParticleHandler<dim>     &particle_handler,
    const typename DEM::dem_data_structures<dim>::cell_index_int_map
      mobility_status)
  {
    // Clear and connect a new cell weight function
    triangulation.signals.weight.disconnect_all_slots();

    triangulation.signals.weight.connect(
      [](const typename Triangulation<dim>::cell_iterator &,
         const CellStatus) -> unsigned int { return 1000; });

    triangulation.signals.weight.connect(
      [&](const typename parallel::distributed::Triangulation<
            dim>::cell_iterator &cell,
          const CellStatus       status) -> unsigned int {
        return this->calculate_total_cell_weight_with_mobility_status(
          cell, status, particle_handler, mobility_status);
      });
  }

private:
  /**
   * @brief Indicates to the triangulation how much
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
   *
   * @param cell The cell for which the load is calculated
   * @param status The status of the cell related to the coarsening level
   */
  unsigned int
  calculate_total_cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                                          &cell,
    const CellStatus                       status,
    const Particles::ParticleHandler<dim> &particle_handler) const;

  /**
   * Similar to the cell_weight() function, this function is used when the cell
   * weight is adapted to the mobility status. For instance, if the
   * cell is inactive, its computational load will be significantly lower than
   * if it is a mobile cell since there is no force calculation and no velocity
   * integration for the particles that lie within it. The weight of the cells
   * must thus be adapted to the status of the cell.
   *
   * cell load = cell weight + load balancing factor * n particles * particle
   * weight
   *
   * @param cell The cell for which the load is calculated
   * @param status The status of the cell related to the coarsening level
   * @param mobility_status The mobility status of the cell
   */

  unsigned int
  calculate_total_cell_weight_with_mobility_status(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                                          &cell,
    const CellStatus                       status,
    const Particles::ParticleHandler<dim> &particle_handler,
    const typename DEM::dem_data_structures<dim>::cell_index_int_map
      mobility_status) const;

  // Weight of the cell and particles in the cells
  const unsigned int cell_weight = 1000;
  unsigned int       particle_weight;

  // Factors with dynamic load balancing with adaptive sparse contacts
  double inactive_load_balancing_factor;
  double active_load_balancing_factor;
};


#endif
