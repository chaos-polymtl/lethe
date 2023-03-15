/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
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

#include <core/dem_properties.h>

#include <dem/data_containers.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_vector.templates.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>

using namespace dealii;

#ifndef lethe_disable_particle_contact_h
#  define lethe_disable_particle_contact_h

// Special template instance for this class.
// Unsigned integer would have been a better choice, but is not working
// with the Vector class (error because of calling of abs() function)
template class LinearAlgebra::distributed::Vector<int>;

template <int dim>
class DisableParticleContact
{
public:
  DisableParticleContact<dim>();

  /**
   * Mobility status flag used to identify the status at nodes and the status
   * of the cell:
   *
   * inactive:
   * the movement of particles is negligible in the cell, particles within this
   * cell are not considered in the contact detection and force calculation
   * procedure.
   * active:
   * movement of particles is negligible, but there's at least one neighbor cell
   * that is flagged as mobile, meaning that particles need to be in contact
   * candidates lists
   * mobile:
   * movement of particles is significant or there is at least one neighbor cell
   * that is mobile by criteria, particles need to be in contact candidates
   * lists
   * empty:
   * cell is empty, only useful for mobility at nodes
   */
  // TODO : remove n_mobility_status
  enum mobility_status : unsigned int
  {
    inactive,
    active,
    mobile,
    empty,
    n_mobility_status
  };

  /**
   * Carries out the identification of the mobility status of each cell through
   * processing at nodes.
   *
   * The following 4 checks (search loops) are done:
   * 1. Check if the cell is empty (n_particle = 0), if so, nodes and cells are
   * flagged as empty mobility status
   * 2. Check if the cell is mobile by criterion (n_particle > 0, average
   * granular temperature > threshold, solid fraction < threshold, has at least
   * one empty node from previous check), if so, nodes and cells are flagged as
   * mobile mobility status
   * 3. Check if the cell is mobile by neighbor (at least a node is flagged as
   * mobile from previous check), if so, cells are flagged as mobile status and
   * nodes that are not mobile are flagged as active
   * 4. Check if the cell is active (at least a node is flagged as active from
   * previous check), if so, cells are flagged as active status
   *
   * The remaining cells are flagged as inactive by default (vector is
   * initialized with 0)
   *
   * @param background_dh The dof handler of the background grid
   * @param particle_handler The particle handler that contains all the particles
   * @param mpi_communicator The MPI communicator
   */
  void
  identify_mobility_status(
    const DoFHandler<dim> &                background_dh,
    const Particles::ParticleHandler<dim> &particle_handler,
    MPI_Comm                               mpi_communicator);

  /**
   * Carries out the calculation of the granular temperature and solid
   * fraction approximation in each local cell. Those values are a criteria
   * for cell mobility
   *
   * @param triangulation The triangulation
   * @param particle_handler The particle handler that contains all the particles
   */
  void
  calculate_cell_granular_temperature(
    const Particles::ParticleHandler<dim> &particle_handler,
    const unsigned int                     n_active_cells);

  void
  update_active_ghost_cell_set(const DoFHandler<dim> &background_dh);

  /**
   * Find the mobility status of a cell
   *
   * @param cell The iterator of the cell that needs mobility evaluation
   */
  inline unsigned int
  check_cell_mobility(
    const typename Triangulation<dim>::active_cell_iterator &cell) const
  {
    // Look for the mobility status from map, if not found, return inactive
    // status since no inactive status is stored in the map
    auto it = cell_mobility_status_map.find(cell->active_cell_index());
    if (it != cell_mobility_status_map.end())
      return it->second;
    else
      return mobility_status::inactive;
  }

  /**
   * Convert the vector of mobility status set to a vector of mobility status
   * with cell id order
   * status_to_cell is a vector of 3 sets of cell iterator (inactive, active,
   * mobile), it can't be used as is in the pvd post-processing or any data out,
   * it needs to be converted to a vector of mobility status by active cell
   * index
   *
   * @param status The initiated vector for the conversion
   */
  void
  get_mobility_status_vector(std::vector<unsigned int> &status)
  {
    for (auto &cell_to_status : cell_mobility_status_map)
      {
        status[cell_to_status.first] = cell_to_status.second;
      }
  };

  void
  get_mobility_status_vector(Vector<float> &status)
  {
    for (auto &cell_to_status : cell_mobility_status_map)
      {
        status[cell_to_status.first] = cell_to_status.second;
      }
  };

  std::unordered_map<types::global_cell_index, unsigned int> &
  get_mobility_status_map()
  {
    return cell_mobility_status_map;
  }

  LinearAlgebra::distributed::Vector<int>
  get_mobility_at_nodes()
  {
    return mobility_at_nodes;
  }

  unsigned int
  get_n_mobility_status()
  {
    return mobility_status::n_mobility_status;
  }

  void
  set_threshold_values(const double granular_temperature,
                       const double solid_fraction)
  {
    granular_temperature_threshold = granular_temperature;
    solid_fraction_threshold       = solid_fraction;
  }

private:
  void
  check_granular_temperature(
    const typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &                                    cells_neighbor_list,
    const Particles::ParticleHandler<dim> &particle_handler);

  std::unordered_map<types::global_cell_index, unsigned int>
                 cell_mobility_status_map;
  Vector<double> solid_fractions;

  Vector<double> granular_temperature_average;
  std::set<typename DoFHandler<dim>::active_cell_iterator> active_ghost_cells;
  LinearAlgebra::distributed::Vector<int>                  mobility_at_nodes;

  double granular_temperature_threshold;
  double solid_fraction_threshold;
};


#endif // lethe_disable_particle_contact_h
