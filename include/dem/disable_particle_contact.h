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
   * of the cell
   * inactive: movement of particles is negligible in the cell, particles won't
   *           be in the contact candidates list
   * active:   movement of particles is negligible, but there's at least one
   * neighbor cell that is flag mobile, meaning that particles need to be in
   * contact candidates lists mobile:   movement of particles is significant or
   * there is at least one neighbor cell that is mobile by criteria, particles
   * need to be in contact candidates lists empty:    cell is empty, only useful
   * for mobility at nodes
   */
  enum mobility_status
  {
    inactive,
    active,
    mobile,
    empty,
    n_mobility_status
  };

  /**
   * Carries out the identification of the mobility status of each cell thought
   * processing at nodes.
   *
   * 4 checks (search loops) are done is order:
   * 1. Check if the cell is empty (n_particle = 0), if so, nodes and cells are
   * flagged as empty mobility status
   * 2. Check if the cell is mobile by criterion (n_particle > 0, average
   * granular temperature > threshold, solid fraction < limit, as one empty node
   * from previous check), if so, nodes and cells are flagged as mobile mobility
   * status
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
  calculate_average_granular_temperature(
    const DoFHandler<dim> &                background_dh,
    const Particles::ParticleHandler<dim> &particle_handler);

  /**
   * Find the mobility status of a cell
   *
   * @param cell The iterator of the cell that needs mobility evaluation
   */
  inline unsigned int
  check_cell_mobility(
    const typename Triangulation<dim>::active_cell_iterator &cell) const
  {
    unsigned int status;
    // Check if the cell is in the mobile status set
    status = mobility_status::mobile;
    auto cell_iterator_from_mobile_container =
      status_to_cell[status].find(cell);

    if (cell_iterator_from_mobile_container != status_to_cell[status].end())
      {
        return status;
      }

    // Check if the cell is in the active status set
    status = mobility_status::active;
    auto cell_iterator_from_active_container =
      status_to_cell[status].find(cell);

    if (cell_iterator_from_active_container != status_to_cell[status].end())
      {
        return status;
      }

    // If not mobile or active, cell is inactive
    return mobility_status::inactive;
  }

  /**
   * Convert the vector of mobility status set to a vector of mobility status
   * with cell id order
   * status_to_cell is a vector of 3 sets of cell iterator (inactive, active,
   * mobile), it can't be used as is in the pvd post-processing, it needs to be
   * converted to a vector of mobility status by active cell index
   *
   * @param n_active_cells The number of active cells in the triangulation
   */
  std::vector<unsigned int>
  get_mobility_status_vector(unsigned int n_active_cells)
  {
    std::vector<unsigned int> status(n_active_cells);

    // Loop over all set with different mobility status (except empty)
    for (unsigned int i_set = 0; i_set < mobility_status::n_mobility_status - 1;
         i_set++)
      {
        for (auto &cell : status_to_cell[i_set])
          status[cell->active_cell_index()] = i_set;
      }

    return status;
  };

  std::vector<typename DEM::dem_data_structures<dim>::cell_set>
  get_mobility_status()
  {
    return status_to_cell;
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
  set_limit_value(const double granular_temperature,
                  const double solid_fraction)
  {
    granular_temperature_limit = granular_temperature;
    solid_fraction_limit       = solid_fraction;
  }

private:
  void
  check_granular_temperature(
    const typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &                                    cells_neighbor_list,
    const Particles::ParticleHandler<dim> &particle_handler);

  std::vector<typename DEM::dem_data_structures<dim>::cell_set> status_to_cell;
  Vector<double>                                                solid_fractions;

  Vector<double> granular_temperature_average;
  std::set<typename DoFHandler<dim>::active_cell_iterator> active_ghost_cells;
  LinearAlgebra::distributed::Vector<int>                  mobility_at_nodes;

  double granular_temperature_limit;
  double solid_fraction_limit;
};



#endif // lethe_disable_particle_contact_h
