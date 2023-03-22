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

/**
 * The general idea behind the algorithm :
 * It uses the granular temperature to determinate if particles in a cell are
 * mobile enough that the contact forces are worth computation. A cell having a
 * granular temperature under a value (default is 1e-4) will have an inactive
 * status which makes them rejected in the broad search step. This results in no
 * contact forces or velocity integration in future steps for those particles,
 * which make less computation cost for the simulation.
 *
 * Cells may be flagged as that called mobility status :
 * mobile (everything is calculated as if feature is not enabled),
 * active (particles with low motion next to mobile particles),
 * inactive (particles are not in a neighbor’s candidate list of cells around)
 *
 * There are some edge cases that need some attention. In these cases, particles
 * in cells can't be deactivated :
 *
 * The solid fraction of the cell is under a value (default = 40%) : particles
 * in the cell may have other forces that can move them and it is supposed that
 * with this fraction, there is not enough particles around to hold them in
 * their position.
 *
 * Cells having empty cell neighbors : it means that particles are next to a
 * floating wall/mesh and they may have a change of contact forces from it if
 * the wall disappears or is moving.
 *
 * The cell has cell neighbors which is flagged as mobile from the 3 criteria
 * mentioned above (granular temperature, solid fractions or floating walls) :
 * this is because motion is badly propagated to particles to the cell around
 * without this additional "layer" of mobile cells.
 *
 * The cell has cell neighbors which is flagged as mobile from the previous
 * criterion (additional layer of mobile cells) : this is again to allow the
 * propagation of motion, but only the particles in contact with particles
 * from the mobile cells are considered for the contact force calculation but
 * their position is not computed at the integration step. Those cells are
 * flagged as active cells. It works with the assignment and verification of the
 * mobility status at nodes to check the status of the neighboring cells.
 *
 */
template <int dim>
class DisableParticleContact
{
public:
  DisableParticleContact<dim>();

  /**
   * Mobility status flag used to identify the status at nodes and the status
   * of the cell:
   *
   * inactive (0)
   * The movement of particles in the cell is considered as negligible,
   * particles within this cell are not considered in the contact detection
   * (rejected at the broad search step), so no force calculation or
   * integration is applied
   *
   * active (1)
   * The movement of particles in the cell is considered as negligible, but
   * there's at least one neighbor cell that is flagged as mobile, meaning that
   * particles need to be in contact candidates lists at the broad search step,
   * particles directly in contact with the mobile cell are also considered in
   * force calculation, but none of the particles in the cell are integrated.
   *
   * mobile (2)
   * The movement of particles in the cell is significant or there is at least
   * one neighbor cell that is mobile by criteria (see the
   * identify_mobility_status() description), particles need to be in contact
   * candidates lists at the broad search step and particles are treated as
   * usual (force calculation and integration)
   *
   * empty (3)
   * This status is only used for the node-based mobility status identification,
   * no cells are flagged as empty, only node can by identify as empty. Without
   * this identification of the empty cells, we can't identify the cell that
   * have a empty neighbor cell, which is critical for simulationd using
   * floating walls or mesh
   */
  enum mobility_status : unsigned int
  {
    inactive = 0, // used for cells and nodes
    active   = 1, // used for cells and nodes
    mobile   = 2, // used for cells and nodes
    empty    = 3  // used for nodes only
  };

  /**
   * @brief Create or update a set of the active and ghost cells so that we don't
   * have to loop over all the cells in the triangulation for the granular
   * temperature and solid fraction calculation, and during the identification
   * of the mobility status. This set prevent 4 iteration steps over all the
   * cells  + the verification if the cell is locally owned, ghost or not.
   * This set is updated at every load balance step since cells are
   * redistributed among processors.
   *
   * @param background_dh The DoFHandler of the background grid
   */
  void
  update_active_ghost_cell_set(const DoFHandler<dim> &background_dh);

  /**
   * @brief Carries out the identification of the mobility status of each cell
   * through a node-based identification and check. Only the active and ghost
   * cells are processed.
   *
   * The following 4 checks (search loops) are done:
   *
   * 1. Check if the cell is empty (n_particle = 0), if so, nodes and cells are
   * flagged as empty mobility status (3)
   *
   * 2. Check if the cell is mobile by criteria (average granular temperature >
   * threshold, solid fraction < threshold or has at least one empty node from
   * previous check), if so, nodes are flagged and cells are stored with mobile
   * mobility status (2)
   *
   * 3. Check if the cell is mobile by neighbor (at least a node is flagged as
   * mobile from previous check), if so, cells are stored in map as mobile
   * status (2) and nodes that are not mobile are flagged as active (1)
   *
   * 4. Check if the cell is active (at least a node is flagged as active from
   * previous check), if so, cells are stored with active status in the map (1)
   *
   * The remaining cells are inactive (0)
   *
   * @param background_dh The dof handler of the background grid
   *
   * @param particle_handler The particle handler that contains all the particles
   *
   * @param mpi_communicator The MPI communicator
   */
  void
  identify_mobility_status(
    const DoFHandler<dim> &                background_dh,
    const Particles::ParticleHandler<dim> &particle_handler,
    const unsigned int                     n_active_cells,
    MPI_Comm                               mpi_communicator);

  /**
   * @brief Find the mobility status of a cell
   *
   * @param cell The iterator of the cell that needs mobility evaluation
   */
  inline unsigned int
  check_cell_mobility(
    const typename Triangulation<dim>::active_cell_iterator &cell) const
  {
    return cell_mobility_status.at(cell->active_cell_index());
  }

  /**
   * @brief Convert the map of mobility status to a vector of mobility status
   * because map can't be used as is in the pvd post-processing or any data out,
   * it needs to be converted to a vector of mobility status by active cell
   * index
   *
   * @param status The initiated vector for the conversion
   */
  void
  get_mobility_status_vector(std::vector<unsigned int> &status)
  {
    for (auto &cell_to_status : cell_mobility_status)
      {
        status[cell_to_status.first] = cell_to_status.second;
      }
  };

  void
  get_mobility_status_vector(Vector<float> &status)
  {
    for (auto &cell_to_status : cell_mobility_status)
      {
        status[cell_to_status.first] = cell_to_status.second;
      }
  };

  typename DEM::dem_data_structures<dim>::cell_index_int_map &
  get_mobility_status()
  {
    return cell_mobility_status;
  }

  void
  set_threshold_values(const double granular_temperature,
                       const double solid_fraction)
  {
    granular_temperature_threshold = granular_temperature;
    solid_fraction_threshold       = solid_fraction;
  }

private:
  /**
   * @brief Carries out the calculation of the granular temperature and solid
   * fraction approximation in each local cell. Those values are criteria for
   * cell mobility
   *
   * @param particle_handler The particle handler that contains all the particles
   *
   * @param cell_granular_temperature The empty vector of granular temperature
   *
   * @param cell_solid_fraction The empty vector of solid fraction
   */
  void
  calculate_granular_temperature_solid_fraction(
    const Particles::ParticleHandler<dim> &particle_handler,
    Vector<double> &                       cell_granular_temperature,
    Vector<double> &                       cell_solid_fraction);

  // Map of cell mobility status, the key is the active cell index and the value
  // is the mobility status
  typename DEM::dem_data_structures<dim>::cell_index_int_map
    cell_mobility_status;

  // Set of active and ghost cells, used to loop over only the active and ghost
  // cells without looping over all the cells in the triangulation many times
  std::set<typename DoFHandler<dim>::active_cell_iterator> active_ghost_cells;

  // Vector of mobility status at nodes, used to check the value at node to
  // determine the mobility status of the cell, this type of vector is used
  // to allow update values in parallel
  LinearAlgebra::distributed::Vector<int> mobility_at_nodes;

  // Threshold values for granular temperature and solid fraction
  double granular_temperature_threshold;
  double solid_fraction_threshold;
};


#endif // lethe_disable_particle_contact_h
