// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Unit test for the neighborhood threshold of the particle-particle
 * broad search.
 *
 * Four particles are inserted in two neighboring cells, with pair distances
 * of 0.25, 0.5 and 0.75. The broad search is called without threshold
 * (infinite threshold, all the pairs are stored), then with thresholds that
 * discard some of the pairs. A threshold equal to the
 * distance of a pair keeps that pair, as the fine search does for a pair that
 * is already in contact.
 *
 * The periodic broad search is then checked on a domain that is periodic in
 * x. The distance of the periodic pairs is computed with the nearest periodic
 * image, so pairs that are far apart in real space but close through the
 * periodic boundary are kept.
 */

// Deal.II
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/periodic_boundary.h>

#include <dem/data_containers.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/periodic_boundaries_manipulator.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <algorithm>
#include <array>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace dealii;

/**
 * @brief Insert a particle at a given location in the particle handler.
 *
 * @param[in] triangulation Triangulation in which the particle is located.
 * @param[in] mapping Mapping of the triangulation.
 * @param[in] location Location of the particle.
 * @param[in] id Id of the particle.
 * @param[in,out] particle_handler Particle handler in which the particle is
 * inserted.
 */
template <int dim>
void
insert_particle(const parallel::distributed::Triangulation<dim> &triangulation,
                const MappingQ1<dim>                            &mapping,
                const Point<dim>                                &location,
                const types::particle_index                      id,
                Particles::ParticleHandler<dim> &particle_handler)
{
  const auto [cell, reference_location] =
    GridTools::find_active_cell_around_point(mapping, triangulation, location);
  Particles::Particle<dim> particle(location, reference_location, id);
  particle_handler.insert_particle(particle, cell);
}

/**
 * @brief Print the candidate pairs, sorted, each pair being written with the
 * lowest id first.
 *
 * @param[in] label Label printed before the pairs.
 * @param[in] contact_pair_candidates Candidates of the broad search.
 */
template <int dim>
void
print_candidates(
  const std::string &label,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &contact_pair_candidates)
{
  std::vector<std::pair<types::particle_index, types::particle_index>> pairs;
  for (const auto &[particle_one_id, particle_two_ids] :
       contact_pair_candidates)
    for (const types::particle_index particle_two_id : particle_two_ids)
      pairs.emplace_back(std::min(particle_one_id, particle_two_id),
                         std::max(particle_one_id, particle_two_id));
  std::sort(pairs.begin(), pairs.end());

  deallog << label << ":";
  for (const auto &[particle_one_id, particle_two_id] : pairs)
    deallog << " (" << particle_one_id << ", " << particle_two_id << ")";
  deallog << std::endl;
}

/**
 * @brief Check the threshold of the local particle-particle broad search.
 */
template <int dim>
void
test_local()
{
  // Cube [-1, 1]^dim refined twice: cells of size 0.5
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(2);

  MappingQ1<dim>                  mapping;
  Particles::ParticleHandler<dim> particle_handler(triangulation, mapping);

  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;
  find_cell_neighbors<dim>(triangulation,
                           cells_local_neighbor_list,
                           cells_ghost_neighbor_list);

  // Particles 0 and 3 are in the cell [-0.5, 0] in x, particles 1 and 2 in
  // the neighboring cell [0, 0.5]. The coordinates are exactly representable,
  // so the distances 0.25 (pairs 0-1, 0-3, 1-2), 0.5 (pairs 0-2, 1-3) and 0.75
  // (pair 2-3) are exact.
  const std::array<double, 4> x_locations = {-0.125, 0.125, 0.375, -0.375};
  for (unsigned int i = 0; i < x_locations.size(); ++i)
    {
      Point<dim> location;
      location[0] = x_locations[i];
      for (int d = 1; d < dim; ++d)
        location[d] = 0.25;
      insert_particle<dim>(
        triangulation, mapping, location, i, particle_handler);
    }

  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_candidates;
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_candidates;

  const std::vector<std::pair<std::string, double>> thresholds = {
    {"no threshold", std::numeric_limits<double>::infinity()},
    {"threshold 0.6", 0.6 * 0.6},
    {"threshold 0.3", 0.3 * 0.3},
    {"threshold 0.25", 0.25 * 0.25}};

  for (const auto &[label, neighborhood_threshold_squared] : thresholds)
    {
      find_particle_particle_contact_pairs<dim>(particle_handler,
                                                cells_local_neighbor_list,
                                                cells_ghost_neighbor_list,
                                                local_contact_pair_candidates,
                                                ghost_contact_pair_candidates,
                                                neighborhood_threshold_squared);
      print_candidates<dim>("dim = " + std::to_string(dim) + ", local, " +
                              label,
                            local_contact_pair_candidates);
    }
}

/**
 * @brief Check the threshold of the periodic particle-particle broad search.
 */
template <int dim>
void
test_periodic()
{
  // Unit cube periodic in x, refined twice: cells of size 0.25. The faces of
  // the periodic boundary in x have the ids 0 and 1.
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, 0, 1, true);

  std::vector<GridTools::PeriodicFacePair<
    typename parallel::distributed::Triangulation<dim>::cell_iterator>>
    matched_pairs;
  GridTools::collect_periodic_faces(triangulation, 0, 1, 0, matched_pairs);
  triangulation.add_periodicity(matched_pairs);
  triangulation.refine_global(2);

  Parameters::PeriodicBoundaries periodic_boundaries;
  periodic_boundaries[0] = {.neighbor_id = 1, .direction = 0};

  PeriodicBoundariesManipulator<dim> periodic_boundaries_object;
  periodic_boundaries_object.set_periodic_boundaries_information(
    periodic_boundaries);
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    periodic_boundaries_cells_information;
  periodic_boundaries_object.map_periodic_cells(
    triangulation, periodic_boundaries_cells_information);

  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    cells_local_periodic_neighbor_list;
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_periodic_neighbor_list;
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_local_periodic_neighbor_list;
  find_cell_periodic_neighbors<dim>(triangulation,
                                    periodic_boundaries_cells_information,
                                    cells_local_periodic_neighbor_list,
                                    cells_ghost_periodic_neighbor_list,
                                    cells_ghost_local_periodic_neighbor_list);

  MappingQ1<dim>                  mapping;
  Particles::ParticleHandler<dim> particle_handler(triangulation, mapping);

  // Particles 0 and 3 are in the cell [0, 0.25] in x, at the periodic
  // boundary 0, and particles 1 and 2 in the cell [0.75, 1], at the periodic
  // boundary 1. Through the periodic boundary, the distances are 0.125 (pair
  // 0-1), 0.25 (pairs 0-2 and 1-3) and 0.375 (pair 2-3), while the real
  // distances are all larger than 0.6.
  const std::array<double, 4> x_locations = {0.0625, 0.9375, 0.8125, 0.1875};
  for (unsigned int i = 0; i < x_locations.size(); ++i)
    {
      Point<dim> location;
      location[0] = x_locations[i];
      for (int d = 1; d < dim; ++d)
        location[d] = 0.375;
      insert_particle<dim>(
        triangulation, mapping, location, i, particle_handler);
    }

  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_periodic_candidates;
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_periodic_candidates;
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    ghost_local_contact_pair_periodic_candidates;

  const std::vector<std::pair<std::string, double>> thresholds = {
    {"no threshold", std::numeric_limits<double>::infinity()},
    {"threshold 0.25", 0.25 * 0.25},
    {"threshold 0.2", 0.2 * 0.2}};

  for (const auto &[label, neighborhood_threshold_squared] : thresholds)
    {
      find_particle_particle_periodic_contact_pairs<dim>(
        particle_handler,
        cells_local_periodic_neighbor_list,
        cells_ghost_periodic_neighbor_list,
        cells_ghost_local_periodic_neighbor_list,
        local_contact_pair_periodic_candidates,
        ghost_contact_pair_periodic_candidates,
        ghost_local_contact_pair_periodic_candidates,
        neighborhood_threshold_squared,
        periodic_boundaries_object.get_periodic_offset_per_direction());
      print_candidates<dim>("dim = " + std::to_string(dim) + ", periodic, " +
                              label,
                            local_contact_pair_periodic_candidates);
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test_local<2>();
      test_local<3>();
      test_periodic<2>();
      test_periodic<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
