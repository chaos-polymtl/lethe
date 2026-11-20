// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_wall_broad_search.h>

#include <numeric>
#include <vector>

using namespace dealii;

namespace
{
  /**
   * @brief Call a function once for every pair of a solid cell and a
   * background cell that neighbors a background cell to which the solid cell
   * is mapped.
   *
   * Only the mapped background cells that are locally owned or ghost are
   * considered. The background cells of a solid cell are visited in the order
   * of the mapped pairs, then in the order of the neighbor lists. A background
   * cell that neighbors several mapped background cells of the same solid cell
   * is visited only once, the first time it is reached.
   *
   * @tparam dim An integer that denotes the number of spatial dimensions.
   * @tparam VisitFunction Type of the function called for every pair. It takes
   * the solid cell and the background cell as arguments.
   *
   * @param[in] mapped_pairs The (background cell, solid cell) pairs of one
   * solid surface mapped in the background triangulation.
   * @param[in] cells_total_neighbor_list A container in which all the neighbor
   * cells of the local cells are stored.
   * @param[in] visit_cell Function called for every pair of a solid cell and a
   * background cell.
   */
  template <int dim, typename VisitFunction>
  void
  for_each_solid_cell_and_neighbor_cell(
    const typename DEM::dem_data_structures<
      dim>::solid_surfaces_mesh_information::value_type &mapped_pairs,
    const typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
                        &cells_total_neighbor_list,
    const VisitFunction &visit_cell)
  {
    if (mapped_pairs.empty())
      return;

    // Group the pairs by solid cell with a counting sort, which keeps the
    // order of the pairs of each solid cell
    const unsigned int n_solid_cells =
      mapped_pairs.front().second->get_triangulation().n_active_cells();
    std::vector<unsigned int> first_pair_of_solid_cell(n_solid_cells + 1, 0);
    for (const auto &mapped_pair : mapped_pairs)
      ++first_pair_of_solid_cell[mapped_pair.second->active_cell_index() + 1];
    std::partial_sum(first_pair_of_solid_cell.begin(),
                     first_pair_of_solid_cell.end(),
                     first_pair_of_solid_cell.begin());

    std::vector<unsigned int> pairs_grouped_by_solid_cell(mapped_pairs.size());
    std::vector<unsigned int> next_position_of_solid_cell(
      first_pair_of_solid_cell.begin(), first_pair_of_solid_cell.end() - 1);
    for (unsigned int pair_index = 0; pair_index < mapped_pairs.size();
         ++pair_index)
      {
        const unsigned int solid_cell_index =
          mapped_pairs[pair_index].second->active_cell_index();
        pairs_grouped_by_solid_cell
          [next_position_of_solid_cell[solid_cell_index]++] = pair_index;
      }

    // Index of the last solid cell for which each background cell was visited
    std::vector<unsigned int> last_visiting_solid_cell(
      mapped_pairs.front().first->get_triangulation().n_active_cells(),
      numbers::invalid_unsigned_int);

    for (unsigned int solid_cell_index = 0; solid_cell_index < n_solid_cells;
         ++solid_cell_index)
      for (unsigned int position = first_pair_of_solid_cell[solid_cell_index];
           position < first_pair_of_solid_cell[solid_cell_index + 1];
           ++position)
        {
          const auto &[background_cell, solid_cell] =
            mapped_pairs[pairs_grouped_by_solid_cell[position]];

          if (!background_cell->is_locally_owned() &&
              !background_cell->is_ghost())
            continue;

          for (const auto &neighbor_cell : cells_total_neighbor_list.at(
                 background_cell->global_active_cell_index()))
            {
              unsigned int &last_solid_cell =
                last_visiting_solid_cell[neighbor_cell->active_cell_index()];
              if (last_solid_cell == solid_cell_index)
                continue;
              last_solid_cell = solid_cell_index;

              visit_cell(solid_cell, neighbor_cell);
            }
        }
  }
} // namespace

template <int dim>
void
find_particle_wall_contact_pairs(
  const std::map<int, boundary_cells_info_struct<dim>>
                                        &boundary_cells_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &particle_wall_contact_candidates)
{
  // Clearing particle_wall_contact_candidates (output of this function)
  particle_wall_contact_candidates.clear();

  // Iterating over the boundary_cells_information, which is the output of
  // the find_boundary_cells_information find_boundary_cells_information class.
  // This map contains all the required information of the system boundary
  // cells and faces. In this loop we find the particles located in each of
  // these boundary cells
  for (auto boundary_cells_information_iterator =
         boundary_cells_information.begin();
       boundary_cells_information_iterator != boundary_cells_information.end();
       ++boundary_cells_information_iterator)
    {
      auto boundary_cells_content = boundary_cells_information_iterator->second;
      auto cell                   = boundary_cells_content.cell;

      // Finding particles located in the corresponding cell
      // (boundary_cells_content.cell)
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(cell);

      const bool particles_exist_in_main_cell = !particles_in_cell.empty();

      // If the main cell is not empty
      if (particles_exist_in_main_cell)
        {
          for (typename Particles::ParticleHandler<dim>::
                 particle_iterator_range::iterator particle_in_cell_iterator =
                   particles_in_cell.begin();
               particle_in_cell_iterator != particles_in_cell.end();
               ++particle_in_cell_iterator)
            {
              // Store particles candidate to wall with boundary cell info
              store_candidates(particle_in_cell_iterator,
                               boundary_cells_content,
                               particle_wall_contact_candidates);
            }
        }
    }
}

template <int dim>
void
find_particle_floating_wall_contact_pairs(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<dim>::active_cell_iterator>>
                                        &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      simulation_time,
  typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
    &particle_floating_wall_candidates)
{
  // Clearing particle_floating_wall_candidates(output of this function)
  particle_floating_wall_candidates.clear();

  // Iterating over the boundary_cells_for_floating_walls, which is the output
  // of the find_boundary_cells_for_floating_walls function in
  // find_boundary_cells_information class. This unordered_map contains all the
  // boundary cells of floating walls. In this loop
  // we find the particles located in boundary cells of floating
  // walls
  for (auto fw_boundary_cells_information_iterator =
         boundary_cells_for_floating_walls.begin();
       fw_boundary_cells_information_iterator !=
       boundary_cells_for_floating_walls.end();
       ++fw_boundary_cells_information_iterator)
    {
      auto floating_wall_id = fw_boundary_cells_information_iterator->first;

      // Checking simulation time for temporary floating walls
      if (simulation_time >=
            floating_wall_properties.time_start[floating_wall_id] &&
          simulation_time <=
            floating_wall_properties.time_end[floating_wall_id])
        {
          auto boundary_cells_content =
            fw_boundary_cells_information_iterator->second;

          for (auto cell = boundary_cells_content.begin();
               cell != boundary_cells_content.end();
               ++cell)
            {
              // Finding particles located in the corresponding cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_cell = particle_handler.particles_in_cell(*cell);

              const bool particles_exist_in_main_cell =
                !particles_in_cell.empty();

              // If the main cell is not empty
              if (particles_exist_in_main_cell)
                {
                  for (typename Particles::ParticleHandler<
                         dim>::particle_iterator_range::iterator
                         particle_in_cell_iterator = particles_in_cell.begin();
                       particle_in_cell_iterator != particles_in_cell.end();
                       ++particle_in_cell_iterator)
                    {
                      store_candidates<dim>(particle_in_cell_iterator,
                                            floating_wall_id,
                                            particle_floating_wall_candidates);
                    }
                }
            }
        }
    }
}

template <int dim>
void
particle_solid_surfaces_contact_search(
  const typename DEM::dem_data_structures<dim>::solid_surfaces_mesh_information
                                        &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
    &cells_total_neighbor_list)
{
  // Clear the candidate container
  particle_floating_mesh_contact_candidates.clear();
  particle_floating_mesh_contact_candidates.resize(
    solid_surfaces_mesh_information.size());

  for (unsigned int solid_counter = 0;
       solid_counter < solid_surfaces_mesh_information.size();
       ++solid_counter)
    {
      auto &candidates =
        particle_floating_mesh_contact_candidates[solid_counter];

      // BB NOTE
      // To look more into this, I am not sure to which extent this is true
      // The code has never really been tested with more than one floating solid

      // "I am assuming that triangles in different solids have different unique
      // global ids. If it's not the case, we have to modify the code" - ???

      // The particles of every background cell that neighbors a background
      // cell to which a solid cell is mapped are candidates of the solid cell
      for_each_solid_cell_and_neighbor_cell<dim>(
        solid_surfaces_mesh_information[solid_counter],
        cells_total_neighbor_list,
        [&](const auto &solid_cell, const auto &neighbor_cell) {
          const auto particles_in_cell =
            particle_handler.particles_in_cell(neighbor_cell);

          if (particles_in_cell.empty())
            return;

          auto &solid_cell_candidates = candidates[solid_cell];
          for (auto particle_in_cell_iterator = particles_in_cell.begin();
               particle_in_cell_iterator != particles_in_cell.end();
               ++particle_in_cell_iterator)
            solid_cell_candidates.insert(
              {particle_in_cell_iterator->get_id(), particle_in_cell_iterator});
        });
    }
}

template <int dim, typename PropertiesIndex>
void
find_particle_wall_contact_pairs(
  const std::map<int, boundary_cells_info_struct<dim>>
                                        &boundary_cells_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  // Clearing particle_wall_contact_candidates (output of this function)
  particle_wall_contact_candidates.clear();

  // Iterating over the boundary_cells_information, which is the output of
  // the find_boundary_cells_information find_boundary_cells_information class.
  // This map contains all the required information of the system boundary
  // cells and faces. In this loop we find the particles located in each of
  // these boundary cells
  for (auto boundary_cells_information_iterator =
         boundary_cells_information.begin();
       boundary_cells_information_iterator != boundary_cells_information.end();
       ++boundary_cells_information_iterator)
    {
      auto boundary_cells_content = boundary_cells_information_iterator->second;
      auto cell                   = boundary_cells_content.cell;

      // If main cell has status other than mobile, skip to next cell
      // No needs to check if the main cell has any particle after this
      // step since mobile cells have particles
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(cell);
      if (main_cell_mobility_status !=
          AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
        continue;

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(cell);

      for (typename Particles::ParticleHandler<
             dim>::particle_iterator_range::iterator particle_in_cell_iterator =
             particles_in_cell.begin();
           particle_in_cell_iterator != particles_in_cell.end();
           ++particle_in_cell_iterator)
        {
          // Store particles candidate to wall with boundary cell info
          store_candidates(particle_in_cell_iterator,
                           boundary_cells_content,
                           particle_wall_contact_candidates);
        }
    }
}

template <int dim, typename PropertiesIndex>
void
find_particle_floating_wall_contact_pairs(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<dim>::active_cell_iterator>>
                                        &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      simulation_time,
  typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  // Clearing particle_floating_wall_candidates(output of this function)
  particle_floating_wall_candidates.clear();

  // Iterating over the boundary_cells_for_floating_walls, which is the output
  // of the find_boundary_cells_for_floating_walls function in
  // find_boundary_cells_information class. This unordered_map contains all the
  // boundary cells of floating walls. In this loop
  // we find the particles located in boundary cells of floating
  // walls
  for (auto fw_boundary_cells_information_iterator =
         boundary_cells_for_floating_walls.begin();
       fw_boundary_cells_information_iterator !=
       boundary_cells_for_floating_walls.end();
       ++fw_boundary_cells_information_iterator)
    {
      auto floating_wall_id = fw_boundary_cells_information_iterator->first;

      // Checking simulation time for temporary floating walls
      if (simulation_time >=
            floating_wall_properties.time_start[floating_wall_id] &&
          simulation_time <=
            floating_wall_properties.time_end[floating_wall_id])
        {
          auto boundary_cells_content =
            fw_boundary_cells_information_iterator->second;

          for (auto cell = boundary_cells_content.begin();
               cell != boundary_cells_content.end();
               ++cell)
            {
              // If main cell has status other than mobile, skip to next cell
              // No needs to check if the main cell has any particle after this
              // step since mobile cells have particles
              unsigned int main_cell_mobility_status =
                sparse_contacts_object.check_cell_mobility(*cell);
              if (main_cell_mobility_status !=
                  AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
                continue;

              // Finding particles located in the corresponding cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_cell = particle_handler.particles_in_cell(*cell);

              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particle_in_cell_iterator = particles_in_cell.begin();
                   particle_in_cell_iterator != particles_in_cell.end();
                   ++particle_in_cell_iterator)
                {
                  store_candidates<dim>(particle_in_cell_iterator,
                                        floating_wall_id,
                                        particle_floating_wall_candidates);
                }
            }
        }
    }
}

template <int dim, typename PropertiesIndex>
void
particle_solid_surfaces_contact_search(
  const typename DEM::dem_data_structures<dim>::solid_surfaces_mesh_information
                                        &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
                                                     &cells_total_neighbor_list,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  // Clear the candidate container
  particle_floating_mesh_contact_candidates.clear();
  particle_floating_mesh_contact_candidates.resize(
    solid_surfaces_mesh_information.size());

  for (unsigned int solid_counter = 0;
       solid_counter < solid_surfaces_mesh_information.size();
       ++solid_counter)
    {
      auto &candidates =
        particle_floating_mesh_contact_candidates[solid_counter];

      // I am assuming that triangles in different solids have different unique
      // global ids. If it's not the case, we have to modify the code

      // The particles of every mobile background cell that neighbors a
      // background cell to which a solid cell is mapped are candidates of the
      // solid cell
      for_each_solid_cell_and_neighbor_cell<dim>(
        solid_surfaces_mesh_information[solid_counter],
        cells_total_neighbor_list,
        [&](const auto &solid_cell, const auto &neighbor_cell) {
          // No need to check if the cell has any particle after this step
          // since mobile cells have particles
          if (sparse_contacts_object.check_cell_mobility(neighbor_cell) !=
              AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
            return;

          const auto particles_in_cell =
            particle_handler.particles_in_cell(neighbor_cell);

          auto &solid_cell_candidates = candidates[solid_cell];
          for (auto particle_in_cell_iterator = particles_in_cell.begin();
               particle_in_cell_iterator != particles_in_cell.end();
               ++particle_in_cell_iterator)
            solid_cell_candidates.insert(
              {particle_in_cell_iterator->get_id(), particle_in_cell_iterator});
        });
    }
}

template <int dim>
void
store_candidates(
  const typename Particles::ParticleHandler<
    dim>::particle_iterator_range::iterator &particle_iterator,
  const DEM::global_face_id                 &floating_wall_id,
  typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
    &contact_pair_candidates)
{
  // Find the contact candidate container of the particle
  const types::particle_index particle_id = particle_iterator->get_id();
  auto candidates_container_it = contact_pair_candidates.find(particle_id);

  // Reserve arbitrary vector capacity and store if the particle does not have
  // contact candidate yet
  if (candidates_container_it == contact_pair_candidates.end())
    {
      auto pair_it_bool = contact_pair_candidates.emplace(
        particle_id,
        std::unordered_map<DEM::global_face_id,
                           Particles::ParticleIterator<dim>>());

      candidates_container_it = pair_it_bool.first;
    }

  // Store particle ids from the selected particle iterator
  candidates_container_it->second.emplace(floating_wall_id, particle_iterator);
}


// Explicit Instantiation
template void
find_particle_wall_contact_pairs<2>(
  const std::map<int, boundary_cells_info_struct<2>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_wall_candidates
    &particle_wall_contact_candidates);

template void
find_particle_wall_contact_pairs<3>(
  const std::map<int, boundary_cells_info_struct<3>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_wall_candidates
    &particle_wall_contact_candidates);

template void
find_particle_wall_contact_pairs<2, DEM::DEMProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<2>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<3, DEM::DEMProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<3>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<2>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<2, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<3>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<3, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<2, DEM::DEMMPProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<2>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<2, DEM::CFDDEMMPProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<2>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<2, DEM::CFDDEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<3, DEM::DEMMPProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<3>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_wall_contact_pairs<3, DEM::CFDDEMMPProperties::PropertiesIndex>(
  const std::map<int, boundary_cells_info_struct<3>>
                                      &boundary_cells_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<3, DEM::CFDDEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<2>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<2>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<2> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &particle_floating_wall_candidates);

template void
find_particle_floating_wall_contact_pairs<3>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<3>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<3> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &particle_floating_wall_candidates);

template void
find_particle_floating_wall_contact_pairs<2,
                                          DEM::DEMProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<2>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<2> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<3,
                                          DEM::DEMProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<3>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<3> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<
  2,
  DEM::CFDDEMProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<2>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<2> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<2, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<
  3,
  DEM::CFDDEMProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<3>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<3> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<3, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<
  2,
  DEM::DEMMPProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<2>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<2> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<
  2,
  DEM::CFDDEMMPProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<2>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<2> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<2, DEM::CFDDEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<
  3,
  DEM::DEMMPProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<3>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<3> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_floating_wall_contact_pairs<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex>(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<3>::active_cell_iterator>>
                                      &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<3> &floating_wall_properties,
  const double                                    simulation_time,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<3, DEM::CFDDEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<2>(
  const DEM::dem_data_structures<2>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<2>::cells_total_neighbor_list
    &cells_total_neighbor_list);

template void
particle_solid_surfaces_contact_search<3>(
  const DEM::dem_data_structures<3>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<3>::cells_total_neighbor_list
    &cells_total_neighbor_list);

template void
particle_solid_surfaces_contact_search<2, DEM::DEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<2>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<3, DEM::DEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<3>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<2,
                                       DEM::CFDDEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<2>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<2, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<3,
                                       DEM::CFDDEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<3>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<3, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<2,
                                       DEM::DEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<2>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<
  2,
  DEM::CFDDEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<2> &particle_handler,
  DEM::dem_data_structures<2>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<2>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<2, DEM::CFDDEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<3,
                                       DEM::DEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<3>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
particle_solid_surfaces_contact_search<
  3,
  DEM::CFDDEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::solid_surfaces_mesh_information
                                      &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<3> &particle_handler,
  DEM::dem_data_structures<3>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  DEM::dem_data_structures<3>::cells_total_neighbor_list
    &cells_total_neighbor_list,
  const AdaptiveSparseContacts<3, DEM::CFDDEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);
