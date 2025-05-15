// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_broad_search.h>

using namespace DEM;

template <int dim>
void
find_particle_particle_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_candidates)
{
  // Clear containers
  local_contact_pair_candidates.clear();
  ghost_contact_pair_candidates.clear();

  // First we handle the local-local candidate pairs
  // Looping over the potential cells which may contain particles.
  // This includes the cell itself as well as the neighbouring cells that
  // were identified.
  // cell_neighbor_list_iterator is [cell_it, neighbor_0_it, neighbor_1_it, ...]
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      // Check to see if the main cell has any particles
      if (particles_exist_in_main_cell)
        {
          // Find local-local collision pairs in the main cell, 1st particle
          // iterator is skipped since the main particle will not be
          // considered as collision partner with itself
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(particle_in_main_cell->get_id(),
                                    std::next(particle_in_main_cell, 1),
                                    particles_in_main_cell,
                                    local_contact_pair_candidates);
            }

          // Going through neighbor cells of the main cell
          ++cell_neighbor_iterator;
          for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
               ++cell_neighbor_iterator)
            {
              // Defining iterator on local particles in the neighbor cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_neighbor_cell =
                  particle_handler.particles_in_cell(*cell_neighbor_iterator);

              // Capturing particle pairs, the first particle in the main
              // cell and the second particle in the neighbor cells
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates<dim>(particle_in_main_cell->get_id(),
                                        particles_in_neighbor_cell.begin(),
                                        particles_in_neighbor_cell,
                                        local_contact_pair_candidates);
                }
            }
        }
    }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)

  // Looping over cells_ghost_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      if (particles_exist_in_main_cell)
        {
          // Going through ghost neighbor cells of the main cell
          ++cell_neighbor_iterator;

          for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
               ++cell_neighbor_iterator)
            {
              // Defining iterator on ghost particles in the neighbor cells
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_neighbor_cell =
                  particle_handler.particles_in_cell(*cell_neighbor_iterator);

              // Capturing particle pairs, the first particle (local) in
              // the main cell and the second particle (ghost) in the
              // neighbor cells
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates<dim>(particle_in_main_cell->get_id(),
                                        particles_in_neighbor_cell.begin(),
                                        particles_in_neighbor_cell,
                                        ghost_contact_pair_candidates);
                }
            }
        }
    }
}

template <int dim, typename PropertiesIndex>
void
find_particle_particle_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  // Clear containers
  local_contact_pair_candidates.clear();
  ghost_contact_pair_candidates.clear();

  // First we handle the local-local candidate pairs

  // Looping over the potential cells which may contain particles.
  // This includes the cell itself as well as the neighbouring cells that
  // were identified.
  // cell_neighbor_list_iterator is [cell_it, neighbor_0_it, neighbor_1_it, ...]
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell & its mobility status
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(*cell_neighbor_iterator);

      // If main cell has status inactive or advected, skip to next main
      // cell No needs to check if the main cell has any particle after this
      // step since empty cells have inactive or advected mobility status
      // and are skipped
      if (main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::inactive ||
          main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
        continue;

      // Get particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);

      // Store other particles in the main cell as contact candidates if
      // main cell is mobile only (this is equivalent to when adaptive sparse
      // contacts feature is not enabled)
      if (main_cell_mobility_status ==
          AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
        {
          // Find local-local collision pairs in the main cell, 1st particle
          // iterator is skipped since the main particle will not be
          // considered as collision partner with itself
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(particle_in_main_cell->get_id(),
                                    std::next(particle_in_main_cell, 1),
                                    particles_in_main_cell,
                                    local_contact_pair_candidates);
            }
        }

      // Going through neighbor cells of the main cell
      ++cell_neighbor_iterator;
      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator)
        {
          // Get mobility status of current neighbor cell
          unsigned int neighbor_cell_mobility_status =
            sparse_contacts_object.check_cell_mobility(*cell_neighbor_iterator);

          // No storing of particles as candidate if main cell has any active
          // status, but neighbor is not mobile (active or inactive). Particle
          // contacts between particles in 2 active cells or with particles in
          // inactive cells are irrelevant. In this case, we skip this neighbor
          // cell.
          if (neighbor_cell_mobility_status !=
              AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
            {
              if (main_cell_mobility_status ==
                  AdaptiveSparseContacts<dim, PropertiesIndex>::static_active)
                continue;
              if (main_cell_mobility_status ==
                  AdaptiveSparseContacts<dim, PropertiesIndex>::advected_active)
                continue;
            }

          // Store particles in the neighbor cell as contact candidates
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_neighbor_cell =
              particle_handler.particles_in_cell(*cell_neighbor_iterator);
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(particle_in_main_cell->get_id(),
                                    particles_in_neighbor_cell.begin(),
                                    particles_in_neighbor_cell,
                                    local_contact_pair_candidates);
            }
        }
    }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)

  // Looping over cells_ghost_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell & its mobility status
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(*cell_neighbor_iterator);

      // If main cell has status inactive or advected, skip to next main
      // cell No needs to check if the main cell has any particle after this
      // step since empty cells have inactive or advected mobility status and
      // are skipped
      if (main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::inactive ||
          main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
        continue;

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);


      // Going through ghost neighbor cells of the main cell
      ++cell_neighbor_iterator;
      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator)
        {
          unsigned int neighbor_cell_mobility_status =
            sparse_contacts_object.check_cell_mobility(*cell_neighbor_iterator);

          // No storing of particles as candidate if main cell is active
          // but neighbor is not mobile
          if (neighbor_cell_mobility_status !=
              AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
            {
              if (main_cell_mobility_status ==
                  AdaptiveSparseContacts<dim, PropertiesIndex>::static_active)
                continue;
              if (main_cell_mobility_status ==
                  AdaptiveSparseContacts<dim, PropertiesIndex>::advected_active)
                continue;
            }


          // Defining iterator on ghost particles in the neighbor cells
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_neighbor_cell =
              particle_handler.particles_in_cell(*cell_neighbor_iterator);

          // Capturing particle pairs, the first particle (local) in
          // the main cell and the second particle (ghost) in the
          // neighbor cells
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(particle_in_main_cell->get_id(),
                                    particles_in_neighbor_cell.begin(),
                                    particles_in_neighbor_cell,
                                    ghost_contact_pair_candidates);
            }
        }
    }
}


template <int dim>
void
find_particle_particle_periodic_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list,
  typename dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<dim>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates)
{
  // Clear containers
  local_contact_pair_periodic_candidates.clear();
  ghost_contact_pair_periodic_candidates.clear();
  ghost_local_contact_pair_periodic_candidates.clear();

  // Looping over the potential periodic cells which may contain particles.
  for (auto cell_periodic_neighbor_list_iterator =
         cells_local_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_local_periodic_neighbor_list.end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell is on periodic boundary 0
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      // Check to see if the main cell has any particles
      if (particles_exist_in_main_cell)
        {
          // Going through periodic neighbor cells on the periodic boundary 1
          // of the main cell. We start at the second iterator since the first
          // one is the main cell.
          ++cell_periodic_neighbor_iterator;
          for (; cell_periodic_neighbor_iterator !=
                 cell_periodic_neighbor_list_iterator->end();
               ++cell_periodic_neighbor_iterator)
            {
              // Defining iterator on particles in the local periodic neighbor
              // cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_neighbor_cell =
                  particle_handler.particles_in_cell(
                    *cell_periodic_neighbor_iterator);

              // Capturing particle pairs, the first particle in the main
              // cell and the second particle in the local periodic neighbor
              // cell
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates<dim>(
                    particle_in_main_cell->get_id(),
                    particles_in_periodic_neighbor_cell.begin(),
                    particles_in_periodic_neighbor_cell,
                    local_contact_pair_periodic_candidates);
                }
            }
        }
    }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)

  // Looping over cells_ghost_neighbor_list
  for (auto cell_periodic_neighbor_list_iterator =
         cells_ghost_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_ghost_periodic_neighbor_list.end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      if (particles_exist_in_main_cell)
        {
          // Going through ghost neighbor cells of the main cell
          ++cell_periodic_neighbor_iterator;

          for (; cell_periodic_neighbor_iterator !=
                 cell_periodic_neighbor_list_iterator->end();
               ++cell_periodic_neighbor_iterator)
            {
              // Defining iterator on ghost particles in the neighbor cells
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_neighbor_cell =
                  particle_handler.particles_in_cell(
                    *cell_periodic_neighbor_iterator);

              // Capturing particle pairs, the first particle (local) in
              // the main cell and the second particle (ghost) in the
              // neighbor cells
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates<dim>(
                    particle_in_main_cell->get_id(),
                    particles_in_periodic_neighbor_cell.begin(),
                    particles_in_periodic_neighbor_cell,
                    ghost_contact_pair_periodic_candidates);
                }
            }
        }
    }

  // Looping over cells_ghost_local_periodic_neighbor_list
  for (auto cell_periodic_neighbor_list_iterator =
         cells_ghost_local_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_ghost_local_periodic_neighbor_list.end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      if (particles_exist_in_main_cell)
        {
          // Going through ghost neighbor cells of the main cell
          ++cell_periodic_neighbor_iterator;

          for (; cell_periodic_neighbor_iterator !=
                 cell_periodic_neighbor_list_iterator->end();
               ++cell_periodic_neighbor_iterator)
            {
              // Defining iterator on local particles in the neighbor cells
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_periodic_neighbor_cell =
                  particle_handler.particles_in_cell(
                    *cell_periodic_neighbor_iterator);

              // Capturing particle pairs, the first particle (ghost) in
              // the main cell and the second particles (local) in the
              // neighbor cells
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates<dim>(
                    particle_in_main_cell->get_id(),
                    particles_in_periodic_neighbor_cell.begin(),
                    particles_in_periodic_neighbor_cell,
                    ghost_local_contact_pair_periodic_candidates);
                }
            }
        }
    }
}

template <int dim, typename PropertiesIndex>
void
find_particle_particle_periodic_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list,
  typename dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<dim>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  // Clear containers
  local_contact_pair_periodic_candidates.clear();
  ghost_contact_pair_periodic_candidates.clear();
  ghost_local_contact_pair_periodic_candidates.clear();

  // Looping over the potential periodic cells which may contain particles.
  for (auto cell_periodic_neighbor_list_iterator =
         cells_local_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_local_periodic_neighbor_list.end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell on periodic boundary 0 & its mobility status
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(
          *cell_periodic_neighbor_iterator);

      // If main cell has status inactive or advected, skip to next main cell
      // No needs to check if the main cell has any particle after this
      // step since empty cells have inactive or advected mobility status
      // and are skipped
      if (main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::inactive ||
          main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
        continue;

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      // Going through periodic neighbor cells on the periodic boundary 1
      // of the main cell
      ++cell_periodic_neighbor_iterator;
      for (; cell_periodic_neighbor_iterator !=
             cell_periodic_neighbor_list_iterator->end();
           ++cell_periodic_neighbor_iterator)
        {
          unsigned int neighbor_cell_mobility_status =
            sparse_contacts_object.check_cell_mobility(
              *cell_periodic_neighbor_iterator);

          // No storing of particles if main cell is active but neighbor is
          // not mobile
          if ((main_cell_mobility_status ==
                 AdaptiveSparseContacts<dim, PropertiesIndex>::static_active ||
               main_cell_mobility_status ==
                 AdaptiveSparseContacts<dim,
                                        PropertiesIndex>::advected_active) &&
              neighbor_cell_mobility_status !=
                AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
            continue;

          // Defining iterator on local particles in the periodic neighbor
          // cell
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_periodic_neighbor_cell =
              particle_handler.particles_in_cell(
                *cell_periodic_neighbor_iterator);

          // Capturing particle pairs, the first particle in the main
          // cell and the second particle in the neighbor cells
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(particle_in_main_cell->get_id(),
                                    particles_in_periodic_neighbor_cell.begin(),
                                    particles_in_periodic_neighbor_cell,
                                    local_contact_pair_periodic_candidates);
            }
        }
    }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)
  for (auto cell_periodic_neighbor_list_iterator =
         cells_ghost_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_ghost_periodic_neighbor_list.end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell & its mobility status
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(
          *cell_periodic_neighbor_iterator);

      // If main cell has status inactive or advected, skip to next main
      // cell No needs to check if the main cell has any particle after this
      // step since empty cells have inactive or advected mobility status
      // and are skipped
      if (main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::inactive ||
          main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
        continue;

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      // Going through ghost neighbor cells of the main cell
      ++cell_periodic_neighbor_iterator;
      for (; cell_periodic_neighbor_iterator !=
             cell_periodic_neighbor_list_iterator->end();
           ++cell_periodic_neighbor_iterator)
        {
          unsigned int neighbor_cell_mobility_status =
            sparse_contacts_object.check_cell_mobility(
              *cell_periodic_neighbor_iterator);

          // No storing of particles if main cell is active but neighbor is
          // not mobile
          if ((main_cell_mobility_status ==
                 AdaptiveSparseContacts<dim, PropertiesIndex>::static_active ||
               main_cell_mobility_status ==
                 AdaptiveSparseContacts<dim,
                                        PropertiesIndex>::advected_active) &&
              neighbor_cell_mobility_status !=
                AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
            continue;

          // Defining iterator on ghost particles in the neighbor cells
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_periodic_neighbor_cell =
              particle_handler.particles_in_cell(
                *cell_periodic_neighbor_iterator);

          // Capturing particle pairs, the first particle (local) in
          // the main cell and the second particle (ghost) in the
          // neighbor cells
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(particle_in_main_cell->get_id(),
                                    particles_in_periodic_neighbor_cell.begin(),
                                    particles_in_periodic_neighbor_cell,
                                    ghost_contact_pair_periodic_candidates);
            }
        }
    }

  // Looping over cells_ghost_local_periodic_neighbor_list
  for (auto cell_periodic_neighbor_list_iterator =
         cells_ghost_local_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_ghost_local_periodic_neighbor_list.end();
       ++cell_periodic_neighbor_list_iterator)
    {
      // The main cell
      auto cell_periodic_neighbor_iterator =
        cell_periodic_neighbor_list_iterator->begin();

      // If main cell has status inactive or advected, skip to next main
      // cell No needs to check if the main cell has any particle after this
      // step since empty cells have inactive or advected mobility status
      // and are skipped
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(
          *cell_periodic_neighbor_iterator);
      if (main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::inactive ||
          main_cell_mobility_status ==
            AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
        continue;

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_periodic_neighbor_iterator);

      // Going through ghost neighbor cells of the main cell
      ++cell_periodic_neighbor_iterator;
      for (; cell_periodic_neighbor_iterator !=
             cell_periodic_neighbor_list_iterator->end();
           ++cell_periodic_neighbor_iterator)
        {
          unsigned int neighbor_cell_mobility_status =
            sparse_contacts_object.check_cell_mobility(
              *cell_periodic_neighbor_iterator);

          // No storing of particles if main cell is active but neighbor is
          // not mobile
          if ((main_cell_mobility_status ==
                 AdaptiveSparseContacts<dim, PropertiesIndex>::static_active ||
               main_cell_mobility_status ==
                 AdaptiveSparseContacts<dim,
                                        PropertiesIndex>::advected_active) &&
              neighbor_cell_mobility_status !=
                AdaptiveSparseContacts<dim, PropertiesIndex>::mobile)
            continue;

          // Defining iterator on local particles in the neighbor cells
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_periodic_neighbor_cell =
              particle_handler.particles_in_cell(
                *cell_periodic_neighbor_iterator);

          // Capturing particle pairs, the first particle (ghost) in
          // the main cell and the second particles (local) in the
          // neighbor cells
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates<dim>(
                particle_in_main_cell->get_id(),
                particles_in_periodic_neighbor_cell.begin(),
                particles_in_periodic_neighbor_cell,
                ghost_local_contact_pair_periodic_candidates);
            }
        }
    }
}

template <int dim>
void
store_candidates(
  const types::particle_index &main_particle_id,
  const typename Particles::ParticleHandler<
    dim>::particle_iterator_range::iterator &particle_begin,
  const typename Particles::ParticleHandler<dim>::particle_iterator_range
    &particles_to_evaluate,
  typename dem_data_structures<dim>::particle_particle_candidates
    &contact_pair_candidates)
{
  // Find the contact candidate container of the main particle
  auto candidates_container_it = contact_pair_candidates.find(main_particle_id);

  // Reserve arbitrary vector capacity and store if the particle does not have
  // contact candidate yet
  if (candidates_container_it == contact_pair_candidates.end())
    {
      std::vector<types::particle_index> candidates_container;
      candidates_container.reserve(40);

      // Insert the empty vector and get the iterator to the inserted element
      // prior storing the particle ids
      auto pair_it_bool =
        contact_pair_candidates.emplace(main_particle_id, candidates_container);
      candidates_container_it = pair_it_bool.first;
    }

  // Store particle ids from the selected particle iterator
  for (auto particle_iterator = particle_begin;
       particle_iterator != particles_to_evaluate.end();
       ++particle_iterator)
    {
      candidates_container_it->second.emplace_back(particle_iterator->get_id());
    }
}

template void
find_particle_particle_contact_pairs<2>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_candidates);

template void
find_particle_particle_contact_pairs<3>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_candidates);

template void
find_particle_particle_contact_pairs<2, DEM::DEMProperties::PropertiesIndex>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_contact_pairs<3, DEM::DEMProperties::PropertiesIndex>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_periodic_contact_pairs<2>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates);

template void
find_particle_particle_periodic_contact_pairs<3>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates);

template void
find_particle_particle_periodic_contact_pairs<2>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_periodic_contact_pairs<3>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_contact_pairs<2, DEM::CFDDEMProperties::PropertiesIndex>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<2, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_contact_pairs<3, DEM::CFDDEMProperties::PropertiesIndex>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<3, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_periodic_contact_pairs<2>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<2, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_periodic_contact_pairs<3>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<3, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_contact_pairs<2, DEM::DEMMPProperties::PropertiesIndex>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_contact_pairs<3, DEM::DEMMPProperties::PropertiesIndex>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_periodic_contact_pairs<2>(
  dealii::Particles::ParticleHandler<2> &particle_handler,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<2>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<2>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
find_particle_particle_periodic_contact_pairs<3>(
  dealii::Particles::ParticleHandler<3> &particle_handler,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_local_neighbor_list,
  typename dem_data_structures<3>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename dem_data_structures<3>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);
