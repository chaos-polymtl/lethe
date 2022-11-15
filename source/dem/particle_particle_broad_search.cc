#include <dem/dem_container_manager.h>
#include <dem/particle_particle_broad_search.h>

using namespace dealii;

template <int dim>
ParticleParticleBroadSearch<dim>::ParticleParticleBroadSearch()
{}

template <int dim>
void
ParticleParticleBroadSearch<dim>::find_particle_particle_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  DEMContainerManager<dim> &               container_manager)
{
  // Pre-fetch and clear containers
  auto &local_contact_pair_candidates =
    container_manager.local_contact_pair_candidates;
  auto &ghost_contact_pair_candidates =
    container_manager.ghost_contact_pair_candidates;
  auto &cells_local_neighbor_list = container_manager.cells_local_neighbor_list;
  auto &cells_ghost_neighbor_list = container_manager.cells_ghost_neighbor_list;
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
              store_candidates(
                std::next(particle_in_main_cell, 1),
                particles_in_main_cell,
                local_contact_pair_candidates[particle_in_main_cell->get_id()]);
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
                  store_candidates(particles_in_neighbor_cell.begin(),
                                   particles_in_neighbor_cell,
                                   local_contact_pair_candidates
                                     [particle_in_main_cell->get_id()]);
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
                  store_candidates(particles_in_neighbor_cell.begin(),
                                   particles_in_neighbor_cell,
                                   ghost_contact_pair_candidates
                                     [particle_in_main_cell->get_id()]);
                }
            }
        }
    }
}

template <int dim>
void
ParticleParticleBroadSearch<dim>::find_particle_particle_periodic_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  DEMContainerManager<dim> &               container_manager)
{
  // Pre-fetch and clear containers
  auto &local_contact_pair_periodic_candidates =
    container_manager.local_contact_pair_periodic_candidates;
  auto &ghost_contact_pair_periodic_candidates =
    container_manager.ghost_contact_pair_periodic_candidates;
  auto &cells_local_periodic_neighbor_list =
    container_manager.cells_local_periodic_neighbor_list;
  auto &cells_ghost_periodic_neighbor_list =
    container_manager.cells_ghost_periodic_neighbor_list;
  local_contact_pair_periodic_candidates.clear();
  ghost_contact_pair_periodic_candidates.clear();

  // Looping over the potential cells which may contain particles.
  for (auto cell_periodic_neighbor_list_iterator =
         cells_local_periodic_neighbor_list.begin();
       cell_periodic_neighbor_list_iterator !=
       cells_local_periodic_neighbor_list.end();
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

      // Check to see if the main cell has any particles
      if (particles_exist_in_main_cell)
        {
          // Going through periodic neighbor cells of the main cell
          ++cell_periodic_neighbor_iterator;
          for (; cell_periodic_neighbor_iterator !=
                 cell_periodic_neighbor_list_iterator->end();
               ++cell_periodic_neighbor_iterator)
            {
              // Defining iterator on local particles in the neighbor cell
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
                  store_candidates(particles_in_periodic_neighbor_cell.begin(),
                                   particles_in_periodic_neighbor_cell,
                                   local_contact_pair_periodic_candidates
                                     [particle_in_main_cell->get_id()]);
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

      std::cout << "main cell: "
                << (*cell_periodic_neighbor_iterator)->active_cell_index()
                << std::endl;

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
              std::cout
                << "  neighbor cell: "
                << (*cell_periodic_neighbor_iterator)->active_cell_index()
                << std::endl;
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
                  std::cout
                    << "    particle: " << particle_in_main_cell->get_id()
                    << std::endl;

                  // Create a arbitrary temporary empty container
                  if (ghost_contact_pair_periodic_candidates
                        [particle_in_main_cell->get_id()]
                          .empty())
                    {
                      ghost_contact_pair_periodic_candidates
                        [particle_in_main_cell->get_id()]
                          .reserve(40);
                    }

                  // Store particle ids from the selected particle iterator
                  for (auto particle_iterator =
                         particles_in_periodic_neighbor_cell.begin();
                       particle_iterator !=
                       particles_in_periodic_neighbor_cell.end();
                       ++particle_iterator)
                    {
                      std::cout
                        << "      candidate: " << particle_iterator->get_id()
                        << std::endl;
                      ghost_contact_pair_periodic_candidates
                        [particle_in_main_cell->get_id()]
                          .emplace_back(particle_iterator->get_id());
                    }
                }
            }
        }
    }
}

template class ParticleParticleBroadSearch<2>;
template class ParticleParticleBroadSearch<3>;
