#include <dem/pp_broad_search.h>

using namespace dealii;

template <int dim>
PPBroadSearch<dim>::PPBroadSearch()
{}

template <int dim>
void
PPBroadSearch<dim>::find_PP_Contact_Pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const std::vector<
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
    *cells_local_neighbor_list,
  const std::vector<
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
    *                                        cells_ghost_neighbor_list,
  std::unordered_map<int, std::vector<int>> &local_contact_pair_candidates,
  std::unordered_map<int, std::vector<int>> &ghost_contact_pair_candidates)
{
  // First we will handle the local-lcoal candidate pairs

  // Clearing local_contact_pair_candidates
  local_contact_pair_candidates.clear();
  ghost_contact_pair_candidates.clear();

  // Looping over cells_local_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list->begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list->end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

      // Check to see if the main cell has any particles
      if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) > 0)
        {
          // Particles in the main cell
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_main_cell =
              particle_handler.particles_in_cell(*cell_neighbor_iterator);

          // finding local-local collision pairs in the main cell, particle
          // counter starts from 1, becasue each particle will not be considered
          // as collision partner with itself
          for (typename Particles::ParticleHandler<
                 dim>::particle_iterator_range::iterator
                 particles_in_main_cell_iterator_one =
                   particles_in_main_cell.begin();
               particles_in_main_cell_iterator_one !=
               particles_in_main_cell.end();
               ++particles_in_main_cell_iterator_one)
            {
              const unsigned int particle_one_id =
                particles_in_main_cell_iterator_one->get_id();

              std::vector<int> &particle_candidate_container =
                local_contact_pair_candidates[particle_one_id];

              if (particle_candidate_container.empty())
                {
                  particle_candidate_container.reserve(40);
                }

              // Advancing the second iterator to capture all the particle pairs
              // in the main cell
              auto particles_in_main_cell_iterator_two =
                std::next(particles_in_main_cell_iterator_one, 1);

              for (; particles_in_main_cell_iterator_two !=
                     particles_in_main_cell.end();
                   ++particles_in_main_cell_iterator_two)
                {
                  // Capturing all the local-local particle pairs in the main
                  // cell
                  particle_candidate_container.emplace_back(
                    particles_in_main_cell_iterator_two->get_id());
                }
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

              // Capturing particle pairs, the first particle in the main cell
              // and the second particle in the neighbor cells
              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_main_cell_iterator =
                       particles_in_main_cell.begin();
                   particles_in_main_cell_iterator !=
                   particles_in_main_cell.end();
                   ++particles_in_main_cell_iterator)
                {
                  const unsigned int particle_one_id =
                    particles_in_main_cell_iterator->get_id();

                  std::vector<int> &particle_candidate_container =
                    local_contact_pair_candidates[particle_one_id];
                  if (particle_candidate_container.empty())
                    particle_candidate_container.reserve(40);

                  for (typename Particles::ParticleHandler<
                         dim>::particle_iterator_range::iterator
                         particles_in_neighbor_cell_iterator =
                           particles_in_neighbor_cell.begin();
                       particles_in_neighbor_cell_iterator !=
                       particles_in_neighbor_cell.end();
                       ++particles_in_neighbor_cell_iterator)
                    {
                      particle_candidate_container.emplace_back(
                        particles_in_neighbor_cell_iterator->get_id());
                    }
                }
            }
        }
    }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)

  // Clearing ghost_contact_pair_candidates
  ghost_contact_pair_candidates.clear();

  // Looping over cells_ghost_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list->begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list->end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);

      // If the main cell is not empty
      if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) > 0)
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

              // Capturing particle pairs, the first particle (local) in the
              // main cell and the second particle (ghost) in the neighbor cells
              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_main_cell_iterator =
                       particles_in_main_cell.begin();
                   particles_in_main_cell_iterator !=
                   particles_in_main_cell.end();
                   ++particles_in_main_cell_iterator)
                {
                  unsigned int particle_one_id =
                    particles_in_main_cell_iterator->get_id();

                  std::vector<int> &particle_candidate_container =
                    ghost_contact_pair_candidates[particle_one_id];
                  particle_candidate_container.reserve(40);

                  for (typename Particles::ParticleHandler<
                         dim>::particle_iterator_range::iterator
                         particles_in_neighbor_cell_iterator =
                           particles_in_neighbor_cell.begin();
                       particles_in_neighbor_cell_iterator !=
                       particles_in_neighbor_cell.end();
                       ++particles_in_neighbor_cell_iterator)
                    {
                      unsigned int particle_two_id =
                        particles_in_neighbor_cell_iterator->get_id();

                      particle_candidate_container.push_back(particle_two_id);
                    }
                }
            }
        }
    }
}

template class PPBroadSearch<2>;
template class PPBroadSearch<3>;
