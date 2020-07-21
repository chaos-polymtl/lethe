#include <dem/pp_broad_search.h>

using namespace dealii;

template <int dim> PPBroadSearch<dim>::PPBroadSearch() {}

template <int dim>
void PPBroadSearch<dim>::find_PP_Contact_Pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    const std::vector<
        std::vector<typename Triangulation<dim>::active_cell_iterator>>
        *cells_local_neighbor_list,
    const std::vector<
        std::vector<typename Triangulation<dim>::active_cell_iterator>>
        *cells_ghost_neighbor_list,
    std::map<std::pair<int, int>,
             std::pair<typename Particles::ParticleIterator<dim>,
                       typename Particles::ParticleIterator<dim>>>
        &local_contact_pair_candidates,
    std::map<std::pair<int, int>,
             std::pair<typename Particles::ParticleIterator<dim>,
                       typename Particles::ParticleIterator<dim>>>
        &ghost_contact_pair_candidates) {

  // First we will handle the local candiated pairs

  // Since the local_contact_pair_candidates (which is the real output of the
  // function) is defined as an input of the function, it should be cleared
  local_contact_pair_candidates.clear();

  // Looping over cell_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list->begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list->end();
       ++cell_neighbor_list_iterator) {

    // The main cell
    auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

    // NOTE - Important, this is very weird
    //      //      std::cout << "looking if main cell has particles " <<
    //      std::endl; std::cout << "number of particles : "
    //                <<
    //                particle_handler.n_particles_in_cell(*cell_neighbor_iterator)
    //                << std::endl;
    //      // Check to see if the main cell has any particles
    if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) > 0) {

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
          particles_in_main_cell =
              particle_handler.particles_in_cell(*cell_neighbor_iterator);

      // finding collision pairs in the main cell, particle counter starts
      // from 1, becasue each particle will not be considered as collision
      // partner with itself
      for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
               iterator particles_in_main_cell_iterator_one =
                   particles_in_main_cell.begin();
           particles_in_main_cell_iterator_one != particles_in_main_cell.end();
           ++particles_in_main_cell_iterator_one) {
        // Advancing the second iterator to capture all the particle pairs
        // in the main cell
        auto particles_in_main_cell_iterator_two =
            std::next(particles_in_main_cell_iterator_one, 1);

        for (; particles_in_main_cell_iterator_two !=
               particles_in_main_cell.end();
             ++particles_in_main_cell_iterator_two) {

          // Capturing all the particle pairs in the main cell
          auto contact_pair =
              std::make_pair(particles_in_main_cell_iterator_one,
                             particles_in_main_cell_iterator_two);
          local_contact_pair_candidates.insert(
              {std::make_pair(particles_in_main_cell_iterator_one->get_id(),
                              particles_in_main_cell_iterator_two->get_id()),
               contact_pair});
        }
      }

      // Going through neighbor cells of the main cell
      ++cell_neighbor_iterator;

      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator) {
        // Defining iterator on particles in the neighbor cell
        typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_neighbor_cell =
                particle_handler.particles_in_cell(*cell_neighbor_iterator);

        // Capturing particle pairs, the first particle in the main cell
        // and the second particle in the neighbor cells
        for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
                 iterator particles_in_main_cell_iterator =
                     particles_in_main_cell.begin();
             particles_in_main_cell_iterator != particles_in_main_cell.end();
             ++particles_in_main_cell_iterator) {
          for (typename Particles::ParticleHandler<
                   dim>::particle_iterator_range::iterator
                   particles_in_neighbor_cell_iterator =
                       particles_in_neighbor_cell.begin();
               particles_in_neighbor_cell_iterator !=
               particles_in_neighbor_cell.end();
               ++particles_in_neighbor_cell_iterator) {
            auto contact_pair =
                std::make_pair(particles_in_main_cell_iterator,
                               particles_in_neighbor_cell_iterator);
            local_contact_pair_candidates.insert(
                {std::make_pair(particles_in_main_cell_iterator->get_id(),
                                particles_in_neighbor_cell_iterator->get_id()),
                 contact_pair});
          }
        }
      }
    }
  }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)

  // Since the ghost_contact_pair_candidates (which is the real output of the
  // function) is defined as an input of the function, it should be cleared
  ghost_contact_pair_candidates.clear();

  // Looping over cell_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list->begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list->end();
       ++cell_neighbor_list_iterator) {

    // The main cell
    auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

    // Particles in the main cell
    typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
            particle_handler.particles_in_cell(*cell_neighbor_iterator);

    // If the main cell is not empty
    if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) > 0) {

      // Going through ghost neighbor cells of the main cell
      ++cell_neighbor_iterator;

      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator) {
        // Defining iterator on particles in the neighbor cell
        typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_neighbor_cell =
                particle_handler.particles_in_cell(*cell_neighbor_iterator);

        // Capturing particle pairs, the first particle in the main cell
        // and the second particle in the neighbor cells
        for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
                 iterator particles_in_main_cell_iterator =
                     particles_in_main_cell.begin();
             particles_in_main_cell_iterator != particles_in_main_cell.end();
             ++particles_in_main_cell_iterator) {
          for (typename Particles::ParticleHandler<
                   dim>::particle_iterator_range::iterator
                   particles_in_neighbor_cell_iterator =
                       particles_in_neighbor_cell.begin();
               particles_in_neighbor_cell_iterator !=
               particles_in_neighbor_cell.end();
               ++particles_in_neighbor_cell_iterator) {
            auto contact_pair =
                std::make_pair(particles_in_main_cell_iterator,
                               particles_in_neighbor_cell_iterator);
            ghost_contact_pair_candidates.insert(
                {std::make_pair(particles_in_main_cell_iterator->get_id(),
                                particles_in_neighbor_cell_iterator->get_id()),
                 contact_pair});
          }
        }
      }
    }
  }
}

template class PPBroadSearch<2>;
template class PPBroadSearch<3>;
