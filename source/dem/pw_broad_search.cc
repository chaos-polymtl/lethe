#include <dem/pw_broad_search.h>

using namespace dealii;

template <int dim> PWBroadSearch<dim>::PWBroadSearch() {}

template <int dim>
void PWBroadSearch<dim>::find_PW_Contact_Pairs(
    std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
    Particles::ParticleHandler<dim> &particle_handler,
    std::map<int, std::tuple<std::pair<Particles::ParticleIterator<dim>, int>,
                             Tensor<1, dim>, Point<dim>>>
        &pw_contact_candidates) {
  // Since the pw_contact_candidates (which is the real output of the
  // function) is defined as an input of the function, it should be cleared
  pw_contact_candidates.clear();

  // Defining and reseting a local particle-particle candidate counter. This is
  // used as a key to the output map
  int contact_candidate_counter = 0;

  // Iterating over the boundary_cells_information vector which is the output of
  // the find_boundary_cells_information find_boundary_cells_information class.
  // This vector contains all the required information of the system boundary
  // cells and faces. In this loop we find the particles located in each of
  // these boundary cells
  for (auto boundary_cells_information_iterator =
           boundary_cells_information.begin();
       boundary_cells_information_iterator != boundary_cells_information.end();
       ++boundary_cells_information_iterator) {
    // Finding particles located in the corresponding cell
    // (boundary_cells_information_iterator.cell)
    typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(
            boundary_cells_information_iterator->cell);

    for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_cell_iterator = particles_in_cell.begin();
         particles_in_cell_iterator != particles_in_cell.end();
         ++particles_in_cell_iterator) {
      // Making the tuple and adding it to the pw_contact_candidates
      // vector. This vector is the output of this function
      pw_contact_candidates.insert(
          {contact_candidate_counter,
           std::make_tuple(
               std::make_pair(
                   particles_in_cell_iterator,
                   boundary_cells_information_iterator->boundary_face_id),
               boundary_cells_information_iterator->normal_vector,
               boundary_cells_information_iterator->point_on_face)});
      ++contact_candidate_counter;
    }
  }
}

template class PWBroadSearch<2>;
template class PWBroadSearch<3>;
