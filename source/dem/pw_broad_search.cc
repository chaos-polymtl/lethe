#include <dem/pw_broad_search.h>

using namespace dealii;

template <int dim> PWBroadSearch<dim>::PWBroadSearch() {}

template <int dim>
std::vector<
    std::tuple<std::pair<typename Particles::ParticleIterator<dim>, int>,
               Tensor<1, dim>, Point<dim>>>
PWBroadSearch<dim>::find_PW_Contact_Pairs(
    std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
    Particles::ParticleHandler<dim> &particle_handler) {

  // A vector of tuples which contains all the canditates for particle-wall
  // collision at each time step. Each tuple contains a collision pair (a
  // particle located near the boundaries, boundary id), the normal vector of
  // the corresponding boundary face and a point on the face
  std::vector<
      std::tuple<std::pair<typename Particles::ParticleIterator<dim>, int>,
                 Tensor<1, dim>, Point<dim>>>
      pw_contact_pair_candidates;

  // Iterating over the boundary_cells_information vector which is the output of
  // the find_boundary_cells_information find_boundary_cells_information class.
  // This vector contains all the required information of the system boundary
  // cells and faces. In this loop we find the particles located in each of
  // these boundary cells
  for (auto boundary_cells_information_iterator =
           boundary_cells_information.begin();
       boundary_cells_information_iterator != boundary_cells_information.end();
       ++boundary_cells_information_iterator) {

    // Fidning particles located in the corresponding cell
    // (boundary_cells_information_iterator.cell)
    typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(
            boundary_cells_information_iterator->cell);

    for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_cell_iterator = particles_in_cell.begin();
         particles_in_cell_iterator != particles_in_cell.end();
         ++particles_in_cell_iterator) {

      // Making the tuple and adding it to the pw_contact_pair_candidates
      // vector. This vector is the output of this function
      pw_contact_pair_candidates.push_back(std::make_tuple(
          std::make_pair(particles_in_cell_iterator,
                         boundary_cells_information_iterator->boundary_id),
          boundary_cells_information_iterator->normal_vector,
          boundary_cells_information_iterator->point_on_face));
    }
  }

  return pw_contact_pair_candidates;
}

template class PWBroadSearch<2>;
template class PWBroadSearch<3>;
