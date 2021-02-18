#include <dem/pw_broad_search.h>

using namespace dealii;

template <int dim>
PWBroadSearch<dim>::PWBroadSearch()
{}

template <int dim>
void
PWBroadSearch<dim>::find_particle_wall_contact_pairs(
  const std::map<int, boundary_cells_info_struct<dim>>
    &                                    boundary_cells_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  unsigned int>>> &pw_contact_candidates)
{
  // Clearing pw_contact_candidates (output of this function)
  pw_contact_candidates.clear();

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
                 particle_iterator_range::iterator particles_in_cell_iterator =
                   particles_in_cell.begin();
               particles_in_cell_iterator != particles_in_cell.end();
               ++particles_in_cell_iterator)
            {
              // Making the tuple and adding it to the pw_contact_candidates
              // vector. This vector is the output of this function

              std::tuple map_content =
                std::make_tuple(particles_in_cell_iterator,
                                boundary_cells_content.normal_vector,
                                boundary_cells_content.point_on_face,
                                boundary_cells_content.boundary_id);

              pw_contact_candidates[particles_in_cell_iterator->get_id()]
                .insert({boundary_cells_content.boundary_face_id, map_content});
            }
        }
    }
}

template <int dim>
void
PWBroadSearch<dim>::find_particle_floating_wall_contact_pairs(
  const std::unordered_map<
    types::particle_index,
    std::set<typename Triangulation<dim>::active_cell_iterator>>
    &                                    boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double &                                    simulation_time,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>>
    &pfw_contact_candidates)
{
  // Clearing pfw_contact_candidates(output of this function)
  pfw_contact_candidates.clear();

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
      unsigned int floating_wall_id =
        fw_boundary_cells_information_iterator->first;

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
                         particles_in_cell_iterator = particles_in_cell.begin();
                       particles_in_cell_iterator != particles_in_cell.end();
                       ++particles_in_cell_iterator)
                    {
                      pfw_contact_candidates[particles_in_cell_iterator
                                               ->get_id()]
                        .insert({floating_wall_id, particles_in_cell_iterator});
                    }
                }
            }
        }
    }
}

template class PWBroadSearch<2>;
template class PWBroadSearch<3>;
