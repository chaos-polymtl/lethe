#include <dem/update_local_particle_particle_contact_container.h>

using namespace dealii;

template <int dim>
void
update_local_particle_particle_contact_container_iterators(
  typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container)
{
  for (auto adjacent_particles_iterator = local_adjacent_particles.begin();
       adjacent_particles_iterator != local_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      unsigned int particle_one_id  = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto particle_particle_map_iterator =
             pairs_in_contant_content->begin();
           particle_particle_map_iterator != pairs_in_contant_content->end();
           ++particle_particle_map_iterator)
        {
          unsigned int particle_two_id = particle_particle_map_iterator->first;

          particle_particle_map_iterator->second.particle_one =
            particle_container[particle_one_id];
          particle_particle_map_iterator->second.particle_two =
            particle_container[particle_two_id];
        }
    }
}

template void update_local_particle_particle_contact_container_iterators<2>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_local_particle_particle_contact_container_iterators<3>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);
