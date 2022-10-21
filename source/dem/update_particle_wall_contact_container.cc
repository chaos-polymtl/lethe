#include <dem/update_particle_wall_contact_container.h>

using namespace dealii;

template <int dim>
void
update_particle_wall_contact_container_iterators(
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container)
{
  for (auto particle_wall_pairs_in_contact_iterator =
         particle_wall_pairs_in_contact.begin();
       particle_wall_pairs_in_contact_iterator !=
       particle_wall_pairs_in_contact.end();
       ++particle_wall_pairs_in_contact_iterator)
    {
      auto particle_id = particle_wall_pairs_in_contact_iterator->first;

      auto pairs_in_contant_content =
        &particle_wall_pairs_in_contact_iterator->second;

      for (auto particle_wall_map_iterator = pairs_in_contant_content->begin();
           particle_wall_map_iterator != pairs_in_contant_content->end();
           ++particle_wall_map_iterator)
        {
          particle_wall_map_iterator->second.particle =
            particle_container[particle_id];
        }
    }
}

template <int dim>
void
update_particle_floating_mesh_contact_container_iterators(
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
    &particle_floating_mesh_in_contact,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container)
{
  for (unsigned int solid_counter = 0;
       solid_counter < particle_floating_mesh_in_contact.size();
       ++solid_counter)
    {
      auto &particle_floating_wall_contact_pair =
        particle_floating_mesh_in_contact[solid_counter];
      for (auto particle_floating_wall_iterator =
             particle_floating_wall_contact_pair.begin();
           particle_floating_wall_iterator !=
           particle_floating_wall_contact_pair.end();
           ++particle_floating_wall_iterator)
        {
          auto pairs_in_contant_content =
            &particle_floating_wall_iterator->second;

          for (auto particle_wall_map_iterator =
                 pairs_in_contant_content->begin();
               particle_wall_map_iterator != pairs_in_contant_content->end();
               ++particle_wall_map_iterator)
            {
              auto particle_id = particle_wall_map_iterator->first;

              particle_wall_map_iterator->second.particle =
                particle_container[particle_id];
            }
        }
    }
}

template void update_particle_wall_contact_container_iterators<2>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_particle_wall_contact_container_iterators<3>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

template void update_particle_floating_mesh_contact_container_iterators<2>(
  typename DEM::dem_data_structures<2>::particle_floating_mesh_in_contact
    &particle_floating_mesh_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_particle_floating_mesh_contact_container_iterators<3>(
  typename DEM::dem_data_structures<3>::particle_floating_mesh_in_contact
    &particle_floating_mesh_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);
