#include <dem/update_particle_point_line_contact_container.h>

using namespace dealii;

template <int dim>
void
update_particle_point_line_contact_container_iterators(
  typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
    &particle_points_in_contact,
  typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
    &particle_lines_in_contact,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container)
{
  for (auto particle_point_pairs_in_contact_iterator =
         particle_points_in_contact.begin();
       particle_point_pairs_in_contact_iterator !=
       particle_points_in_contact.end();
       ++particle_point_pairs_in_contact_iterator)
    {
      unsigned int particle_id =
        particle_point_pairs_in_contact_iterator->first;
      auto pairs_in_contact_content =
        &particle_point_pairs_in_contact_iterator->second;
      pairs_in_contact_content->particle = particle_container[particle_id];
    }

  for (auto particle_line_pairs_in_contact_iterator =
         particle_lines_in_contact.begin();
       particle_line_pairs_in_contact_iterator !=
       particle_lines_in_contact.end();
       ++particle_line_pairs_in_contact_iterator)
    {
      unsigned int particle_id = particle_line_pairs_in_contact_iterator->first;
      auto         pairs_in_contant_content =
        &particle_line_pairs_in_contact_iterator->second;
      pairs_in_contant_content->particle = particle_container[particle_id];
    }
}

template void update_particle_point_line_contact_container_iterators<2>(
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &particle_points_in_contact,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &particle_lines_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_particle_point_line_contact_container_iterators<3>(
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &particle_points_in_contact,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &particle_lines_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);
