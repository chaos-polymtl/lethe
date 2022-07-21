#include <dem/update_particle_wall_contact_container.h>

using namespace dealii;

template <int dim>
void
update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<unsigned int, particle_wall_contact_info_struct<dim>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &particle_container)
{
  for (auto particle_wall_pairs_in_contact_iterator =
         particle_wall_pairs_in_contact.begin();
       particle_wall_pairs_in_contact_iterator !=
       particle_wall_pairs_in_contact.end();
       ++particle_wall_pairs_in_contact_iterator)
    {
      unsigned int particle_id = particle_wall_pairs_in_contact_iterator->first;

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
void update_particle_moving_wall_contact_container_iterators(
  std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<dim>>,
           dem_data_containers::cut_cell_comparison<dim>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &particle_container)
{
  for (auto particle_moving_wall_iterator =
         particle_wall_pairs_in_contact.begin();
       particle_moving_wall_iterator != particle_wall_pairs_in_contact.end();
       ++particle_moving_wall_iterator)
    {
      auto pairs_in_contant_content = &particle_moving_wall_iterator->second;

      for (auto particle_wall_map_iterator = pairs_in_contant_content->begin();
           particle_wall_map_iterator != pairs_in_contant_content->end();
           ++particle_wall_map_iterator)
        {
          unsigned int particle_id = particle_wall_map_iterator->first;

          particle_wall_map_iterator->second.particle =
            particle_container[particle_id];
        }
    }
}

template void
update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<unsigned int, particle_wall_contact_info_struct<2>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &particle_container);

template void
update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<unsigned int, particle_wall_contact_info_struct<3>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &particle_container);

template void update_particle_moving_wall_contact_container_iterators(
  std::map<typename Triangulation<1, 2>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<2>>,
           dem_data_containers::cut_cell_comparison<2>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &particle_container);

template void update_particle_moving_wall_contact_container_iterators(
  std::map<typename Triangulation<2, 3>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<3>>,
           dem_data_containers::cut_cell_comparison<3>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &particle_container);
