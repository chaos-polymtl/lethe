#include <dem/update_particle_wall_contact_container.h>

using namespace dealii;

template <int dim>
void
update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<types::boundary_id, particle_wall_contact_info_struct<dim>>>
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
void update_particle_floating_wall_contact_container_iterators(
  std::vector<
    std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
             std::unordered_map<types::particle_index,
                                particle_wall_contact_info_struct<dim>>,
             dem_data_containers::cut_cell_comparison<dim>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &particle_container)
{
  for (unsigned int solid_counter = 0;
       solid_counter < particle_wall_pairs_in_contact.size();
       ++solid_counter)
    {
      auto &particle_floating_wall_contact_pair =
        particle_wall_pairs_in_contact[solid_counter];
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

template void update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<types::boundary_id, particle_wall_contact_info_struct<2>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &particle_container);

template void update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<types::boundary_id, particle_wall_contact_info_struct<3>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &particle_container);

template void update_particle_floating_wall_contact_container_iterators(
  std::vector<std::map<typename Triangulation<1, 2>::active_cell_iterator,
                       std::unordered_map<types::particle_index,
                                          particle_wall_contact_info_struct<2>>,
                       dem_data_containers::cut_cell_comparison<2>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &particle_container);

template void update_particle_floating_wall_contact_container_iterators(
  std::vector<std::map<typename Triangulation<2, 3>::active_cell_iterator,
                       std::unordered_map<types::particle_index,
                                          particle_wall_contact_info_struct<3>>,
                       dem_data_containers::cut_cell_comparison<3>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &particle_container);
