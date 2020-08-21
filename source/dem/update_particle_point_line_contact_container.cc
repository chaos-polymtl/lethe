#include <dem/update_particle_point_line_contact_container.h>

using namespace dealii;

template <int dim>
void
update_particle_point_line_contact_container_iterators(
  std::map<int, particle_point_line_contact_info_struct<dim>>
    &particle_points_in_contact,
  std::map<int, particle_point_line_contact_info_struct<dim>>
    &particle_lines_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &particle_container)
{
  for (auto particle_point_pairs_in_contact_iterator =
         particle_points_in_contact.begin();
       particle_point_pairs_in_contact_iterator !=
       particle_points_in_contact.end();
       ++particle_point_pairs_in_contact_iterator)
    {
      int  particle_id = particle_point_pairs_in_contact_iterator->first;
      auto pairs_in_contact_content =
        &particle_point_pairs_in_contact_iterator->second;
      pairs_in_contact_content->particle = particle_container.at(particle_id);
    }

  for (auto particle_line_pairs_in_contact_iterator =
         particle_lines_in_contact.begin();
       particle_line_pairs_in_contact_iterator !=
       particle_lines_in_contact.end();
       ++particle_line_pairs_in_contact_iterator)
    {
      int  particle_id = particle_line_pairs_in_contact_iterator->first;
      auto pairs_in_contant_content =
        &particle_line_pairs_in_contact_iterator->second;
      pairs_in_contant_content->particle = particle_container.at(particle_id);
    }
}

template void
update_particle_point_line_contact_container_iterators(
  std::map<int, particle_point_line_contact_info_struct<2>>
    &particle_points_in_contact,
  std::map<int, particle_point_line_contact_info_struct<2>>
    &particle_lines_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<2>>
    &particle_container);

template void
update_particle_point_line_contact_container_iterators(
  std::map<int, particle_point_line_contact_info_struct<3>>
    &particle_points_in_contact,
  std::map<int, particle_point_line_contact_info_struct<3>>
    &particle_lines_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<3>>
    &particle_container);
