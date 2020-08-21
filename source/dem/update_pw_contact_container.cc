#include <dem/update_pw_contact_container.h>

using namespace dealii;

template <int dim>
void
update_pw_contact_container_iterators(
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    &pw_pairs_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &particle_container)
{
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact.begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact.end();
       ++pw_pairs_in_contact_iterator)
    {
      int particle_id = pw_pairs_in_contact_iterator->first;

      auto pairs_in_contant_content = &pw_pairs_in_contact_iterator->second;

      for (auto pw_map_iterator = pairs_in_contant_content->begin();
           pw_map_iterator != pairs_in_contant_content->end();
           ++pw_map_iterator)
        {
          pw_map_iterator->second.particle = particle_container.at(particle_id);
        }
    }
}

template void
update_pw_contact_container_iterators(
  std::map<int, std::map<int, pw_contact_info_struct<2>>> &pw_pairs_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<2>>
    &particle_container);

template void
update_pw_contact_container_iterators(
  std::map<int, std::map<int, pw_contact_info_struct<3>>> &pw_pairs_in_contact,
  const std::unordered_map<int, Particles::ParticleIterator<3>>
    &particle_container);
