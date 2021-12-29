#include <dem/update_local_particle_particle_contact_container.h>

using namespace dealii;

template <int dim>
void
update_local_particle_particle_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
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

template void update_local_particle_particle_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<2>>>
    &local_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &particle_container);

template void update_local_particle_particle_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<3>>>
    &local_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &particle_container);
