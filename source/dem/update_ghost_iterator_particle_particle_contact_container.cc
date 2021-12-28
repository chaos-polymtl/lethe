#include <dem/update_ghost_iterator_particle_particle_contact_container.h>

using namespace dealii;

template <int dim>
void
update_ghost_iterator_particle_particle_contact_container(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &ghost_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &ghost_particle_container)
{
  for (auto adjacent_particles_iterator = ghost_adjacent_particles.begin();
       adjacent_particles_iterator != ghost_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto particle_particle_map_iterator =
             pairs_in_contant_content->begin();
           particle_particle_map_iterator != pairs_in_contant_content->end();
           ++particle_particle_map_iterator)
        {
          unsigned int particle_two_id = particle_particle_map_iterator->first;
          particle_particle_map_iterator->second.particle_two =
            ghost_particle_container[particle_two_id];
        }
    }
}

template void update_ghost_iterator_particle_particle_contact_container(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<2>>>
    &ghost_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &ghost_particle_container);

template void update_ghost_iterator_particle_particle_contact_container(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<3>>>
    &ghost_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &ghost_particle_container);
