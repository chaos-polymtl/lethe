#include <dem/update_particle_container.h>

using namespace dealii;

template <int dim>
void
update_particle_container(
  std::unordered_map<int, Particles::ParticleIterator<dim>> &particle_container,
  const Particles::ParticleHandler<dim> *                    particle_handler)
{
  particle_container.clear();

  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }

  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }
}

template void
update_particle_container(
  std::unordered_map<int, Particles::ParticleIterator<2>> &particle_container,
  const Particles::ParticleHandler<2> *                    particle_handler);

template void
update_particle_container(
  std::unordered_map<int, Particles::ParticleIterator<3>> &particle_container,
  const Particles::ParticleHandler<3> *                    particle_handler);
