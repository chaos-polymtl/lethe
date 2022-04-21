#include <dem/update_particle_container.h>

using namespace dealii;

template <int dim>
void
update_particle_container(
  std::vector<Particles::ParticleIterator<dim>> &particle_container,
  const Particles::ParticleHandler<dim> *        particle_handler)
{
  particle_container.clear();
  particle_container.resize(particle_handler->n_global_particles());

  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
#if DEAL_II_VERSION_GTE(10, 0, 0)
      particle_container[particle_iterator->get_local_index()] =
#else
      particle_container[particle_iterator->get_id()] =
#endif
        particle_iterator;
    }

  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
#if DEAL_II_VERSION_GTE(10, 0, 0)
      particle_container[particle_iterator->get_local_index()] =
#else
      particle_container[particle_iterator->get_id()] =
#endif
        particle_iterator;
    }
}

template void update_particle_container(
  std::vector<Particles::ParticleIterator<2>> &particle_container,
  const Particles::ParticleHandler<2> *        particle_handler);

template void update_particle_container(
  std::vector<Particles::ParticleIterator<3>> &particle_container,
  const Particles::ParticleHandler<3> *        particle_handler);
