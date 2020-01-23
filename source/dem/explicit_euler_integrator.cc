#include <dem/explicit_euler_integrator.h>

template <int dim, int spacedim>
void
ExplicitEulerIntegrator<dim, spacedim>::integrate(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  Point<dim>                                 g,
  double                                     dt)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Acceleration calculation:
      particle->get_properties()[10] =
        g[0] +
        (particle->get_properties()[13]) / (particle->get_properties()[19]);
      particle->get_properties()[11] =
        g[1] +
        (particle->get_properties()[14]) / (particle->get_properties()[19]);
      particle->get_properties()[12] =
        g[2] +
        (particle->get_properties()[15]) / (particle->get_properties()[19]);

      // Velocity integration:
      particle->get_properties()[7] =
        particle->get_properties()[7] + dt * particle->get_properties()[10];
      particle->get_properties()[8] =
        particle->get_properties()[8] + dt * particle->get_properties()[11];
      particle->get_properties()[9] =
        particle->get_properties()[9] + dt * particle->get_properties()[12];
      // Position integration:
      particle->get_properties()[4] =
        particle->get_properties()[4] + dt * particle->get_properties()[7];
      particle->get_properties()[5] =
        particle->get_properties()[5] + dt * particle->get_properties()[8];
      particle->get_properties()[6] =
        particle->get_properties()[6] + dt * particle->get_properties()[9];

      particle->set_location({particle->get_properties()[4],
                              particle->get_properties()[5],
                              particle->get_properties()[6]});
      // Angular velocity:
      /*
      particle->get_properties()[16] = particle->get_properties()[16] +
        (particle->get_properties()[21]) / (particle->get_properties()[20]);
      particle->get_properties()[17] =particle->get_properties()[17] +
        (particle->get_properties()[22]) / (particle->get_properties()[20]);
      particle->get_properties()[18] = particle->get_properties()[18] +
        (particle->get_properties()[23]) / (particle->get_properties()[20]);
        */
    }
}

template class ExplicitEulerIntegrator<3, 3>;
