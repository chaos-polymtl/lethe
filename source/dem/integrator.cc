/*
 * Integration.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */

#include "dem/integrator.h"

#include <deal.II/particles/particle_handler.h>
using namespace dealii;

template <int dim, int spacedim>
Integrator<dim, spacedim>::Integrator()
{}


// Why are all the numbers for the properties hard-coded? This should
// Come out of an enum class or something like this...

template <int dim, int spacedim>
void
Integrator<dim, spacedim>::eulerIntegration(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  Point<dim>                                 g,
  float                                      dt)
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

template <int dim, int spacedim>
void
Integrator<dim, spacedim>::rk2Integration(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  Point<dim>                                 g,
  float                                      dt)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Acceleration calculation:
      double axStar = particle->get_properties()[10];
      double ayStar = particle->get_properties()[11];
      double azStar = particle->get_properties()[12];

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
      double vxStar = particle->get_properties()[7];
      double vyStar = particle->get_properties()[8];
      double vzStar = particle->get_properties()[9];

      particle->get_properties()[7] =
        particle->get_properties()[7] +
        (dt / 2.0) * (axStar + particle->get_properties()[10]);

      particle->get_properties()[8] =
        particle->get_properties()[8] +
        (dt / 2.0) * (ayStar + particle->get_properties()[11]);
      particle->get_properties()[9] =
        particle->get_properties()[9] +
        (dt / 2.0) * (azStar + particle->get_properties()[12]);

      // Position integration:
      particle->get_properties()[4] =
        particle->get_properties()[4] +
        (dt / 2.0) * (vxStar + particle->get_properties()[7]);
      particle->get_properties()[5] =
        particle->get_properties()[5] +
        (dt / 2.0) * (vyStar + particle->get_properties()[8]);
      particle->get_properties()[6] =
        particle->get_properties()[6] +
        (dt / 2.0) * (vzStar + particle->get_properties()[9]);

      particle->set_location({particle->get_properties()[4],
                              particle->get_properties()[5],
                              particle->get_properties()[6]});

      // Angular velocity using Euler:
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


template <int dim, int spacedim>
void
Integrator<dim, spacedim>::velocity_verlet_integration(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  Point<dim>                                 g,
  float                                      dt)
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

      Point<dim> vel   = {particle->get_properties()[7],
                        particle->get_properties()[8],
                        particle->get_properties()[9]};
      Point<dim> acc   = {particle->get_properties()[10],
                        particle->get_properties()[11],
                        particle->get_properties()[12]};
      Point<dim> vStar = vel + (acc * dt / 2);

      particle->get_properties()[4] =
        particle->get_properties()[4] + (vStar[0] * dt);
      particle->get_properties()[5] =
        particle->get_properties()[5] + (vStar[1] * dt);
      particle->get_properties()[6] =
        particle->get_properties()[6] + (vStar[2] * dt);
      particle->set_location({particle->get_properties()[4],
                              particle->get_properties()[5],
                              particle->get_properties()[6]});

      particle->get_properties()[7] = vStar[0] + (acc[0] * dt / 2);
      particle->get_properties()[8] = vStar[1] + (acc[1] * dt / 2);
      particle->get_properties()[9] = vStar[2] + (acc[2] * dt / 2);

      // Angular velocity using Euler:
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


template class Integrator<3, 3>;
