/*
 * Integration.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */

#include "dem/integration.h"

#include <deal.II/particles/particle_handler.h>
using namespace dealii;

template <int dim, int spacedim>
Integration<dim, spacedim>::Integration()
{}

template <int dim, int spacedim>
void
Integration<dim, spacedim>::eulerIntegration(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  ParametersDEM<dim>                         DEMparam)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Acceleration calculation:
      particle->get_properties()[10] =
        DEMparam.physicalProperties.gx +
        (particle->get_properties()[13]) / (particle->get_properties()[19]);
      particle->get_properties()[11] =
        DEMparam.physicalProperties.gy +
        (particle->get_properties()[14]) / (particle->get_properties()[19]);
      particle->get_properties()[12] =
        DEMparam.physicalProperties.gz +
        (particle->get_properties()[15]) / (particle->get_properties()[19]);

      // Velocity integration:
      particle->get_properties()[7] =
        particle->get_properties()[7] +
        DEMparam.simulationControl.dt * particle->get_properties()[10];
      particle->get_properties()[8] =
        particle->get_properties()[8] +
        DEMparam.simulationControl.dt * particle->get_properties()[11];
      particle->get_properties()[9] =
        particle->get_properties()[9] +
        DEMparam.simulationControl.dt * particle->get_properties()[12];
      // Position integration:
      particle->get_properties()[4] =
        particle->get_properties()[4] +
        DEMparam.simulationControl.dt * particle->get_properties()[7];
      particle->get_properties()[5] =
        particle->get_properties()[5] +
        DEMparam.simulationControl.dt * particle->get_properties()[8];
      particle->get_properties()[6] =
        particle->get_properties()[6] +
        DEMparam.simulationControl.dt * particle->get_properties()[9];

      particle->set_location({particle->get_properties()[4],
                              particle->get_properties()[5],
                              particle->get_properties()[6]});
    }
}

template <int dim, int spacedim>
void
Integration<dim, spacedim>::rk2Integration(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  ParametersDEM<dim>                         DEMparam)
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
        DEMparam.physicalProperties.gx +
        (particle->get_properties()[13]) / (particle->get_properties()[19]);
      particle->get_properties()[11] =
        DEMparam.physicalProperties.gy +
        (particle->get_properties()[14]) / (particle->get_properties()[19]);
      particle->get_properties()[12] =
        DEMparam.physicalProperties.gz +
        (particle->get_properties()[15]) / (particle->get_properties()[19]);


      // Velocity integration:
      double vxStar = particle->get_properties()[7];
      double vyStar = particle->get_properties()[8];
      double vzStar = particle->get_properties()[9];

      particle->get_properties()[7] =
        particle->get_properties()[7] +
        (DEMparam.simulationControl.dt / 2.0) *
          (axStar + particle->get_properties()[10]);

      particle->get_properties()[8] =
        particle->get_properties()[8] +
        (DEMparam.simulationControl.dt / 2.0) *
          (ayStar + particle->get_properties()[11]);
      particle->get_properties()[9] =
        particle->get_properties()[9] +
        (DEMparam.simulationControl.dt / 2.0) *
          (azStar + particle->get_properties()[12]);

      // Position integration:
      particle->get_properties()[4] =
        particle->get_properties()[4] +
        (DEMparam.simulationControl.dt / 2.0) *
          (vxStar + particle->get_properties()[7]);
      particle->get_properties()[5] =
        particle->get_properties()[5] +
        (DEMparam.simulationControl.dt / 2.0) *
          (vyStar + particle->get_properties()[8]);
      particle->get_properties()[6] =
        particle->get_properties()[6] +
        (DEMparam.simulationControl.dt / 2.0) *
          (vzStar + particle->get_properties()[9]);

      particle->set_location({particle->get_properties()[4],
                              particle->get_properties()[5],
                              particle->get_properties()[6]});
    }
}



/*
void Integration::velVerIntegration(Particles::ParticleHandler<3,3>
&particle_handler, float dt)
{
  for (auto particle = particle_handler.begin(); particle !=
particle_handler.end(); ++particle)
  {
    //Position integration:
    particle->get_properties()[4] = particle->get_properties()[4] + dt *
(particle->get_properties()[7] + 0.5 * dt * (particle->get_properties()[10]));
    particle->get_properties()[5] = particle->get_properties()[5] + dt *
(particle->get_properties()[8] + 0.5 * dt * (particle->get_properties()[11]));
    particle->get_properties()[6] = particle->get_properties()[6] + dt *
(particle->get_properties()[9] + 0.5 * dt * (particle->get_properties()[12]));
    /Velocity integration:
    /here it need the new acceleration which should be updated after the
contact force calculation: particle->get_properties()[7] = () + 0.5 * dt * ();


    particle->set_location({particle->get_properties()[4],
particle->get_properties()[5], particle->get_properties()[6]});
  }

}
*/

template class Integration<3, 3>;
