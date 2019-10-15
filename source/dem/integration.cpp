/*
 * Integration.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */

#include "integration.h"
#include <deal.II/particles/particle_handler.h>
using namespace dealii;

Integration::Integration() {
}


void Integration::eulerIntegration(Particles::ParticleHandler<3,3> &particle_handler, float dt)
{
	for (auto particle = particle_handler.begin(); particle != particle_handler.end(); ++particle)
	{
		//Velocity integration:
		particle->get_properties()[7] = particle->get_properties()[7] + dt * particle->get_properties()[10];
		particle->get_properties()[8] = particle->get_properties()[8] + dt * particle->get_properties()[11];
		particle->get_properties()[9] = particle->get_properties()[9] + dt * particle->get_properties()[12];
		//Position integration:
		particle->get_properties()[4] = particle->get_properties()[4] + dt * particle->get_properties()[7];
		particle->get_properties()[5] = particle->get_properties()[5] + dt * particle->get_properties()[8];
		particle->get_properties()[6] = particle->get_properties()[6] + dt * particle->get_properties()[9];

		particle->set_location({particle->get_properties()[4], particle->get_properties()[5], particle->get_properties()[6]});
	}
}

/*
void Integration::velVerIntegration(Particles::ParticleHandler<3,3> &particle_handler, float dt)
{
	for (auto particle = particle_handler.begin(); particle != particle_handler.end(); ++particle)
	{
		//Position integration:
		particle->get_properties()[4] = particle->get_properties()[4] + dt * (particle->get_properties()[7] + 0.5 * dt * (particle->get_properties()[10]));
		particle->get_properties()[5] = particle->get_properties()[5] + dt * (particle->get_properties()[8] + 0.5 * dt * (particle->get_properties()[11]));
		particle->get_properties()[6] = particle->get_properties()[6] + dt * (particle->get_properties()[9] + 0.5 * dt * (particle->get_properties()[12]));
		//Velocity integration:
		//******************************************************************
		//here it need the new acceleration which should be updated after the contact force calculation:
		particle->get_properties()[7] = () + 0.5 * dt * ();
		//******************************************************************





		particle->set_location({particle->get_properties()[4], particle->get_properties()[5], particle->get_properties()[6]});
	}

}
*/
