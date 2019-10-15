/*
 * Integration.h
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include <deal.II/particles/particle_handler.h>
using namespace dealii;

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

class Integration {
public:
	Integration();
	void eulerIntegration(Particles::ParticleHandler<3,3> &, float);
	void velVerIntegration(Particles::ParticleHandler<3,3> &, float);
	void gearIntegration();
};

#endif /* INTEGRATION_H_ */
