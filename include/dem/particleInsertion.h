/*
 * ParticleInsertion.h
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */
#include <deal.II/distributed/tria.h>
#include <deal.II/particles/particle_handler.h>

#include "readInputScript.h"


#ifndef PARTICLEINSERTION_H_
#define PARTICLEINSERTION_H_


class ParticleInsertion {
public:

	ParticleInsertion(ReadInputScript readInput);
	void uniformInsertion(dealii::Particles::ParticleHandler<3, 3>&, const dealii::Triangulation<3, 3>&, ReadInputScript readInput, int&, dealii::Particles::PropertyPool&, dealii::Particles::Particle<3>&);


};

#endif /* PARTICLEINSERTION_H_ */
