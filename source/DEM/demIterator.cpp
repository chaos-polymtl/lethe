/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */

#include "demIterator.h"

#include "integration.h"
#include "readInputScript.h"
#include "integration.h"
using namespace dealii;

DEM_iterator::DEM_iterator() {
}



void DEM_iterator::engine(Particles::ParticleHandler<3,3> &particle_handler, int &step, float &time, ReadInputScript readInput, Integration  Integ1)
{
	// contact search

	//Integration
	Integ1.eulerIntegration(particle_handler, readInput.dt);

	//***********************************************************************
	//Verlet should be updated after writing the contact force
	//Integ1.velVerIntegration(particle_handler, readInput.dt);
	//***********************************************************************


	step = step + 1;
	time = step * readInput.dt;

}
