/*
 * DEMiterator.h
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include "readInputScript.h"
#include "integration.h"

#ifndef DEMITERATOR_H_
#define DEMITERATOR_H_

class DEM_iterator {
public:
	DEM_iterator();
	void engine(dealii::Particles::ParticleHandler<3, 3>&, int&, float&, ReadInputScript, Integration);
};

#endif /* DEMITERATOR_H_ */
