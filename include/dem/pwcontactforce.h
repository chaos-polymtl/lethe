/*
 * pwcontactforce.h
 *
 *  Created on: Dec 5, 2019
 *      Author: shahab
 */
#include <iostream>
#include <vector>
#include <tuple>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include "readinputscript.h"

using namespace dealii;

#ifndef PWCONTACTFORCE_H_
#define PWCONTACTFORCE_H_

class PWContactForce {
public:
	PWContactForce();
	void pwLinearCF(std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, int>, Point<3>, Point<3>, double, double, double, Point<3>, double >>, Particles::ParticleHandler<3,3>&, ReadInputScript);

private:
	double vecValue(Point<3>);
	int sgn(float);
};

#endif /* PWCONTACTFORCE_H_ */
