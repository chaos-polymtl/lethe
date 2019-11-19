/*
 * contactForce.h
 *
 *  Created on: Oct 31, 2019
 *      Author: shahab
 */
#include <iostream>
#include <vector>
#include <tuple>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>
#include "readInputScript.h"

using namespace dealii;

#ifndef CONTACTFORCE_H_
#define CONTACTFORCE_H_

class ContactForce {
public:
	ContactForce();
	void linearCF(std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double>> , Particles::ParticleHandler<3,3>&, ReadInputScript);
private:
	double dotProduct(std::vector<double>, std::vector<double>);
	std::vector<double> crossProduct(std::vector<double>, std::vector<double>);
	std::vector<double> vecSubtract(std::vector<double>, std::vector<double>);
	std::vector<double> vecAdd(std::vector<double>, std::vector<double>);
	std::vector<double> numVecProd(double, std::vector<double>);
	double vecValue(std::vector<double>);
	int sgn(float);
};

#endif /* CONTACTFORCE_H_ */
