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
	void engine(int&, dealii::Particles::ParticleHandler<3, 3>&, const dealii::Triangulation<3, 3>&, int&, float&, ReadInputScript, std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,std::vector<Triangulation<3>::active_cell_iterator>>, std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>, std::vector<double>, double, std::vector<double>,double, std::vector<double>, std::vector<double>, double, double>> &);

private:
	void forceReinit(Particles::ParticleHandler<3,3> &);
	void checkSimBound(Particles::ParticleHandler<3,3> &, ReadInputScript);
};

#endif /* DEMITERATOR_H_ */
