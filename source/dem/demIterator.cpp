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
#include "contactSearch.h"
#include "contactForce.h"
using namespace dealii;

DEM_iterator::DEM_iterator() {
}



void DEM_iterator::engine(int &nPart, Particles::ParticleHandler<3,3> &particle_handler, const Triangulation<3,3> &tr, int &step, float &time, ReadInputScript readInput, std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,std::vector<Triangulation<3>::active_cell_iterator>> cellNeighbor, std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double>> &contactInfo)
{

	// contact search
	std::vector<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3, 3>>> contactPairs;
	//std::vector<std::tuple<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3,3>,double,std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> >> contactInfo;

	ContactSearch cs;
	contactPairs = cs.findContactPairs(nPart, particle_handler,tr, cellNeighbor.second, cellNeighbor.first);
	cs.fineSearch(contactPairs, particle_handler, contactInfo, readInput.dt);

	// contact force
	ContactForce cf;
	cf.linearCF(contactInfo, particle_handler, readInput);





	//Integration
	Integration	 Integ1;
	Integ1.eulerIntegration(particle_handler, readInput.dt);

	//***********************************************************************
	//Verlet should be updated after writing the contact force
	//Integ1.velVerIntegration(particle_handler, readInput.dt);
	//***********************************************************************


	step = step + 1;
	time = step * readInput.dt;
}
