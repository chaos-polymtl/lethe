/*
 * contactSearch.h
 *
 *  Created on: Oct 29, 2019
 *      Author: shahab
 */
#include <iostream>
#include <vector>
#include <deal.II/particles/particle.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/particles/particle_handler.h>


using namespace dealii;

#ifndef CONTACTSEARCH_H_
#define CONTACTSEARCH_H_

class ContactSearch {
public:
	ContactSearch();
	std::vector<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3, 3>>> findContactPairs(int, dealii::Particles::ParticleHandler<3, 3>&, const Triangulation<3,3>&, std::vector<Triangulation<3>::active_cell_iterator>, std::vector<std::set<Triangulation<3>::active_cell_iterator>>);
	//question? why doesnt work by passing by value?
	std::pair <std::vector<std::set<Triangulation<3>::active_cell_iterator>>,std::vector<Triangulation<3>::active_cell_iterator>> findCellNeighbors(int, const Triangulation<3,3>&);
	void fineSearch(std::vector<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3, 3>>>, dealii::Particles::ParticleHandler<3, 3>&, std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double>> & , float);

	//std::directSearch();

private:
	std::vector<double> normVector(Point<3>, Point<3>);
	double dotProduct(std::vector<double>, std::vector<double>);
	std::vector<double> crossProduct(std::vector<double>, std::vector<double>);
	std::vector<double> vecSubtract(std::vector<double>, std::vector<double>);
	std::vector<double> vecAdd(std::vector<double>, std::vector<double>);
	std::vector<double> numVecProd(double, std::vector<double>);
	double vecValue(std::vector<double>);


};

#endif /* CONTACTSEARCH_H_ */
