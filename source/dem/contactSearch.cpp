/*
 * contactSearch.cpp
 *
 *  Created on: Oct 29, 2019
 *      Author: shahab
 */

#include "contactSearch.h"
#include <iostream>
#include <vector>
#include <deal.II/particles/particle.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/particle_handler.h>

using namespace dealii;

ContactSearch::ContactSearch() {
}



std::pair <std::vector<std::set<Triangulation<3>::active_cell_iterator>> , std::vector<Triangulation<3>::active_cell_iterator>>   ContactSearch::findCellNeighbors(int cellNum, const Triangulation<3,3> &tr)
{
		std::vector<std::set<Triangulation<3>::active_cell_iterator>> cellNeighborList(cellNum);
		std::vector<Triangulation<3>::active_cell_iterator> totallCellList;

		int iter = 0;
		auto v_to_c   = GridTools::vertex_to_cell_map(tr);
		for (Triangulation<3>::active_cell_iterator cell = tr.begin_active(); cell != tr.end(); ++cell)
		{
			cellNeighborList[iter].insert(cell);
			totallCellList.push_back(cell);

			for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_cell; ++v)
			{
				for (const auto &neighbor: v_to_c[cell->vertex_index(v)])
				{
					auto it = std::find( totallCellList.begin(), totallCellList.end() , neighbor);
					if(it == totallCellList.end())
						cellNeighborList[iter].insert(neighbor);
				}
			}
			iter++;
		}
		std::pair <std::vector<std::set<Triangulation<3>::active_cell_iterator>> , std::vector<Triangulation<3>::active_cell_iterator>> out = std::make_pair(cellNeighborList , totallCellList);
		return out;
	}





std::vector<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3, 3>>> ContactSearch::findContactPairs(int nPart, Particles::ParticleHandler<3,3> &particle_handler, const Triangulation<3,3> &tr, std::vector<Triangulation<3>::active_cell_iterator> totallCellList, std::vector<std::set<Triangulation<3>::active_cell_iterator>> cellNeighborList)
{
	std::vector<std::vector<Particles::ParticleIterator<3, 3>>> contactList(nPart+1, std::vector<Particles::ParticleIterator<3, 3>> ());
	std::vector<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3, 3>>> contactPairs;
	if (!contactPairs.empty())
	{
		contactPairs.clear();
		std::cout<<"wtf!"<<std::endl;
	}

	Triangulation<3>::active_cell_iterator currrentCell;
	for (auto particleIter=particle_handler.begin(); particleIter != particle_handler.end(); ++particleIter)
	{

		currrentCell = GridTools::find_active_cell_around_point(tr, particleIter->get_location());
		auto it1 = std::find(totallCellList.begin(), totallCellList.end(), currrentCell );
		int index = std::distance(totallCellList.begin(), it1);

		for(auto cellIt=cellNeighborList[index].begin(); cellIt!=cellNeighborList[index].end(); cellIt++)
		{
			const Particles::ParticleHandler<3,3>::particle_iterator_range particle_range = particle_handler.particles_in_cell(*cellIt);
			for (typename Particles::ParticleHandler<3,3>::particle_iterator_range::iterator partIter =particle_range.begin(); partIter != particle_range.end();++partIter)
			{
				contactList[particleIter->get_id()].push_back(partIter);
				auto cPair = std::make_pair(particleIter,partIter);
				auto cPair2 = std::make_pair(partIter,particleIter);
				auto it2 = std::find( contactPairs.begin(), contactPairs.end() , cPair);
				auto it3 = std::find( contactPairs.begin(), contactPairs.end() , cPair2);
				if(it2 == contactPairs.end())
					if(it3 == contactPairs.end())
						if(particleIter->get_id() != partIter->get_id())
							contactPairs.push_back(cPair);
			}
		}
	}
	return contactPairs;
	}


void ContactSearch::fineSearch(std::vector<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>> contactPairs, dealii::Particles::ParticleHandler<3,3> &particle_handler, std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3, 3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double>> &contactInfo, float dt)
		{
		Point<3, double> v1, v2;
		std::tuple<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3,3>>, std::vector<double>, double, std::vector<double>, double, std::vector<double>, std::vector<double>, double, double > infoTuple;
		std::vector<std::pair<Particles::ParticleIterator<3,3>, Particles::ParticleIterator<3,3>>> searchPair;
		//std::vector<std::tuple<std::pair<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3, 3>>, double >> lastStepInfo;
		//std::vector<std::tuple<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3,3>,double, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double> >> contactInfo;
		double distance;
		if (!searchPair.empty())
		{
			searchPair.clear();
			std::cout<<"wtf!"<<std::endl;
		}



		for (unsigned int i = 0; i < contactInfo.size(); i++)
		{
			v1 = {std::get<0>(contactInfo[i]).first->get_properties()[4] ,std::get<0>(contactInfo[i]).first->get_properties()[5] ,std::get<0>(contactInfo[i]).first->get_properties()[6]};
			v2 = {std::get<0>(contactInfo[i]).second->get_properties()[4],std::get<0>(contactInfo[i]).second->get_properties()[5],std::get<0>(contactInfo[i]).second->get_properties()[6]};

			distance = ((std::get<0>(contactInfo[i]).first->get_properties()[2] + std::get<0>(contactInfo[i]).first->get_properties()[2])/2) - v1.distance(v2);
				if(distance > 0)
				{


					std::vector<double> normVec = normVector({contactPairs[i].second->get_properties()[4] ,contactPairs[i].second->get_properties()[5] ,contactPairs[i].second->get_properties()[6]} , {contactPairs[i].first->get_properties()[4],contactPairs[i].first->get_properties()[5],contactPairs[i].first->get_properties()[6]});
					std::vector<double> relVel = vecAdd(vecSubtract({contactPairs[i].first->get_properties()[7], contactPairs[i].first->get_properties()[8], contactPairs[i].first->get_properties()[9]} , {contactPairs[i].second->get_properties()[7], contactPairs[i].second->get_properties()[8], contactPairs[i].second->get_properties()[9]}) , (crossProduct(vecAdd(numVecProd(contactPairs[i].first->get_properties()[2] ,{contactPairs[i].first->get_properties()[16], contactPairs[i].first->get_properties()[17], contactPairs[i].first->get_properties()[18]} ) , numVecProd( contactPairs[i].second->get_properties()[2],{contactPairs[i].second->get_properties()[16], contactPairs[i].second->get_properties()[17], contactPairs[i].second->get_properties()[18]} )) , normVec)));
					double normRelVel = dotProduct(relVel,normVec);
					std::vector<double> relNormVel = numVecProd(normRelVel , normVec);
					std::vector<double> relTangVel = vecSubtract(relVel , relNormVel);

					std::vector<double> tangVec = {0, 0, 0};
					if(vecValue(relTangVel) != 0)
					{
					std::vector<double> tangVec = { relTangVel[0]/vecValue(relTangVel) , relTangVel[1]/vecValue(relTangVel), relTangVel[2]/vecValue(relTangVel)};
					}

					double tangRelVel = dotProduct(relVel , tangVec);
					double tangOverlap = std::get<8>(contactInfo[i]) + ((dotProduct(std::get<1>(contactInfo[i]) , std::get<6>(contactInfo[i]))) * dt);

					infoTuple = std::make_tuple(std::make_pair(contactPairs[i].first, contactPairs[i].second), relVel, distance, normVec, normRelVel, relTangVel, tangVec, tangRelVel, tangOverlap);
					contactInfo[i] = infoTuple;



					//std::cout<< "the first particle is: " << (contactInfo[i]).first->get_properties()[0] << " and second particle: " << (contactInfo[i]).second->get_properties()[0]<<std::endl;
					std :: cout << " contact information2: " << "rel vel: " << relVel[0] << " " <<  relVel[1] << " " << relVel[2] << std::endl;
					std :: cout << " contact information2: " << "distance: " << distance << std::endl;
					std :: cout << " contact information2: " << "normVec: " << normVec[0] << " " <<  normVec[1] << " " << normVec[2] << std::endl;
					std :: cout << " contact information2: " << "normRelVel: " << normRelVel << std::endl;

					std :: cout << " contact information2: " << "relTangVel: " << relTangVel[0] << " " <<  relTangVel[1] << " " << relTangVel[2] << std::endl;
					std :: cout << " contact information2: " << "tangRelVel: " << tangRelVel << std::endl;
					std :: cout << " contact information2: " << "tangVec: " << tangVec[0] << " " <<  tangVec[1] << " " << tangVec[2] << std::endl;
					std :: cout << " contact information2: " << "tangOverlap: " << tangOverlap << std::endl;


				}
				else
				{

					std::cout<< "heyyyy erasing!!!!!!"<<std::endl;
					contactInfo.erase(contactInfo.begin()+i);

				}
				searchPair.push_back(std::get<0>(contactInfo[i]));


		}




		for (unsigned int i = 0; i < contactPairs.size(); i++)
			{
				v1 = {contactPairs[i].first->get_properties()[4] ,contactPairs[i].first->get_properties()[5] ,contactPairs[i].first->get_properties()[6]};
				v2 = {contactPairs[i].second->get_properties()[4],contactPairs[i].second->get_properties()[5],contactPairs[i].second->get_properties()[6]};
				distance = ((contactPairs[i].first->get_properties()[2] + contactPairs[i].second->get_properties()[2])/2) - v1.distance(v2);


				auto it4 = std::find(searchPair.begin(), searchPair.end(), contactPairs[i]);
				if(it4 == searchPair.end())
				{
				if(distance > 0)
					{
						std::vector<double> normVec = normVector({contactPairs[i].second->get_properties()[4] ,contactPairs[i].second->get_properties()[5] ,contactPairs[i].second->get_properties()[6]} , {contactPairs[i].first->get_properties()[4],contactPairs[i].first->get_properties()[5],contactPairs[i].first->get_properties()[6]});
						std::vector<double> relVel = vecAdd(vecSubtract({contactPairs[i].first->get_properties()[7], contactPairs[i].first->get_properties()[8], contactPairs[i].first->get_properties()[9]} , {contactPairs[i].second->get_properties()[7], contactPairs[i].second->get_properties()[8], contactPairs[i].second->get_properties()[9]}) , (crossProduct(vecAdd(numVecProd(contactPairs[i].first->get_properties()[2] ,{contactPairs[i].first->get_properties()[16], contactPairs[i].first->get_properties()[17], contactPairs[i].first->get_properties()[18]} ) , numVecProd( contactPairs[i].second->get_properties()[2],{contactPairs[i].second->get_properties()[16], contactPairs[i].second->get_properties()[17], contactPairs[i].second->get_properties()[18]} )) , normVec)));
						double normRelVel = dotProduct(relVel,normVec);
						std::vector<double> relNormVel = numVecProd(normRelVel , normVec);
						std::vector<double> relTangVel = vecSubtract(relVel , relNormVel);
						std::vector<double> tangVec = {0, 0, 0};
							if(vecValue(relTangVel) != 0)
							{
							tangVec = { relTangVel[0]/vecValue(relTangVel) , relTangVel[1]/vecValue(relTangVel), relTangVel[2]/vecValue(relTangVel)};
							}
						double tangRelVel = dotProduct(relVel , tangVec);
						double tangOverlap = 0;

						infoTuple = std::make_tuple(std::make_pair(contactPairs[i].first, contactPairs[i].second), relVel, distance, normVec, normRelVel, relTangVel, tangVec, tangRelVel, tangOverlap);
						contactInfo.push_back(infoTuple);

						std::cout<< "the first particle is: " << contactPairs[i].first->get_properties()[0] << " and second particle: " << contactPairs[i].second->get_properties()[0]<<std::endl;
						std :: cout << " contact information: " << "rel vel: " << relVel[0] << " " <<  relVel[1] << " " << relVel[2] << std::endl;
						std :: cout << " contact information: " << "distance: " << distance << std::endl;
						std :: cout << " contact information: " << "normVec: " << normVec[0] << " " <<  normVec[1] << " " << normVec[2] << std::endl;
						std :: cout << " contact information: " << "normRelVel: " << normRelVel << std::endl;

						std :: cout << " contact information: " << "relTangVel: " << relTangVel[0] << " " <<  relTangVel[1] << " " << relTangVel[2] << std::endl;
						std :: cout << " contact information: " << "tangRelVel: " << tangRelVel << std::endl;
						std :: cout << " contact information: " << "tangVec: " << tangVec[0] << " " <<  tangVec[1] << " " << tangVec[2] << std::endl;
						std :: cout << " contact information: " << "tangOverlap: " << tangOverlap << std::endl;




					}
				}



			}


		}




std::vector<double> ContactSearch::normVector(Point<3> point1, Point<3> point2)
		{
			Point<3> point3;
			point3 = (point1 - point2) / (point1.distance(point2));
			return {point3[0], point3[1], point3[2]};
		}


double ContactSearch::dotProduct(std::vector<double> A, std::vector<double> B)
{
	return (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]);
	}

std::vector<double> ContactSearch::crossProduct(std::vector<double> A, std::vector<double> B)
		{
			return {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]};
		}

std::vector<double> ContactSearch::vecSubtract(std::vector<double> A, std::vector<double> B)
		{
			return {A[0] - B[0], A[1] - B[1], A[2] - B[2]};
		}

std::vector<double> ContactSearch::vecAdd(std::vector<double> A, std::vector<double> B)
{
	return {A[0] + B[0], A[1] + B[1], A[2] + B[2]};
	}

std::vector<double> ContactSearch::numVecProd(double A, std::vector<double> B)
		{
	return {A * B[0], A * B[1], A * B[2]};
		}

double ContactSearch::vecValue(std::vector<double> A)
{
	return (sqrt(pow(A[0],2)+pow(A[1],2)+pow(A[2],2)));
	}



