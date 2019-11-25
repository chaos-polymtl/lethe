/*
 * DEMiterator.cpp
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include "dem/demIterator.h"
#include "dem/contactForce.h"
#include "dem/contactSearch.h"
#include "dem/integration.h"
#include "dem/readInputScript.h"
using namespace dealii;

DEM_iterator::DEM_iterator() {}

void DEM_iterator::forceReinit(
    Particles::ParticleHandler<3, 3> &particle_handler) {
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {
    particle->get_properties()[13] = 0;
    particle->get_properties()[14] = 0;
    particle->get_properties()[15] = 0;
  }
}

void DEM_iterator::checkSimBound(
    Particles::ParticleHandler<3, 3> &particle_handler,
    ReadInputScript readInput) {
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {
    if (particle->get_properties()[4] < readInput.x_min ||
        particle->get_properties()[4] > readInput.x_max ||
        particle->get_properties()[5] < readInput.y_min ||
        particle->get_properties()[5] > readInput.y_max ||
        particle->get_properties()[6] < readInput.z_min ||
        particle->get_properties()[6] > readInput.z_max) {
      particle_handler.remove_particle(particle);
    }
  }
}

void DEM_iterator::engine(
    int &nPart, Particles::ParticleHandler<3, 3> &particle_handler,
    const Triangulation<3, 3> &tr, int &step, float &time,
    ReadInputScript readInput,
    std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,
              std::vector<Triangulation<3>::active_cell_iterator>>
        cellNeighbor,
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                                     Particles::ParticleIterator<3, 3>>,
                           std::vector<double>, double, std::vector<double>,
                           double, std::vector<double>, std::vector<double>,
                           double, double>> &contactInfo) {

  // check simulation boundaries
  checkSimBound(particle_handler, readInput);

  // contact search
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>>
      contactPairs;
  // std::vector<std::tuple<Particles::ParticleIterator<3,3>,Particles::ParticleIterator<3,3>,double,std::vector<double>,
  // std::vector<double>, std::vector<double>, std::vector<double> >>
  // contactInfo;

  // force reinitilization
  forceReinit(particle_handler);

  ContactSearch cs;
  contactPairs = cs.findContactPairs(nPart, particle_handler, tr,
                                     cellNeighbor.second, cellNeighbor.first);

  cs.fineSearch(contactPairs, particle_handler, contactInfo, readInput.dt);

  // contact force
  ContactForce cf;
  cf.linearCF(contactInfo, particle_handler, readInput);

  // for (auto particle = particle_handler.begin(); particle !=
  // particle_handler.end(); ++particle)
  //	{
  //	std:: cout << particle->get_properties()[13]<< " " <<
  //particle->get_properties()[14] << " " << particle->get_properties()[15] <<
  //std::endl;
  //	}

  // Integration
  Integration Integ1;
  Integ1.eulerIntegration(particle_handler, readInput);

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {
    std::cout << "force: " << particle->get_properties()[13] << " "
              << particle->get_properties()[14] << " "
              << particle->get_properties()[15] << std::endl;

    std::cout << "acceleration: " << particle->get_properties()[10] << " "
              << particle->get_properties()[11] << " "
              << particle->get_properties()[12] << std::endl;
    std::cout << "velocity: " << particle->get_properties()[7] << " "
              << particle->get_properties()[8] << " "
              << particle->get_properties()[9] << std::endl;
  }

  //***********************************************************************
  // Verlet should be updated after writing the contact force
  // Integ1.velVerIntegration(particle_handler, readInput.dt);
  //***********************************************************************

  step = step + 1;
  time = step * readInput.dt;
}
