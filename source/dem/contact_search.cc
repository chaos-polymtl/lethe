/*
 * contactSearch.cpp
 *
 *  Created on: Oct 29, 2019
 *      Author: shahab
 */

#include "dem/contact_search.h"

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

using namespace dealii;

template <int dim, int spacedim>
ContactSearch<dim, spacedim>::ContactSearch() {}

template <int dim, int spacedim>
void ContactSearch<dim, spacedim>::fineSearch(
    std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                          Particles::ParticleIterator<dim, spacedim>>>
        contactPairs,
    std::vector<std::map<int, Particles::ParticleIterator<dim, spacedim>>>
        &inContactPairs,
    std::vector<std::map<int, ContactInfoStruct<dim, spacedim>>> &inContactInfo,
    float dt, Particles::ParticleHandler<dim, spacedim> &particle_handler) {
  Point<dim, double> loc1, loc2;
  double distance;

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end(); ++particle) {
    typename std::map<int, Particles::ParticleIterator<dim, spacedim>>::iterator
        map_it;

    auto temporary_map = inContactPairs[particle->get_id()];
    //    std::cout << "The temporary map size is : " << temporary_map.size()
    //              << std::endl;

    for (map_it = temporary_map.begin(); map_it != temporary_map.end();
         map_it++)
    // while (map_it != inContactPairs[particle->get_id()].end())
    {
      auto particleTwo = map_it->second;

      loc1 = particle->get_location();
      loc2 = particleTwo->get_location();

      distance =
          ((particle->get_properties()[2] + particleTwo->get_properties()[2]) /
           2) -
          loc1.distance(loc2);
      if (distance > 0) {
        Point<dim> contactVector;
        contactVector = (loc2 - loc1);
        Point<dim> normVec = contactVector / contactVector.norm();

        Point<dim> part1Vel = {particle->get_properties()[7],
                               particle->get_properties()[8],
                               particle->get_properties()[9]};
        Point<dim> part2Vel = {particleTwo->get_properties()[7],
                               particleTwo->get_properties()[8],
                               particleTwo->get_properties()[9]};

        Point<dim> part1AngVel = {particle->get_properties()[16],
                                  particle->get_properties()[17],
                                  particle->get_properties()[18]};

        Point<dim> part2AngVel = {particleTwo->get_properties()[16],
                                  particleTwo->get_properties()[17],
                                  particleTwo->get_properties()[18]};

        // ************************************************
        // in tempelate <3> <2> take care
        Point<dim> relVel;
        relVel = (part1Vel - part2Vel) +
                 (cross_product_3d(
                     (((particle->get_properties()[2] / 2.0) * part1AngVel) +
                      ((particleTwo->get_properties()[2] / 2.0) * part2AngVel)),
                     normVec));

        double normRelVel = relVel * normVec;
        Point<dim> relNormVel = normRelVel * normVec;
        Point<dim> relTangVel;
        relTangVel = relVel - relNormVel;
        Point<dim> tangVec = {0, 0, 0};

        double relTangVelVal = relTangVel.norm();
        if (relTangVelVal != 0) {
          tangVec = relTangVel / relTangVelVal;
        }

        double tangRelVel = relVel * tangVec;
        // add this to property of particles // ****************
        double tangOverlap =
            inContactInfo[particle->get_id()][particleTwo->get_id()]
                .tangOverlap +
            (tangRelVel * dt);

        (inContactInfo[particle->get_id()])[particleTwo->get_id()].normOverlap =
            distance;
        (inContactInfo[particle->get_id()])[particleTwo->get_id()].normVec =
            normVec;
        (inContactInfo[particle->get_id()])[particleTwo->get_id()].normRelVel =
            normRelVel;
        (inContactInfo[particle->get_id()])[particleTwo->get_id()].tangVec =
            tangVec;
        (inContactInfo[particle->get_id()])[particleTwo->get_id()].tangRelVel =
            tangRelVel;
        (inContactInfo[particle->get_id()])[particleTwo->get_id()].tangOverlap =
            tangOverlap;

        (inContactInfo[particle->get_id()])[particleTwo->get_id()].particleI =
            particle;
        (inContactInfo[particle->get_id()])[particleTwo->get_id()].particleJ =
            particleTwo;
      } else {
        //(inContactPairs[particle->get_id()]).erase(particleTwo->get_id());auto
        // map_it3 =
        // inContactPairs[(contactPairs[i].second)->get_id()].find(
        //(inContactInfo[particle->get_id()]).erase(particleTwo->get_id());
        //(inContactPairs[particle->get_id()]).erase(map_it);
        (inContactPairs[particle->get_id()]).erase(particleTwo->get_id());
        (inContactInfo[particle->get_id()]).erase(particleTwo->get_id());
      }
      //     searchPair.push_back(std::get<0>(contactInfo[i]));

      // map_it++;
      //  j++;
    }
  }

  typename std::map<int, Particles::ParticleIterator<dim, spacedim>>::iterator
      map_it;
  for (unsigned int i = 0; i < contactPairs.size(); i++) {
    // which one is faster? finding in a map or calculation of distance?

    loc1 = (contactPairs[i].first)->get_location();
    loc2 = (contactPairs[i].second)->get_location();

    distance = ((contactPairs[i].first->get_properties()[2] +
                 contactPairs[i].second->get_properties()[2]) /
                2) -
               loc1.distance(loc2);

    if (distance > 0) {
      auto map_it2 = inContactPairs[(contactPairs[i].first)->get_id()].find(
          (contactPairs[i].second)->get_id());
      auto map_it3 = inContactPairs[(contactPairs[i].second)->get_id()].find(
          (contactPairs[i].first)->get_id());
      if (map_it2 == inContactPairs[(contactPairs[i].first)->get_id()].end() &&
          map_it3 == inContactPairs[(contactPairs[i].second)->get_id()].end()) {
        inContactPairs[(contactPairs[i].first)->get_id()].insert(
            {(contactPairs[i].second)->get_id(), contactPairs[i].second});
        //    std::cout << (contactPairs[i].first)->get_id() << " "
        //            << (contactPairs[i].second)->get_id() << std::endl;

        Point<dim> contactVector;
        contactVector = (loc2 - loc1);
        Point<dim> normVec = contactVector / contactVector.norm();

        Point<dim> part1Vel = {contactPairs[i].first->get_properties()[7],
                               contactPairs[i].first->get_properties()[8],
                               contactPairs[i].first->get_properties()[9]};
        Point<dim> part2Vel = {contactPairs[i].second->get_properties()[7],
                               contactPairs[i].second->get_properties()[8],
                               contactPairs[i].second->get_properties()[9]};

        Point<dim> part1AngVel = {contactPairs[i].first->get_properties()[16],
                                  contactPairs[i].first->get_properties()[17],
                                  contactPairs[i].first->get_properties()[18]};

        Point<dim> part2AngVel = {contactPairs[i].second->get_properties()[16],
                                  contactPairs[i].second->get_properties()[17],
                                  contactPairs[i].second->get_properties()[18]};

        // ************************************************
        // in tempelate <3> <2> take care
        Point<dim> relVel;
        relVel = (part1Vel - part2Vel) +
                 (cross_product_3d(
                     (((contactPairs[i].first->get_properties()[2] / 2.0) *
                       part1AngVel) +
                      ((contactPairs[i].second->get_properties()[2] / 2.0) *
                       part2AngVel)),
                     normVec));

        double normRelVel = relVel * normVec;
        Point<dim> relNormVel = normRelVel * normVec;
        Point<dim> relTangVel;
        relTangVel = relVel - relNormVel;

        Point<dim> tangVec = {0, 0, 0};
        double relTangVelVal = relTangVel.norm();
        if (relTangVelVal != 0) {
          tangVec = relTangVel / relTangVelVal;
        }
        double tangRelVel = relVel * tangVec;
        double tangOverlap = 0;

        ContactInfoStruct<dim, spacedim> contactInfo;
        contactInfo.normOverlap = distance;
        contactInfo.normVec = normVec;
        contactInfo.normRelVel = normRelVel;
        contactInfo.tangVec = tangVec;
        contactInfo.tangRelVel = tangRelVel;
        contactInfo.tangOverlap = tangOverlap;
        contactInfo.particleI = contactPairs[i].first;
        contactInfo.particleJ = contactPairs[i].second;
        inContactInfo[(contactPairs[i].first)->get_id()].insert(
            {(contactPairs[i].second)->get_id(), contactInfo});

        //     std::cout << contactInfo.particleI->get_id() << " "
        //               << contactInfo.particleJ->get_id() << std::endl;
        //     std::cout << "---------------------------------------" <<
        //     std::endl;
      }
    }
  }
}

template class ContactSearch<3, 3>;
