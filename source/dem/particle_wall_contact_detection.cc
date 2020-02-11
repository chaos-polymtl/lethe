/*
 * pwcontactdetection.cpp
 *
 *  Created on: Nov 26, 2019
 *      Author: shahab
 */
#include "dem/particle_wall_contact_detection.h"

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>

using namespace dealii;

template <int dim, int spacedim>
ParticleWallContactDetection<dim, spacedim>::ParticleWallContactDetection() {}

template <int dim, int spacedim>
std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>, Point<dim>>>
ParticleWallContactDetection<dim, spacedim>::pwcontactlist(
    std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
    Particles::ParticleHandler<dim, spacedim> &particle_handler) {
  std::vector<std::tuple<
      std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
      Point<dim>, Point<dim>>>
      pwContactList;
  std::tuple<
      std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
      Point<dim>, Point<dim>>
      pwInfoTuple;
  // std::vector<
  //  std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>>
  //  searchPair;

  for (unsigned int i = 0; i < boundary_cells_information.size(); ++i) {
    typename Triangulation<dim>::active_cell_iterator workingCell =
        (boundary_cells_information[i]).cell;
    typename Particles::ParticleHandler<dim, spacedim>::particle_iterator_range
        particle_range = particle_handler.particles_in_cell(workingCell);

    for (typename Particles::ParticleHandler<
             dim, spacedim>::particle_iterator_range::iterator partIter =
             particle_range.begin();
         partIter != particle_range.end(); ++partIter) {
      auto pwPair =
          std::make_pair(partIter, (boundary_cells_information[i]).boundary_id);

      pwInfoTuple =
          std::make_tuple(pwPair, (boundary_cells_information[i]).normal_vector,
                          (boundary_cells_information[i]).point_on_face);
      pwContactList.push_back(pwInfoTuple);
    }
  }

  return pwContactList;
}

template <int dim, int spacedim>
void ParticleWallContactDetection<dim, spacedim>::pwFineSearch(
    std::vector<std::tuple<
        std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
        Point<dim>, Point<dim>>>
        pwContactList,
    std::vector<std::tuple<
        std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
        Point<dim>, Point<dim>, double, double, double, Point<dim>, double>>
        &pwContactInfo,
    float dt) {
  std::tuple<
      std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
      Point<dim>, Point<dim>, double, double, double, Point<dim>, double>
      pwInfoTuple;
  std::vector<
      std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>>
      pwSearchPair;
  if (!pwSearchPair.empty()) {
    pwSearchPair.clear();
  }
  // reread this part:
  for (unsigned int i = 0; i < pwContactInfo.size(); i++) {
    Point<dim> partCenter = std::get<0>(pwContactInfo[i]).first->get_location();
    Point<dim> partVector;
    partVector = partCenter - std::get<2>(pwContactInfo[i]);
    Point<dim> projection =
        findProjection(partVector, std::get<1>(pwContactInfo[i]));

    double distance =
        ((std::get<0>(pwContactInfo[i]).first->get_properties()[2]) / 2) -
        (sqrt(projection.square()));
    if (distance > 0) {
      Point<dim> normVec = std::get<1>(pwContactInfo[i]);
      Point<dim> pVel = {
          std::get<0>(pwContactInfo[i]).first->get_properties()[7],
          std::get<0>(pwContactInfo[i]).first->get_properties()[8],
          std::get<0>(pwContactInfo[i]).first->get_properties()[9]};
      Point<dim> pOmega = {
          std::get<0>(pwContactInfo[i]).first->get_properties()[16],
          std::get<0>(pwContactInfo[i]).first->get_properties()[17],
          std::get<0>(pwContactInfo[i]).first->get_properties()[18]};
      Point<dim> relVel =
          pVel +
          cross_product_3d(
              (((std::get<0>(pwContactInfo[i]).first->get_properties()[2]) /
                2) *
               pOmega),
              normVec);
      // find how to perform vec = std::get<0>(anothervec)

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
      double tangOverlap = std::get<5>(pwContactInfo[i]) + (tangRelVel * dt);

      pwInfoTuple = std::make_tuple(
          std::get<0>(pwContactInfo[i]), normVec, std::get<2>(pwContactInfo[i]),
          distance, normRelVel, tangOverlap, tangVec, tangRelVel);
      pwContactInfo[i] = pwInfoTuple;
    }

    else {
      pwContactInfo.erase(pwContactInfo.begin() + i);
    }
    // find how to perform vec = std::get<0>(anothervec)
    pwSearchPair.push_back(std::get<0>(pwContactInfo[i]));
  }

  for (unsigned int i = 0; i < pwContactList.size(); i++) {
    Point<dim> partCenter = std::get<0>(pwContactList[i]).first->get_location();
    Point<dim> partVector;
    partVector = partCenter - std::get<2>(pwContactList[i]);
    Point<dim> projection =
        findProjection(partVector, std::get<1>(pwContactList[i]));
    double distance =
        ((std::get<0>(pwContactList[i]).first->get_properties()[2]) / 2) -
        (sqrt(projection.square()));

    auto it4 = std::find(pwSearchPair.begin(), pwSearchPair.end(),
                         std::get<0>(pwContactList[i]));
    if (it4 == pwSearchPair.end()) {
      if (distance > 0) {
        Point<dim> normVec = std::get<1>(pwContactList[i]);

        Point<dim> pVel = {
            std::get<0>(pwContactList[i]).first->get_properties()[7],
            std::get<0>(pwContactList[i]).first->get_properties()[8],
            std::get<0>(pwContactList[i]).first->get_properties()[9]};
        Point<dim> pOmega = {
            std::get<0>(pwContactList[i]).first->get_properties()[16],
            std::get<0>(pwContactList[i]).first->get_properties()[17],
            std::get<0>(pwContactList[i]).first->get_properties()[18]};
        Point<dim> relVel =
            pVel +
            cross_product_3d(
                (((std::get<0>(pwContactList[i]).first->get_properties()[2]) /
                  2) *
                 pOmega),
                normVec);
        // cross_prod_3d should change for 2d simulations
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

        pwInfoTuple =
            std::make_tuple(std::get<0>(pwContactList[i]), normVec,
                            std::get<2>(pwContactList[i]), distance, normRelVel,
                            tangOverlap, tangVec, tangRelVel);

        pwContactInfo.push_back(pwInfoTuple);
      }
    }
  }
}

template <int dim, int spacedim>
Point<dim>
ParticleWallContactDetection<dim, spacedim>::findProjection(Point<dim> pointA,
                                                            Point<dim> pointB) {
  Point<dim> pointC;
  pointC = ((pointA * pointB) / (pointB.square())) * pointB;

  return pointC;
}

template class ParticleWallContactDetection<3, 3>;
