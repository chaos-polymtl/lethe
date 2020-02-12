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
void ParticleWallContactDetection<dim, spacedim>::pwFineSearch(
    std::vector<std::tuple<typename Particles::ParticleIterator<dim, spacedim>,
                           Tensor<1, dim>, Point<dim>>>
        pwContactList,
    std::vector<std::tuple<typename Particles::ParticleIterator<dim, spacedim>,
                           Tensor<1, dim>, Point<dim>, double, double, double,
                           Point<dim>, double>> &pwContactInfo,
    float dt) {
  std::tuple<typename Particles::ParticleIterator<dim, spacedim>,
             Tensor<1, dim>, Point<dim>, double, double, double, Point<dim>,
             double>
      pwInfoTuple;
  std::vector<typename Particles::ParticleIterator<dim, spacedim>> pwSearchPair;
  if (!pwSearchPair.empty()) {
    pwSearchPair.clear();
  }
  // reread this part:
  for (unsigned int i = 0; i < pwContactInfo.size(); i++) {
    Point<dim> partCenter = std::get<0>(pwContactInfo[i])->get_location();
    Point<dim> partVector;
    partVector = partCenter - std::get<2>(pwContactInfo[i]);
    Point<dim> projection =
        findProjection(partVector, std::get<1>(pwContactInfo[i]));

    double distance =
        ((std::get<0>(pwContactInfo[i])->get_properties()[2]) / 2) -
        (sqrt(projection.square()));
    if (distance > 0) {
      Tensor<1, dim> normVec = std::get<1>(pwContactInfo[i]);
      Point<dim> pVel = {std::get<0>(pwContactInfo[i])->get_properties()[7],
                         std::get<0>(pwContactInfo[i])->get_properties()[8],
                         std::get<0>(pwContactInfo[i])->get_properties()[9]};
      Point<dim> pOmega = {std::get<0>(pwContactInfo[i])->get_properties()[16],
                           std::get<0>(pwContactInfo[i])->get_properties()[17],
                           std::get<0>(pwContactInfo[i])->get_properties()[18]};
      Point<dim> relVel =
          pVel +
          cross_product_3d(
              (((std::get<0>(pwContactInfo[i])->get_properties()[2]) / 2) *
               pOmega),
              normVec);
      // find how to perform vec = std::get<0>(anothervec)

      double normRelVel = relVel * normVec;
      Tensor<1, dim> relNormVel = normRelVel * normVec;
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
    Point<dim> partCenter = std::get<0>(pwContactList[i])->get_location();
    Point<dim> partVector;
    partVector = partCenter - std::get<2>(pwContactList[i]);
    Point<dim> projection =
        findProjection(partVector, std::get<1>(pwContactList[i]));
    double distance =
        ((std::get<0>(pwContactList[i])->get_properties()[2]) / 2) -
        (sqrt(projection.square()));

    auto it4 = std::find(pwSearchPair.begin(), pwSearchPair.end(),
                         std::get<0>(pwContactList[i]));
    if (it4 == pwSearchPair.end()) {
      if (distance > 0) {
        Tensor<1, dim> normVec = std::get<1>(pwContactList[i]);

        Point<dim> pVel = {std::get<0>(pwContactList[i])->get_properties()[7],
                           std::get<0>(pwContactList[i])->get_properties()[8],
                           std::get<0>(pwContactList[i])->get_properties()[9]};
        Point<dim> pOmega = {
            std::get<0>(pwContactList[i])->get_properties()[16],
            std::get<0>(pwContactList[i])->get_properties()[17],
            std::get<0>(pwContactList[i])->get_properties()[18]};
        Point<dim> relVel =
            pVel +
            cross_product_3d(
                (((std::get<0>(pwContactList[i])->get_properties()[2]) / 2) *
                 pOmega),
                normVec);
        // cross_prod_3d should change for 2d simulations
        double normRelVel = relVel * normVec;
        Tensor<1, dim> relNormVel = normRelVel * normVec;
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

// fix tensor and point
template <int dim, int spacedim>
Point<dim> ParticleWallContactDetection<dim, spacedim>::findProjection(
    Point<dim> pointA, Tensor<1, dim> pointB) {
  Point<dim> PointBprime;
  PointBprime = pointB;
  Point<dim> pointC;
  pointC = ((pointA * PointBprime) / (PointBprime.square())) * PointBprime;

  return pointC;
}

template class ParticleWallContactDetection<3, 3>;
