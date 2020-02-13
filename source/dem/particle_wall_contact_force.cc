/*
 * pwcontactforce.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: shahab
 */

#include "dem/particle_wall_contact_force.h"

#include <deal.II/base/point.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

#include "dem/dem_iterator.h"

template <int dim, int spacedim>
ParticleWallContactForce<dim, spacedim>::ParticleWallContactForce() {}

template <int dim, int spacedim>
void ParticleWallContactForce<dim, spacedim>::pwLinearCF(
    std::vector<
        std::map<int, std::tuple<Particles::ParticleIterator<dim, spacedim>,
                                 Tensor<1, dim>, Point<dim>, double, double,
                                 double, Tensor<1, dim>, double>>>
        &pwContactInfo,
    physical_info_struct<dim> &physical_info_struct) {

  for (auto vector_iterator = pwContactInfo.begin();
       vector_iterator != pwContactInfo.end(); ++vector_iterator) {
    // for (unsigned int i = 0; i < pwContactInfo.size(); i++) {
    for (auto map_iterator = vector_iterator->begin();
         map_iterator != vector_iterator->end(); ++map_iterator) {
      Point<dim> totalForce;
      double yEff =
          pow((((1.0 - pow(physical_info_struct.Poisson_ratio_particle, 2.0)) /
                physical_info_struct.Young_modulus_particle) +
               ((1.0 - pow(physical_info_struct.Poisson_ratio_wall, 2.0)) /
                physical_info_struct.Young_modulus_wall)),
              -1.0);
      double kn =
          1.2024 *
          pow((pow(std::get<0>(map_iterator->second)->get_properties()[19],
                   0.5) *
               pow(yEff, 2.0) *
               (std::get<0>(map_iterator->second)->get_properties()[2] / 2.0) *
               abs(std::get<4>(map_iterator->second))),
              0.4);
      double kt =
          1.2024 *
          pow((pow(std::get<0>(map_iterator->second)->get_properties()[19],
                   0.5) *
               pow(yEff, 2.0) *
               (std::get<0>(map_iterator->second)->get_properties()[2] / 2.0) *
               abs(std::get<7>(map_iterator->second))),
              0.4);
      double ethan =
          (-2.0 * log(physical_info_struct.restitution_coefficient_wall) *
           sqrt(std::get<0>(map_iterator->second)->get_properties()[19] * kn)) /
          (sqrt((pow(log(physical_info_struct.restitution_coefficient_wall),
                     2.0)) +
                pow(3.1415, 2.0)));
      double ethat = 0;
      if (physical_info_struct.restitution_coefficient_wall == 0) {
        ethat =
            2.0 *
            sqrt(2.0 / 7.0 *
                 std::get<0>(map_iterator->second)->get_properties()[19] * kt);
      } else {
        ethat =
            (-2.0 * log(physical_info_struct.restitution_coefficient_wall) *
             sqrt(2.0 / 7.0 *
                  std::get<0>(map_iterator->second)->get_properties()[19] *
                  kt)) /
            (sqrt(pow(3.1415, 2.0) +
                  pow(log(physical_info_struct.restitution_coefficient_wall),
                      2.0)));
      }

      Point<dim> springNormForce;
      springNormForce = (kn * std::get<3>(map_iterator->second)) *
                        std::get<1>(map_iterator->second);
      Point<dim> dashpotNormForce;
      dashpotNormForce = (ethan * std::get<4>(map_iterator->second)) *
                         std::get<1>(map_iterator->second);

      Point<dim> normalForce;
      normalForce = springNormForce - dashpotNormForce;

      Point<dim> temp_point;
      temp_point = std::get<6>(map_iterator->second);
      Point<dim> springTangForce =
          (kt * std::get<5>(map_iterator->second)) * temp_point;
      Point<dim> dashpotTangForce =
          (ethat * std::get<7>(map_iterator->second)) * temp_point;

      Point<dim> tangForce;
      tangForce = springTangForce - dashpotTangForce;

      if (tangForce.norm() < (physical_info_struct.friction_coefficient_wall *
                              normalForce.norm())) {
        totalForce = normalForce + tangForce;
      } else {
        Point<dim> coulumbTangForce =
            (-1.0 * physical_info_struct.friction_coefficient_wall *
             normalForce.norm() * sgn(std::get<5>(map_iterator->second))) *
            temp_point;

        totalForce = normalForce + coulumbTangForce;
      }

      std::get<0>(map_iterator->second)->get_properties()[13] =
          std::get<0>(map_iterator->second)->get_properties()[13] +
          totalForce[0];
      std::get<0>(map_iterator->second)->get_properties()[14] =
          std::get<0>(map_iterator->second)->get_properties()[14] +
          totalForce[1];
      std::get<0>(map_iterator->second)->get_properties()[15] =
          std::get<0>(map_iterator->second)->get_properties()[15] +
          totalForce[2];

      // calculation of torque
      /*
       Point<dim> torqueTi;
       torqueTi =
  ((std::get<0>(map_iterator->second)->get_properties()[2])/2.0) *
  cross_product_3d( std::get<1>(map_iterator->second) , totalForce); Point<dim>
  omegai = {std::get<0>(map_iterator->second)->get_properties()[16] ,
  std::get<0>(map_iterator->second)->get_properties()[17] ,
  std::get<0>(map_iterator->second)->get_properties()[18]};

      Point<dim> omegaiw = {0.0, 0.0, 0.0};
      double omegaNorm = omegai.norm();
      if(omegaNorm != 0)
      {omegaiw = omegai / omegaNorm ;}
      Point<dim> torquer;
     torquer = -1.0 * physical_info_struct.rolling_friction_coefficient_wall *
  ((std::get<0>(map_iterator->second)->get_properties()[2])/2.0) *
  normalForce.norm() * omegaiw;

     std::get<0>(map_iterator->second)->get_properties()[21] =
  std::get<0>(map_iterator->second)->get_properties()[21] + torqueTi[0] +
  torquer[0]; std::get<0>(map_iterator->second)->get_properties()[22] =
  std::get<0>(map_iterator->second)->get_properties()[22] + torqueTi[1] +
  torquer[1]; std::get<0>(map_iterator->second)->get_properties()[23] =
  std::get<0>(map_iterator->second)->get_properties()[23] + torqueTi[2] +
  torquer[2];
  */
    }
  }
}

template <int dim, int spacedim>
void ParticleWallContactForce<dim, spacedim>::pwNonLinearCF(
    std::vector<
        std::map<int, std::tuple<Particles::ParticleIterator<dim, spacedim>,
                                 Tensor<1, dim>, Point<dim>, double, double,
                                 double, Tensor<1, dim>, double>>>
        &pwContactInfo,
    physical_info_struct<dim> &physical_info_struct) {
  for (auto vector_iterator = pwContactInfo.begin();
       vector_iterator != pwContactInfo.end(); ++vector_iterator) {
    // for (unsigned int i = 0; i < pwContactInfo.size(); i++) {
    for (auto map_iterator = vector_iterator->begin();
         map_iterator != vector_iterator->end(); ++map_iterator) {
      Point<dim> totalForce;
      double yEff =
          pow((((1.0 - pow(physical_info_struct.Poisson_ratio_particle, 2.0)) /
                physical_info_struct.Young_modulus_particle) +
               ((1.0 - pow(physical_info_struct.Poisson_ratio_wall, 2.0)) /
                physical_info_struct.Young_modulus_wall)),
              -1.0);
      double gEff =
          pow(((2.0 * (2.0 - physical_info_struct.Poisson_ratio_particle) *
                (1.0 + physical_info_struct.Poisson_ratio_particle)) /
               (physical_info_struct.Young_modulus_particle)) +
                  ((2.0 * (2.0 - physical_info_struct.Poisson_ratio_wall) *
                    (1.0 + physical_info_struct.Poisson_ratio_wall)) /
                   (physical_info_struct.Young_modulus_wall)),
              -1.0);

      double betha =
          log(physical_info_struct.Poisson_ratio_wall) /
          sqrt(pow(log(physical_info_struct.Poisson_ratio_wall), 2.0) + 9.8696);
      double sn =
          2.0 * yEff *
          sqrt((std::get<0>(map_iterator->second)->get_properties()[2] / 2.0) *
               std::get<3>(map_iterator->second));
      double st =
          8.0 * gEff *
          sqrt((std::get<0>(map_iterator->second)->get_properties()[2] / 2.0) *
               std::get<3>(map_iterator->second));
      double kn =
          1.3333 * yEff *
          sqrt((std::get<0>(map_iterator->second)->get_properties()[2] / 2.0) *
               std::get<3>(map_iterator->second));
      double ethan =
          -1.8257 * betha *
          sqrt(sn * std::get<0>(map_iterator->second)->get_properties()[19]);
      double kt =
          8.0 * gEff *
          sqrt((std::get<0>(map_iterator->second)->get_properties()[2] / 2.0) *
               std::get<3>(map_iterator->second));
      double ethat =
          -1.8257 * betha *
          sqrt(st * std::get<0>(map_iterator->second)->get_properties()[19]);
      Point<dim> springNormForce;
      springNormForce = (kn * std::get<3>(map_iterator->second)) *
                        std::get<1>(map_iterator->second);
      Point<dim> dashpotNormForce;
      dashpotNormForce = (ethan * std::get<4>(map_iterator->second)) *
                         std::get<1>(map_iterator->second);

      Point<dim> normalForce;
      normalForce = springNormForce - dashpotNormForce;

      Point<dim> temp_point;
      temp_point = std::get<6>(map_iterator->second);
      Point<dim> springTangForce =
          (kt * std::get<5>(map_iterator->second)) * temp_point;
      Point<dim> dashpotTangForce =
          (ethat * std::get<7>(map_iterator->second)) * temp_point;

      Point<dim> tangForce;
      tangForce = springTangForce - dashpotTangForce;

      double coulumbLimit =
          physical_info_struct.friction_coefficient_wall * normalForce.norm();
      if (tangForce.norm() < coulumbLimit) {
        totalForce = normalForce + tangForce;
      } else {
        Point<dim> coulumbTangForce =
            (-1.0 * coulumbLimit * sgn(std::get<5>(map_iterator->second))) *
            temp_point;

        totalForce = normalForce + coulumbTangForce;
      }

      std::get<0>(map_iterator->second)->get_properties()[13] =
          std::get<0>(map_iterator->second)->get_properties()[13] +
          totalForce[0];
      std::get<0>(map_iterator->second)->get_properties()[14] =
          std::get<0>(map_iterator->second)->get_properties()[14] +
          totalForce[1];
      std::get<0>(map_iterator->second)->get_properties()[15] =
          std::get<0>(map_iterator->second)->get_properties()[15] +
          totalForce[2];

      // calculation of torque
      /*
       Point<dim> torqueTi;
       torqueTi =
  ((std::get<0>(map_iterator->second)->get_properties()[2])/2.0) *
  cross_product_3d( std::get<1>(map_iterator->second) , totalForce); Point<dim>
  omegai = {std::get<0>(map_iterator->second)->get_properties()[16] ,
  std::get<0>(map_iterator->second)->get_properties()[17] ,
  std::get<0>(map_iterator->second)->get_properties()[18]};

      Point<dim> omegaiw = {0.0, 0.0, 0.0};
      double omegaNorm = omegai.norm();
      if(omegaNorm != 0)
      {omegaiw = omegai / omegaNorm ;}
      Point<dim> torquer;
     torquer = -1.0 * physical_info_struct.rolling_friction_coefficient_wall *
  ((std::get<0>(map_iterator->second)->get_properties()[2])/2.0) *
  normalForce.norm() * omegaiw;

     std::get<0>(map_iterator->second)->get_properties()[21] =
  std::get<0>(map_iterator->second)->get_properties()[21] + torqueTi[0] +
  torquer[0]; std::get<0>(map_iterator->second)->get_properties()[22] =
  std::get<0>(map_iterator->second)->get_properties()[22] + torqueTi[1] +
  torquer[1]; std::get<0>(map_iterator->second)->get_properties()[23] =
  std::get<0>(map_iterator->second)->get_properties()[23] + torqueTi[2] +
  torquer[2];
  */
    }
  }
}

template <int dim, int spacedim>
int ParticleWallContactForce<dim, spacedim>::sgn(float a) {
  int b = 0;
  if (a > 0) {
    b = 1;
  } else if (a < 0) {
    b = -1;
  } else if (a == 0) {
    b = 0;
  }
  return b;
}

template class ParticleWallContactForce<3, 3>;
