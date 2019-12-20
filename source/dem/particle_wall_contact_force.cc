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

ParticleWallContactForce::ParticleWallContactForce()
{}

void ParticleWallContactForce::pwLinearCF(
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                         Point<3>,
                         Point<3>,
                         double,
                         double,
                         double,
                         Point<3>,
                         double>>   pwContactInfo,
  Particles::ParticleHandler<3, 3> &particle_handler,
  ParametersDEM<3>                  DEMparam)
{
  for (unsigned int i = 0; i < pwContactInfo.size(); i++)
    {
      {
        Point<3> springNormForce =
          (DEMparam.physicalProperties.kn * std::get<3>(pwContactInfo[i])) *
          std::get<1>(pwContactInfo[i]);
        Point<3> dashpotNormForce =
          (DEMparam.physicalProperties.ethan * std::get<4>(pwContactInfo[i])) *
          std::get<1>(pwContactInfo[i]);
        Point<3> normalForce = springNormForce + dashpotNormForce;

        std::get<0>(pwContactInfo[i]).first->get_properties()[13] =
          normalForce[0];
        std::get<0>(pwContactInfo[i]).first->get_properties()[14] =
          normalForce[1];
        std::get<0>(pwContactInfo[i]).first->get_properties()[15] =
          normalForce[2];

        Point<3> springTangForce =
          (DEMparam.physicalProperties.kt * std::get<5>(pwContactInfo[i])) *
          std::get<6>(pwContactInfo[i]);
        Point<3> dashpotTangForce =
          (DEMparam.physicalProperties.ethat * std::get<7>(pwContactInfo[i])) *
          std::get<6>(pwContactInfo[i]);
        Point<3> tangForce = springTangForce + dashpotTangForce;

        if (vecValue(tangForce) <
            (DEMparam.physicalProperties.mu * vecValue(normalForce)))
          {
            std::get<0>(pwContactInfo[i]).first->get_properties()[13] =
              std::get<0>(pwContactInfo[i]).first->get_properties()[13] +
              tangForce[0];
            std::get<0>(pwContactInfo[i]).first->get_properties()[14] =
              std::get<0>(pwContactInfo[i]).first->get_properties()[14] +
              tangForce[1];
            std::get<0>(pwContactInfo[i]).first->get_properties()[15] =
              std::get<0>(pwContactInfo[i]).first->get_properties()[15] +
              tangForce[2];
          }
        else
          {
            Point<3> coulumbTangForce =
              (-1.0 * DEMparam.physicalProperties.mu * vecValue(normalForce) *
               sgn(std::get<5>(pwContactInfo[i]))) *
              std::get<6>(pwContactInfo[i]);
            std::get<0>(pwContactInfo[i]).first->get_properties()[13] =
              std::get<0>(pwContactInfo[i]).first->get_properties()[13] +
              coulumbTangForce[0];
            std::get<0>(pwContactInfo[i]).first->get_properties()[14] =
              std::get<0>(pwContactInfo[i]).first->get_properties()[14] +
              coulumbTangForce[1];
            std::get<0>(pwContactInfo[i]).first->get_properties()[15] =
              std::get<0>(pwContactInfo[i]).first->get_properties()[15] +
              coulumbTangForce[2];
          }
      }
    }
}


double ParticleWallContactForce::vecValue(Point<3> A)
{
  return (sqrt(pow(A[0], 2) + pow(A[1], 2) + pow(A[2], 2)));
}

int
ParticleWallContactForce::sgn(float a)
{
  int b;
  if (a > 0)
    {
      b = 1;
    }
  else if (a < 0)
    {
      b = -1;
    }
  else if (a == 0)
    {
      b = 0;
    }
  return b;
}
