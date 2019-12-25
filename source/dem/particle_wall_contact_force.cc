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
                         double>> pwContactInfo,
  ParametersDEM<3>                DEMparam)
{
  for (unsigned int i = 0; i < pwContactInfo.size(); i++)
    {
      Point<3> totalForce;
      Point<3> springNormForce =
        (DEMparam.physicalProperties.kn * std::get<3>(pwContactInfo[i])) *
        std::get<1>(pwContactInfo[i]);
      Point<3> dashpotNormForce =
        (DEMparam.physicalProperties.ethan * std::get<4>(pwContactInfo[i])) *
        std::get<1>(pwContactInfo[i]);
      Point<3> normalForce = springNormForce + dashpotNormForce;

      Point<3> springTangForce =
        (DEMparam.physicalProperties.kt * std::get<5>(pwContactInfo[i])) *
        std::get<6>(pwContactInfo[i]);
      Point<3> dashpotTangForce =
        (DEMparam.physicalProperties.ethat * std::get<7>(pwContactInfo[i])) *
        std::get<6>(pwContactInfo[i]);
      Point<3> tangForce = springTangForce + dashpotTangForce;

      if (tangForce.norm() <
          (DEMparam.physicalProperties.mu * normalForce.norm()))
        {
          totalForce = normalForce + tangForce;
        }
      else
        {
          Point<3> coulumbTangForce =
            (-1.0 * DEMparam.physicalProperties.mu * normalForce.norm() *
             sgn(std::get<5>(pwContactInfo[i]))) *
            std::get<6>(pwContactInfo[i]);

          totalForce = normalForce + coulumbTangForce;
        }


      std::get<0>(pwContactInfo[i]).first->get_properties()[13] = totalForce[0];
      std::get<0>(pwContactInfo[i]).first->get_properties()[14] = totalForce[1];
      std::get<0>(pwContactInfo[i]).first->get_properties()[15] = totalForce[2];
    }
}



int
ParticleWallContactForce::sgn(float a)
{
  int b=0;
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
