/*
 * contactForce.cpp
 *
 *  Created on: Oct 31, 2019
 *      Author: shahab
 */

#include "dem/contact_force.h"

#include <deal.II/base/point.h>

#include "dem/dem_iterator.h"

using namespace dealii;

template <int dim, int spacedim>
ContactForce<dim,spacedim>::ContactForce()
{}

template <int dim, int spacedim>
void ContactForce<dim,spacedim>::linearCF(
  std::vector<std::tuple<std::pair<typename Particles::ParticleIterator<dim,spacedim>,
                                   typename Particles::ParticleIterator<dim,spacedim>>,
                         double,
                         Point<dim>,
                         double,
                         Point<dim>,
                         double,
                         double>>   contactInfo,
  ParametersDEM<dim>                  DEMparam)
{
  for (unsigned int i = 0; i < contactInfo.size(); i++)
    {
      {
        Point<dim> totalForce;
        Point<dim> normalForce =
          ((-1.0 * DEMparam.physicalProperties.kn *
            std::get<1>(contactInfo[i])) -
           (DEMparam.physicalProperties.ethan * std::get<3>(contactInfo[i]))) *
          std::get<2>(contactInfo[i]);


        Point<dim> tangForce =
          (-1.0 * DEMparam.physicalProperties.kt * std::get<6>(contactInfo[i]) -
           DEMparam.physicalProperties.ethat * std::get<5>(contactInfo[i])) *
          std::get<4>(contactInfo[i]);


        if (tangForce.norm() <
            (DEMparam.physicalProperties.mu * normalForce.norm()))
          {
            totalForce = normalForce + tangForce;
          }
        else
          {
            Point<dim> coulumbTangForce =
              DEMparam.physicalProperties.mu * normalForce.norm() *
              sgn(std::get<6>(contactInfo[i])) * std::get<4>(contactInfo[i]);

            totalForce = normalForce + coulumbTangForce;
          }

        std::get<0>(contactInfo[i]).first->get_properties()[13] = totalForce[0];
        std::get<0>(contactInfo[i]).first->get_properties()[14] = totalForce[1];
        std::get<0>(contactInfo[i]).first->get_properties()[15] = totalForce[2];

        std::get<0>(contactInfo[i]).second->get_properties()[13] =
          -1 * totalForce[0];
        std::get<0>(contactInfo[i]).second->get_properties()[14] =
          -1 * totalForce[1];
        std::get<0>(contactInfo[i]).second->get_properties()[15] =
          -1 * totalForce[2];

        // calculation of torque
        // Point<dim> Torquei =
        // ((std::get<0>(contactInfo[i]).first->get_properties()[2])/2.0) *
        // crossProduct( , )
      }
    }
}


template <int dim, int spacedim>
int
ContactForce<dim,spacedim>::sgn(float a)
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

template class ContactForce<3,3>;
