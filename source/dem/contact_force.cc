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
ContactForce<dim, spacedim>::ContactForce()
{}

template <int dim, int spacedim>
void
ContactForce<dim, spacedim>::linearCF(
  std::vector<
    std::tuple<std::pair<typename Particles::ParticleIterator<dim, spacedim>,
                         typename Particles::ParticleIterator<dim, spacedim>>,
               double,
               Point<dim>,
               double,
               Point<dim>,
               double,
               double>> contactInfo,
  ParametersDEM<dim>    DEMparam)
{
  for (unsigned int i = 0; i < contactInfo.size(); i++)
    {
      {
        Point<dim> totalForce;
        double     mEff =
          (std::get<0>(contactInfo[i]).first->get_properties()[19] *
           std::get<0>(contactInfo[i]).second->get_properties()[19]) /
          (std::get<0>(contactInfo[i]).first->get_properties()[19] +
           std::get<0>(contactInfo[i]).second->get_properties()[19]);
        double rEff =
          (std::get<0>(contactInfo[i]).first->get_properties()[2] *
           std::get<0>(contactInfo[i]).second->get_properties()[2]) /
          (2 * (std::get<0>(contactInfo[i]).first->get_properties()[2] +
                std::get<0>(contactInfo[i]).second->get_properties()[2]));
        double yEff = DEMparam.physicalProperties.Yp /
                      (2 * (1 - pow(DEMparam.physicalProperties.vp, 2)));
        // double gEff = (DEMparam.physicalProperties.Yp) / (4 *
        // (2-DEMparam.physicalProperties.vp) *
        // (1+DEMparam.physicalProperties.vp));
        double kn = 1.2024 * pow((pow(mEff, 0.5) * pow(yEff, 2) * rEff *
                                  abs(std::get<3>(contactInfo[i]))),
                                 0.4);
        double kt = 1.2024 * pow((pow(mEff, 0.5) * pow(yEff, 2) * rEff *
                                  abs(std::get<5>(contactInfo[i]))),
                                 0.4);
        double ethan =
          (-2 * log(DEMparam.physicalProperties.ep) * sqrt(mEff * kn)) /
          (sqrt((pow(log(DEMparam.physicalProperties.ep), 2)) +
                pow(3.1415, 2)));
        double ethat = 0;
        if (DEMparam.physicalProperties.ep == 0)
          {
            ethat = 2 * sqrt(2 / 7 * mEff * kt);
          }
        else
          {
            ethat = (-2 * log(DEMparam.physicalProperties.ep) *
                     sqrt(2 / 7 * mEff * kt)) /
                    (sqrt(pow(3.1415, 2) +
                          pow(log(DEMparam.physicalProperties.ep), 2)));
          }


        Point<dim> normalForce = ((-1.0 * kn * std::get<1>(contactInfo[i])) -
                                  (ethan * std::get<3>(contactInfo[i]))) *
                                 std::get<2>(contactInfo[i]);


        Point<dim> tangForce = (-1.0 * kt * std::get<6>(contactInfo[i]) -
                                ethat * std::get<5>(contactInfo[i])) *
                               std::get<4>(contactInfo[i]);


        if (tangForce.norm() <
            (DEMparam.physicalProperties.mup * normalForce.norm()))
          {
            totalForce = normalForce + tangForce;
          }
        else
          {
            Point<dim> coulumbTangForce =
              DEMparam.physicalProperties.mup * normalForce.norm() *
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
ContactForce<dim, spacedim>::sgn(float a)
{
  int b = 0;
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

template class ContactForce<3, 3>;
