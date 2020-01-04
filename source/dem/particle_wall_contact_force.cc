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
ParticleWallContactForce<dim, spacedim>::ParticleWallContactForce()
{}

template <int dim, int spacedim>
void
ParticleWallContactForce<dim, spacedim>::pwLinearCF(
  std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>,
    double,
    double,
    double,
    Point<dim>,
    double>>         pwContactInfo,
  ParametersDEM<dim> DEMparam)
{
  for (unsigned int i = 0; i < pwContactInfo.size(); i++)
    {
      Point<dim> totalForce;
      double     yEff = pow((((1.0 - pow(DEMparam.physicalProperties.vp, 2.0)) /
                          DEMparam.physicalProperties.Yp) +
                         ((1.0 - pow(DEMparam.physicalProperties.vw, 2.0)) /
                          DEMparam.physicalProperties.Yw)),
                        -1.0);
      double     kn =
        1.2024 *
        pow((pow(std::get<0>(pwContactInfo[i]).first->get_properties()[19],
                 0.5) *
             pow(yEff, 2.0) *
             (std::get<0>(pwContactInfo[i]).first->get_properties()[2] / 2.0) *
             abs(std::get<4>(pwContactInfo[i]))),
            0.4);
      double kt =
        1.2024 *
        pow((pow(std::get<0>(pwContactInfo[i]).first->get_properties()[19],
                 0.5) *
             pow(yEff, 2.0) *
             (std::get<0>(pwContactInfo[i]).first->get_properties()[2] / 2.0) *
             abs(std::get<7>(pwContactInfo[i]))),
            0.4);
      double ethan =
        (-2.0 * log(DEMparam.physicalProperties.ew) *
         sqrt(std::get<0>(pwContactInfo[i]).first->get_properties()[19] * kn)) /
        (sqrt((pow(log(DEMparam.physicalProperties.ew), 2.0)) + pow(3.1415, 2.0)));
      double ethat = 0;
      if (DEMparam.physicalProperties.ew == 0)
        {
          ethat =
            2.0 * sqrt(2.0 / 7.0 *
                     std::get<0>(pwContactInfo[i]).first->get_properties()[19] *
                     kt);
        }
      else
        {
          ethat =
            (-2.0 * log(DEMparam.physicalProperties.ew) *
             sqrt(2.0 / 7.0 *
                  std::get<0>(pwContactInfo[i]).first->get_properties()[19] *
                  kt)) /
            (sqrt(pow(3.1415, 2.0) +
                  pow(log(DEMparam.physicalProperties.ew), 2.0)));
        }


      Point<dim> springNormForce =
        (kn * std::get<3>(pwContactInfo[i])) * std::get<1>(pwContactInfo[i]);
      Point<dim> dashpotNormForce =
        (ethan * std::get<4>(pwContactInfo[i])) * std::get<1>(pwContactInfo[i]);

      Point<dim> normalForce;
      normalForce = springNormForce - dashpotNormForce;


      Point<dim> springTangForce =
        (kt * std::get<5>(pwContactInfo[i])) * std::get<6>(pwContactInfo[i]);
      Point<dim> dashpotTangForce =
        (ethat * std::get<7>(pwContactInfo[i])) * std::get<6>(pwContactInfo[i]);

      Point<dim> tangForce;
      tangForce= springTangForce - dashpotTangForce;

      if (tangForce.norm() <
          (DEMparam.physicalProperties.muw * normalForce.norm()))
        {
          totalForce = normalForce + tangForce;
        }
      else
        {
          Point<dim> coulumbTangForce =
            (-1.0 * DEMparam.physicalProperties.muw * normalForce.norm() *
             sgn(std::get<5>(pwContactInfo[i]))) *
            std::get<6>(pwContactInfo[i]);

          totalForce = normalForce + coulumbTangForce;
        }


      std::get<0>(pwContactInfo[i]).first->get_properties()[13] = std::get<0>(pwContactInfo[i]).first->get_properties()[13] + totalForce[0];
      std::get<0>(pwContactInfo[i]).first->get_properties()[14] = std::get<0>(pwContactInfo[i]).first->get_properties()[14] + totalForce[1];
      std::get<0>(pwContactInfo[i]).first->get_properties()[15] = std::get<0>(pwContactInfo[i]).first->get_properties()[15] + totalForce[2];

      // calculation of torque
      /*
       Point<dim> torqueTi;
       torqueTi = ((std::get<0>(pwContactInfo[i]).first->get_properties()[2])/2.0) * cross_product_3d( std::get<1>(pwContactInfo[i]) , totalForce);
 Point<dim> omegai = {std::get<0>(pwContactInfo[i]).first->get_properties()[16] ,   std::get<0>(pwContactInfo[i]).first->get_properties()[17] ,   std::get<0>(pwContactInfo[i]).first->get_properties()[18]};

      Point<dim> omegaiw = {0.0, 0.0, 0.0};
      double omegaNorm = omegai.norm();
      if(omegaNorm != 0)
      {omegaiw = omegai / omegaNorm ;}
      Point<dim> torquer;
     torquer = -1.0 * DEMparam.physicalProperties.murw * ((std::get<0>(pwContactInfo[i]).first->get_properties()[2])/2.0) * normalForce.norm() * omegaiw;

     std::get<0>(pwContactInfo[i]).first->get_properties()[21] = std::get<0>(pwContactInfo[i]).first->get_properties()[21] + torqueTi[0] + torquer[0];
     std::get<0>(pwContactInfo[i]).first->get_properties()[22] = std::get<0>(pwContactInfo[i]).first->get_properties()[22] + torqueTi[1] + torquer[1];
     std::get<0>(pwContactInfo[i]).first->get_properties()[23] = std::get<0>(pwContactInfo[i]).first->get_properties()[23] + torqueTi[2] + torquer[2];
  */
}
}


template <int dim, int spacedim>
void
ParticleWallContactForce<dim, spacedim>::pwNonLinearCF(
  std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>,
    double,
    double,
    double,
    Point<dim>,
    double>>         pwContactInfo,
  ParametersDEM<dim> DEMparam)
{
  for (unsigned int i = 0; i < pwContactInfo.size(); i++)
    {
      Point<dim> totalForce;
      double     yEff = pow((((1.0 - pow(DEMparam.physicalProperties.vp, 2.0)) /
                          DEMparam.physicalProperties.Yp) +
                         ((1.0 - pow(DEMparam.physicalProperties.vw, 2.0)) /
                          DEMparam.physicalProperties.Yw)),
                        -1.0);
      double gEff = pow(((2.0*(2.0-DEMparam.physicalProperties.vp)*(1.0+DEMparam.physicalProperties.vp)) / (DEMparam.physicalProperties.Yp)) + ((2.0*(2.0-DEMparam.physicalProperties.vw)*(1.0+DEMparam.physicalProperties.vw)) / (DEMparam.physicalProperties.Yw)),-1.0);

              double betha = log(DEMparam.physicalProperties.vw)/sqrt(pow(log(DEMparam.physicalProperties.vw),2.0) + 9.8696) ;
              double sn = 2.0 * yEff * sqrt((std::get<0>(pwContactInfo[i]).first->get_properties()[2] / 2.0) * std::get<3>(pwContactInfo[i]));
              double st = 8.0 * gEff * sqrt((std::get<0>(pwContactInfo[i]).first->get_properties()[2] / 2.0) * std::get<3>(pwContactInfo[i]));
              double kn = 1.3333 * yEff * sqrt((std::get<0>(pwContactInfo[i]).first->get_properties()[2] / 2.0) * std::get<3>(pwContactInfo[i]));
              double ethan = -1.8257 * betha * sqrt(sn * std::get<0>(pwContactInfo[i]).first->get_properties()[19]);
              double kt = 8.0 * gEff * sqrt((std::get<0>(pwContactInfo[i]).first->get_properties()[2] / 2.0) * std::get<3>(pwContactInfo[i]));
              double ethat = -1.8257 * betha * sqrt(st * std::get<0>(pwContactInfo[i]).first->get_properties()[19]);
      Point<dim> springNormForce =
        (kn * std::get<3>(pwContactInfo[i])) * std::get<1>(pwContactInfo[i]);
      Point<dim> dashpotNormForce =
        (ethan * std::get<4>(pwContactInfo[i])) * std::get<1>(pwContactInfo[i]);

      Point<dim> normalForce;
      normalForce = springNormForce - dashpotNormForce;


      Point<dim> springTangForce =
        (kt * std::get<5>(pwContactInfo[i])) * std::get<6>(pwContactInfo[i]);
      Point<dim> dashpotTangForce =
        (ethat * std::get<7>(pwContactInfo[i])) * std::get<6>(pwContactInfo[i]);

      Point<dim> tangForce;
      tangForce= springTangForce - dashpotTangForce;

      double coulumbLimit = DEMparam.physicalProperties.muw * normalForce.norm();
      if (tangForce.norm() < coulumbLimit)
        {
          totalForce = normalForce + tangForce;
        }
      else
        {
          Point<dim> coulumbTangForce =
            (-1.0 * coulumbLimit *
             sgn(std::get<5>(pwContactInfo[i]))) *
            std::get<6>(pwContactInfo[i]);

          totalForce = normalForce + coulumbTangForce;
        }


      std::get<0>(pwContactInfo[i]).first->get_properties()[13] = std::get<0>(pwContactInfo[i]).first->get_properties()[13] + totalForce[0];
      std::get<0>(pwContactInfo[i]).first->get_properties()[14] = std::get<0>(pwContactInfo[i]).first->get_properties()[14] + totalForce[1];
      std::get<0>(pwContactInfo[i]).first->get_properties()[15] = std::get<0>(pwContactInfo[i]).first->get_properties()[15] + totalForce[2];

      // calculation of torque
      /*
       Point<dim> torqueTi;
       torqueTi = ((std::get<0>(pwContactInfo[i]).first->get_properties()[2])/2.0) * cross_product_3d( std::get<1>(pwContactInfo[i]) , totalForce);
 Point<dim> omegai = {std::get<0>(pwContactInfo[i]).first->get_properties()[16] ,   std::get<0>(pwContactInfo[i]).first->get_properties()[17] ,   std::get<0>(pwContactInfo[i]).first->get_properties()[18]};

      Point<dim> omegaiw = {0.0, 0.0, 0.0};
      double omegaNorm = omegai.norm();
      if(omegaNorm != 0)
      {omegaiw = omegai / omegaNorm ;}
      Point<dim> torquer;
     torquer = -1.0 * DEMparam.physicalProperties.murw * ((std::get<0>(pwContactInfo[i]).first->get_properties()[2])/2.0) * normalForce.norm() * omegaiw;

     std::get<0>(pwContactInfo[i]).first->get_properties()[21] = std::get<0>(pwContactInfo[i]).first->get_properties()[21] + torqueTi[0] + torquer[0];
     std::get<0>(pwContactInfo[i]).first->get_properties()[22] = std::get<0>(pwContactInfo[i]).first->get_properties()[22] + torqueTi[1] + torquer[1];
     std::get<0>(pwContactInfo[i]).first->get_properties()[23] = std::get<0>(pwContactInfo[i]).first->get_properties()[23] + torqueTi[2] + torquer[2];
  */
}
}

template <int dim, int spacedim>
int
ParticleWallContactForce<dim, spacedim>::sgn(float a)
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

template class ParticleWallContactForce<3, 3>;
