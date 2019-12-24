/*
 * ParticleInsertion.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */

#include "dem/particle_insertion.h"

#include <deal.II/base/array_view.h>
#include <deal.II/base/data_out_base.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <math.h>

#include <fstream>
#include <string>

#include "iostream"



using namespace dealii;

ParticleInsertion::ParticleInsertion(ParametersDEM<3> DEMparam)
{
  int n_exp;

  n_exp = int((DEMparam.insertionInfo.x_max - DEMparam.insertionInfo.x_min) /
              (2 * DEMparam.physicalProperties.diameter)) *
          int((DEMparam.insertionInfo.y_max - DEMparam.insertionInfo.y_min) /
              (2 * DEMparam.physicalProperties.diameter)) *
          int((DEMparam.insertionInfo.z_max - DEMparam.insertionInfo.z_min) /
              (2 * DEMparam.physicalProperties.diameter));
  if (DEMparam.insertionInfo.nInsert > n_exp)
    std::cout << "The inserted number of particles ("
              << DEMparam.insertionInfo.nInsert
              << ") is higher than maximum expected number of particles ("
              << n_exp << ")" << std::endl;
} // add error here


void ParticleInsertion::uniformInsertion(
  Particles::ParticleHandler<3, 3> &particle_handler,
  const Triangulation<3, 3> &       tr,
  ParametersDEM<3>                  DEMparam,
  int &                             nPart,
  Particles::PropertyPool &         pool)
{

  int nx = int((DEMparam.insertionInfo.x_max - DEMparam.insertionInfo.x_min) /
               (2 * DEMparam.physicalProperties.diameter));
  int ny = int((DEMparam.insertionInfo.y_max - DEMparam.insertionInfo.y_min) /
               (2 * DEMparam.physicalProperties.diameter));
  int nz = int((DEMparam.insertionInfo.z_max - DEMparam.insertionInfo.z_min) /
               (2 * DEMparam.physicalProperties.diameter));
  int nP = 0;


  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k)
        if (nP < DEMparam.insertionInfo.nInsert)
          {
            Point<3>     position;
            Point<3>     reference_position;
            unsigned int id;

            position[0] = DEMparam.insertionInfo.x_min +
                          (DEMparam.physicalProperties.diameter / 2) +
                          (i * 1.1 * DEMparam.physicalProperties.diameter);
            position[1] = DEMparam.insertionInfo.y_min +
                          (DEMparam.physicalProperties.diameter / 2) +
                          (j * 1.1 * DEMparam.physicalProperties.diameter);
            position[2] = DEMparam.insertionInfo.z_min +
                          (DEMparam.physicalProperties.diameter / 2) +
                          (k * 1.1 * DEMparam.physicalProperties.diameter);
            id = i * ny * nz + j * nz + k + nPart + 1;
            Particles::Particle<3> particle(position, reference_position, id);
            Triangulation<3, 3>::active_cell_iterator cell =
              GridTools::find_active_cell_around_point(tr,
                                                       particle.get_location());


            Particles::ParticleIterator<3, 3> pit =
              particle_handler.insert_particle(particle, cell);

            particle.set_property_pool(pool);

            pit->get_properties()[0] = id;
            pit->get_properties()[1] = 1;
            pit->get_properties()[2] = DEMparam.physicalProperties.diameter;
            pit->get_properties()[3] = DEMparam.physicalProperties.density;
            // Position
            pit->get_properties()[4] = position[0];
            pit->get_properties()[5] = position[1];
            pit->get_properties()[6] = position[2];
            // Velocity
            pit->get_properties()[7] = 0;
            pit->get_properties()[8] = 0;
            pit->get_properties()[9] = 0;
            // Acceleration
            pit->get_properties()[10] = 0 + DEMparam.physicalProperties.gx;
            pit->get_properties()[11] = 0 + DEMparam.physicalProperties.gy;
            pit->get_properties()[12] = 0 + DEMparam.physicalProperties.gz;
            // Force
            pit->get_properties()[13] = 0;
            pit->get_properties()[14] = 0;
            pit->get_properties()[15] = 0;
            // w
            pit->get_properties()[16] = 0;
            pit->get_properties()[17] = 0;
            pit->get_properties()[18] = 0;
            // mass and moi
            pit->get_properties()[19] =
              DEMparam.physicalProperties.density *
              ((4.0 / 3.0) * 3.1415 *
               pow((pit->get_properties()[2] / 2.0), 3.0));
            pit->get_properties()[20] =
              (2.0 / 5.0) * (pit->get_properties()[19]) *
              pow((pit->get_properties()[2] / 2.0), 2.0);



            ++nP;
          }

  nPart = nPart + DEMparam.insertionInfo.nInsert;
}
