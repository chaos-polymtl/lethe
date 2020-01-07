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

template <int dim, int spacedim>
ParticleInsertion<dim, spacedim>::ParticleInsertion(float x_min, float y_min, float z_min, float x_max, float y_max, float z_max, double dp, int nInsert)
{
  int n_exp;

  n_exp = int((x_max - x_min) /
              (2 * dp)) *
          int((y_max - y_min) /
              (2 * dp)) *
          int((z_max - z_min) /
              (2 * dp));
  if (nInsert > n_exp)
    std::cout << "The inserted number of particles ("
              << nInsert
              << ") is higher than maximum expected number of particles ("
              << n_exp << ")" << std::endl;
} // add error here

template <int dim, int spacedim>
void
ParticleInsertion<dim, spacedim>::uniformInsertion(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  const Triangulation<dim, spacedim> &       tr,
  int &                                      nPart,
  Particles::PropertyPool &                  pool, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max, double dp, int nInsert, int rhop, Point<dim> g)
{
  int nx = int((x_max - x_min) /
               (2 * dp));
  int ny = int((y_max - y_min) /
               (2 * dp));
  int nz = int((z_max - z_min) /
               (2 * dp));
  int nP = 0;


  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k)
        if (nP < nInsert)
          {
            Point<dim>   position;
            Point<dim>   reference_position;
            unsigned int id;

            position[0] = x_min +
                          (dp / 2) +
                          (i * 1.5 * dp);
            position[1] = y_min +
                          (dp / 2) +
                          (j * 1.5 * dp);
            position[2] = z_min +
                          (dp / 2) +
                          (k * 1.5 * dp);
            id = i * ny * nz + j * nz + k + nPart + 1;
            Particles::Particle<dim> particle(position, reference_position, id);
            typename Triangulation<dim, spacedim>::active_cell_iterator cell =
              GridTools::find_active_cell_around_point(tr,
                                                       particle.get_location());


            Particles::ParticleIterator<dim, spacedim> pit =
              particle_handler.insert_particle(particle, cell);

            particle.set_property_pool(pool);

            pit->get_properties()[0] = id;
            pit->get_properties()[1] = 1;
            pit->get_properties()[2] = dp;
            pit->get_properties()[3] = rhop;
            // Position
            pit->get_properties()[4] = position[0];
            pit->get_properties()[5] = position[1];
            pit->get_properties()[6] = position[2];
            // Velocity
            pit->get_properties()[7] = 0;
            pit->get_properties()[8] = 0;
            pit->get_properties()[9] = 0;
            // Acceleration
            pit->get_properties()[10] = 0 + g[0];
            pit->get_properties()[11] = 0 + g[1];
            pit->get_properties()[12] = 0 + g[2];
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
              rhop *
              ((4.0 / 3.0) * 3.1415 *
               pow((pit->get_properties()[2] / 2.0), 3.0));
            pit->get_properties()[20] =
              (2.0 / 5.0) * (pit->get_properties()[19]) *
              pow((pit->get_properties()[2] / 2.0), 2.0);
            // Torque
            pit->get_properties()[21] = 0;
            pit->get_properties()[22] = 0;
            pit->get_properties()[23] = 0;

              ++nP;
          }

  nPart = nPart + nInsert;
}


template <int dim, int spacedim>
void
ParticleInsertion<dim, spacedim>::nonUniformInsertion(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  const Triangulation<dim, spacedim> &       tr,
  int &                                      nPart,
  Particles::PropertyPool &                  pool, float x_min, float y_min, float z_min, float x_max, float y_max, float z_max, double dp, int nInsert, int rhop, Point<dim> g)
{
  int nx = int((x_max - x_min) /
               (2 * dp));
  int ny = int((y_max - y_min) /
               (2 * dp));
  int nz = int((z_max - z_min) /
               (2 * dp));
  int nP = 0;


  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k)
        if (nP < nInsert)
          {
            Point<dim>   position;
            Point<dim>   reference_position;
            unsigned int id;

            position[0] = x_min +
                          (dp / 2) +
                          (i * 1.5 * dp);
            position[1] = y_min +
                          (dp / 2) +
                          (j * 1.5 * dp);
            position[2] = z_min +
                          (dp / 2) +
                          (k * 1.5 * dp);
            id = i * ny * nz + j * nz + k + nPart + 1;
            Particles::Particle<dim> particle(position, reference_position, id);
            typename Triangulation<dim, spacedim>::active_cell_iterator cell =
              GridTools::find_active_cell_around_point(tr,
                                                       particle.get_location());


            Particles::ParticleIterator<dim, spacedim> pit =
              particle_handler.insert_particle(particle, cell);

            particle.set_property_pool(pool);

            pit->get_properties()[0] = id;
            pit->get_properties()[1] = 1;
            pit->get_properties()[2] = dp;
            pit->get_properties()[3] = rhop;
            // Position
            int randNum1 = rand() % 101;
            int randNum2 = rand() % 101;
            pit->get_properties()[4] = position[0] + randNum1 * (dp/400.0);
            pit->get_properties()[5] = position[1] + randNum2 * (dp/400.0);
            pit->get_properties()[6] = position[2];
            // Velocity
            pit->get_properties()[7] = 0;
            pit->get_properties()[8] = 0;
            pit->get_properties()[9] = 0;
            // Acceleration
            pit->get_properties()[10] = 0 + g[0];
            pit->get_properties()[11] = 0 + g[1];
            pit->get_properties()[12] = 0 + g[2];
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
              rhop *
              ((4.0 / 3.0) * 3.1415 *
               pow((pit->get_properties()[2] / 2.0), 3.0));
            pit->get_properties()[20] =
              (2.0 / 5.0) * (pit->get_properties()[19]) *
              pow((pit->get_properties()[2] / 2.0), 2.0);
            // Torque
            pit->get_properties()[21] = 0;
            pit->get_properties()[22] = 0;
            pit->get_properties()[23] = 0;

              ++nP;

          }

  nPart = nPart + nInsert;
}




template class ParticleInsertion<3, 3>;
