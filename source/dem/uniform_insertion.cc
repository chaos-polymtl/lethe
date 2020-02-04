/*
 * ParticleInsertion.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */

#include "dem/uniform_insertion.h"
#include <dem/dem_properties.h>

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

using namespace DEM;

template <int dim, int spacedim>
UniformInsertion<dim, spacedim>::UniformInsertion(double x_min, double y_min,
                                                  double z_min, double x_max,
                                                  double y_max, double z_max,
                                                  double dp,
                                                  int inserted_number_at_step,
                                                  double distance_threshold) {
  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number;

  // distance_threshold shows the ratio of the distance between the centers of
  // two adjacent particles to the diameter of particles
  maximum_particle_number = int((x_max - x_min) / (distance_threshold * dp)) *
                            int((y_max - y_min) / (distance_threshold * dp)) *
                            int((z_max - z_min) / (distance_threshold * dp));

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (inserted_number_at_step > maximum_particle_number)
    std::cout << "The inserted number of particles (" << inserted_number_at_step
              << ") is higher than maximum expected number of particles ("
              << maximum_particle_number << ")" << std::endl;
} // add error here

template <int dim, int spacedim>
void UniformInsertion<dim, spacedim>::insert(
    Particles::ParticleHandler<dim, spacedim> &particle_handler,
    const Triangulation<dim, spacedim> &tr, int &total_particle_number_insystem,
    Particles::PropertyPool &pool, double x_min, double y_min, double z_min,
    double x_max, double y_max, double z_max, double dp,
    int inserted_number_at_step, int rhop, Tensor<1, dim> g,
    double distance_threshold) {

  // nx, ny and nz are the results of discretization of the insertion domain in
  // x, y and z directions
  int nx = int((x_max - x_min) / (distance_threshold * dp));
  int ny = int((y_max - y_min) / (distance_threshold * dp));
  int nz = int((z_max - z_min) / (distance_threshold * dp));

  // inserted_sofar_step shows the number of inserted particles, so far, at this
  // step
  int inserted_sofar_step = 0;

  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k)

        // We need to check if the number of inserted particles so far at this
        // step reached the total desired number of inserted particles at this
        // step
        if (inserted_sofar_step < inserted_number_at_step) {
          Point<dim> position;
          Point<dim> reference_position;
          unsigned int id;

          // Obtaning position of the inserted particle
          position[0] = x_min + (dp / 2) + (i * distance_threshold * dp);
          position[1] = y_min + (dp / 2) + (j * distance_threshold * dp);
          position[2] = z_min + (dp / 2) + (k * distance_threshold * dp);

          // Since the id of each particle should be unique, we need to use the
          // total number of particles in the system
          // (total_particle_number_insystem) to calculate the ids of new
          // particles
          id = i * ny * nz + j * nz + k + total_particle_number_insystem + 1;

          // Inserting the new particle using its location, id and containing
          // cell
          Particles::Particle<dim> particle(position, reference_position, id);
          typename Triangulation<dim, spacedim>::active_cell_iterator cell =
              GridTools::find_active_cell_around_point(tr,
                                                       particle.get_location());
          Particles::ParticleIterator<dim, spacedim> pit =
              particle_handler.insert_particle(particle, cell);

          // Setting property pool of inserted particle
          particle.set_property_pool(pool);

          // Initialization of the properties of the new particle
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
              rhop * ((4.0 / 3.0) * 3.1415 *
                      pow((pit->get_properties()[2] / 2.0), 3.0));
          pit->get_properties()[20] =
              (2.0 / 5.0) * (pit->get_properties()[19]) *
              pow((pit->get_properties()[2] / 2.0), 2.0);
          // Torque
          pit->get_properties()[21] = 0;
          pit->get_properties()[22] = 0;
          pit->get_properties()[23] = 0;

          ++inserted_sofar_step;
        }

  // Updating the total number of particles in the system
  total_particle_number_insystem =
      total_particle_number_insystem + inserted_number_at_step;
}

template class UniformInsertion<3, 3>;
