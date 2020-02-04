/*
 * ParticleInsertion.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: meteor
 */

#include "dem/nonuniform_insertion.h"

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

#include <dem/dem_properties.h>
#include <math.h>

#include <fstream>
#include <string>

#include "iostream"

using namespace DEM;

template <int dim, int spacedim>
NonUniformInsertion<dim, spacedim>::NonUniformInsertion(
  double                             dp,
  InsertionInfoStruct<dim, spacedim> insertion_info_struct)
{
  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number;

  // distance_threshold shows the ratio of the distance between the centers of
  // two adjacent particles to the diameter of particles
  maximum_particle_number =
    int((insertion_info_struct.x_max - insertion_info_struct.x_min) /
        (insertion_info_struct.distance_threshold * dp)) *
    int((insertion_info_struct.y_max - insertion_info_struct.y_min) /
        (insertion_info_struct.distance_threshold * dp)) *
    int((insertion_info_struct.z_max - insertion_info_struct.z_min) /
        (insertion_info_struct.distance_threshold * dp));

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_info_struct.inserted_number_at_step > maximum_particle_number)
    std::cout << "The inserted number of particles ("
              << insertion_info_struct.inserted_number_at_step
              << ") is higher than maximum expected number of particles ("
              << maximum_particle_number << ")" << std::endl;
  std::cout << "Inserting " << maximum_particle_number << " at this step"
            << std::endl;
} // add error here

template <int dim, int spacedim>
void
NonUniformInsertion<dim, spacedim>::insert(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  const Triangulation<dim, spacedim> &       tr,
  int &                                      total_particle_number_insystem,
  Particles::PropertyPool &                  pool,
  double                                     dp,
  int                                        rhop,
  Tensor<1, dim>                             g,
  InsertionInfoStruct<dim, spacedim>         insertion_info_struct)
{
  // nx, ny and nz are the results of discretization of the insertion domain in
  // x, y and z directions
  int nx = int((insertion_info_struct.x_max - insertion_info_struct.x_min) /
               (insertion_info_struct.distance_threshold * dp));
  int ny = int((insertion_info_struct.y_max - insertion_info_struct.y_min) /
               (insertion_info_struct.distance_threshold * dp));
  int nz = int((insertion_info_struct.z_max - insertion_info_struct.z_min) /
               (insertion_info_struct.distance_threshold * dp));

  // inserted_sofar_step shows the number of inserted particles, so far, at
  // this step
  int inserted_sofar_step = 0;

  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      if (dim == 3)
        for (int k = 0; k < nz; ++k)

          // We need to check if the number of inserted particles so far at this
          // step reached the total desired number of inserted particles at this
          // step
          if (inserted_sofar_step <
              insertion_info_struct.inserted_number_at_step)
            {
              Point<dim>   position;
              Point<dim>   reference_position;
              unsigned int id;

              // Obtaning position of the inserted particle
              // In non-uniform insertion, two random numbers are created and
              // added to the position of particles
              int randNum1 = rand() % 101;
              int randNum2 = rand() % 101;
              position[0] =
                insertion_info_struct.x_min + (dp / 2) +
                (i * insertion_info_struct.distance_threshold * dp) +
                randNum1 * (dp / 400.0);
              position[1] =
                insertion_info_struct.y_min + (dp / 2) +
                (j * insertion_info_struct.distance_threshold * dp) +
                randNum2 * (dp / 400.0);
              if (dim == 3)
                position[2] =
                  insertion_info_struct.z_min + (dp / 2) +
                  (k * insertion_info_struct.distance_threshold * dp);

              // Since the id of each particle should be unique, we need to use
              // the total number of particles in the system
              // (total_particle_number_insystem) to calculate the ids of new
              // particles
              if (dim == 3)
                id =
                  i * ny * nz + j * nz + k + total_particle_number_insystem + 1;
              if (dim == 2)
                id = i * ny + j + total_particle_number_insystem;

              // Inserting the new particle using its location, id and
              // containing cell
              Particles::Particle<dim> particle(position,
                                                reference_position,
                                                id);
              typename Triangulation<dim, spacedim>::active_cell_iterator cell =
                GridTools::find_active_cell_around_point(
                  tr, particle.get_location());
              Particles::ParticleIterator<dim, spacedim> pit =
                particle_handler.insert_particle(particle, cell);

              // Setting property pool of inserted particle
              particle.set_property_pool(pool);

              // Initialization of the properties of the new particle
              pit->get_properties()[PropertiesIndex::id]   = id;
              pit->get_properties()[PropertiesIndex::type] = 1;
              pit->get_properties()[PropertiesIndex::dp]   = dp;
              pit->get_properties()[PropertiesIndex::rho]  = rhop;
              // Velocity
              pit->get_properties()[PropertiesIndex::v_x] = 0;
              pit->get_properties()[PropertiesIndex::v_y] = 0;
              if (dim == 3)
                pit->get_properties()[PropertiesIndex::v_z] = 0;
              // Acceleration
              pit->get_properties()[PropertiesIndex::acc_x] = 0 + g[0];
              pit->get_properties()[PropertiesIndex::acc_y] = 0 + g[1];
              if (dim == 3)
                pit->get_properties()[PropertiesIndex::acc_z] = 0 + g[2];
              // Force
              pit->get_properties()[PropertiesIndex::force_x] = 0;
              pit->get_properties()[PropertiesIndex::force_y] = 0;
              if (dim == 3)
                pit->get_properties()[PropertiesIndex::force_z] = 0;
              // w
              pit->get_properties()[PropertiesIndex::omega_x] = 0;
              pit->get_properties()[PropertiesIndex::omega_y] = 0;
              if (dim == 3)
                pit->get_properties()[PropertiesIndex::omega_z] = 0;
              // mass and moi
              pit->get_properties()[PropertiesIndex::mass] =
                rhop *
                ((4.0 / 3.0) * 3.1415 *
                 pow((pit->get_properties()[PropertiesIndex::dp] / 2.0), 3.0));

              pit->get_properties()[PropertiesIndex::mom_inertia] =
                (2.0 / 5.0) * (pit->get_properties()[PropertiesIndex::mass]) *
                pow((pit->get_properties()[PropertiesIndex::dp] / 2.0), 2.0);
              // Torque
              pit->get_properties()[PropertiesIndex::M_x] = 0;
              pit->get_properties()[PropertiesIndex::M_y] = 0;
              if (dim == 3)
                pit->get_properties()[PropertiesIndex::M_z] = 0;

              ++inserted_sofar_step;
            }

  // Updating the total number of particles in the system
  total_particle_number_insystem =
    total_particle_number_insystem +
    insertion_info_struct.inserted_number_at_step;
}

template class NonUniformInsertion<3, 3>;
