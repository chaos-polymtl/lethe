/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <dem/insertion.h>

// Prints the insertion information
template <int dim>
void
Insertion<dim>::print_insertion_info(const unsigned int &inserted_this_step,
                                     const unsigned int &remained_particles,
                                     const unsigned int &particle_type,
                                     const ConditionalOStream &pcout)
{
  pcout << "***************************************************************** "
           "\n";
  pcout << inserted_this_step << " particles of type " << particle_type
        << " were inserted, " << remained_particles << " particles of type "
        << particle_type << " remaining" << std::endl;
  pcout << "***************************************************************** "
           "\n";
}

// Carries out assigning the properties of inserted particles. The output vector
// is used in insert_global_particles as input argument
template <int dim>
void
Insertion<dim>::assign_particle_properties(
  const DEMSolverParameters<dim> &  dem_parameters,
  const unsigned int &              inserted_this_step_this_proc,
  const unsigned int &              current_inserting_particle_type,
  std::vector<std::vector<double>> &particle_properties)
{
  // Clearing and resizing particle_properties
  particle_properties.clear();
  particle_properties.reserve(inserted_this_step_this_proc);

  // Getting properties as local parameters
  // TODO: MAYBE CHANGE THE INPUT TO PHYSICAL PROPERTIES DIRECTLY
  auto physical_properties = dem_parameters.lagrangian_physical_properties;

  particle_size_sampling(
    particle_sizes,
    physical_properties
      .particle_average_diameter[current_inserting_particle_type],
    physical_properties.particle_size_std[current_inserting_particle_type],
    inserted_this_step_this_proc);

  // A loop is defined over the number of particles which are going to be
  // inserted at this step
  for (unsigned int particle_counter = 0;
       particle_counter < inserted_this_step_this_proc;
       ++particle_counter)
    {
      double type     = current_inserting_particle_type;
      double diameter = 0.;
      (particle_sizes[particle_counter] >= 0) ?
        diameter = particle_sizes[particle_counter] :
        -particle_sizes[particle_counter];
      double density =
        physical_properties.density_particle[current_inserting_particle_type];
      double vel_x                   = 0.;
      double vel_y                   = 0.;
      double vel_z                   = 0.;
      double omega_x                 = 0.;
      double omega_y                 = 0.;
      double omega_z                 = 0.;
      double fem_force_x             = 0.;
      double fem_force_y             = 0.;
      double fem_force_z             = 0.;
      double volumetric_contribution = 0.;
      double mass                    = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);


      std::vector<double> properties_of_one_particle{type,
                                                     diameter,
                                                     vel_x,
                                                     vel_y,
                                                     vel_z,
                                                     omega_x,
                                                     omega_y,
                                                     omega_z,
                                                     fem_force_x,
                                                     fem_force_y,
                                                     fem_force_z,
                                                     volumetric_contribution,
                                                     mass};

      particle_properties.push_back(properties_of_one_particle);
      properties_of_one_particle.clear();
    }
}

template <int dim>
void
Insertion<dim>::particle_size_sampling(std::vector<double> &particle_sizes,
                                       const double         average,
                                       const double         standard_deviation,
                                       const double         particle_number)
{
  particle_sizes.clear();
  particle_sizes.reserve(particle_number);

  std::random_device         rd{};
  std::mt19937               gen{rd()};
  std::normal_distribution<> distribution{average, standard_deviation};

  for (unsigned int n = 0; n < particle_number; ++n)
    particle_sizes.push_back(distribution(gen));
}

template <int dim>
void
Insertion<dim>::calculate_insertion_domain_maximum_particle_number(
  const DEMSolverParameters<dim> &dem_parameters,
  const ConditionalOStream &      pcout)
{
  // Getting properties as local parameters
  const auto insertion_information = dem_parameters.insertion_info;

  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number;

  // distance_threshold shows the ratio of the distance between the centers of
  // two adjacent particles to the diameter of particles
  maximum_particle_number =
    int((insertion_information.x_max - insertion_information.x_min) /
        (insertion_information.distance_threshold * this->maximum_diameter)) *
    int((insertion_information.y_max - insertion_information.y_min) /
        (insertion_information.distance_threshold * this->maximum_diameter));
  if (dim == 3)
    {
      maximum_particle_number =
        maximum_particle_number *
        int(
          (insertion_information.z_max - insertion_information.z_min) /
          (insertion_information.distance_threshold * this->maximum_diameter));
    }

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step > maximum_particle_number)
    {
      pcout << "Warning: the requested number of particles for insertion ("
            << insertion_information.inserted_this_step
            << ") is higher than maximum expected number of particles ("
            << maximum_particle_number << ")" << std::endl;

      // Updating the number of inserted particles at each step
      inserted_this_step = maximum_particle_number;
    }
  else
    {
      inserted_this_step = insertion_information.inserted_this_step;
    }

  // number_of_particles_x_direction, number_of_particles_y_direction and
  // number_of_particles_z_direction are the results of discretization of the
  // insertion domain in x, y and z directions
  number_of_particles_x_direction =
    int((insertion_information.x_max - insertion_information.x_min) /
        (insertion_information.distance_threshold * this->maximum_diameter));
  number_of_particles_y_direction =
    int((insertion_information.y_max - insertion_information.y_min) /
        (insertion_information.distance_threshold * this->maximum_diameter));
  if (dim == 3)
    {
      number_of_particles_z_direction = int(
        (insertion_information.z_max - insertion_information.z_min) /
        (insertion_information.distance_threshold * this->maximum_diameter));
    }
}

template class Insertion<2>;
template class Insertion<3>;
