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

#include <sstream>

// Prints the insertion information
template <int dim>
void
Insertion<dim>::print_insertion_info(const unsigned int &inserted_this_step,
                                     const unsigned int &remained_particles,
                                     const unsigned int &particle_type,
                                     const ConditionalOStream &pcout)
{
  std::stringstream ss;

  ss << inserted_this_step << " particles of type " << particle_type
     << " were inserted, " << remained_particles << " particles of type "
     << particle_type << " remaining";

  announce_string(pcout, ss.str(), '*');
}

// Carries out assigning the properties of inserted particles. The output vector
// is used in insert_global_particles as input argument
template <int dim>
void
Insertion<dim>::assign_particle_properties(
  const DEMSolverParameters<dim>   &dem_parameters,
  const unsigned int               &inserted_this_step_this_proc,
  const unsigned int               &current_inserting_particle_type,
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
      double vel_x        = dem_parameters.insertion_info.vel_x;
      double vel_y        = dem_parameters.insertion_info.vel_y;
      double vel_z        = dem_parameters.insertion_info.vel_z;
      double omega_x      = dem_parameters.insertion_info.omega_x;
      double omega_y      = dem_parameters.insertion_info.omega_y;
      double omega_z      = dem_parameters.insertion_info.omega_z;
      double fem_force_x  = 0.;
      double fem_force_y  = 0.;
      double fem_force_z  = 0.;
      double fem_torque_x = 0.;
      double fem_torque_y = 0.;
      double fem_torque_z = 0.;
      double mass         = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);
      double volumetric_contribution = 0.;

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
                                                     fem_torque_x,
                                                     fem_torque_y,
                                                     fem_torque_z,
                                                     mass,
                                                     volumetric_contribution};

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
  const ConditionalOStream       &pcout)
{
  // Getting properties as local parameters
  const auto insertion_information = dem_parameters.insertion_info;

  // Checking if the insertion directions are valid (no repetition)
  int axis_sum = 0;
  if constexpr (dim == 2)
    {
      axis_sum = insertion_information.axis_0 + insertion_information.axis_1;

      AssertThrow(
        axis_sum == 1,
        ExcMessage("Invalid insertion directions: 2 directions are the same "));
    }
  if constexpr (dim == 3)
    {
      axis_sum = insertion_information.axis_0 + insertion_information.axis_1 +
                 insertion_information.axis_2;

      AssertThrow(
        axis_sum == 3,
        ExcMessage(
          "Invalid insertion directions: at least 2 directions are the same "));
    }


  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number = 1;

  number_of_particles_directions.resize(dim);
  axis_min.resize(dim);
  axis_max.resize(dim);

  std::vector<unsigned int> axis_list = {insertion_information.axis_0,
                                         insertion_information.axis_1};

  if constexpr (dim == 3)
    {
      axis_list.push_back(insertion_information.axis_2);
    }

  // Assigning the minimum and maximum positions of the insertion box in respect
  // to the axis order
  for (unsigned int axis : axis_list)
    {
      switch (axis)
        {
          case 0:
            axis_min[0] = insertion_information.x_min;
            axis_max[0] = insertion_information.x_max;
            break;
          case 1:
            axis_min[1] = insertion_information.y_min;
            axis_max[1] = insertion_information.y_max;
            break;
          case 2:
            axis_min[2] = insertion_information.z_min;
            axis_max[2] = insertion_information.z_max;
            break;
          default:
            AssertThrow(false,
                        ExcMessage("Insertion direction must be 0, 1 or 2"));
        }

      // Assign max number of particles according to the direction and calculate
      // the total max number (maximum_particle_number = max_x * max_y * max_z)
      int number_of_particles = static_cast<int>(
        (axis_max[axis] - axis_min[axis]) /
        (insertion_information.distance_threshold * this->maximum_diameter));
      number_of_particles_directions[axis] = number_of_particles;
      maximum_particle_number *= number_of_particles;
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
}

template class Insertion<2>;
template class Insertion<3>;
