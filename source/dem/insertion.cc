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

// Carries out assigning the properties of inserted particles. The output vector
// is used in insert_global_particles as input argument
template <int dim>
std::vector<std::vector<double>> Insertion<dim>::assign_particle_properties(
    const DEMSolverParameters<dim> &dem_parameters,
    const unsigned int &inserted_this_step,
    const unsigned int &inserted_so_far) {
  // Defining output vector
  std::vector<std::vector<double>> properties;

  // Getting properties as local parameters
  // TODO: MAYBE CHANGE THE INPUT TO PHYSICAL PROPERTIES DIRECTLY
  const auto physical_properties = dem_parameters.physical_properties;

  // A loop is defined over the number of particles which are going to be
  // inserted at this step
  for (unsigned int particle_counter = 0; particle_counter < inserted_this_step;
       ++particle_counter) {
    double id = (inserted_so_far + particle_counter);
    double type = 1.0;
    double diameter = physical_properties.diameter;
    double density = physical_properties.density;
    double vel_x = 0.0;
    double vel_y = 0.0;
    double vel_z = 0.0;
    double acc_x = 0.0;
    double acc_y = 0.0;
    double acc_z = 0.0;
    double f_x = 0.0;
    double f_y = 0.0;
    double f_z = 0.0;
    double w_x = 0.0;
    double w_y = 0.0;
    double w_z = 0.0;
    double mass = physical_properties.density *
                  ((4.0 / 3.0) * 3.1415 * pow((diameter / 2.0), 3.0));
    double MOI = (2.0 / 5.0) * (mass)*pow((diameter / 2.0), 2.0);
    double T_x = 0;
    double T_y = 0;
    double T_z = 0;

    std::vector<double> properties_of_one_particle{
        id,    type,  diameter, density, vel_x, vel_y, vel_z,
        acc_x, acc_y, acc_z,    f_x,     f_y,   f_z,   w_x,
        w_y,   w_z,   mass,     MOI,     T_x,   T_y,   T_z};

    properties.push_back(properties_of_one_particle);
    properties_of_one_particle.clear();
  }
  return properties;
}

template class Insertion<2>;
template class Insertion<3>;
