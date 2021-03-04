/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Lucka Barbeau, Bruno Blais, Polytechnique Montreal, 2019 -
 */

#include <deal.II/base/point.h>


#ifndef lethe_ib_particle_h
#  define lethe_ib_particle_h

using namespace dealii;

template <int dim>
class IBParticle
{
public:
    Point<dim> position;
    Point<dim> last_position;
    Tensor<1, dim> forces;
    Tensor<1, dim> last_forces;

    Tensor<1, 3> torques;
    double masses;
    Tensor<2, 3>  inertia;
    // Translational velocity
    Tensor<1, dim> velocity;
    Tensor<1, dim> last_velocity;
    Tensor<1, dim> velocity_iter;
    // Angular velocity
    Tensor<1, 3> angular_position;
    // by default the angular position is always 0 use to integrate motion
    Tensor<1, 3> last_angular_position;

    Tensor<1, 3> omega;
    Tensor<1, 3> last_omega;
    Tensor<1, 3> omega_iter;
    Tensor<1, 3> last_torques;
    double local_alpha_torque ;
    double local_alpha_force;


    double radius;

    // Pressure imposition location
    Point<dim> pressure_location;
};

#endif
