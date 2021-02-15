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



// This class defined values related to a particle used for in the sharp interface immersed boundary.
    // Each particle defined will have these value used in the solver
    void initialise_all();
    void initialise_last();

    Point<dim> position;
    Point<dim> last_position;
    Tensor<1, dim> forces;
    Tensor<1, dim> last_forces;

    Tensor<1, 3> torques;
    double mass;
    Tensor<2, 3>  inertia;
    // Translational velocity
    Tensor<1, dim> velocity;
    // Store the last velocity of the fix point iteration.
    Tensor<1, dim> last_velocity;
    Tensor<1, dim> velocity_iter;
    // Angular velocity

    // By default the angular position is always 0 on every axis.
    Tensor<1, 3> angular_position;
    // Store the last angular position of the particle for integration.
    Tensor<1, 3> last_angular_position;

    // Angular velocity
    Tensor<1, 3> omega;
    // Store the last angular velocity of the particle for integration.
    Tensor<1, 3> last_omega;
    // Store the last  angular velocity of the fix point iteration.
    Tensor<1, 3> omega_iter;

    // Allow the definition of a local relaxation parameter for each particle in the integration process.
    double local_alpha_torque ;
    double local_alpha_force;


    double radius;


    // Pressure imposition location
    Point<dim> pressure_location;
};

#endif
