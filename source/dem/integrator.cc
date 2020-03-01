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

#include <dem/integrator.h>

using namespace dealii;

// Why are all the numbers for the properties hard-coded? This should
// Come out of an enum class or something like this...

// Not sure if this implementation of the RK2 integrator is valid and there are
// currently no unit tests for it

// template <int dim>
// void
// Integrator<dim>::rk2Integration(
//  Particles::ParticleHandler<dim> &particle_handler,
//  Point<dim>                                 g,
//  float                                      dt)
//{
//  for (auto particle = particle_handler.begin();
//       particle != particle_handler.end();
//       ++particle)
//    {
//      // Acceleration calculation:
//      double axStar = particle->get_properties()[10];
//      double ayStar = particle->get_properties()[11];
//      double azStar = particle->get_properties()[12];

//      particle->get_properties()[10] =
//        g[0] +
//        (particle->get_properties()[13]) / (particle->get_properties()[19]);
//      particle->get_properties()[11] =
//        g[1] +
//        (particle->get_properties()[14]) / (particle->get_properties()[19]);
//      particle->get_properties()[12] =
//        g[2] +
//        (particle->get_properties()[15]) / (particle->get_properties()[19]);

//      // Velocity integration:
//      double vxStar = particle->get_properties()[7];
//      double vyStar = particle->get_properties()[8];
//      double vzStar = particle->get_properties()[9];

//      particle->get_properties()[7] =
//        particle->get_properties()[7] +
//        (dt / 2.0) * (axStar + particle->get_properties()[10]);

//      particle->get_properties()[8] =
//        particle->get_properties()[8] +
//        (dt / 2.0) * (ayStar + particle->get_properties()[11]);
//      particle->get_properties()[9] =
//        particle->get_properties()[9] +
//        (dt / 2.0) * (azStar + particle->get_properties()[12]);

//      // Position integration:
//      particle->get_properties()[4] =
//        particle->get_properties()[4] +
//        (dt / 2.0) * (vxStar + particle->get_properties()[7]);
//      particle->get_properties()[5] =
//        particle->get_properties()[5] +
//        (dt / 2.0) * (vyStar + particle->get_properties()[8]);
//      particle->get_properties()[6] =
//        particle->get_properties()[6] +
//        (dt / 2.0) * (vzStar + particle->get_properties()[9]);

//      particle->set_location({particle->get_properties()[4],
//                              particle->get_properties()[5],
//                              particle->get_properties()[6]});

//      // Angular velocity using Euler:
//      /*
//      particle->get_properties()[16] = particle->get_properties()[16] +
//        (particle->get_properties()[21]) / (particle->get_properties()[20]);
//      particle->get_properties()[17] =particle->get_properties()[17] +
//        (particle->get_properties()[22]) / (particle->get_properties()[20]);
//      particle->get_properties()[18] = particle->get_properties()[18] +
//        (particle->get_properties()[23]) / (particle->get_properties()[20]);
//        */
//    }
//}
