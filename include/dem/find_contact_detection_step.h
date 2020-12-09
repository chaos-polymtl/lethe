/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#ifndef find_contact_detection_step_h
#define find_contact_detection_step_h

#include <deal.II/base/config.h>

#include <deal.II/particles/particle_handler.h>


using namespace dealii;


/**
 * Carries out finding steps for dynamic contact search
 *
 * @param particle_handler
 * @param dt DEM time step
 * @param smallest_contact_search_criterion A criterion for finding
 * dynamic contact search steps. This value is defined as the minimum of
 * particle-particle and particle-wall displacement threshold values
 * @param mpi_communicator
 * @param contact_detection_step Returns 1 if the maximum cumulative
 * displacement of particles exceeds the threshold and 0 otherwise
 *
 */

template <int dim>
void
find_contact_detection_step(Particles::ParticleHandler<dim> &particle_handler,
                            const double &                   dt,
                            const double &smallest_contact_search_criterion,
                            MPI_Comm &    mpi_communicator,
                            unsigned int &contact_detection_step);

#endif
