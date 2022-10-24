/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2021 by the Lethe authors
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
 */

#ifndef dem_post_processing_h
#define dem_post_processing_h

#include <core/utilities.h>
#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Calculate total kinetic energy of particles
 *
 * @param particle_handler A reference to the particle handler being used for DEM
 * @param mpi_communicator
 */

namespace DEM
{
  template <int dim>
  statistics
  calculate_granular_kinetic_energy(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm &                       mpi_communicator);
}
#endif
