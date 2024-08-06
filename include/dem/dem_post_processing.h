/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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

#ifndef lethe_dem_post_processing_h
#define lethe_dem_post_processing_h

#include <core/utilities.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

namespace DEM
{
  /**
   * @brief enum that is used to identify which variables are calculated in the granular statistics
   *
   */
  enum class dem_statistic_variable
  {
    translational_kinetic_energy,
    rotational_kinetic_energy,
    velocity,
    omega,
  };

  /**
   * @brief Calculate statistics on a DEM ParticleHandler. At the moment, the
   * following statistics are supported:
   * - Translational kinetic energy
   * - Rotational (angular) kinetic energy
   * - Translational velocity
   * - Rotational (angular) velocity
   *
   * @tparam dim Dimensionality of the problem (2D or 3D)
   * @tparam dem_statistics_variable Enum variable used to identify which
   * granular statistics is being calculated
   *
   * @param particle_handler A reference to the particle handler being used for
   * DEM
   * @param mpi_communicator The MPI communicator
   */
  template <int dim, dem_statistic_variable var>
  statistics
  calculate_granular_statistics(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm                        &mpi_communicator);
} // namespace DEM

#endif
