// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
   * @tparam dim An integer that denotes the number of spatial dimensions.
   * @tparam solve_type Type of solver used for the DEM.
   * @tparam dem_statistics_variable Enum variable used to identify which
   * granular statistics is being calculated
   *
   * @param particle_handler A reference to the particle handler being used for
   * DEM
   * @param mpi_communicator The MPI communicator
   */
  template <int dim, DEM::SolverType solver_type, dem_statistic_variable var>
  statistics
  calculate_granular_statistics(
    const Particles::ParticleHandler<dim> &particle_handler,
    const MPI_Comm                        &mpi_communicator);
} // namespace DEM

#endif
