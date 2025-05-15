// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_stabilization_h
#define lethe_stabilization_h

#include <deal.II/base/utilities.h>

#include <cmath>

using namespace dealii;

/**
 * @brief Calculate the stabilization parameter for the Navier-Stokes equations
 * in steady-state
 * @return Value of the stabilization parameter - tau
 *
 * @param u_mag Magnitude of the velocity
 *
 * @param kinematic_viscosity Kinematic viscosity
 *
 * @param h Cell size; it should be calculated using the diameter of a sphere of
 * equal volume to that of the cell.
 */
inline double
calculate_navier_stokes_gls_tau_steady(const double u_mag,
                                       const double kinematic_viscosity,
                                       const double h)
{
  return 1. / std::sqrt(Utilities::fixed_power<2>(2. * u_mag / h) +
                        9 * Utilities::fixed_power<2>(4 * kinematic_viscosity /
                                                      (h * h)));
}

/**
 * @brief Calculate the stabilization parameter for the transient Navier-Stokes
 * equations
 * @return Value of the stabilization parameter - tau
 *
 * @param u_mag Magnitude of the velocity
 *
 * @param kinematic_viscosity Kinematic viscosity
 *
 * @param h Cell size; it should be calculated using the diameter of a sphere of
 * equal volume to that of the cell.
 *
 * @param sdt Inverse of the time-step (1/dt)
 */

inline double
calculate_navier_stokes_gls_tau_transient(const double u_mag,
                                          const double kinematic_viscosity,
                                          const double h,
                                          const double sdt)
{
  return 1. / std::sqrt(Utilities::fixed_power<2>(sdt) +
                        Utilities::fixed_power<2>(2. * u_mag / h) +
                        9 * Utilities::fixed_power<2>(4 * kinematic_viscosity /
                                                      (h * h)));
}


#endif
