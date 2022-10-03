/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2022 -  by the Lethe authors
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
 */

#include <deal.II/base/utilities.h>

#include <cmath>

using namespace dealii;

/**
 * @brief Calculate the stabilization parameter for the Navier-Stokes equations in steady-state
 * @return Value of the stabilization parameter - tau
 *
 * @param u_mag Magnitude of the velocity
 *
 * @param viscosity Kinematic viscosity
 *
 * @param h Cell size. Should be calculated using the diameter of a sphere of equal volume to that of the cell
 *
 * @param density Density of the fluid, assumed to be 1 by default.
 *
 */
inline double
calculate_navier_stokes_gls_tau_steady(const double u_mag,
                                       const double viscosity,
                                       const double h,
                                       const double density = 1)
{
  return 1. / std::sqrt(Utilities::fixed_power<2>(2. * density * u_mag / h) +
                        9 * Utilities::fixed_power<2>(4 * density * viscosity /
                                                      (h * h)));
}

/**
 * @brief Calculate the stabilization parameter for the transient Navier-Stokes equations
 * @return Value of the stabilization parameter - tau
 *
 * @param u_mag Magnitude of the velocity
 *
 * @param viscosity Kinematic viscosity
 *
 * @param h Cell size. Should be calculated using the diameter of a sphere of equal volume to that of the cell
 *
 * @param sdt Inverse of the time-step (1/dt)
 *
 * @param density Density of the fluid, assumed to be 1 by default.
 */

inline double
calculate_navier_stokes_gls_tau_transient(const double u_mag,
                                          const double viscosity,
                                          const double h,
                                          const double sdt,
                                          const double density = 1)
{
  return 1. / std::sqrt(Utilities::fixed_power<2>(density * sdt) +
                        Utilities::fixed_power<2>(2. * density * u_mag / h) +
                        9 * Utilities::fixed_power<2>(4 * density * viscosity /
                                                      (h * h)));
}
