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

#include <cmath>

/**
 * @brief Calculate the stabilization parameter for the Navier-Stokes equations in steady-state
 * @return Value of the stabilization parameter - tau
 *
 * @param u_mag Magnitude of the velocity
 *
 * @param viscosity Kinematic viscosity
 *
 * @param h Cell size. Should be calculated using the diameter of a sphere of equal volume to that of the cell
 */
inline double
calculate_navier_stokes_tau_steady(const double u_mag,
                                   const double viscosity,
                                   const double h)
{
  return 1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                        9 * std::pow(4 * viscosity / (h * h), 2));
}

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
 * \param sdt Inverse of the time-step (1/dt)
 */

inline double
calculate_navier_stokes_tau_transient(const double u_mag,
                                      const double viscosity,
                                      const double h,
                                      const double sdt)
{
  return 1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                        9 * std::pow(4 * viscosity / (h * h), 2));
}
