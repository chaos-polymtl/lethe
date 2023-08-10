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
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */


#ifndef lethe_sdirk_h
#define lethe_sdirk_h

#include <core/parameters.h>

#include <deal.II/lac/full_matrix.h>

#include <vector>

using namespace dealii;

/**
 * @brief Calculates the SDIRK coefficient
 *
 * @param order Order of the SDIRK method to be used
 *
 * @param time_step Value of the time step used for time integration.
 *
 * SDIRK 22 -  Coefficients from Kennedy and Carpenter 2016
 * Intermediary step is a t+ (1-sqrt(2)/2) * dt
 *
 * #SDIRK 33 - Coefficients from Kennedy and Carpenter 2016
 * gamma = 0.435866521508458999416019
 * b = 1.20849664917601007033648
 * c =0.717933260754229499708010
 *
 * Butcher's Tableau for SDIRK33
 * --------------------------
 * | gamma      0        0  |
 * |c-gamma   gamma      0  |
 * |   b    1-b-gamma  gamma|
 * --------------------------
 *
 * Intermediary 1 step is at
 * Intermediary 2 step is at
 *
 */
FullMatrix<double>
sdirk_coefficients(const unsigned int order, const double time_step);



#endif
