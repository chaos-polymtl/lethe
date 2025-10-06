// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_sdirk_h
#define lethe_sdirk_h

#include <core/parameters.h>

#include <deal.II/lac/full_matrix.h>

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
