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

#include <deal.II/lac/full_matrix.h>

#include <vector>

using namespace dealii;

/**
 * Returns the SDIRK coefficient for a given order of the scheme
 */

FullMatrix<double>
sdirk_coefficients(const unsigned int order, const double time_step);

#endif
