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

#ifndef LETHE_BDF_H
#define LETHE_BDF_H

#include <deal.II/lac/vector.h>
#include <vector>

using namespace dealii;

/**
 * Calculate the coefficients required for BDF integration of order n
 * from order $n=1$ to order $n=5$.
 * The formulas are derived analytically, but the coefficients
 * Could also be determined through recursion on the fly.
 */
Vector<double> bdf_coefficients(unsigned int order, std::vector<double> dt);
Vector<double> delta(unsigned int order, unsigned int n, unsigned int j,
                     Vector<double> times);

#endif
