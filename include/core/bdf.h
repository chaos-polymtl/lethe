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

#ifndef lethe_bdf_h
#define lethe_bdf_h

#include <deal.II/lac/vector.h>

#include <vector>

using namespace dealii;

/**
 * @brief Calculate the coefficients required for BDF integration of order n.
 * The coefficients are determined through a recursion algorithm.
 * The algorithm is taken from : Hay, Alexander, et al. "hp-Adaptive time
 * integration based on the BDF for viscous flows." Journal of Computational
 * Physics 291 (2015): 151-176.
 *
 * @param order The order of the BDF method. The BDF method of order n requires n+1 arrays
 *
 * @param time_steps a vector containing all the time steps. The time steps should be in reverse order.
 * For example, if the method is a BDF2, it uses three values for the time : (t)
 * (t-dt_1) and (t-dt_2). Thus the time step vector should contain dt_1 and
 * dt_2.
 */
Vector<double>
bdf_coefficients(unsigned int order, const std::vector<double> time_steps);


Vector<double>
delta(unsigned int order, unsigned int n, unsigned int j, Vector<double> times);

#endif
