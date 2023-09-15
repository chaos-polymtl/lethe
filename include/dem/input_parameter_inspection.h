/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <dem/dem_solver_parameters.h>

using namespace std;

#ifndef input_parameter_inspection_h
#  define input_parameter_inspection_h

/**
 * Manages checking input parameters in the parameter handler to be in the
 * correct range. Warnings or error would appear if the parameters are not in
 * their acceptable range.
 *
 * @param dem_parameters Input DEM parameters in the parameter handler file
 * @param pcout Printing in parallel
 * @param standard_deviation_multiplier Constant standard deviation multiplier.
 * It is defined in the dem constructor and it is equal to 2.5 to cover 99% of
 * the particle size distribution
 *
 */

template <int dim>
void
input_parameter_inspection(const DEMSolverParameters<dim> &dem_parameters,
                           const ConditionalOStream       &pcout,
                           const double standard_deviation_multiplier);

#endif /* input_parameter_inspection_h */
