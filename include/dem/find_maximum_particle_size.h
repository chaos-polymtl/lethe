/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#ifndef find_maximum_particle_size_h
#  define find_maximum_particle_size_h

/**
 * Find maximum particle size based on defined size values in the parameter
 * handler file. This value is extremely important in polydisperse systems. Note
 * that this value is calculated here using average value + 5 * standard
 * deviation of the diameter, since more than 99.999% of values in a normal
 * distribution lie within +-5 standard deviations. This value can be modified
 * using standard_deviation_multiplier variable.
 *
 * @param lagrangian_physical_properties DEM physical properties declared in the
 * .prm file
 * @param standard_deviation_multiplier Constant standard deviation multiplier.
 * It is defined in the dem constructor and it is equal to 2.5 to cover 99% of
 * the particle size distribution
 * @return Maximum particle size
 *
 */

double
find_maximum_particle_size(
  const Parameters::Lagrangian::LagrangianPhysicalProperties
    &           lagrangian_physical_properties,
  const double &standard_deviation_multiplier);

#endif
