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
#include <dem/distributions.h>

#ifndef find_maximum_particle_size_h
#  define find_maximum_particle_size_h

/**
 * @brief Find maximum particle size based on defined size values in the parameter
 * handler file. This value is extremely important in polydisperse systems. Loop
 * over all the distribution and return the maximum diameter value.
 *
 * @param lagrangian_physical_properties DEM physical properties declared in the
 * .prm file
 * @param distribution_object_container Contains all of the distribution per particle type.
 *
 */
double
find_maximum_particle_size(const std::vector<std::shared_ptr<Distribution>>
                             &distribution_object_container);

#endif
