/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_input_parameter_inspection_h
#define lethe_input_parameter_inspection_h

#include <dem/dem_solver_parameters.h>
#include <dem/distributions.h>

/**
 * @brief Check input parameters in the parameter handler to be in the
 * correct range. Warnings or error would appear if the parameters are not in
 * their acceptable range.
 *
 * @param dem_parameters Input DEM parameters in the parameter handler file
 * @param pcout Printing in parallel
 * @param size_distribution_object_container Contain all the types of distribution
 * being used for each type of particle.
 */

template <int dim>
void
input_parameter_inspection(const DEMSolverParameters<dim> &dem_parameters,
                           const ConditionalOStream       &pcout,
                           const std::vector<std::shared_ptr<Distribution>>
                             &size_distribution_object_container);

#endif
