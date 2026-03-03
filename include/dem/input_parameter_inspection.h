// SPDX-FileCopyrightText: Copyright (c) 2020, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
