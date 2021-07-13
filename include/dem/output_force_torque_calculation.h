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
 * Author: Lethe's community, 2021
 */

#include <dem/dem_solver_parameters.h>

#ifndef lethe_output_force_torque_calculation_h
#  define lethe_output_force_torque_calculation_h

/**
 * @brief write_forces_torques_output_results
 * Writes the results of force and torque calculations in a file, and depending
 * on the verbosity, in the terminal
 */
template <int dim>
void
write_forces_torques_output_results(
  std::map<unsigned int, Tensor<1, dim>> force_on_walls,
  std::map<unsigned int, Tensor<1, dim>> torque_on_walls,
  const ConditionalOStream &             pcout);

#endif // lethe_output_force_torque_calculation_h
