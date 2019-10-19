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
 * Author: Simon Gauvin, Polytechnique Montreal, 2019
 */

#ifndef LETHE_PHYSICSSOLVER
#define LETHE_PHYSICSSOLVER

#include "parameters.h"

/**
 * An interface for all physics solver classes to derive from.
 */

class PhysicsSolver
{
public:
  virtual void
  assemble_matrix_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method) = 0;

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) = 0;

  virtual void
  solve_linear_system(const bool       initial_step,
                      const double     absolute_residual,
                      const double     relative_residual) = 0;
};

#endif