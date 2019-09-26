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

#include "solvers/boundaryconditions.h"

namespace BoundaryConditions
{
  extern template class NSBoundaryConditions<2>;
  extern template class NSBoundaryConditions<3>;
} // namespace BoundaryConditions

// extern template class PeriodicBoundaryValues<2>;
// extern template class PeriodicBoundaryValues<3>;
