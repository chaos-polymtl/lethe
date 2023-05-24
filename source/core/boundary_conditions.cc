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

#include "core/boundary_conditions.h"

namespace BoundaryConditions
{
  extern template class BoundaryConditions<2>;
  extern template class BoundaryConditions<3>;
  extern template class NSBoundaryConditions<2>;
  extern template class NSBoundaryConditions<3>;
  extern template class HTBoundaryConditions<2>;
  extern template class HTBoundaryConditions<3>;
  extern template class TracerBoundaryConditions<2>;
  extern template class TracerBoundaryConditions<3>;
  extern template class VOFBoundaryConditions<2>;
  extern template class VOFBoundaryConditions<3>;
  extern template class CahnHilliardBoundaryConditions<2>;
  extern template class CahnHilliardBoundaryConditions<3>;
} // namespace BoundaryConditions
