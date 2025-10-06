// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/boundary_conditions.h>

namespace BoundaryConditions
{
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
