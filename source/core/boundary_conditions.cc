// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/cahn_hilliard_boundary_conditions.h>
#include <core/cls_boundary_conditions.h>
#include <core/fluid_dynamics_boundary_conditions.h>
#include <core/heat_transfer_boundary_conditions.h>
#include <core/time_harmonic_maxwell_boundary_conditions.h>
#include <core/tracer_boundary_conditions.h>

namespace BoundaryConditions
{
  extern template class FluidDynamicsBoundaryConditions<2>;
  extern template class FluidDynamicsBoundaryConditions<3>;
  extern template class HeatTransferBoundaryConditions<2>;
  extern template class HeatTransferBoundaryConditions<3>;
  extern template class TracerBoundaryConditions<2>;
  extern template class TracerBoundaryConditions<3>;
  extern template class CLSBoundaryConditions<2>;
  extern template class CLSBoundaryConditions<3>;
  extern template class CahnHilliardBoundaryConditions<2>;
  extern template class CahnHilliardBoundaryConditions<3>;
  extern template class TimeHarmonicMaxwellBoundaryConditions<2>;
  extern template class TimeHarmonicMaxwellBoundaryConditions<3>;
} // namespace BoundaryConditions
