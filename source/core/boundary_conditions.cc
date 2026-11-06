// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/boundary_conditions_cahn_hilliard.h>
#include <core/boundary_conditions_cls.h>
#include <core/boundary_conditions_heat_transfer.h>
#include <core/boundary_conditions_navier_stokes.h>
#include <core/boundary_conditions_time_harmonic_maxwell.h>
#include <core/boundary_conditions_tracer.h>

namespace BoundaryConditions
{
  extern template class NSBoundaryConditions<2>;
  extern template class NSBoundaryConditions<3>;
  extern template class HTBoundaryConditions<2>;
  extern template class HTBoundaryConditions<3>;
  extern template class TracerBoundaryConditions<2>;
  extern template class TracerBoundaryConditions<3>;
  extern template class CLSBoundaryConditions<2>;
  extern template class CLSBoundaryConditions<3>;
  extern template class CahnHilliardBoundaryConditions<2>;
  extern template class CahnHilliardBoundaryConditions<3>;
  extern template class TimeHarmonicMaxwellBoundaryConditions<2>;
  extern template class TimeHarmonicMaxwellBoundaryConditions<3>;
} // namespace BoundaryConditions
