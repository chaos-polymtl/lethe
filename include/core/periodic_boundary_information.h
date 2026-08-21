// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_periodic_boundary_information_h
#define lethe_periodic_boundary_information_h

#include <deal.II/base/types.h>

#include <map>

namespace Parameters
{
  /**
   * @brief The subset of boundary condition information required to set up the
   * periodicity of a triangulation.
   *
   * CFD and DEM simulations declare their boundary conditions through different
   * parameter classes (BoundaryConditions::BoundaryConditions and
   * Parameters::Lagrangian::BCDEM respectively), but both describe periodicity
   * with the same two maps. Both classes derive from this structure so that a
   * single routine, setup_periodic_boundary_conditions(), can set up the
   * periodicity of the triangulation for either solver.
   *
   * A periodic boundary pair is stored once, keyed by the principal (first)
   * boundary id of the pair.
   */
  struct PeriodicBoundaryInformation
  {
    /// Neighbor periodic boundary id of each periodic boundary pair, keyed by
    /// the principal periodic boundary id.
    std::map<dealii::types::boundary_id, dealii::types::boundary_id>
      periodic_neighbor_id;

    /// Direction of periodicity (0, 1 or 2) of each periodic boundary pair,
    /// keyed by the principal periodic boundary id.
    std::map<dealii::types::boundary_id, unsigned int> periodic_direction;
  };
} // namespace Parameters

#endif
