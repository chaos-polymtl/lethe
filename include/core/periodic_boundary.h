// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_periodic_boundary_h
#define lethe_periodic_boundary_h

#include <deal.II/base/types.h>

#include <map>

namespace Parameters
{
  /**
   * @brief A periodic boundary pair, keyed elsewhere by its principal boundary id.
   */
  struct PeriodicBoundary
  {
    /// Neighbor periodic boundary id of the pair.
    dealii::types::boundary_id neighbor_id;

    /// Direction of periodicity (0, 1 or 2) of the pair.
    unsigned int direction;
  };

  /**
   * @brief The periodic boundary pairs of a simulation, keyed by the principal
   * (first) boundary id of each pair.
   *
   * CFD and DEM simulations declare their boundary conditions through different
   * parameter classes (BoundaryConditions::BoundaryConditions and
   * Parameters::Lagrangian::BCDEM respectively), but both describe periodicity
   * with this container, so that a single routine,
   * setup_periodic_boundary_conditions(), sets up the periodicity of the
   * triangulation for either solver.
   */
  using PeriodicBoundaries =
    std::map<dealii::types::boundary_id, PeriodicBoundary>;
} // namespace Parameters

#endif
