// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_boundary_conditions_h
#define lethe_boundary_conditions_h

#include <core/periodic_boundary.h>

#include <deal.II/base/types.h>

#include <map>
#include <type_traits>

using namespace dealii;

namespace BoundaryConditions
{
  /**
   * @brief This class is the base class for all boundary conditions. It stores
   * the general information that all boundary condition share.
   * In Lethe, boundary conditions are identified with an id, a type and, in the
   * special case of periodic boundary condition, a periodic matching id
   * (periodic_id) and a periodic direction (0, 1 or 2).
   *
   * @tparam BoundaryType Enumeration of the boundary condition types of the
   * physics (e.g. FluidDynamicsBoundaryType).
   */
  template <typename BoundaryType>
    requires std::is_enum_v<BoundaryType>
  class BoundaryConditions
  {
  public:
    /// Map containing the boundary id and the boundary type
    std::map<types::boundary_id, BoundaryType> type;

    /// Number of boundary conditions
    unsigned int number_of_boundary_conditions;

    /// indicator for transient BCs
    bool time_dependent;

    /// Periodic boundary pairs, keyed by the principal periodic boundary id
    Parameters::PeriodicBoundaries periodic_boundaries;
  };

} // namespace BoundaryConditions

#endif
