// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_boundary_conditions_h
#define lethe_boundary_conditions_h

#include <core/periodic_boundary.h>
#include <core/utilities.h>

#include <map>

using namespace dealii;

namespace BoundaryConditions
{
  enum class BoundaryType
  {
    // common
    none,
    outlet, // outlet is used for fluid dynamics, tracer and eventually other
            // physics
    periodic,
    periodic_neighbor, // The periodic neighbour is used to indicate a boundary
                       // which matches with a main periodic boundary

    // for fluid
    noslip,
    slip,
    function,
    function_weak,
    partial_slip,
    pressure,
    neumann_traction,
    // for heat transfer
    noflux,
    temperature,
    convection_radiation,
    // for tracer
    tracer_dirichlet,
    // for cls
    cls_dirichlet,
    // for cahn hilliard
    cahn_hilliard_noflux,
    cahn_hilliard_dirichlet_phase_order,
    cahn_hilliard_angle_of_contact,
    cahn_hilliard_free_angle,
    // for time harmonic maxwell
    pec,
    pmc,
    silver_muller,
    electric_field,
    magnetic_field,
    impedance_boundary,
    waveguide_port
  };

  /// Deprecated strings of the boundary condition "type" parameters, which
  /// remain accepted in parameter files
  inline const DeprecatedEnumNames<BoundaryType>
    deprecated_boundary_type_names = {
      {"function weak", BoundaryType::function_weak},
      {"partial slip", BoundaryType::partial_slip},
      {"neumann traction", BoundaryType::neumann_traction},
      {"convection-radiation-flux", BoundaryType::convection_radiation},
      {"silver muller", BoundaryType::silver_muller},
      {"electric field", BoundaryType::electric_field},
      {"magnetic field", BoundaryType::magnetic_field},
      {"impedance boundary", BoundaryType::impedance_boundary},
      {"waveguide port", BoundaryType::waveguide_port}};

  /**
   * @brief This class is the base class for all boundary conditions. It stores
   * the general information that all boundary condition share.
   * In Lethe, boundary conditions are identified with an id, a type and, in the
   * special case of periodic boundary condition, a periodic matching id
   * (periodic_id) and a periodic direction (0, 1 or 2).
   */
  class BoundaryConditions
  {
  public:
    /// Map containing the boundary id and the boundary type
    std::map<types::boundary_id, BoundaryType> type;

    /// Map containing the penalization parameter for weak dirichlet BCs and
    /// outlets for a corresponding id.
    std::map<types::boundary_id, double> beta;

    /// Map containing the boundary layer thickness tangent component parameter
    /// for partial slip for a corresponding id.
    std::map<types::boundary_id, double> boundary_layer_thickness;

    /// Number of boundary conditions
    unsigned int number_of_boundary_conditions;

    /// indicator for transient BCs
    bool time_dependent;

    /// Periodic boundary pairs, keyed by the principal periodic boundary id
    Parameters::PeriodicBoundaries periodic_boundaries;
  };

} // namespace BoundaryConditions

#endif
