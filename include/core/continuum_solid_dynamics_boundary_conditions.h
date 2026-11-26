// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_continuum_solid_dynamics_boundary_conditions_h
#define lethe_continuum_solid_dynamics_boundary_conditions_h

#include <core/boundary_conditions.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/types.h>

#include <map>
#include <memory>
#include <vector>

using namespace dealii;

namespace BoundaryConditions
{
  /// Boundary condition types of the continuum solid dynamics. "none" imposes
  /// no solid velocity constraint (the boundary stays impermeable to the solid
  /// phase, as the continuity equation has no boundary flux there).
  /// periodic_neighbor is not a user choice: it is assigned to the boundary
  /// matching a periodic boundary.
  enum class ContinuumSolidDynamicsBoundaryType
  {
    none,
    noslip,
    slip,
    function,
    outlet,
    periodic,
    periodic_neighbor
  };

  /// Boundary condition types of the continuum solid dynamics which can be
  /// selected in the parameter file
  inline const std::vector<ContinuumSolidDynamicsBoundaryType>
    continuum_solid_dynamics_boundary_types = {
      ContinuumSolidDynamicsBoundaryType::none,
      ContinuumSolidDynamicsBoundaryType::noslip,
      ContinuumSolidDynamicsBoundaryType::slip,
      ContinuumSolidDynamicsBoundaryType::function,
      ContinuumSolidDynamicsBoundaryType::outlet,
      ContinuumSolidDynamicsBoundaryType::periodic};

  /**
   * @brief Boundary conditions of the continuum solid dynamics (solid phase of
   * the Euler-Euler model).
   *
   * On "function" boundaries, the solid velocity and the solid volume fraction
   * are imposed through one function with dim + 1 components: the dim
   * components of the solid velocity, followed by the solid volume fraction.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class ContinuumSolidDynamicsBoundaryConditions
    : public BoundaryConditions<ContinuumSolidDynamicsBoundaryType>
  {
  public:
    /// Function imposed on "function" boundaries (solid velocity, solid
    /// volume fraction)
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      function;

    /**
     * @brief Declare the entries of one boundary condition.
     *
     * @param[in,out] prm Parameter handler.
     * @param[in] default_boundary_id Default boundary id of this condition.
     */
    void
    declare_default_entry(ParameterHandler        &prm,
                          const types::boundary_id default_boundary_id);

    /**
     * @brief Declare the "boundary conditions continuum solid dynamics"
     * subsection.
     *
     * @param[in,out] prm Parameter handler.
     * @param[in] number_of_boundary_conditions Number of boundary conditions
     * to declare.
     */
    void
    declare_parameters(ParameterHandler  &prm,
                       const unsigned int number_of_boundary_conditions);

    /**
     * @brief Parse one boundary condition.
     *
     * @param[in] prm Parameter handler.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Parse the "boundary conditions continuum solid dynamics"
     * subsection.
     *
     * @param[in] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace BoundaryConditions

#endif