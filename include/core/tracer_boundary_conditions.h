// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tracer_boundary_conditions_h
#define lethe_tracer_boundary_conditions_h

#include <core/boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <map>
#include <memory>
#include <vector>

using namespace dealii;

namespace BoundaryConditions
{
  /// Boundary condition types of the tracer. periodic_neighbor is not a user
  /// choice: it is assigned to the boundary matching a periodic boundary.
  enum class TracerBoundaryType
  {
    none,
    dirichlet,
    outlet,
    periodic,
    periodic_neighbor
  };

  /// Boundary condition types of the tracer which can be selected in the
  /// parameter file
  inline const std::vector<TracerBoundaryType> tracer_boundary_types = {
    TracerBoundaryType::none,
    TracerBoundaryType::dirichlet,
    TracerBoundaryType::outlet,
    TracerBoundaryType::periodic};

  /**
   * @brief This class manages the boundary conditions for Tracer solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */
  template <int dim>
  class TracerBoundaryConditions : public BoundaryConditions<TracerBoundaryType>
  {
  public:
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      tracer;

    /**
     * @brief Declare the default parameters for one tracer boundary condition.
     *
     * @param[in, out] prm The parameter handler object.
     * @param[in] default_boundary_id The default mesh boundary id.
     */
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declare the tracer boundary condition parameters.
     *
     * @param[in, out] prm The parameter handler object.
     * @param[in] number_of_boundary_conditions The number of boundary
     * conditions to declare.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parse the information for one tracer boundary condition.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Parse the tracer boundary condition parameters.
     *
     * @param[in] prm The parameter handler object.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace BoundaryConditions

#endif
