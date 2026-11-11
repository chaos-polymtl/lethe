// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cls_boundary_conditions_h
#define lethe_cls_boundary_conditions_h

#include <core/boundary_conditions.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <map>
#include <memory>

using namespace dealii;

namespace BoundaryConditions
{
  /**
   * @brief This class manages the boundary conditions for CLS solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */

  template <int dim>
  class CLSBoundaryConditions : public BoundaryConditions
  {
  public:
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      phase_indicator;

    /**
     * @brief Declares the default parameters for one CLS boundary condition.
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] default_boundary_id The default mesh boundary id.
     */
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declares the CLS boundary condition parameters.
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] number_of_boundary_conditions The number of boundary
     * conditions to declare.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the information for one CLS boundary condition.
     *
     * @param[in] prm The parameter handler containing the boundary condition.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Parses the CLS boundary condition parameters.
     *
     * @param[in] prm The parameter handler containing the boundary conditions.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };



} // namespace BoundaryConditions

#endif
