// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_qmom_boundary_conditions_h
#define lethe_qmom_boundary_conditions_h

#include <core/boundary_conditions.h>
#include <core/parameters_multiphysics.h>

#include <deal.II/base/parameter_handler.h>

#include <map>

using namespace dealii;

namespace BoundaryConditions
{
  /**
   * @brief Boundary conditions of the QMOM physics. They apply to the advective
   * flux of all the QMOM fields (scalars and moments), which is meant to be
   * computed with an upwind numerical flux:
   *  - dirichlet: the prescribed values enter the domain where the flow enters
   *    (u.n < 0) and the interior values leave it where the flow leaves;
   *  - outlet: the interior values leave the domain where the flow leaves and
   *    nothing enters where the flow enters (backflow);
   *  - none: the advective flux is zero, which is the natural choice for
   *    walls.
   * Boundaries without a QMOM boundary condition are treated as none.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class QMOMBoundaryConditions : public BoundaryConditions
  {
  public:
    /// Values imposed by the Dirichlet boundary conditions, for each boundary
    /// id of type qmom_dirichlet
    std::map<types::boundary_id, Parameters::QMOMFunctions<dim>> dirichlet;

    /**
     * @brief Declare the parameters of one boundary condition.
     *
     * @param[in,out] prm Parameter handler.
     * @param[in] default_boundary_id Default boundary id of the boundary
     * condition.
     */
    void
    declare_default_entry(ParameterHandler        &prm,
                          const types::boundary_id default_boundary_id);

    /**
     * @brief Declare the "boundary conditions qmom" subsection.
     *
     * @param[in,out] prm Parameter handler.
     * @param[in] number_of_boundary_conditions Maximal number of boundary
     * conditions, used to declare the "bc i" subsections.
     */
    void
    declare_parameters(ParameterHandler  &prm,
                       const unsigned int number_of_boundary_conditions);

    /**
     * @brief Parse one boundary condition.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Parse the "boundary conditions qmom" subsection.
     *
     * @param[in,out] prm Parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace BoundaryConditions

#endif
