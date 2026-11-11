// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_cahn_hilliard_boundary_conditions_h
#define lethe_cahn_hilliard_boundary_conditions_h

#include <core/boundary_conditions.h>

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <map>
#include <memory>

using namespace dealii;

namespace BoundaryConditions
{
  /**
   * @brief This class manages the functions associated with function-type boundary conditions
   * of the Cahn-Hilliard equations
   */
  template <int dim>
  class CahnHilliardBoundaryFunctions
  {
  public:
    // Cahn-Hilliard variables
    Functions::ParsedFunction<dim> phi;
    Functions::ParsedFunction<dim> eta;
  };

  /**
   * @brief This class manages the boundary conditions for the Cahn-Hilliard solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */

  template <int dim>
  class CahnHilliardBoundaryConditions : public BoundaryConditions
  {
  public:
    std::map<types::boundary_id, double> angle_of_contact;
    std::map<types::boundary_id, double> phase_dirichlet_value;
    std::map<types::boundary_id,
             std::shared_ptr<CahnHilliardBoundaryFunctions<dim>>>
      bcFunctions;

    /**
     * @brief Declares the default parameters for one Cahn-Hilliard boundary condition.
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] default_boundary_id The default mesh boundary id.
     */
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declares the Cahn-Hilliard boundary condition parameters.
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] number_of_boundary_conditions The number of boundary
     * conditions to declare.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the information for one Cahn-Hilliard boundary condition.
     *
     * @param[in] prm The parameter handler containing the boundary condition.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Parses the Cahn-Hilliard boundary condition parameters.
     *
     * @param[in] prm The parameter handler containing the boundary conditions.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
} // namespace BoundaryConditions


/**
 * @brief Class that implements a boundary conditions for the Cahn-Hilliard equation
 * where the phase and chemical potential are defined using individual functions
 */
template <int dim>
class CahnHilliardFunctionDefined : public Function<dim>
{
private:
  Functions::ParsedFunction<dim> *phi;

public:
  CahnHilliardFunctionDefined(Functions::ParsedFunction<dim> *p_phi)
    : Function<dim>(2)
    , phi(p_phi)
  {}

  virtual double
  value(const Point<dim> &p, unsigned int component) const override;
};


/**
 * @brief Calculates the value of a Function-type for the Cahn-Hilliard equations
 *
 * @param[in] p A point (generally a gauss point)
 * @param[in] component The vector component of the boundary condition (0-x, 1-y
 * and 2-z)
 */
template <int dim>
double
CahnHilliardFunctionDefined<dim>::value(const Point<dim>  &p,
                                        const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  if (component == 0)
    {
      return phi->value(p);
    }
  if (component == 1)
    {
      return 0.;
    }
  return 0.;
}

#endif
