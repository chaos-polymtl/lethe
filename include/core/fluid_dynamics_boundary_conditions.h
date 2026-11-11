// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_boundary_conditions_h
#define lethe_fluid_dynamics_boundary_conditions_h

#include <core/boundary_conditions.h>

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <map>
#include <memory>
#include <utility>
#include <vector>

using namespace dealii;

namespace BoundaryConditions
{

  /// Boundary condition types of the fluid dynamics
  inline const std::vector<BoundaryType> fluid_dynamics_boundary_types = {
    BoundaryType::none,
    BoundaryType::noslip,
    BoundaryType::slip,
    BoundaryType::function,
    BoundaryType::periodic,
    BoundaryType::pressure,
    BoundaryType::neumann_traction,
    BoundaryType::function_weak,
    BoundaryType::partial_slip,
    BoundaryType::outlet};

  /**
   * @brief This class manages the functions associated with function-type boundary conditions
   * of the Navier-Stokes equations
   */
  template <int dim>
  class NavierStokesBoundaryFunctions
  {
  public:
    // Velocity components
    Functions::ParsedFunction<dim> u;
    Functions::ParsedFunction<dim> v;
    Functions::ParsedFunction<dim> w;

    /// Pressure
    Functions::ParsedFunction<dim> p;

    // Neumann Traction σ = 2με(u)− pI
    // σ.n = g
    Functions::ParsedFunction<dim> traction_fn;

    /// Point for the center of rotation
    Point<dim> center_of_rotation;

    /*
      Constructor required to avoid vector-valued function t of ParsedFunction
      type which in-turn initializes the function_parser object with
      n_components = 1 by default if n_components is not passed.

      Moreover, the test and trial function is dim+1 components including
      pressure. Hence, we have the traction function as a vector-valued function
      with dim components for the traction vector and 1 component for the
      pressure value at the boundary which is dummy and is needed to properly
      assemble the traction term as FEFaceIntegrator uses MatrixFree which in
      turn is based on the DoFHandler which in turn is based on the FE system
      that has dim velocity components and 1 pressure component.

        t = (t₁, t₂, …, t_dim, 0.0 (dummy pressure component))

        φ = (φ₁, φ₂, …, φ_dim, φ_p)

        ∫_{Γ_t} = t · φ ds

    */

    NavierStokesBoundaryFunctions()
      : traction_fn(dim + 1) // The traction function has dim components for the
                             // traction vector plus 1 for the pressure.
    {}
  };

  /**
   * @brief This class manages the boundary conditions forhe fluid dynamics solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */
  template <int dim>
  class FluidDynamicsBoundaryConditions : public BoundaryConditions
  {
  public:
    /// Functions for (u : x component of the velocity, v : y component of the
    /// velocity
    ///, w : z-component of the velocity, p : pressure, Neumann traction:
    /// traction_fn) for all boundaries
    std::map<types::boundary_id,
             std::shared_ptr<NavierStokesBoundaryFunctions<dim>>>
      navier_stokes_functions;

    /**
     * @brief Parses the information for one fluid dynamics boundary condition.
     *
     * @param[in] prm The parameter handler containing the boundary condition.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Declares the default parameters for one fluid dynamics boundary condition.
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] default_boundary_id The default mesh boundary id.
     */
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declares the Navier-Stokes boundary conditions
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] number_of_boundary_conditions The number of boundary
     * conditions to declare.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the Navier-Stokes boundary conditions
     *
     * @param[in] prm The parameter handler containing the boundary conditions.
     */
    void
    parse_parameters(ParameterHandler &prm);

    /**
     * @brief Creates a noslip boundary condition for id=0
     * Used in tests only
     */
    void
    createNoSlip();

    /// Fix pressure constant using a node
    bool fix_pressure_constant;
  };


  /**
   * @brief Creates a noslip boundary condition for id=0
   * Used in tests only
   */
  template <int dim>
  void
  FluidDynamicsBoundaryConditions<dim>::createNoSlip()
  {
    this->type[0]                       = BoundaryType::noslip;
    this->number_of_boundary_conditions = 1;
    this->beta[0]                       = 0;
  }



} // namespace BoundaryConditions


/**
 * @brief This class implements a boundary conditions for the Navier-Stokes equation
 * where the velocity component are defined using individual functions
 */
template <int dim>
class NavierStokesFunctionDefined : public Function<dim>
{
  Functions::ParsedFunction<dim> *u;
  Functions::ParsedFunction<dim> *v;
  Functions::ParsedFunction<dim> *w;

public:
  NavierStokesFunctionDefined(Functions::ParsedFunction<dim> *p_u,
                              Functions::ParsedFunction<dim> *p_v,
                              Functions::ParsedFunction<dim> *p_w)
    : Function<dim>(dim + 1)
    , u(p_u)
    , v(p_v)
    , w(p_w)
  {}

  double
  value(const Point<dim> &point, const unsigned int component) const override;
};


/**
 * @brief Calculates the value of a Function-type Navier-Stokes equations
 *
 * @param[in] point A point at which the function will be evaluated
 * @param[in] component The vector component of the boundary condition (0-x, 1-y
 * and 2-z)
 */
template <int dim>
double
NavierStokesFunctionDefined<dim>::value(const Point<dim>  &point,
                                        const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  if (component == 0)
    {
      return u->value(point);
    }
  if (component == 1)
    {
      return v->value(point);
    }
  if (component == 2)
    {
      return w->value(point);
    }
  return 0.;
}

/**
 * @brief This class implements a pressure boundary condition for the Navier-Stokes equations.
 */
template <int dim>
class NavierStokesPressureFunctionDefined : public Function<dim>
{
  Functions::ParsedFunction<dim> *p;

public:
  NavierStokesPressureFunctionDefined(Functions::ParsedFunction<dim> *p_p)
    : Function<dim>(dim + 1)
    , p(p_p)
  {}

  double
  value(const Point<dim> &point, unsigned int component) const override;
};


/**
 * @brief Calculates the value of a Function-type Navier-Stokes equations
 *
 * @param[in] point A point (generally a gauss point)
 * @param[in] component The vector component of the boundary condition (0-x, 1-y
 * and 2-z)
 */
template <int dim>
double
NavierStokesPressureFunctionDefined<dim>::value(
  const Point<dim>  &point,
  const unsigned int component) const
{
  if (std::cmp_equal(component, dim))
    {
      return p->value(point);
    }
  return 0.;
}

#endif
