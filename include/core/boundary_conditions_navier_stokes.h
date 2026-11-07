// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_boundary_conditions_navier_stokes_h
#define lethe_boundary_conditions_navier_stokes_h

#include <core/boundary_conditions.h>
#include <core/periodic_boundary.h>
#include <core/utilities.h>

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace dealii;

DeclException1(NavierStokesBoundaryDuplicated,
               types::boundary_id,
               << "Fluid dynamics boundary id: " << arg1
               << " has already been declared as a boundary condition");

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
   *
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
   * @brief This class manages the boundary conditions for Navier-Stokes solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */
  template <int dim>
  class NavierStokesBoundaryConditions : public BoundaryConditions
  {
  public:
    /// Functions for (u : x component of the velocity, v : y component of the
    /// velocity
    ///, w : z-component of the velocity, p : pressure, Neumann traction:
    /// traction_fn) for all boundaries
    std::map<types::boundary_id,
             std::shared_ptr<NavierStokesBoundaryFunctions<dim>>>
      navier_stokes_functions;

    void
    parse_boundary(ParameterHandler &prm);
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declares the Navier-Stokes boundary conditions
     *
     * @param prm The parameter file
     * @param number_of_boundary_conditions The number of boundary conditions to be declared. This parameter is generally pre-parsed from a first read of the prm file.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the Navier-Stokes boundary conditions
     *
     * @param prm The parameter file
     */
    void
    parse_parameters(ParameterHandler &prm);
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
  NavierStokesBoundaryConditions<dim>::createNoSlip()
  {
    this->type[0]                       = BoundaryType::noslip;
    this->number_of_boundary_conditions = 1;
    this->beta[0]                       = 0;
  }


  /**
   * @brief Declares the default parameter for a boundary condition id i_bc
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param default_boundary_id Default value of the boundary id. This corresponds to the number of the boundary condition subsection.
   */
  template <int dim>
  void
  NavierStokesBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    declare_enum_entry(prm,
                       "type",
                       BoundaryType::none,
                       "Type of boundary condition.",
                       deprecated_boundary_type_names,
                       fluid_dynamics_boundary_types);


    prm.declare_entry("id",
                      Utilities::to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions.");

    prm.declare_entry(
      "periodic id",
      "-1",
      Patterns::Integer(),
      "Mesh id for periodic face matching. Default entry is -1 to ensure that the periodic id is set by the user");

    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    // Create a dummy NavierStokesBoundaryFunctions object to declare the
    // appropriate parameters for this boundary condition.
    NavierStokesBoundaryFunctions<dim> temporary_fluid_dynamics_functions;


    prm.enter_subsection("u");
    temporary_fluid_dynamics_functions.u.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("v");
    temporary_fluid_dynamics_functions.v.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("w");
    temporary_fluid_dynamics_functions.w.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("p");
    temporary_fluid_dynamics_functions.p.declare_parameters(prm);
    prm.leave_subsection();

    /*
      Vector-valued function for the traction with dim components for the
      traction vector and 1 component for the pressure value at the boundary
      which is dummy and is needed to properly assemble the traction term as
      FEFaceIntegrator uses MatrixFree which in turn is based on the DoFHandler
      which in turn is based on the FE system that has dim velocity components
      and 1 pressure component.
    */

    prm.enter_subsection("traction fn");
    temporary_fluid_dynamics_functions.traction_fn.declare_parameters(prm,
                                                                      dim + 1);
    prm.leave_subsection();

    // Center of rotation of the boundary condition for torque calculation
    prm.enter_subsection("center of rotation");
    prm.declare_entry("x", "0", Patterns::Double(), "X COR");
    prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
    prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
    prm.leave_subsection();

    prm.declare_entry(
      "beta",
      "1",
      Patterns::Double(),
      "Penalty parameter for weak boundary condition imposed through Nitsche's method or outlets");

    prm.declare_entry(
      "boundary layer thickness",
      "0",
      Patterns::Double(),
      "Thickness of the boundary layer used to calculate the penalty parameter for partial slip boundary condition in tangent direction imposed through Nitsche's method or outlets");
  }


  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  NavierStokesBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the fluid dynamics boundary conditions. Please ensure that the id is set for every boundary condition."));

    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    NavierStokesBoundaryDuplicated(boundary_id));

        // Allocate the navier_stokes_functions object for every boundary
        // condition to ensure that they have a defined function and a center of
        // rotation.
        navier_stokes_functions[boundary_id] =
          std::make_shared<NavierStokesBoundaryFunctions<dim>>();

        prm.enter_subsection("u");
        navier_stokes_functions[boundary_id]->u.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("v");
        navier_stokes_functions[boundary_id]->v.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("w");
        navier_stokes_functions[boundary_id]->w.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("p");
        navier_stokes_functions[boundary_id]->p.parse_parameters(prm);
        prm.leave_subsection();


        prm.enter_subsection("traction fn");
        navier_stokes_functions[boundary_id]->traction_fn.parse_parameters(prm);
        prm.leave_subsection();


        prm.enter_subsection("center of rotation");
        navier_stokes_functions[boundary_id]->center_of_rotation[0] =
          prm.get_double("x");
        navier_stokes_functions[boundary_id]->center_of_rotation[1] =
          prm.get_double("y");

        if (dim == 3)
          navier_stokes_functions[boundary_id]->center_of_rotation[2] =
            prm.get_double("z");
        prm.leave_subsection();

        // Establish the type of boundary condition
        const BoundaryType boundary_type =
          string_to_enum<BoundaryType>(prm.get("type"),
                                       deprecated_boundary_type_names,
                                       "type");
        this->type[boundary_id] = boundary_type;
        if (boundary_type == BoundaryType::periodic)
          {
            types::boundary_id periodic_boundary_id =
              prm.get_integer("periodic id");

            this->type[boundary_id] = BoundaryType::periodic;

            // We attribute a periodic neighbor boundary type to the neighbor
            // boundary to ensure that all boundaries have a defined type
            this->type[periodic_boundary_id] = BoundaryType::periodic_neighbor;

            // We store the periodic id and direction
            this->periodic_boundaries[boundary_id] = {
              .neighbor_id = periodic_boundary_id,
              .direction   = static_cast<unsigned int>(
                prm.get_integer("periodic direction"))};

            // Allocate the navier_stokes_functions object for the periodic
            // neighbor condition to ensure that they have a defined function
            // and a center of rotation. Parse the information of the boundary
            // id to duplicate it.
            navier_stokes_functions[periodic_boundary_id] =
              std::make_shared<NavierStokesBoundaryFunctions<dim>>();

            prm.enter_subsection("u");
            navier_stokes_functions[periodic_boundary_id]->u.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("v");
            navier_stokes_functions[periodic_boundary_id]->v.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("w");
            navier_stokes_functions[periodic_boundary_id]->w.parse_parameters(
              prm);
            prm.leave_subsection();


            prm.enter_subsection("traction fn");
            navier_stokes_functions[periodic_boundary_id]
              ->traction_fn.parse_parameters(prm);
            prm.leave_subsection();


            prm.enter_subsection("p");
            navier_stokes_functions[periodic_boundary_id]->p.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("center of rotation");
            navier_stokes_functions[periodic_boundary_id]
              ->center_of_rotation[0] = prm.get_double("x");
            navier_stokes_functions[periodic_boundary_id]
              ->center_of_rotation[1] = prm.get_double("y");

            if (dim == 3)
              navier_stokes_functions[periodic_boundary_id]
                ->center_of_rotation[2] = prm.get_double("z");
            prm.leave_subsection();
          }
        this->beta[boundary_id] = prm.get_double("beta");
        this->boundary_layer_thickness[boundary_id] =
          prm.get_double("boundary layer thickness");
      }
  }


  template <int dim>
  void
  NavierStokesBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions");
    {
      prm.declare_entry("number",
                        Utilities::int_to_string(number_of_boundary_conditions),
                        Patterns::Integer(),
                        "Number of boundary conditions");
      prm.declare_entry(
        "time dependent",
        "false",
        Patterns::Bool(),
        "Bool to define if the boundary condition is time dependent");

      prm.declare_entry(
        "fix pressure constant",
        "false",
        Patterns::Bool(),
        "Bool to define if zero pressure is used as a constraint one node");

      for (unsigned int n = 0; n < number_of_boundary_conditions; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            declare_default_entry(prm, n);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }


  /**
   * @brief Parse the boundary conditions
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  NavierStokesBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions");
    {
      this->number_of_boundary_conditions = prm.get_integer("number");
      this->time_dependent                = prm.get_bool("time dependent");
      this->fix_pressure_constant = prm.get_bool("fix pressure constant");

      for (unsigned int n = 0; n < this->number_of_boundary_conditions; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            parse_boundary(prm);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
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
 * @param point A point at which the function will be evaluated
 *
 * @param component The vector component of the boundary condition (0-x, 1-y and 2-z)
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
 * @param point A point (generally a gauss point)
 *
 * @param component The vector component of the boundary condition (0-x, 1-y and 2-z)
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
