// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/fluid_dynamics_boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

DeclException1(NavierStokesBoundaryDuplicated,
               types::boundary_id,
               << "Fluid dynamics boundary id: " << arg1
               << " has already been declared as a boundary condition");

namespace BoundaryConditions
{
  template <int dim>
  void
  FluidDynamicsBoundaryConditions<dim>::declare_default_entry(
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

  template <int dim>
  void
  FluidDynamicsBoundaryConditions<dim>::declare_parameters(
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

  template <int dim>
  void
  FluidDynamicsBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
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
  FluidDynamicsBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
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

  template class FluidDynamicsBoundaryConditions<2>;
  template class FluidDynamicsBoundaryConditions<3>;
} // namespace BoundaryConditions

template <int dim>
double
FluidDynamicsVelocityFunctionDefined<dim>::value(
  const Point<dim>  &point,
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

template <int dim>
double
FluidDynamicsPressureFunctionDefined<dim>::value(
  const Point<dim>  &point,
  const unsigned int component) const
{
  if (std::cmp_equal(component, dim))
    {
      return p->value(point);
    }
  return 0.;
}

template class FluidDynamicsVelocityFunctionDefined<2>;
template class FluidDynamicsVelocityFunctionDefined<3>;
template class FluidDynamicsPressureFunctionDefined<2>;
template class FluidDynamicsPressureFunctionDefined<3>;
