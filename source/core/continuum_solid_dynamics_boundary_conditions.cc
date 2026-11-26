// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/continuum_solid_dynamics_boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

DeclException1(ContinuumSolidDynamicsBoundaryDuplicated,
               types::boundary_id,
               << "Continuum solid dynamics boundary id: " << arg1
               << " has already been declared as a boundary condition");

namespace BoundaryConditions
{
  template <int dim>
  void
  ContinuumSolidDynamicsBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    declare_enum_entry(
      prm,
      "type",
      ContinuumSolidDynamicsBoundaryType::noslip,
      "Type of boundary condition for the continuum solid dynamics.",
      {},
      continuum_solid_dynamics_boundary_types);

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    // Solid velocity (dim components) followed by the solid volume fraction
    prm.enter_subsection("function");
    Functions::ParsedFunction<dim>::declare_parameters(prm, dim + 1);
    prm.leave_subsection();

    prm.declare_entry(
      "periodic id",
      "-1",
      Patterns::Integer(),
      "Mesh id for periodic face matching. Default entry is -1 to ensure that the periodic id is set by the user");

    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");
  }

  template <int dim>
  void
  ContinuumSolidDynamicsBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions continuum solid dynamics");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");

      prm.declare_entry(
        "time dependent",
        "false",
        Patterns::Bool(),
        "Bool to define if the boundary condition is time-dependent");

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
  ContinuumSolidDynamicsBoundaryConditions<dim>::parse_boundary(
    ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the continuum solid dynamics boundary conditions. Please ensure that the id is set for every boundary condition."));

    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    ContinuumSolidDynamicsBoundaryDuplicated(boundary_id));

        // Allocate the function (solid velocity, solid volume fraction)
        function[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>(dim + 1);

        const ContinuumSolidDynamicsBoundaryType boundary_type =
          string_to_enum<ContinuumSolidDynamicsBoundaryType>(prm.get("type"));
        this->type[boundary_id] = boundary_type;

        if (boundary_type == ContinuumSolidDynamicsBoundaryType::function)
          {
            prm.enter_subsection("function");
            function[boundary_id]->parse_parameters(prm);
            prm.leave_subsection();
          }
        else if (boundary_type == ContinuumSolidDynamicsBoundaryType::periodic)
          {
            types::boundary_id periodic_boundary_id =
              prm.get_integer("periodic id");

            // The neighbor boundary gets a periodic neighbor type so that all
            // boundaries have a defined type
            this->type[periodic_boundary_id] =
              ContinuumSolidDynamicsBoundaryType::periodic_neighbor;

            this->periodic_boundaries[boundary_id] = {
              .neighbor_id = periodic_boundary_id,
              .direction   = static_cast<unsigned int>(
                prm.get_integer("periodic direction"))};
          }
      }
  }

  template <int dim>
  void
  ContinuumSolidDynamicsBoundaryConditions<dim>::parse_parameters(
    ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions continuum solid dynamics");
    {
      this->number_of_boundary_conditions = prm.get_integer("number");
      this->time_dependent                = prm.get_bool("time dependent");

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

  template class ContinuumSolidDynamicsBoundaryConditions<2>;
  template class ContinuumSolidDynamicsBoundaryConditions<3>;
} // namespace BoundaryConditions