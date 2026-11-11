// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/heat_transfer_boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

DeclException1(HeatTransferBoundaryDuplicated,
               types::boundary_id,
               << "Heat transfer boundary id: " << arg1
               << " has already been declared as a boundary condition");

namespace BoundaryConditions
{

  template <int dim>
  void
  HeatTransferBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    declare_enum_entry(prm,
                       "type",
                       BoundaryType::noflux,
                       "Type of boundary condition for heat transfer.",
                       deprecated_boundary_type_names,
                       heat_transfer_boundary_types);

    prm.declare_entry("id",
                      Utilities::to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    Functions::ParsedFunction<dim> temporary_function;

    // Expression for the temperature for an imposed temperature at bc
    prm.enter_subsection("value");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Expression for the h coefficient of convection-radiation-flux bc
    prm.enter_subsection("h");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Temperature of environment for convection-radiation-flux bc
    prm.enter_subsection("Tinf");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Emissivity of the boundary for convection-radiation-flux bc
    prm.enter_subsection("emissivity");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Heat flux (Neumann) at the boundary for convection-radiation-flux bc
    prm.enter_subsection("heat_flux");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Periodic boundary condition parameters for HT physics
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
  HeatTransferBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions heat transfer");
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
      prm.declare_entry("Stefan-Boltzmann constant",
                        "0.000000056703",
                        Patterns::Double(),
                        "Stefan-Boltzmann constant");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  HeatTransferBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the heat transfer boundary conditions. Please ensure that the id is set for every boundary condition."));

    // Loop through all boundary ids to ensure that they are all non-negative
    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    HeatTransferBoundaryDuplicated(boundary_id));


        const BoundaryType boundary_type =
          string_to_enum<BoundaryType>(prm.get("type"),
                                       deprecated_boundary_type_names,
                                       "type");
        this->type[boundary_id] = boundary_type;
        if (boundary_type == BoundaryType::convection_radiation)
          {
            this->has_convection_radiation_bc = true;

            // Emissivity validity (0 <= emissivity <= 1) will be checked at
            // evaluation.
          }
        else if (boundary_type == BoundaryType::periodic)
          {
            types::boundary_id periodic_boundary_id =
              prm.get_integer("periodic id");

            // We attribute a periodic neighbor boundary type to the neighbor
            // boundary to ensure that all boundaries have a defined type
            this->type[periodic_boundary_id] = BoundaryType::periodic_neighbor;

            // We store the periodic id and direction
            this->periodic_boundaries[boundary_id] = {
              .neighbor_id = periodic_boundary_id,
              .direction   = static_cast<unsigned int>(
                prm.get_integer("periodic direction"))};
          }

        // All the functions are parsed since they might be used for
        // post-processing
        prm.enter_subsection("value");
        this->dirichlet_value[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->dirichlet_value[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();
        prm.enter_subsection("h");
        this->h[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->h[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();
        prm.enter_subsection("Tinf");
        this->Tinf[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->Tinf[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();
        prm.enter_subsection("emissivity");
        this->emissivity[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->emissivity[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();
        prm.enter_subsection("heat_flux");
        this->heat_flux_bc[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->heat_flux_bc[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();
      }
  }

  template <int dim>
  void
  HeatTransferBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions heat transfer");
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
      this->Stefan_Boltzmann_constant =
        prm.get_double("Stefan-Boltzmann constant");
    }
    prm.leave_subsection();
  }

  template class HeatTransferBoundaryConditions<2>;
  template class HeatTransferBoundaryConditions<3>;
} // namespace BoundaryConditions
