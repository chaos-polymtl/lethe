// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/cahn_hilliard_boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <string>
#include <vector>

using namespace dealii;

DeclException1(CahnHilliardBoundaryDuplicated,
               types::boundary_id,
               << "Cahn Hilliard boundary id: " << arg1
               << " has already been declared as a boundary condition");

namespace BoundaryConditions
{
  template <int dim>
  void
  CahnHilliardBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection(
        "none|noflux|dirichlet|angle_of_contact|free_angle|periodic"),
      "Type of boundary condition for the Cahn-Hilliard equations"
      "Choices are <none|noflux|dirichlet|angle_of_contact|free_angle|periodic>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    prm.declare_entry(
      "angle value",
      "0",
      Patterns::Double(),
      "Inner angle of contact between the fluid 1 and the boundary (in degrees)");

    Functions::ParsedFunction<dim> temporary_function;
    prm.enter_subsection("phi");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Periodic boundary condition parameters for Cahn-Hilliards physics
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
  CahnHilliardBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions cahn hilliard");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      prm.declare_entry(
        "time dependent",
        "false",
        Patterns::Bool(),
        "Bool to define if the boundary condition is time dependent");

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
  CahnHilliardBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the Cahn Hilliard boundary conditions. Please ensure that the id is set for every boundary condition."));
    // Loop through all boundary ids to ensure that they are all non-negative
    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    CahnHilliardBoundaryDuplicated(boundary_id));

        // Create and parse the phase order boundary condition for all cases.
        prm.enter_subsection("phi");
        this->bcFunctions[boundary_id] =
          std::make_shared<CahnHilliardBoundaryFunctions<dim>>();
        bcFunctions[boundary_id]->phi.parse_parameters(prm);
        prm.leave_subsection();

        // Do the same with the angle of contact
        this->angle_of_contact[boundary_id] = prm.get_double("angle value");

        const std::string op = prm.get("type");
        if (op == "none")
          {
            this->type[boundary_id] = BoundaryType::none;
          }
        else if (op == "noflux")
          {
            this->type[boundary_id] = BoundaryType::cahn_hilliard_noflux;
          }
        else if (op == "dirichlet")
          {
            this->type[boundary_id] =
              BoundaryType::cahn_hilliard_dirichlet_phase_order;
          }
        else if (op == "angle_of_contact")
          {
            this->type[boundary_id] =
              BoundaryType::cahn_hilliard_angle_of_contact;
          }
        else if (op == "free_angle")
          {
            this->type[boundary_id] = BoundaryType::cahn_hilliard_free_angle;
          }
        else if (op == "periodic")
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
          }
        else
          {
            AssertThrow(
              false,
              ExcMessage("Unknown boundary condition type for Cahn-Hilliard."));
          }
      }
  }

  template <int dim>
  void
  CahnHilliardBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions cahn hilliard");
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

  template class CahnHilliardBoundaryConditions<2>;
  template class CahnHilliardBoundaryConditions<3>;
} // namespace BoundaryConditions
