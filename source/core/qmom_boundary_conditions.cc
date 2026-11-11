// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/qmom_boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

#include <string>
#include <vector>

using namespace dealii;

DeclException1(QMOMBoundaryDuplicated,
               types::boundary_id,
               << "QMOM boundary id: " << arg1
               << " has already been declared as a boundary condition");

namespace BoundaryConditions
{
  template <int dim>
  void
  QMOMBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("dirichlet|outlet|none"),
                      "Type of boundary condition for the QMOM fields. "
                      "Choices are <dirichlet|outlet|none>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    prm.enter_subsection("dirichlet");
    Parameters::QMOMFunctions<dim>::declare_parameters(prm);
    prm.leave_subsection();
  }

  template <int dim>
  void
  QMOMBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions qmom");
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

      for (unsigned int i_bc = 0; i_bc < number_of_boundary_conditions; ++i_bc)
        {
          prm.enter_subsection("bc " + std::to_string(i_bc));
          {
            declare_default_entry(prm, i_bc);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  QMOMBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    const std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the QMOM boundary conditions. Please ensure that the id is set for every boundary condition."));

    const std::string boundary_type = prm.get("type");
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    QMOMBoundaryDuplicated(boundary_id));

        if (boundary_type == "dirichlet")
          {
            this->type[boundary_id] = BoundaryType::qmom_dirichlet;
            prm.enter_subsection("dirichlet");
            dirichlet[boundary_id].parse_parameters(prm);
            prm.leave_subsection();
          }
        else if (boundary_type == "outlet")
          {
            this->type[boundary_id] = BoundaryType::outlet;
          }
        else if (boundary_type == "none")
          {
            this->type[boundary_id] = BoundaryType::none;
          }
        else
          {
            AssertThrow(false,
                        ExcMessage(
                          "Unknown boundary condition type for QMOM: " +
                          boundary_type));
          }
      }
  }

  template <int dim>
  void
  QMOMBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions qmom");
    {
      this->number_of_boundary_conditions = prm.get_integer("number");
      this->time_dependent                = prm.get_bool("time dependent");

      for (unsigned int i_bc = 0; i_bc < this->number_of_boundary_conditions;
           ++i_bc)
        {
          prm.enter_subsection("bc " + std::to_string(i_bc));
          {
            parse_boundary(prm);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  template class QMOMBoundaryConditions<2>;
  template class QMOMBoundaryConditions<3>;
} // namespace BoundaryConditions
