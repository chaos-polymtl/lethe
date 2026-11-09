// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_boundary_conditions_cls_h
#define lethe_boundary_conditions_cls_h

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

DeclException1(CLSBoundaryDuplicated,
               types::boundary_id,
               << "CLS boundary id: " << arg1
               << " has already been declared as a boundary condition");


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

    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);
    void
    parse_boundary(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Declares the default parameters for a boundary condition id i_bc
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param default_boundary_id Default value given to the boundary id.
   */
  template <int dim>
  void
  CLSBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("none|dirichlet|periodic"),
                      "Type of boundary condition for CLS"
                      "Choices are <none|dirichlet|periodic>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    Functions::ParsedFunction<dim> temporary_function;
    prm.enter_subsection("dirichlet");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Periodic boundary condition parameters for CLS physics
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

  /**
   * @brief Declare the boundary conditions default parameters
   * Calls declareDefaultEntry method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param number_of_boundary_conditions Number of boundary conditions
   */
  template <int dim>
  void
  CLSBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions CLS");
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

  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

  template <int dim>
  void
  CLSBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the CLS boundary conditions. Please ensure that the id is set for every boundary condition."));

    // Loop through all boundary ids to ensure that they are all non-negative
    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    CLSBoundaryDuplicated(boundary_id));

        if (auto const option = prm.get("type"); option == "none")
          {
            this->type[boundary_id] = BoundaryType::none;
          }

        if (auto const option = prm.get("type"); option == "dirichlet")
          {
            this->type[boundary_id] = BoundaryType::cls_dirichlet;
            prm.enter_subsection("dirichlet");
            phase_indicator[boundary_id] =
              std::make_shared<Functions::ParsedFunction<dim>>();
            phase_indicator[boundary_id]->parse_parameters(prm);
            prm.leave_subsection();
          }
        if (auto const option = prm.get("type"); option == "periodic")
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
      }
  }


  /**
   * @brief Parse the boundary conditions
   * Calls parse_boundary method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

  template <int dim>
  void
  CLSBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions CLS");
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

} // namespace BoundaryConditions

#endif
