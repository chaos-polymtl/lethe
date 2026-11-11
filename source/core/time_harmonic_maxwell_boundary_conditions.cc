// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/time_harmonic_maxwell_boundary_conditions.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <memory>
#include <string>
#include <vector>

using namespace dealii;

DeclException1(TimeHarmonicMaxwellDuplicated,
               types::boundary_id,
               << "Time Harmonic Maxwell boundary id: " << arg1
               << " has already been declared as a boundary condition");

namespace BoundaryConditions
{

  template <int dim>
  void
  TimeHarmonicMaxwellBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    declare_enum_entry(
      prm,
      "type",
      BoundaryType::silver_muller,
      "Type of boundary condition for Time Harmonic Maxwell equations.",
      deprecated_boundary_type_names,
      time_harmonic_maxwell_boundary_types);

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    // Create a dummy TimeHarmonicMaxwellBoundaryFunctions object to
    // declare the appropriate parameters for the relevant boundary conditions.
    TimeHarmonicMaxwellBoundaryFunctions<dim> temporary_em_functions;

    prm.enter_subsection("E x real part");
    temporary_em_functions.e_x_real.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("E y real part");
    temporary_em_functions.e_y_real.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("E z real part");
    temporary_em_functions.e_z_real.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("E x imag part");
    temporary_em_functions.e_x_imag.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("E y imag part");
    temporary_em_functions.e_y_imag.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("E z imag part");
    temporary_em_functions.e_z_imag.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("H x real part");
    temporary_em_functions.h_x_real.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("H y real part");
    temporary_em_functions.h_y_real.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("H z real part");
    temporary_em_functions.h_z_real.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("H x imag part");
    temporary_em_functions.h_x_imag.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("H y imag part");
    temporary_em_functions.h_y_imag.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("H z imag part");
    temporary_em_functions.h_z_imag.declare_parameters(prm);
    prm.leave_subsection();

    // Create a dummy object to declare the appropriate parameters for the
    // relevant boundary conditions.
    Functions::ParsedFunction<dim> temporary_function;

    prm.enter_subsection("surface admittance real part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("surface admittance imag part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("excitation x real part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("excitation x imag part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("excitation y real part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("excitation y imag part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("excitation z real part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("excitation z imag part");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();
  }

  template <int dim>
  void
  TimeHarmonicMaxwellBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions time harmonic maxwell");
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
  TimeHarmonicMaxwellBoundaryConditions<dim>::parse_boundary(
    ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the time harmonic maxwell boundary conditions. Please ensure that the id is set for every boundary condition."));

    // Loop through all boundary ids to ensure that they are all non-negative
    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    TimeHarmonicMaxwellDuplicated(boundary_id));

        // Allocate the imposed electromagnetic field functions object for every
        // boundary condition
        imposed_electromagnetic_fields[boundary_id] =
          std::make_shared<TimeHarmonicMaxwellBoundaryFunctions<dim>>();

        prm.enter_subsection("E x real part");
        imposed_electromagnetic_fields[boundary_id]->e_x_real.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("E y real part");
        imposed_electromagnetic_fields[boundary_id]->e_y_real.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("E z real part");
        imposed_electromagnetic_fields[boundary_id]->e_z_real.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("E x imag part");
        imposed_electromagnetic_fields[boundary_id]->e_x_imag.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("E y imag part");
        imposed_electromagnetic_fields[boundary_id]->e_y_imag.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("E z imag part");
        imposed_electromagnetic_fields[boundary_id]->e_z_imag.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("H x real part");
        imposed_electromagnetic_fields[boundary_id]->h_x_real.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("H y real part");
        imposed_electromagnetic_fields[boundary_id]->h_y_real.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("H z real part");
        imposed_electromagnetic_fields[boundary_id]->h_z_real.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("H x imag part");
        imposed_electromagnetic_fields[boundary_id]->h_x_imag.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("H y imag part");
        imposed_electromagnetic_fields[boundary_id]->h_y_imag.parse_parameters(
          prm);
        prm.leave_subsection();

        prm.enter_subsection("H z imag part");
        imposed_electromagnetic_fields[boundary_id]->h_z_imag.parse_parameters(
          prm);
        prm.leave_subsection();

        /// The following functions are parsed individually because they will be
        /// used for Robin boundary conditions and will not be applied on the
        /// global solution vector.
        prm.enter_subsection("surface admittance real part");
        this->surface_admittance_real[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->surface_admittance_real[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("surface admittance imag part");
        this->surface_admittance_imag[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->surface_admittance_imag[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("excitation x real part");
        this->excitation_x_real[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->excitation_x_real[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("excitation y real part");
        this->excitation_y_real[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->excitation_y_real[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("excitation z real part");
        this->excitation_z_real[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->excitation_z_real[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("excitation x imag part");
        this->excitation_x_imag[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->excitation_x_imag[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("excitation y imag part");
        this->excitation_y_imag[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->excitation_y_imag[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("excitation z imag part");
        this->excitation_z_imag[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        this->excitation_z_imag[boundary_id]->parse_parameters(prm);
        prm.leave_subsection();

        // Establish the type of boundary condition
        this->type[boundary_id] =
          string_to_enum<BoundaryType>(prm.get("type"),
                                       deprecated_boundary_type_names,
                                       "type");
      }
  }

  template <int dim>
  void
  TimeHarmonicMaxwellBoundaryConditions<dim>::parse_parameters(
    ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions time harmonic maxwell");
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

  template class TimeHarmonicMaxwellBoundaryConditions<2>;
  template class TimeHarmonicMaxwellBoundaryConditions<3>;
} // namespace BoundaryConditions
