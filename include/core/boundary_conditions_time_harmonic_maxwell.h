// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_boundary_conditions_time_harmonic_maxwell_h
#define lethe_boundary_conditions_time_harmonic_maxwell_h

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

DeclException1(TimeHarmonicMaxwellDuplicated,
               types::boundary_id,
               << "Time Harmonic Maxwell boundary id: " << arg1
               << " has already been declared as a boundary condition");


namespace BoundaryConditions
{
  /// Boundary condition types of the time-harmonic Maxwell equations
  inline const std::vector<BoundaryType> time_harmonic_maxwell_boundary_types =
    {BoundaryType::pec,
     BoundaryType::pmc,
     BoundaryType::silver_muller,
     BoundaryType::electric_field,
     BoundaryType::magnetic_field,
     BoundaryType::impedance_boundary,
     BoundaryType::waveguide_port};

  /**
   * @brief This class manages the functions associated with both imposed electric and magnetic fields boundary conditions of the Time Harmonic Maxwell equations
   *
   */
  template <int dim>
  class TimeHarmonicMaxwellBoundaryFunctions
  {
  public:
    // Electric field components
    Functions::ParsedFunction<dim> e_x_real;
    Functions::ParsedFunction<dim> e_y_real;
    Functions::ParsedFunction<dim> e_z_real;
    Functions::ParsedFunction<dim> e_x_imag;
    Functions::ParsedFunction<dim> e_y_imag;
    Functions::ParsedFunction<dim> e_z_imag;

    /// Magnetic field components
    Functions::ParsedFunction<dim> h_x_real;
    Functions::ParsedFunction<dim> h_y_real;
    Functions::ParsedFunction<dim> h_z_real;
    Functions::ParsedFunction<dim> h_x_imag;
    Functions::ParsedFunction<dim> h_y_imag;
    Functions::ParsedFunction<dim> h_z_imag;
  };


  /**
   * @brief This class manages the boundary conditions for Time Harmonic Maxwell solver.
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */

  template <int dim>
  class TimeHarmonicMaxwellBoundaryConditions : public BoundaryConditions
  {
  public:
    /// Functions for (e_x_real, e_y_real, e_z_real, e_x_imag, e_y_imag,
    /// e_z_imag, h_x_real, h_y_real, h_z_real, h_x_imag, h_y_imag, h_z_imag)
    std::map<types::boundary_id,
             std::shared_ptr<TimeHarmonicMaxwellBoundaryFunctions<dim>>>
      imposed_electromagnetic_fields;

    /// The following functions are associated with both imperfect conductor
    /// and imposed electromagnetic excitation or absorption boundary conditions
    /// of the Time Harmonic Maxwell equations, but are parsed individually
    /// because they will not be applied to the vector solution and will
    /// therefore be handled manually in the solver.

    // Impedance
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      surface_admittance_real;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      surface_admittance_imag;

    /// Excitation
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      excitation_x_real;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      excitation_y_real;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      excitation_z_real;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      excitation_x_imag;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      excitation_y_imag;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      excitation_z_imag;


    void
    parse_boundary(ParameterHandler &prm);
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);


    /**
     * @brief Declares the Time Harmonic Maxwell boundary conditions
     *
     * @param prm The parameter file
     * @param number_of_boundary_conditions The number of boundary conditions to be declared. This parameter is generally pre-parsed from a first read of the prm file.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the Time Harmonic Maxwell boundary conditions
     *
     * @param prm The parameter file
     */
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Declares the default parameters for a boundary condition id i_bc
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   * @param default_boundary_id Default value given to the boundary id.
   */
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

  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

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


  /**
   * @brief Parse the boundary conditions
   * Calls parse_boundary method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

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

} // namespace BoundaryConditions


/**
 * @brief This class implements a boundary condition type for the Time Harmonic Maxwell equation
 * where the electric field components are defined using individual functions
 */
template <int dim>
class TimeHarmonicMaxwellElectricFieldDefined : public Function<dim>
{
  Functions::ParsedFunction<dim> *e_x_real;
  Functions::ParsedFunction<dim> *e_y_real;
  Functions::ParsedFunction<dim> *e_z_real;
  Functions::ParsedFunction<dim> *e_x_imag;
  Functions::ParsedFunction<dim> *e_y_imag;
  Functions::ParsedFunction<dim> *e_z_imag;

public:
  TimeHarmonicMaxwellElectricFieldDefined(
    Functions::ParsedFunction<dim> *p_e_x_real,
    Functions::ParsedFunction<dim> *p_e_y_real,
    Functions::ParsedFunction<dim> *p_e_z_real,
    Functions::ParsedFunction<dim> *p_e_x_imag,
    Functions::ParsedFunction<dim> *p_e_y_imag,
    Functions::ParsedFunction<dim> *p_e_z_imag)
    : Function<dim>(4 * dim)
    , e_x_real(p_e_x_real)
    , e_y_real(p_e_y_real)
    , e_z_real(p_e_z_real)
    , e_x_imag(p_e_x_imag)
    , e_y_imag(p_e_y_imag)
    , e_z_imag(p_e_z_imag)
  {}

  double
  value(const Point<dim> &point, const unsigned int component) const override;
};


/**
 * @brief Calculates the value of an imposed electric field equations.
 *
 * @param point A point at which the function will be evaluated.
 * @param component The vector component of the boundary condition (0-> E x real, 1-> E y real, 2-> E z real, 3-> E x imag, 4-> E y imag, 5-> E z imag)
 */
template <int dim>
double
TimeHarmonicMaxwellElectricFieldDefined<dim>::value(
  const Point<dim>  &point,
  const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  if (component == 0)
    {
      return e_x_real->value(point);
    }
  if (component == 1)
    {
      return e_y_real->value(point);
    }
  if (component == 2)
    {
      return e_z_real->value(point);
    }
  if (component == 3)
    {
      return e_x_imag->value(point);
    }
  if (component == 4)
    {
      return e_y_imag->value(point);
    }
  if (component == 5)
    {
      return e_z_imag->value(point);
    }
  return 0.;
}

/**
 * @brief This class implements a boundary condition type for the Time Harmonic Maxwell equation
 * where the magnetic field components are defined using individual functions
 */
template <int dim>
class TimeHarmonicMaxwellMagneticFieldDefined : public Function<dim>
{
  Functions::ParsedFunction<dim> *h_x_real;
  Functions::ParsedFunction<dim> *h_y_real;
  Functions::ParsedFunction<dim> *h_z_real;
  Functions::ParsedFunction<dim> *h_x_imag;
  Functions::ParsedFunction<dim> *h_y_imag;
  Functions::ParsedFunction<dim> *h_z_imag;

public:
  TimeHarmonicMaxwellMagneticFieldDefined(
    Functions::ParsedFunction<dim> *p_h_x_real,
    Functions::ParsedFunction<dim> *p_h_y_real,
    Functions::ParsedFunction<dim> *p_h_z_real,
    Functions::ParsedFunction<dim> *p_h_x_imag,
    Functions::ParsedFunction<dim> *p_h_y_imag,
    Functions::ParsedFunction<dim> *p_h_z_imag)
    : Function<dim>(4 * dim)
    , h_x_real(p_h_x_real)
    , h_y_real(p_h_y_real)
    , h_z_real(p_h_z_real)
    , h_x_imag(p_h_x_imag)
    , h_y_imag(p_h_y_imag)
    , h_z_imag(p_h_z_imag)
  {}

  double
  value(const Point<dim> &point, const unsigned int component) const override;
};


/**
 * @brief Calculates the value of an imposed magnetic field equations.
 *
 * @param point A point at which the function will be evaluated.
 * @param component The vector component of the boundary condition (6-> H x real, 7-> H y real, 8-> H z real, 9-> H x imag, 10-> H y imag, 11-> H z imag)
 */
template <int dim>
double
TimeHarmonicMaxwellMagneticFieldDefined<dim>::value(
  const Point<dim>  &point,
  const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  if (component == 6)
    {
      return h_x_real->value(point);
    }
  if (component == 7)
    {
      return h_y_real->value(point);
    }
  if (component == 8)
    {
      return h_z_real->value(point);
    }
  if (component == 9)
    {
      return h_x_imag->value(point);
    }
  if (component == 10)
    {
      return h_y_imag->value(point);
    }
  if (component == 11)
    {
      return h_z_imag->value(point);
    }
  return 0.;
}

#endif
