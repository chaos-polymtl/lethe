// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_time_harmonic_maxwell_boundary_conditions_h
#define lethe_time_harmonic_maxwell_boundary_conditions_h

#include <core/boundary_conditions.h>

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <map>
#include <memory>
#include <vector>

using namespace dealii;

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

    /**
     * @brief Parses the information for one time-harmonic Maxwell boundary condition.
     *
     * @param[in] prm The parameter handler containing the boundary condition.
     */
    void
    parse_boundary(ParameterHandler &prm);

    /**
     * @brief Declares the default parameters for one time-harmonic Maxwell boundary condition.
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] default_boundary_id The default mesh boundary id.
     */
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declares the Time Harmonic Maxwell boundary conditions
     *
     * @param[in] prm The parameter handler in which to declare the parameters.
     * @param[in] number_of_boundary_conditions The number of boundary
     * conditions to declare.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the Time Harmonic Maxwell boundary conditions
     *
     * @param[in] prm The parameter handler containing the boundary conditions.
     */
    void
    parse_parameters(ParameterHandler &prm);
  };
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
 * @param[in] point A point at which the function will be evaluated.
 * @param[in] component The vector component of the boundary condition (0-> E x
 * real, 1-> E y real, 2-> E z real, 3-> E x imag, 4-> E y imag, 5-> E z imag)
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
 * @param[in] point A point at which the function will be evaluated.
 * @param[in] component The vector component of the boundary condition (6-> H x
 * real, 7-> H y real, 8-> H z real, 9-> H x imag, 10-> H y imag, 11-> H z imag)
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
