/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef LETHE_BOUNDARYCONDITIONS_H
#define LETHE_BOUNDARYCONDITIONS_H

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace BoundaryConditions
{
  enum class BoundaryType
  {
    noslip,
    slip,
    function,
    periodic
  };

  /**
   * @brief This class is the base class for all boundary conditions. It stores
   * the general information that all boundary condition share.
   * In Lethe, boundary conditions are identified with an id, a type and, in the
   * special case of periodic boundary conditoin, a periodic matching id
   * (periodic_id) and a periodic direction
   */
  template <int dim>
  class BoundaryConditions
  {
  public:
    // ID of boundary condition
    std::vector<unsigned int> id;

    // List of boundary type for each number
    std::vector<BoundaryType> type;

    // Number of boundary conditions
    unsigned int size;
    unsigned int max_size;

    // Periodic boundary condition matching
    std::vector<unsigned int> periodic_id;
    std::vector<unsigned int> periodic_direction;
  };

  /**
   * @brief This class managed the functions associated with function boundary conditions
   * of the Navier-Stokes equations
   *
   */
  template <int dim>
  class NSBoundaryFunctions
  {
  public:
    // Velocity components
    Functions::ParsedFunction<dim> u;
    Functions::ParsedFunction<dim> v;
    Functions::ParsedFunction<dim> w;

    // Point for the center of rotation
    Point<dim> cor;
  };


  /**
   * @brief This class manages the boundary conditions for Navier-Strokes solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently
   *
   */
  template <int dim>
  class NSBoundaryConditions : public BoundaryConditions<dim>
  {
  public:
    // Functions for (u,v,w) for all boundaries
    NSBoundaryFunctions<dim> *bcFunctions;

    void
    parse_boundary(ParameterHandler &prm, unsigned int i_bc);
    void
    declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc);
    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
    void
    createDefaultNoSlip();
  };

  template <int dim>
  void
  NSBoundaryConditions<dim>::createDefaultNoSlip()
  {
    this->id.resize(1);
    this->id[0] = 0;
    this->type.resize(1);
    this->type[0] = BoundaryType::noslip;
    this->size    = 1;
  }

  template <int dim>
  void
  NSBoundaryConditions<dim>::declareDefaultEntry(ParameterHandler &prm,
                                                 unsigned int      i_bc)
  {
    prm.declare_entry("type",
                      "noslip",
                      Patterns::Selection("noslip|slip|function|periodic"),
                      "Type of boundary conditoin"
                      "Choices are <noslip|slip|function>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry("periodic_id",
                      "0",
                      Patterns::Integer(),
                      "Mesh id for periodic face matching");

    prm.declare_entry("periodic_direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    prm.enter_subsection("u");
    bcFunctions[i_bc].u.declare_parameters(prm, 1);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    prm.enter_subsection("v");
    bcFunctions[i_bc].v.declare_parameters(prm, 1);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    prm.enter_subsection("w");
    bcFunctions[i_bc].w.declare_parameters(prm, 1);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    prm.enter_subsection("cor");
    prm.declare_entry("x", "0", Patterns::Double(), "X COR");
    prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
    prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
    prm.leave_subsection();
  }

  template <int dim>
  void
  NSBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm,
                                            unsigned int      i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "noslip")
      this->type[i_bc] = BoundaryType::noslip;
    if (op == "slip")
      this->type[i_bc] = BoundaryType::slip;
    if (op == "function")
      {
        this->type[i_bc] = BoundaryType::function;
        prm.enter_subsection("u");
        bcFunctions[i_bc].u.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("v");
        bcFunctions[i_bc].v.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("w");
        bcFunctions[i_bc].w.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("cor");
        bcFunctions[i_bc].cor[0] = prm.get_double("x");
        bcFunctions[i_bc].cor[1] = prm.get_double("y");
        if (dim == 3)
          bcFunctions[i_bc].cor[2] = prm.get_double("z");
        prm.leave_subsection();
      }
    if (op == "periodic")
      {
        this->type[i_bc]               = BoundaryType::periodic;
        this->periodic_id[i_bc]        = prm.get_integer("periodic_id");
        this->periodic_direction[i_bc] = prm.get_integer("periodic_direction");
      }

    this->id[i_bc] = prm.get_integer("id");
  }

  template <int dim>
  void
  NSBoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->max_size = 7;

    prm.enter_subsection("boundary conditions");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      this->id.resize(this->max_size);
      this->periodic_id.resize(this->max_size);
      this->periodic_direction.resize(this->max_size);
      this->type.resize(this->max_size);
      bcFunctions = new NSBoundaryFunctions<dim>[this->max_size];

      prm.enter_subsection("bc 0");
      {
        declareDefaultEntry(prm, 0);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 1");
      {
        declareDefaultEntry(prm, 1);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 2");
      {
        declareDefaultEntry(prm, 2);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 3");
      {
        declareDefaultEntry(prm, 3);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 4");
      {
        declareDefaultEntry(prm, 4);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 5");
      {
        declareDefaultEntry(prm, 5);
      }
      prm.leave_subsection();

      prm.enter_subsection("bc 6");
      {
        declareDefaultEntry(prm, 6);
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  NSBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions");
    {
      this->size = prm.get_integer("number");
      this->type.resize(this->size);
      this->id.resize(this->size);
      this->periodic_direction.resize(this->size);
      this->periodic_id.resize(this->size);

      if (this->size >= 1)
        {
          prm.enter_subsection("bc 0");
          {
            parse_boundary(prm, 0);
          }
          prm.leave_subsection();
        }
      if (this->size >= 2)
        {
          prm.enter_subsection("bc 1");
          {
            parse_boundary(prm, 1);
          }
          prm.leave_subsection();
        }
      if (this->size >= 3)
        {
          prm.enter_subsection("bc 2");
          {
            parse_boundary(prm, 2);
          }
          prm.leave_subsection();
        }
      if (this->size >= 4)
        {
          prm.enter_subsection("bc 3");
          {
            parse_boundary(prm, 3);
          }
          prm.leave_subsection();
        }

      if (this->size >= 5)
        {
          prm.enter_subsection("bc 4");
          {
            parse_boundary(prm, 4);
          }
          prm.leave_subsection();
        }

      if (this->size >= 6)
        {
          prm.enter_subsection("bc 5");
          {
            parse_boundary(prm, 5);
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
private:
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

  virtual double
  value(const Point<dim> &p, const unsigned int component) const override;
};

template <int dim>
double
NavierStokesFunctionDefined<dim>::value(const Point<dim> & p,
                                        const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));
  if (component == 0)
    {
      return u->value(p);
    }
  else if (component == 1)
    {
      return v->value(p);
    }
  else if (component == 2)
    {
      return w->value(p);
    }
  return 0.;
}


#endif
