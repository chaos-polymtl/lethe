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

#ifndef lethe_boundary_conditions_h
#define lethe_boundary_conditions_h

#include <deal.II/base/function.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace BoundaryConditions
{
  enum class BoundaryType
  {
    // for fluid
    noslip,
    slip,
    function,
    periodic,
    // for heat transfer
    temperature, // Dirichlet
    convection,  // Robin

    // for tracer
    tracer_dirichlet, // Dirichlet tracer
  };

  /**
   * @brief This class is the base class for all boundary conditions. It stores
   * the general information that all boundary condition share.
   * In Lethe, boundary conditions are identified with an id, a type and, in the
   * special case of periodic boundary condition, a periodic matching id
   * (periodic_id) and a periodic direction (0, 1 or 2).
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
   * @brief This class manages the functions associated with function-type boundary conditions
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
    Point<dim> center_of_rotation;
  };


  /**
   * @brief This class manages the boundary conditions for Navier-Stokes solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
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


  /**
   * @brief Creates a default no-slip boundary condition for id=0
   */
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


  /**
   * @brief Declares the default parameter for a boundary condition id i_bc
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param i_bc The boundary condition id.
   */
  template <int dim>
  void
  NSBoundaryConditions<dim>::declareDefaultEntry(ParameterHandler &prm,
                                                 unsigned int      i_bc)
  {
    prm.declare_entry("type",
                      "noslip",
                      Patterns::Selection("noslip|slip|function|periodic"),
                      "Type of boundary condition"
                      "Choices are <noslip|slip|function|periodic>.");

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

    // Center of rotation of the boundary condition for torque calculation
    prm.enter_subsection("center of rotation");
    prm.declare_entry("x", "0", Patterns::Double(), "X COR");
    prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
    prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
    prm.leave_subsection();
  }


  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param i_bc The boundary condition number (and not necessarily it's id).
   */
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

        prm.enter_subsection("center of rotation");
        bcFunctions[i_bc].center_of_rotation[0] = prm.get_double("x");
        bcFunctions[i_bc].center_of_rotation[1] = prm.get_double("y");
        if (dim == 3)
          bcFunctions[i_bc].center_of_rotation[2] = prm.get_double("z");
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


  /**
   * @brief Declare the boundary conditions default parameters
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  NSBoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->max_size = 14;

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

      for (unsigned int n = 0; n < this->max_size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            declareDefaultEntry(prm, n);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }


  /**
   * @brief Parse the boundary conditions
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
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

      for (unsigned int n = 0; n < this->max_size; n++)
        {
          if (this->size >= n + 1)
            {
              prm.enter_subsection("bc " + std::to_string(n));
              {
                parse_boundary(prm, n);
              }
              prm.leave_subsection();
            }
        }
    }
    prm.leave_subsection();
  }

  /**
   * @brief This class manages the boundary conditions for Heat-Transfer solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   * The members "value" and "Tenv" contain double used for bc calculation :
   *  - if bc type is "temperature" (Dirichlet condition), "value" is the
   * double passed to the deal.ii ConstantFunction
   *  - if bc type is "convection" (Robin condition), "value" is the
   * convective heat transfer coefficient (h) and "Tenv" is the
   * environment temperature at the boundary
   */

  template <int dim>
  class HTBoundaryConditions : public BoundaryConditions<dim>
  {
  public:
    std::vector<double> value;
    std::vector<double> h;
    std::vector<double> Tinf;

    void
    declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc);
    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_boundary(ParameterHandler &prm, unsigned int i_bc);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Declares the default parameters for a boundary condition id i_bc
   * i.e. Dirichlet condition (ConstantFunction) with value 0
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param i_bc The boundary condition id.
   */
  template <int dim>
  void
  HTBoundaryConditions<dim>::declareDefaultEntry(ParameterHandler &prm,
                                                 unsigned int      i_bc)
  {
    prm.declare_entry("type",
                      "temperature",
                      Patterns::Selection("temperature|convection"),
                      "Type of boundary condition for heat transfer"
                      "Choices are <temperature|convection>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry("value",
                      "0",
                      Patterns::Double(),
                      "Value (Double) for constant temperature at bc");

    prm.declare_entry("h",
                      "0",
                      Patterns::Double(),
                      "Value (Double) for the h coefficient of convection bc");

    prm.declare_entry("Tinf",
                      "0",
                      Patterns::Double(),
                      "Temperature (Double) of environment for convection bc");
  }

  /**
   * @brief Declare the boundary conditions default parameters
   * Calls declareDefaultEntry method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  HTBoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->max_size = 14;

    prm.enter_subsection("boundary conditions heat transfer");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      this->id.resize(this->max_size);
      this->type.resize(this->max_size);

      for (unsigned int n = 0; n < this->max_size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            declareDefaultEntry(prm, n);
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
   *
   * @param i_bc The boundary condition number (and not necessarily it's id).
   */

  template <int dim>
  void
  HTBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm,
                                            unsigned int      i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "temperature")
      {
        this->type[i_bc]  = BoundaryType::temperature;
        this->value[i_bc] = prm.get_double("value");
      }
    else if (op == "convection")
      {
        this->type[i_bc] = BoundaryType::convection;
        this->h[i_bc]    = prm.get_double("h");
        this->Tinf[i_bc] = prm.get_double("Tinf");
      }

    this->id[i_bc] = prm.get_integer("id");
  }

  /**
   * @brief Parse the boundary conditions
   * Calls parse_boundary method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

  template <int dim>
  void
  HTBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions heat transfer");
    {
      this->size = prm.get_integer("number");

      this->type.resize(this->size);
      this->id.resize(this->size);
      this->value.resize(this->size);
      this->h.resize(this->size);

      this->Tinf.resize(this->size);

      for (unsigned int n = 0; n < this->max_size; n++)
        {
          if (this->size >= n + 1)
            {
              prm.enter_subsection("bc " + std::to_string(n));
              {
                parse_boundary(prm, n);
              }
              prm.leave_subsection();
            }
        }
    }
    prm.leave_subsection();
  }


  /**
   * @brief This class manages the boundary conditions for Tracer solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   *  - if bc type is "dirichlet" (Dirichlet condition), "value" is the
   * double passed to the deal.ii ConstantFunction

   */

  template <int dim>
  class TracerBoundaryConditions : public BoundaryConditions<dim>
  {
  public:
    Functions::ParsedFunction<dim> u;


    void
    declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc);
    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_boundary(ParameterHandler &prm, unsigned int i_bc);
    void
    parse_parameters(ParameterHandler &prm);
  };

  /**
   * @brief Declares the default parameters for a boundary condition id i_bc
   * i.e. Dirichlet condition (ConstantFunction) with value 0
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param i_bc The boundary condition id.
   */
  template <int dim>
  void
  TracerBoundaryConditions<dim>::declareDefaultEntry(ParameterHandler &prm,
                                                     unsigned int      i_bc)
  {
    prm.declare_entry("type",
                      "temperature",
                      Patterns::Selection("dirichlet"),
                      "Type of boundary condition for tracer"
                      "Choices are <dirichlet>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");
  }

  /**
   * @brief Declare the boundary conditions default parameters
   * Calls declareDefaultEntry method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  TracerBoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->max_size = 14;

    prm.enter_subsection("boundary conditions tracer");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      this->id.resize(this->max_size);
      this->type.resize(this->max_size);

      for (unsigned int n = 0; n < this->max_size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            declareDefaultEntry(prm, n);
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
   *
   * @param i_bc The boundary condition number (and not necessarily it's id).
   */

  template <int dim>
  void
  TracerBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm,
                                                unsigned int      i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "temperature")
      {
        this->type[i_bc]  = BoundaryType::tracer_dirichlet;
        this->value[i_bc] = prm.get_double("value");
      }

    this->id[i_bc] = prm.get_integer("id");
  }

  /**
   * @brief Parse the boundary conditions
   * Calls parse_boundary method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

  template <int dim>
  void
  TracerBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions heat transfer");
    {
      this->size = prm.get_integer("number");

      this->type.resize(this->size);

      for (unsigned int n = 0; n < this->max_size; n++)
        {
          if (this->size >= n + 1)
            {
              prm.enter_subsection("bc " + std::to_string(n));
              {
                parse_boundary(prm, n);
              }
              prm.leave_subsection();
            }
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


/**
 * @brief Calculates the value of a Function-type Navier-Stokes equations
 *
 * @param p A point (generally a gauss point)
 *
 * @param component The vector component of the boundary condition (0-x, 1-y and 2-z)
 */
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
