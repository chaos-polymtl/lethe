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
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

using namespace dealii;

DeclException1(
  EmissivityError,
  double,
  << "Emissivity : " << arg1
  << " cannot be larger than 1.0 (black body) or smaller than 0.0");

namespace BoundaryConditions
{
  enum class BoundaryType
  {
    // common
    none,
    // for fluid
    noslip,
    slip,
    function,
    function_weak,
    partial_slip,
    periodic,
    pressure,
    outlet,
    // for heat transfer
    noflux,
    temperature,
    convection_radiation,
    // for tracer
    tracer_dirichlet,
    // for vof
    pw,
    vof_dirichlet,
    // for cahn hilliard
    ch_noflux,
    ch_dirichlet_phase_order,
    ch_angle_of_contact,
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

    // Penalization parameter for weak dirichlet BCs and outlets
    std::vector<double> beta;

    // Boundary layer size tangent component parameter for partial slip
    // dirichlet BCs
    std::vector<double> boundary_layer_thickness;

    // Number of boundary conditions
    unsigned int size;
    unsigned int max_size;
    bool         time_dependent;

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
   * @brief This class manages the functions associated with function-type boundary conditions
   * of the Cahn-Hilliard equations
   *
   */
  template <int dim>
  class CahnHilliardBoundaryFunctions
  {
  public:
    // Cahn-Hilliard variables
    Functions::ParsedFunction<dim> phi;
    Functions::ParsedFunction<dim> eta;
  };

  /**
   * @brief This class manages the functions associated with function-type  pressure boundary conditions
   * for the Navier-Stokes equations
   *
   */
  template <int dim>
  class NSPressureBoundaryFunctions
  {
  public:
    // Velocity components
    Functions::ParsedFunction<dim> p;

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
    NSBoundaryFunctions<dim> *        bcFunctions;
    NSPressureBoundaryFunctions<dim> *bcPressureFunction;

    void
    parse_boundary(ParameterHandler &prm, unsigned int i_bc);
    void
    declareDefaultEntry(ParameterHandler &prm, unsigned int i_bc);
    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
    void
    createNoSlip();
  };


  /**
   * @brief Creates a noslip boundary condition for id=0
   * Used in tests only
   */
  template <int dim>
  void
  NSBoundaryConditions<dim>::createNoSlip()
  {
    this->id.resize(1);
    this->id[0] = 0;
    this->type.resize(1);
    this->type[0] = BoundaryType::noslip;
    this->size    = 1;
    this->beta.resize(1);
    this->beta[0] = 0;
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
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection(
        "none|noslip|slip|function|periodic|pressure|function weak|partial slip|outlet"),
      "Type of boundary condition"
      "Choices are <noslip|slip|function|periodic|pressure|function weak|partial slip|outlet>.");


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

    prm.enter_subsection("p");
    bcPressureFunction[i_bc].p.declare_parameters(prm, 1);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    // Center of rotation of the boundary condition for torque calculation
    prm.enter_subsection("center of rotation");
    prm.declare_entry("x", "0", Patterns::Double(), "X COR");
    prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
    prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
    prm.leave_subsection();

    // Penalization parameter for weakly imposed dirichlet BCs and outlets
    prm.declare_entry(
      "beta",
      "0",
      Patterns::Double(),
      "penalty parameter for weak boundary condition imposed through Nitsche's method or outlets");

    // Penalization parameter for weakly imposed dirichlet BCs and outlets
    prm.declare_entry(
      "boundary layer thickness",
      "0",
      Patterns::Double(),
      "thickness of the boundary layer used to calculate the penalty parameter for partial slip boundary condition in tangent direction imposed through Nitsche's method or outlets");
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
    if (op == "none")
      this->type[i_bc] = BoundaryType::none;
    if (op == "noslip")
      this->type[i_bc] = BoundaryType::noslip;
    if (op == "slip")
      this->type[i_bc] = BoundaryType::slip;
    if (op == "function" || op == "function weak" || op == "partial slip")
      {
        if (op == "function")
          this->type[i_bc] = BoundaryType::function;
        else if (op == "partial slip")
          this->type[i_bc] = BoundaryType::partial_slip;
        else
          this->type[i_bc] = BoundaryType::function_weak;
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
    if (op == "pressure")
      {
        this->type[i_bc] = BoundaryType::pressure;
        prm.enter_subsection("p");
        bcPressureFunction[i_bc].p.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("center of rotation");
        bcPressureFunction[i_bc].center_of_rotation[0] = prm.get_double("x");
        bcPressureFunction[i_bc].center_of_rotation[1] = prm.get_double("y");
        if (dim == 3)
          bcPressureFunction[i_bc].center_of_rotation[2] = prm.get_double("z");
        prm.leave_subsection();
      }
    if (op == "periodic")
      {
        this->type[i_bc]               = BoundaryType::periodic;
        this->periodic_id[i_bc]        = prm.get_integer("periodic_id");
        this->periodic_direction[i_bc] = prm.get_integer("periodic_direction");
      }
    if (op == "outlet")
      {
        this->type[i_bc] = BoundaryType::outlet;
      }

    this->id[i_bc]   = prm.get_integer("id");
    this->beta[i_bc] = prm.get_double("beta");
    this->boundary_layer_thickness[i_bc] =
      prm.get_double("boundary layer thickness");
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
      prm.declare_entry(
        "time dependent",
        "false",
        Patterns::Bool(),
        "Bool to define if the boundary condition is time dependent");

      this->id.resize(this->max_size);
      this->beta.resize(this->max_size);
      this->boundary_layer_thickness.resize(this->max_size);
      this->periodic_id.resize(this->max_size);
      this->periodic_direction.resize(this->max_size);
      this->type.resize(this->max_size);
      bcFunctions        = new NSBoundaryFunctions<dim>[this->max_size];
      bcPressureFunction = new NSPressureBoundaryFunctions<dim>[this->max_size];
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
      this->size           = prm.get_integer("number");
      this->time_dependent = prm.get_bool("time dependent");
      this->type.resize(this->size);
      this->id.resize(this->size);
      this->periodic_direction.resize(this->size);
      this->periodic_id.resize(this->size);

      for (unsigned int n = 0; n < this->size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            parse_boundary(prm, n);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  /**
   * @brief This class manages the boundary conditions for Heat-Transfer solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   * The members "value", "h" and "Tinf" contain double used for bc calculation:
   *
   *  - if bc type is "noflux", no flux is assembled at the boundary
   *
   *  - if bc type is "temperature" (Dirichlet condition), "value" is the
   * double passed to the deal.ii ConstantFunction
   *
   *  - if bc type is "convection-radiation" (Robin condition), "h" is the
   * convective heat transfer coefficient and "Tinf" is the
   * environment temperature at the boundary, "emissivity" is the emissivity
   * coefficient, and "Stefan-Boltzmann constant" is the Stefan-Boltzmann
   * constant = 5.6703*10-8 (W.m-2.K-4)
   */

  template <int dim>
  class HTBoundaryConditions : public BoundaryConditions<dim>
  {
  public:
    std::vector<double> value;
    std::vector<double> h;
    std::vector<double> Tinf;
    std::vector<double> emissivity;
    double              Stefan_Boltzmann_constant;

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
                      "noflux",
                      Patterns::Selection(
                        "noflux|temperature|convection-radiation"),
                      "Type of boundary condition for heat transfer"
                      "Choices are <noflux|temperature|convection-radiation>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry("value",
                      "0",
                      Patterns::Double(),
                      "Value (Double) for constant temperature at bc");

    prm.declare_entry(
      "h",
      "0",
      Patterns::Double(),
      "Value (Double) for the h coefficient of convection-radiation bc");

    prm.declare_entry(
      "Tinf",
      "0",
      Patterns::Double(),
      "Temperature (Double) of environment for convection-radiation bc");

    prm.declare_entry("emissivity",
                      "0.0",
                      Patterns::Double(),
                      "Emissivity of the boundary for convection-radiation bc");
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
      prm.declare_entry("Stefan-Boltzmann constant",
                        "0.000000056703",
                        Patterns::Double(),
                        "Stefan-Boltzmann constant");
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
    if (op == "noflux")
      {
        this->type[i_bc] = BoundaryType::noflux;
      }
    if (op == "temperature")
      {
        this->type[i_bc]  = BoundaryType::temperature;
        this->value[i_bc] = prm.get_double("value");
      }
    else if (op == "convection-radiation")
      {
        this->type[i_bc]       = BoundaryType::convection_radiation;
        this->h[i_bc]          = prm.get_double("h");
        this->Tinf[i_bc]       = prm.get_double("Tinf");
        this->emissivity[i_bc] = prm.get_double("emissivity");

        Assert(this->emissivity[i_bc] <= 1.0 && this->emissivity[i_bc] >= 0.0,
               EmissivityError(this->emissivity[i_bc]));
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
      this->emissivity.resize(this->size);

      for (unsigned int n = 0; n < this->size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            parse_boundary(prm, n);
          }
          prm.leave_subsection();
        }
      this->Stefan_Boltzmann_constant =
        prm.get_double("Stefan-Boltzmann constant");
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
    std::vector<std::shared_ptr<Functions::ParsedFunction<dim>>> tracer;


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
                      "dirichlet",
                      Patterns::Selection("dirichlet"),
                      "Type of boundary condition for tracer"
                      "Choices are <function>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.enter_subsection("dirichlet");
    tracer[i_bc] = std::make_shared<Functions::ParsedFunction<dim>>();
    tracer[i_bc]->declare_parameters(prm);
    prm.set("Function expression", "0");
    prm.leave_subsection();
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
      tracer.resize(this->max_size);

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
    if (op == "dirichlet")
      {
        this->type[i_bc] = BoundaryType::tracer_dirichlet;
        prm.enter_subsection("dirichlet");
        tracer[i_bc]->parse_parameters(prm);
        prm.leave_subsection();
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
    prm.enter_subsection("boundary conditions tracer");
    {
      this->size = prm.get_integer("number");

      this->type.resize(this->size);

      for (unsigned int n = 0; n < this->size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            parse_boundary(prm, n);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  /**
   * @brief This class manages the boundary conditions for the Cahn-Hilliard solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   *  - if bc type is "dirichlet" (Dirichlet condition), "value" is the
   * double passed to the deal.ii ConstantFunction
   */

  template <int dim>
  class CahnHilliardBoundaryConditions : public BoundaryConditions<dim>
  {
  public:
    std::vector<double>                 angle_of_contact;
    std::vector<double>                 phase_dirichlet_value;
    CahnHilliardBoundaryFunctions<dim> *bcFunctions;

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
  CahnHilliardBoundaryConditions<dim>::declareDefaultEntry(
    ParameterHandler &prm,
    unsigned int      i_bc)
  {
    prm.declare_entry(
      "type",
      "noflux",
      Patterns::Selection("noflux|dirichlet|angle_of_contact"),
      "Type of boundary condition for the Cahn-Hilliard equations"
      "Choices are <noflux|dirichlet|angle_of_contact>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.declare_entry(
      "angle value",
      "0",
      Patterns::Double(),
      "Inner angle of contact between the fluid 1 and the boundary (in degrees)");

    prm.enter_subsection("phi");
    bcFunctions[i_bc].phi.declare_parameters(prm, 1);
    prm.set("Function expression", "0");
    prm.leave_subsection();

    return;
  }

  /**
   * @brief Declare the boundary conditions default parameters
   * Calls declareDefaultEntry method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  CahnHilliardBoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->max_size = 14;

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

      this->id.resize(this->max_size);
      this->type.resize(this->max_size);
      bcFunctions = new CahnHilliardBoundaryFunctions<dim>[this->max_size];

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
   * @param i_bc The boundary condition number (and not necessarily its id).
   */

  template <int dim>
  void
  CahnHilliardBoundaryConditions<dim>::parse_boundary(ParameterHandler & prm,
                                                      const unsigned int i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "noflux")
      {
        this->type[i_bc] = BoundaryType::ch_noflux;
      }
    if (op == "dirichlet")
      {
        this->type[i_bc] = BoundaryType::ch_dirichlet_phase_order;
        prm.enter_subsection("phi");
        bcFunctions[i_bc].phi.parse_parameters(prm);
        prm.leave_subsection();
      }
    if (op == "angle_of_contact")
      {
        this->type[i_bc]             = BoundaryType::ch_angle_of_contact;
        this->angle_of_contact[i_bc] = prm.get_double("angle value");
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
  CahnHilliardBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions cahn hilliard");
    {
      this->size           = prm.get_integer("number");
      this->time_dependent = prm.get_bool("time dependent");
      this->type.resize(this->size);
      this->id.resize(this->size);
      this->phase_dirichlet_value.resize(this->size);
      this->angle_of_contact.resize(this->size);

      for (unsigned int n = 0; n < this->size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            parse_boundary(prm, n);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  /**
   * @brief This class manages the boundary conditions for VOF solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   *
   *  - if bc type is "peeling/wetting", peeling/wetting of the free surface
   * will be applied. See vof.cc for further implementation details.
   *
   * - if bc type is "dirichlet", the function is applied on the selected
   * boundary
   *
   * - if bc type is "none", nothing happens
   */

  template <int dim>
  class VOFBoundaryConditions : public BoundaryConditions<dim>
  {
  public:
    std::vector<std::shared_ptr<Functions::ParsedFunction<dim>>> phase_fraction;

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
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param i_bc The boundary condition id.
   */
  template <int dim>
  void
  VOFBoundaryConditions<dim>::declareDefaultEntry(ParameterHandler &prm,
                                                  unsigned int      i_bc)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("none|dirichlet|peeling/wetting"),
                      "Type of boundary condition for VOF"
                      "Choices are <none|dirichlet|peeling/wetting>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(i_bc, 2),
                      Patterns::Integer(),
                      "Mesh id for boundary conditions");

    prm.enter_subsection("dirichlet");
    phase_fraction[i_bc] = std::make_shared<Functions::ParsedFunction<dim>>();
    phase_fraction[i_bc]->declare_parameters(prm);
    prm.set("Function expression", "0");
    prm.leave_subsection();
  }

  /**
   * @brief Declare the boundary conditions default parameters
   * Calls declareDefaultEntry method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  VOFBoundaryConditions<dim>::declare_parameters(ParameterHandler &prm)
  {
    this->max_size = 14;

    prm.enter_subsection("boundary conditions VOF");
    {
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of boundary conditions");
      this->id.resize(this->max_size);
      this->type.resize(this->max_size);
      phase_fraction.resize(this->max_size);

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
  VOFBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm,
                                             unsigned int      i_bc)
  {
    const std::string op = prm.get("type");
    if (op == "none")
      {
        this->type[i_bc] = BoundaryType::none;
      }
    else if (op == "dirichlet")
      {
        this->type[i_bc] = BoundaryType::vof_dirichlet;
        prm.enter_subsection("dirichlet");
        phase_fraction[i_bc]->parse_parameters(prm);
        prm.leave_subsection();
      }
    else if (op == "peeling/wetting")
      {
        this->type[i_bc] = BoundaryType::pw;
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
  VOFBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions VOF");
    {
      this->size = prm.get_integer("number");

      this->type.resize(this->size);
      this->id.resize(this->size);

      for (unsigned int n = 0; n < this->size; n++)
        {
          prm.enter_subsection("bc " + std::to_string(n));
          {
            parse_boundary(prm, n);
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
         ExcIndexRange(component, 0, this->n_components)) if (component == 0)
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

/**
 * @brief This class implements a pressure boundary condition for the Navier-Stokes equations.
 */
template <int dim>
class NavierStokesPressureFunctionDefined : public Function<dim>
{
private:
  Functions::ParsedFunction<dim> *p;

public:
  NavierStokesPressureFunctionDefined(Functions::ParsedFunction<dim> *p_p)
    : Function<dim>(dim + 1)
    , p(p_p)
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
NavierStokesPressureFunctionDefined<dim>::value(
  const Point<dim> & point,
  const unsigned int component) const
{
  if (component == dim)
    {
      return p->value(point);
    }
  return 0.;
}

/**
 * @brief This class implements a boundary conditions for the Cahn-Hilliard equation
 * where the phase and chemical potential are defined using individual functions
 */
template <int dim>
class CahnHilliardFunctionDefined : public Function<dim>
{
private:
  Functions::ParsedFunction<dim> *phi;

public:
  CahnHilliardFunctionDefined(Functions::ParsedFunction<dim> *p_phi)
    : Function<dim>(2)
    , phi(p_phi)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const override;
};


/**
 * @brief Calculates the value of a Function-type for the Cahn-Hilliard equations
 *
 * @param p A point (generally a gauss point)
 *
 * @param component The vector component of the boundary condition (0-x, 1-y and 2-z)
 */
template <int dim>
double
CahnHilliardFunctionDefined<dim>::value(const Point<dim> & p,
                                        const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components)) if (component == 0)
  {
    return phi->value(p);
  }
  else if (component == 1)
  {
    return 0.;
  }
  return 0.;
}

#endif
