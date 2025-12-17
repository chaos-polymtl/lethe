// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_boundary_conditions_h
#define lethe_boundary_conditions_h

#include <core/utilities.h>

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/patterns.h>

#include <map>

using namespace dealii;

DeclException1(
  EmissivityError,
  double,
  << "Emissivity : " << arg1
  << " cannot be larger than 1.0 (black body) or smaller than 0.0");

DeclException1(NavierStokesBoundaryDuplicated,
               types::boundary_id,
               << "Fluid dynamics boundary id: " << arg1
               << " has already been declared as a boundary condition");

DeclException1(HeatTransferBoundaryDuplicated,
               types::boundary_id,
               << "Heat transfer boundary id: " << arg1
               << " has already been declared as a boundary condition");

DeclException1(TracerBoundaryDuplicated,
               types::boundary_id,
               << "Tracer boundary id: " << arg1
               << " has already been declared as a boundary condition");

DeclException1(CahnHilliardBoundaryDuplicated,
               types::boundary_id,
               << "Cahn Hilliard boundary id: " << arg1
               << " has already been declared as a boundary condition");

DeclException1(VOFBoundaryDuplicated,
               types::boundary_id,
               << "VOF boundary id: " << arg1
               << " has already been declared as a boundary condition");


namespace BoundaryConditions
{

  enum class BoundaryType
  {
    // common
    none,
    outlet, // outlet is used for fluid dynamics, tracer and eventually other
            // physics
    periodic,
    periodic_neighbor, // The periodic neighbour is used to indicate a boundary
                       // which matches with a main periodic boundary

    // for fluid
    noslip,
    slip,
    function,
    function_weak,
    partial_slip,
    pressure,
    //  for heat transfer
    noflux,
    temperature,
    convection_radiation,
    // for tracer
    tracer_dirichlet,
    // for vof
    vof_dirichlet,
    // for cahn hilliard
    cahn_hilliard_noflux,
    cahn_hilliard_dirichlet_phase_order,
    cahn_hilliard_angle_of_contact,
    cahn_hilliard_free_angle
  };

  /**
   * @brief This class is the base class for all boundary conditions. It stores
   * the general information that all boundary condition share.
   * In Lethe, boundary conditions are identified with an id, a type and, in the
   * special case of periodic boundary condition, a periodic matching id
   * (periodic_id) and a periodic direction (0, 1 or 2).
   */
  class BoundaryConditions
  {
  public:
    /// Map containing the boundary id and the boundary type
    std::map<types::boundary_id, BoundaryType> type;

    ///  Map containing the penalization parameter for weak dirichlet BCs and
    ///  outlets for a corresponding id
    std::map<types::boundary_id, double> beta;

    /// Map containing the boundary layer thickness tangent component parameter
    /// for partial slip for a corresponding id Map containing Dirichlet BCs for
    /// a corresponding boundary id
    std::map<types::boundary_id, double> boundary_layer_thickness;

    /// Number of boundary conditions
    unsigned int number_of_boundary_conditions;

    /// indicator for transient BCs
    bool time_dependent;

    /// Map containing the boundary id and its corresponding periodic boundary
    /// condition match
    std::map<types::boundary_id, types::boundary_id> periodic_neighbor_id;

    /// Map containing the boundary id and its periodic direction
    std::map<types::boundary_id, unsigned int> periodic_direction;
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

    /// Pressure
    Functions::ParsedFunction<dim> p;

    /// Point for the center of rotation
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
   * @brief This class manages the boundary conditions for Navier-Stokes solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */
  template <int dim>
  class NSBoundaryConditions : public BoundaryConditions
  {
  public:
    /// Functions for (u, v, w, p) for all boundaries
    std::map<types::boundary_id, std::shared_ptr<NSBoundaryFunctions<dim>>>
      navier_stokes_functions;

    void
    parse_boundary(ParameterHandler &prm);
    void
    declare_default_entry(ParameterHandler  &prm,
                          types::boundary_id default_boundary_id);

    /**
     * @brief Declares the Navier-Stokes boundary conditions
     *
     * @param prm The parameter file
     * @param number_of_boundary_conditions The number of boundary conditions to be declared. This parameter is generally pre-parsed from a first read of the prm file.
     */
    void
    declare_parameters(ParameterHandler &prm,
                       unsigned int      number_of_boundary_conditions);

    /**
     * @brief Parses the Navier-Stokes boundary conditions
     *
     * @param prm The parameter file
     */
    void
    parse_parameters(ParameterHandler &prm);
    void
    createNoSlip();

    /// Fix pressure constant using a node
    bool fix_pressure_constant;
  };


  /**
   * @brief Creates a noslip boundary condition for id=0
   * Used in tests only
   */
  template <int dim>
  void
  NSBoundaryConditions<dim>::createNoSlip()
  {
    this->type[0]                       = BoundaryType::noslip;
    this->number_of_boundary_conditions = 1;
    this->beta[0]                       = 0;
  }


  /**
   * @brief Declares the default parameter for a boundary condition id i_bc
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param default_boundary_id Default value of the boundary id. This corresponds to the number of the boundary condition subsection.
   */
  template <int dim>
  void
  NSBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection(
        "none|noslip|slip|function|periodic|pressure|function weak|partial slip|outlet"),
      "Type of boundary condition"
      "Choices are <noslip|slip|function|periodic|pressure|function weak|partial slip|outlet>.");


    prm.declare_entry("id",
                      Utilities::to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions.");

    prm.declare_entry(
      "periodic id",
      "-1",
      Patterns::Integer(),
      "Mesh id for periodic face matching. Default entry is -1 to ensure that the periodic id is set by the user");

    prm.declare_alias("periodic id", "periodic_id", true);

    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    prm.declare_alias("periodic direction", "periodic_direction", true);


    // Create a dummy NSBoundaryFunctions object to declare the appropriate
    // parameters for this boundary condition.
    NSBoundaryFunctions<dim> temporary_fluid_dynamics_functions;


    prm.enter_subsection("u");
    temporary_fluid_dynamics_functions.u.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("v");
    temporary_fluid_dynamics_functions.v.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("w");
    temporary_fluid_dynamics_functions.w.declare_parameters(prm);
    prm.leave_subsection();

    prm.enter_subsection("p");
    temporary_fluid_dynamics_functions.p.declare_parameters(prm);
    prm.leave_subsection();

    // Center of rotation of the boundary condition for torque calculation
    prm.enter_subsection("center of rotation");
    prm.declare_entry("x", "0", Patterns::Double(), "X COR");
    prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
    prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
    prm.leave_subsection();

    prm.declare_entry(
      "beta",
      "1",
      Patterns::Double(),
      "Penalty parameter for weak boundary condition imposed through Nitsche's method or outlets");

    prm.declare_entry(
      "boundary layer thickness",
      "0",
      Patterns::Double(),
      "Thickness of the boundary layer used to calculate the penalty parameter for partial slip boundary condition in tangent direction imposed through Nitsche's method or outlets");
  }


  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */
  template <int dim>
  void
  NSBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the fluid dynamics boundary conditions. Please ensure that the id is set for every boundary condition."));

    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    NavierStokesBoundaryDuplicated(boundary_id));

        // Allocate the navier_stokes_functions object for every boundary
        // condition to ensure that they have a defined function and a center of
        // rotation.
        navier_stokes_functions[boundary_id] =
          std::make_shared<NSBoundaryFunctions<dim>>();

        prm.enter_subsection("u");
        navier_stokes_functions[boundary_id]->u.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("v");
        navier_stokes_functions[boundary_id]->v.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("w");
        navier_stokes_functions[boundary_id]->w.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("p");
        navier_stokes_functions[boundary_id]->p.parse_parameters(prm);
        prm.leave_subsection();

        prm.enter_subsection("center of rotation");
        navier_stokes_functions[boundary_id]->center_of_rotation[0] =
          prm.get_double("x");
        navier_stokes_functions[boundary_id]->center_of_rotation[1] =
          prm.get_double("y");

        if (dim == 3)
          navier_stokes_functions[boundary_id]->center_of_rotation[2] =
            prm.get_double("z");
        prm.leave_subsection();

        // Establish the type of boundary condition
        const std::string op = prm.get("type");
        if (op == "none")
          this->type[boundary_id] = BoundaryType::none;
        if (op == "noslip")
          this->type[boundary_id] = BoundaryType::noslip;
        if (op == "slip")
          this->type[boundary_id] = BoundaryType::slip;
        if (op == "function" || op == "function weak" || op == "partial slip")
          {
            if (op == "function")
              this->type[boundary_id] = BoundaryType::function;
            else if (op == "partial slip")
              this->type[boundary_id] = BoundaryType::partial_slip;
            else
              this->type[boundary_id] = BoundaryType::function_weak;
          }
        if (op == "pressure")
          {
            this->type[boundary_id] = BoundaryType::pressure;
          }
        if (op == "periodic")
          {
            types::boundary_id periodic_boundary_id =
              prm.get_integer("periodic id");

            this->type[boundary_id] = BoundaryType::periodic;

            // We attribute a periodic neighbor boundary type to the neighbor
            // boundary to ensure that all boundaries have a defined type
            this->type[periodic_boundary_id] = BoundaryType::periodic_neighbor;

            // We store the periodic id and direction
            this->periodic_neighbor_id[boundary_id] = periodic_boundary_id;
            this->periodic_direction[boundary_id] =
              prm.get_integer("periodic direction");

            // Allocate the navier_stokes_functions object for the periodic
            // neighbor condition to ensure that they have a defined function
            // and a center of rotation. Parse the information of the boundary
            // id to duplicate it.
            navier_stokes_functions[periodic_boundary_id] =
              std::make_shared<NSBoundaryFunctions<dim>>();

            prm.enter_subsection("u");
            navier_stokes_functions[periodic_boundary_id]->u.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("v");
            navier_stokes_functions[periodic_boundary_id]->v.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("w");
            navier_stokes_functions[periodic_boundary_id]->w.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("p");
            navier_stokes_functions[periodic_boundary_id]->p.parse_parameters(
              prm);
            prm.leave_subsection();

            prm.enter_subsection("center of rotation");
            navier_stokes_functions[periodic_boundary_id]
              ->center_of_rotation[0] = prm.get_double("x");
            navier_stokes_functions[periodic_boundary_id]
              ->center_of_rotation[1] = prm.get_double("y");

            if (dim == 3)
              navier_stokes_functions[periodic_boundary_id]
                ->center_of_rotation[2] = prm.get_double("z");
            prm.leave_subsection();
          }
        if (op == "outlet")
          {
            this->type[boundary_id] = BoundaryType::outlet;
          }

        this->beta[boundary_id] = prm.get_double("beta");
        this->boundary_layer_thickness[boundary_id] =
          prm.get_double("boundary layer thickness");
      }
  }


  template <int dim>
  void
  NSBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions");
    {
      prm.declare_entry("number",
                        Utilities::int_to_string(number_of_boundary_conditions),
                        Patterns::Integer(),
                        "Number of boundary conditions");
      prm.declare_entry(
        "time dependent",
        "false",
        Patterns::Bool(),
        "Bool to define if the boundary condition is time dependent");

      prm.declare_entry(
        "fix pressure constant",
        "false",
        Patterns::Bool(),
        "Bool to define if zero pressure is used as a constraint one node");

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
      this->number_of_boundary_conditions = prm.get_integer("number");
      this->time_dependent                = prm.get_bool("time dependent");
      this->fix_pressure_constant = prm.get_bool("fix pressure constant");

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
   *  - if bc type is "convection-radiation-flux" (Robin condition), "h" is the
   * convective heat transfer coefficient and "Tinf" is the
   * environment temperature at the boundary, "emissivity" is the emissivity
   * coefficient, and "Stefan-Boltzmann constant" is the Stefan-Boltzmann
   * constant = 5.6703*10-8 \f$(W.m^{-2}.K^{-4})\f$. It is also possible to
   * impose a heat flux using "heat_flux"
   *
   */

  template <int dim>
  class HTBoundaryConditions : public BoundaryConditions
  {
  public:
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      dirichlet_value;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      h;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      Tinf;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      emissivity;
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
           heat_flux_bc;
    double Stefan_Boltzmann_constant;

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


    // TODO REMOVE THIS FLAG IT IS NOT USEFUL
    bool has_convection_radiation_bc = false;
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
  HTBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry(
      "type",
      "noflux",
      Patterns::Selection(
        "noflux|temperature|convection-radiation-flux|periodic"),
      "Type of boundary condition for heat transfer"
      "Choices are <noflux|temperature|convection-radiation-flux|periodic>.");

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

    prm.declare_alias("periodic id", "periodic_id", true);

    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    prm.declare_alias("periodic direction", "periodic_direction", true);
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
  HTBoundaryConditions<dim>::declare_parameters(
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

  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   */

  template <int dim>
  void
  HTBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
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


        const std::string op = prm.get("type");
        if (op == "noflux")
          {
            this->type[boundary_id] = BoundaryType::noflux;
          }
        else if (op == "temperature")
          {
            this->type[boundary_id] = BoundaryType::temperature;
          }
        else if (op == "convection-radiation-flux")
          {
            this->type[boundary_id] = BoundaryType::convection_radiation;
            this->has_convection_radiation_bc = true;

            // Emissivity validity (0 <= emissivity <= 1) will be checked at
            // evaluation.
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
            this->periodic_neighbor_id[boundary_id] = periodic_boundary_id;
            this->periodic_direction[boundary_id] =
              prm.get_integer("periodic direction");
          }
        else
          {
            AssertThrow(
              false,
              ExcMessage("Unknown boundary condition type for heat transfer."));
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


  /**
   * @brief This class manages the boundary conditions for Tracer solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */

  template <int dim>
  class TracerBoundaryConditions : public BoundaryConditions
  {
  public:
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      tracer;


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
   * i.e. Dirichlet condition (ConstantFunction) with value 0
   *
   * @param prm A parameter handler which is currently used to parse the simulation information.
   *
   * @param default_boundary_id Default value given to the boundary id.
   */
  template <int dim>
  void
  TracerBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry("type",
                      "outlet",
                      Patterns::Selection("dirichlet|outlet|periodic"),
                      "Type of boundary condition for tracer"
                      "Choices are <dirichlet|outlet|periodic>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    Functions::ParsedFunction<dim> temporary_function;

    prm.enter_subsection("dirichlet");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Periodic boundary condition parameters for Tracer physics
    prm.declare_entry(
      "periodic id",
      "-1",
      Patterns::Integer(),
      "Mesh id for periodic face matching. Default entry is -1 to ensure that the periodic id is set by the user");

    prm.declare_alias("periodic id", "periodic_id", true);


    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    prm.declare_alias("periodic direction", "periodic_direction", true);
  }

  /**
   * @brief Declare the boundary conditions default parameters
   * Calls declareDefaultEntry method for each boundary (max 14 boundaries)
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   * @param number_of_boundary_conditions Number of tracer boundary conditions
   */
  template <int dim>
  void
  TracerBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions tracer");
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
   *
   */

  template <int dim>
  void
  TracerBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the tracer boundary conditions. Please ensure that the id is set for every boundary condition."));


    // Loop through all boundary ids to ensure that they are all non-negative
    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    TracerBoundaryDuplicated(boundary_id));

        // Allocate function for tracer
        tracer[boundary_id] =
          std::make_shared<Functions::ParsedFunction<dim>>();
        const std::string op = prm.get("type");
        if (op == "dirichlet")
          {
            this->type[boundary_id] = BoundaryType::tracer_dirichlet;
            prm.enter_subsection("dirichlet");

            tracer[boundary_id]->parse_parameters(prm);
            prm.leave_subsection();
          }
        else if (op == "outlet")
          {
            this->type[boundary_id] = BoundaryType::outlet;
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
            this->periodic_neighbor_id[boundary_id] = periodic_boundary_id;
            this->periodic_direction[boundary_id] =
              prm.get_integer("periodic direction");
          }
        else
          {
            AssertThrow(false,
                        ExcMessage(
                          "Unknown boundary condition type for tracers."));
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
  TracerBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions tracer");
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

  /**
   * @brief This class manages the boundary conditions for the Cahn-Hilliard solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */

  template <int dim>
  class CahnHilliardBoundaryConditions : public BoundaryConditions
  {
  public:
    std::map<types::boundary_id, double> angle_of_contact;
    std::map<types::boundary_id, double> phase_dirichlet_value;
    std::map<types::boundary_id,
             std::shared_ptr<CahnHilliardBoundaryFunctions<dim>>>
      bcFunctions;

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
   * @brief Declares the default parameters for the Cahn-Hilliard boundary conditions
   *
   * @param prm A parameter handler which is currently used to parse the simulation information.
   *
   * @param default_boundary_id Default value given to the boundary id.
   */
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

    prm.declare_alias("periodic id", "periodic_id", true);

    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    prm.declare_alias("periodic direction", "periodic_direction", true);
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

  /**
   * @brief Parse the information for a boundary condition
   *
   * @param prm A parameter handler which is currently used to parse the simulation information
   *
   */

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
            this->periodic_neighbor_id[boundary_id] = periodic_boundary_id;
            this->periodic_direction[boundary_id] =
              prm.get_integer("periodic direction");
          }
        else
          {
            AssertThrow(
              false,
              ExcMessage("Unknown boundary condition type for Cahn-Hilliard."));
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

  /**
   * @brief This class manages the boundary conditions for VOF solver
   * It introduces the boundary functions and declares the boundary conditions
   * coherently.
   */

  template <int dim>
  class VOFBoundaryConditions : public BoundaryConditions
  {
  public:
    std::map<types::boundary_id,
             std::shared_ptr<Functions::ParsedFunction<dim>>>
      phase_fraction;

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
  VOFBoundaryConditions<dim>::declare_default_entry(
    ParameterHandler        &prm,
    const types::boundary_id default_boundary_id)
  {
    prm.declare_entry("type",
                      "none",
                      Patterns::Selection("none|dirichlet|periodic"),
                      "Type of boundary condition for VOF"
                      "Choices are <none|dirichlet|periodic>.");

    prm.declare_entry("id",
                      Utilities::int_to_string(default_boundary_id, 2),
                      Patterns::List(Patterns::Integer()),
                      "Mesh id for boundary conditions");

    Functions::ParsedFunction<dim> temporary_function;
    prm.enter_subsection("dirichlet");
    temporary_function.declare_parameters(prm);
    prm.leave_subsection();

    // Periodic boundary condition parameters for VOF physics
    prm.declare_entry(
      "periodic id",
      "-1",
      Patterns::Integer(),
      "Mesh id for periodic face matching. Default entry is -1 to ensure that the periodic id is set by the user");

    prm.declare_alias("periodic id", "periodic_id", true);

    prm.declare_entry("periodic direction",
                      "0",
                      Patterns::Integer(),
                      "Direction for periodic boundary condition");

    prm.declare_alias("periodic direction", "periodic_direction", true);
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
  VOFBoundaryConditions<dim>::declare_parameters(
    ParameterHandler  &prm,
    const unsigned int number_of_boundary_conditions)
  {
    prm.enter_subsection("boundary conditions VOF");
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
  VOFBoundaryConditions<dim>::parse_boundary(ParameterHandler &prm)
  {
    // Parse the list of boundary ids
    std::vector<types::boundary_id> boundary_ids =
      convert_string_to_vector<types::boundary_id>(prm, "id");

    // Check if the list contains at least one boundary id
    AssertThrow(
      boundary_ids.size() > 0,
      ExcMessage(
        "A boundary id has not been set for one of the VOF boundary conditions. Please ensure that the id is set for every boundary condition."));

    // Loop through all boundary ids to ensure that they are all non-negative
    // and unique
    for (const auto boundary_id : boundary_ids)
      {
        AssertThrow(this->type.find(boundary_id) == this->type.end(),
                    VOFBoundaryDuplicated(boundary_id));

        if (auto const option = prm.get("type"); option == "none")
          {
            this->type[boundary_id] = BoundaryType::none;
          }

        if (auto const option = prm.get("type"); option == "dirichlet")
          {
            this->type[boundary_id] = BoundaryType::vof_dirichlet;
            prm.enter_subsection("dirichlet");
            phase_fraction[boundary_id] =
              std::make_shared<Functions::ParsedFunction<dim>>();
            phase_fraction[boundary_id]->parse_parameters(prm);
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
            this->periodic_neighbor_id[boundary_id] = periodic_boundary_id;
            this->periodic_direction[boundary_id] =
              prm.get_integer("periodic direction");
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
  VOFBoundaryConditions<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("boundary conditions VOF");
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
 * @brief This class implements a boundary conditions for the Navier-Stokes equation
 * where the velocity component are defined using individual functions
 */
template <int dim>
class NavierStokesFunctionDefined : public Function<dim>
{
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

  double
  value(const Point<dim> &point, const unsigned int component) const override;
};


/**
 * @brief Calculates the value of a Function-type Navier-Stokes equations
 *
 * @param point A point at which the function will be evaluated
 *
 * @param component The vector component of the boundary condition (0-x, 1-y and 2-z)
 */
template <int dim>
double
NavierStokesFunctionDefined<dim>::value(const Point<dim>  &point,
                                        const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  if (component == 0)
    {
      return u->value(point);
    }
  if (component == 1)
    {
      return v->value(point);
    }
  if (component == 2)
    {
      return w->value(point);
    }
  return 0.;
}

/**
 * @brief This class implements a pressure boundary condition for the Navier-Stokes equations.
 */
template <int dim>
class NavierStokesPressureFunctionDefined : public Function<dim>
{
  Functions::ParsedFunction<dim> *p;

public:
  NavierStokesPressureFunctionDefined(Functions::ParsedFunction<dim> *p_p)
    : Function<dim>(dim + 1)
    , p(p_p)
  {}

  double
  value(const Point<dim> &point, unsigned int component) const override;
};


/**
 * @brief Calculates the value of a Function-type Navier-Stokes equations
 *
 * @param point A point (generally a gauss point)
 *
 * @param component The vector component of the boundary condition (0-x, 1-y and 2-z)
 */
template <int dim>
double
NavierStokesPressureFunctionDefined<dim>::value(
  const Point<dim>  &point,
  const unsigned int component) const
{
  if (component == dim)
    {
      return p->value(point);
    }
  return 0.;
}

/**
 * @brief Class that implements a boundary conditions for the Cahn-Hilliard equation
 * where the phase and chemical potential are defined using individual
 * functions
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
  value(const Point<dim> &p, unsigned int component) const override;
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
CahnHilliardFunctionDefined<dim>::value(const Point<dim>  &p,
                                        const unsigned int component) const
{
  Assert(component < this->n_components,
         ExcIndexRange(component, 0, this->n_components));

  if (component == 0)
    {
      return phi->value(p);
    }
  if (component == 1)
    {
      return 0.;
    }
  return 0.;
}

#endif
