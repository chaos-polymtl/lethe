// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters.h>
#include <core/utilities.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>

DeclException2(
  PhaseChangeIntervalError,
  double,
  double,
  << "Liquidus temperature : " << arg1
  << " is not strictly superior to Solidus temperature: " << arg2
  << " The liquidus temperature specific is below or equal to the solidus temperature."
  << " The phase change specific heat model requires that T_liquidus>T_solidus.");

DeclException1(
  NumberOfFluidsError,
  int,
  << "Number of fluids: " << arg1
  << " is not 1 (single phase simulation) or 2 (VOF simulation). This is currently not supported.");

DeclException1(NumberOfSolidsError,
               int,
               << "Number of solids: " << arg1
               << " is larger than 1. This is currently not supported.");

DeclException1(NumberOfMaterialInteractionsError,
               int,
               << "Number of material interactions: " << arg1
               << " is larger than 3. This is currently not supported.");

DeclException2(
  OrderOfFluidIDsError,
  int,
  int,
  << "The first fluid's id is " << arg1 << " and the id of the second fluid is "
  << arg2 << ". The first fluid's id should be lower than the second fluid's.");

DeclException1(
  TwoDimensionalLaserError,
  unsigned int,
  << "Laser beam orientation in : " << arg1
  << "-dimensional simulations cannot be defined in the z direction");

DeclException3(MultipleAdaptationSizeError,
               std::string,
               unsigned int,
               unsigned int,
               << "Error in 'mesh adaptation' : number of '" << arg1 << "' ("
               << arg2 << ") does not correspond to the number of 'variables' ("
               << arg3 << ")");

DeclException3(ParameterStrictlyGreaterThanError,
               std::string,
               double,
               double,
               << "The parameter '" << arg1 << "' is set to: " << arg2
               << ". However, it should be strictly greater than " << arg3
               << ".");

// DeclException

namespace Parameters
{
  SizeOfSubsections
  get_size_of_subsections(const std::string &file_name)
  {
    SizeOfSubsections sizes;
    sizes.boundary_conditions = get_max_subsection_size(file_name);
    sizes.manifolds           = get_max_subsection_size(file_name);
    return sizes;
  }

  void
  SimulationControl::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("simulation control");
    {
      prm.declare_entry(
        "method",
        "steady",
        Patterns::Selection(
          "steady|steady_bdf|bdf1|bdf2|bdf3|sdirk22|sdirk33|sdirk43"),
        "The time integration scheme. "
        "Choices are <steady|steady_bdf|bdf1|bdf2|bdf3|sdirk22|sdirk33|sdirk43>.");

      prm.declare_entry(
        "bdf startup method",
        "multiple step bdf",
        Patterns::Selection("multiple step bdf|initial solution"),
        "The kind of method used to startup high order bdf methods "
        "Choices are <multiple step bdf|initial solution>.");

      prm.declare_entry("time step",
                        "1.",
                        Patterns::Double(),
                        "Time step value");
      prm.declare_entry("time end", "1", Patterns::Double(), "Time step value");
      prm.declare_entry("startup time scaling",
                        "0.4",
                        Patterns::Double(),
                        "Scaling factor used in the iterations necessary to "
                        "start-up the BDF schemes.");
      prm.declare_entry("adapt",
                        "false",
                        Patterns::Bool(),
                        "Adaptative time-stepping <true|false>");
      prm.declare_entry(
        "override time step on restart",
        "false",
        Patterns::Bool(),
        "Override checkpointed time-step upon restart <true|false>");
      prm.declare_entry(
        "time step independent of end time",
        "true",
        Patterns::Bool(),
        "Ensures that the correct time step is kept when using adaptive time step simulations");
      prm.declare_entry("number mesh adapt",
                        "0",
                        Patterns::Integer(),
                        "Number of mesh adaptation (for steady simulations)");
      prm.declare_entry("max cfl",
                        "1",
                        Patterns::Double(),
                        "Maximum CFL value");
      prm.declare_entry("max time step",
                        "1e6",
                        Patterns::Double(),
                        "Maximum time step value");
      prm.declare_entry("stop tolerance",
                        "1e-10",
                        Patterns::Double(),
                        "Tolerance at which the simulation is stopped");

      prm.declare_entry("adaptative time step scaling",
                        "1.1",
                        Patterns::Double(),
                        "Adaptative time step scaling");
      prm.declare_entry("output path",
                        "./",
                        Patterns::FileName(),
                        "File output prefix");

      prm.declare_entry("output name",
                        "out",
                        Patterns::FileName(),
                        "File output prefix");


      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output iteration frequency");

      prm.declare_entry("output time frequency",
                        "-1",
                        Patterns::Double(),
                        "Output time frequency");

      prm.declare_entry(
        "output boundaries",
        "false",
        Patterns::Bool(),
        "Output the boundaries of the domain along with their ID");

      prm.declare_entry("log frequency",
                        "1",
                        Patterns::Integer(),
                        "log frequency");

      prm.declare_entry("log precision",
                        "6",
                        Patterns::Integer(),
                        "Display precision when writing to log",
                        "This setting percolates to all output to the log");

      prm.declare_entry("output times",
                        "-1",
                        Patterns::List(Patterns::Double()),
                        "List of specific output times separated with a comma");

      prm.declare_entry(
        "output control",
        "iteration",
        Patterns::Selection("iteration|time"),
        "The control for the output of the simulation results"
        "Results can be either outputted at constant iteration frequency or at constant time");

      prm.declare_entry("output time interval",
                        "0, 1.7976931348623157e308", // 0, Maximum double
                        Patterns::List(Patterns::Double()),
                        "Output files for a desired time interval");

      prm.declare_entry("subdivision",
                        "1",
                        Patterns::Integer(),
                        "Subdivision of mesh cell in postprocessing");

      prm.declare_entry("group files",
                        "1",
                        Patterns::Integer(),
                        "Maximal number of vtu output files");
    }
    prm.leave_subsection();
  }

  void
  SimulationControl::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("simulation control");
    {
      const std::string sv = prm.get("method");
      if (sv == "steady")
        method = TimeSteppingMethod::steady;
      else if (sv == "steady_bdf")
        method = TimeSteppingMethod::steady_bdf;
      else if (sv == "bdf1")
        method = TimeSteppingMethod::bdf1;
      else if (sv == "bdf2")
        method = TimeSteppingMethod::bdf2;
      else if (sv == "bdf3")
        method = TimeSteppingMethod::bdf3;
      else if (sv == "sdirk22")
        method = TimeSteppingMethod::sdirk22;
      else if (sv == "sdirk33")
        method = TimeSteppingMethod::sdirk33;
      else if (sv == "sdirk43")
        method = TimeSteppingMethod::sdirk43;
      else
        {
          AssertThrow(false, ExcMessage("Invalid time stepping scheme"));
        }
      const std::string bdf_startup_string = prm.get("bdf startup method");
      if (bdf_startup_string == "multiple step bdf")
        bdf_startup_method = BDFStartupMethods::multiple_step_bdf;
      else if (bdf_startup_string == "initial solution")
        bdf_startup_method = BDFStartupMethods::initial_solution;
      else
        {
          AssertThrow(false, ExcMessage("Invalid bdf startup scheme"));
        }

      const std::string osv = prm.get("output control");
      if (osv == "iteration")
        output_control = OutputControl::iteration;
      else if (osv == "time")
        output_control = OutputControl::time;
      else
        {
          AssertThrow(false, ExcMessage("Invalid output control scheme"));
        }
      dt       = prm.get_double("time step");
      time_end = prm.get_double("time end");
      adapt    = prm.get_bool("adapt");
      time_step_independent_of_end_time =
        prm.get_bool("time step independent of end time");
      maxCFL         = prm.get_double("max cfl");
      max_dt         = prm.get_double("max time step");
      stop_tolerance = prm.get_double("stop tolerance");
      adaptative_time_step_scaling =
        prm.get_double("adaptative time step scaling");
      startup_timestep_scaling = prm.get_double("startup time scaling");
      number_mesh_adaptation   = prm.get_integer("number mesh adapt");
      override_time_step_on_restart =
        prm.get_bool("override time step on restart");

      output_folder = prm.get("output path");
      output_name   = prm.get("output name");
      std::erase(output_name, '/');
      output_iteration_frequency = prm.get_integer("output frequency");
      output_time_frequency      = prm.get_double("output time frequency");
      output_times_vector =
        convert_string_to_vector<double>(prm, "output times");
      output_time_interval =
        convert_string_to_vector<double>(prm, "output time interval");
      output_boundaries = prm.get_bool("output boundaries");

      subdivision   = prm.get_integer("subdivision");
      group_files   = prm.get_integer("group files");
      log_frequency = prm.get_integer("log frequency");
      log_precision = prm.get_integer("log precision");
    }
    prm.leave_subsection();
  } // namespace Parameters

  void
  Timer::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("timer");
    {
      prm.declare_entry("type",
                        "none",
                        Patterns::Selection("none|iteration|end"),
                        "Clock monitoring methods "
                        "Choices are <none|iteration|end>.");
      prm.declare_entry(
        "write time in error table",
        "false",
        Patterns::Bool(),
        "Boolean to define if the time is written in the error table");
    }
    prm.leave_subsection();
  }

  void
  Timer::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("timer");
    {
      const std::string cl = prm.get("type");
      if (cl == "none")
        type = Type::none;
      else if (cl == "iteration")
        type = Type::iteration;
      else if (cl == "end")
        type = Type::end;
      write_time_in_error_table = prm.get_bool("write time in error table");
    }
    prm.leave_subsection();
  }

  void
  PowerLawParameters::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("power-law");
    {
      prm.declare_entry("K",
                        "1.0",
                        Patterns::Double(),
                        "Fluid consistency index");
      prm.declare_entry("n", "0.5", Patterns::Double(), "Flow behavior index");
      prm.declare_entry("shear rate min",
                        "0.001",
                        Patterns::Double(),
                        "Minimal shear rate magnitude");
    }
    prm.leave_subsection();
  }

  void
  PowerLawParameters::parse_parameters(ParameterHandler     &prm,
                                       const Dimensionality &dimensions)
  {
    prm.enter_subsection("power-law");
    {
      K = prm.get_double("K");
      // K is in L^2 T^-1
      K *= dimensions.viscosity_scaling;

      // n is dimensionless
      n = prm.get_double("n");

      // The shear rate min is in T^-1
      shear_rate_min = prm.get_double("shear rate min");
      shear_rate_min *= dimensions.time;
    }
    prm.leave_subsection();
  }

  void
  CarreauParameters::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("carreau");
    {
      prm.declare_entry("viscosity_0",
                        "1.0",
                        Patterns::Double(),
                        "Kinematic viscosity at rest");
      prm.declare_entry("viscosity_inf",
                        "1.0",
                        Patterns::Double(),
                        "Kinematic viscosity for an infinite shear rate");
      prm.declare_entry("lambda", "1.0", Patterns::Double(), "Relaxation time");
      prm.declare_entry("a", "2.0", Patterns::Double(), "Carreau parameter");
      prm.declare_entry("n", "0.5", Patterns::Double(), "Power parameter");
    }
    prm.leave_subsection();
  }

  void
  CarreauParameters::parse_parameters(ParameterHandler     &prm,
                                      const Dimensionality &dimensions)
  {
    prm.enter_subsection("carreau");
    {
      kinematic_viscosity_0   = prm.get_double("viscosity_0");
      kinematic_viscosity_inf = prm.get_double("viscosity_inf");

      // Both kinematic viscosities are in L^2 T^-1
      kinematic_viscosity_0 *= dimensions.viscosity_scaling;
      kinematic_viscosity_inf *= dimensions.viscosity_scaling;

      lambda = prm.get_double("lambda");

      // lambda is in T
      lambda *= 1. / dimensions.time;

      // a and n are dimensionless
      a = prm.get_double("a");
      n = prm.get_double("n");
    }
    prm.leave_subsection();
  }

  void
  NonNewtonian::declare_parameters(ParameterHandler &prm) const
  {
    prm.enter_subsection("non newtonian");
    {
      powerlaw_parameters.declare_parameters(prm);
      carreau_parameters.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  void
  NonNewtonian::parse_parameters(ParameterHandler     &prm,
                                 const Dimensionality &dimensions)
  {
    prm.enter_subsection("non newtonian");
    {
      powerlaw_parameters.parse_parameters(prm, dimensions);
      carreau_parameters.parse_parameters(prm, dimensions);
    }
    prm.leave_subsection();
  }


  void
  ImmersedSolidTanhParameters::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("immersed solid tanh");
    {
      prm.declare_entry("tracer diffusivity inside",
                        "1.0",
                        Patterns::Double(),
                        "Tracer diffusivity inside the immersed solid");
      prm.declare_entry("tracer diffusivity outside",
                        "1.0",
                        Patterns::Double(),
                        "Tracer diffusivity outside the immersed solid");
      prm.declare_entry("tracer reaction constant inside",
                        "0.0",
                        Patterns::Double(),
                        "Tracer reaction constant inside the immersed solid");
      prm.declare_entry("tracer reaction constant outside",
                        "0.0",
                        Patterns::Double(),
                        "Tracer reaction constant outside the immersed solid");
      prm.declare_entry("thickness",
                        "1.0",
                        Patterns::Double(),
                        "Thickness to be used with the tanh function");
    }
    prm.leave_subsection();
  }

  void
  ImmersedSolidTanhParameters::parse_parameters(
    ParameterHandler     &prm,
    const Dimensionality &dimensions)
  {
    prm.enter_subsection("immersed solid tanh");
    {
      tracer_diffusivity_inside  = prm.get_double("tracer diffusivity inside");
      tracer_diffusivity_outside = prm.get_double("tracer diffusivity outside");
      tracer_reaction_constant_inside =
        prm.get_double("tracer reaction constant inside");
      tracer_reaction_constant_outside =
        prm.get_double("tracer reaction constant outside");
      thickness = prm.get_double("thickness");

      // Diffusivity is in L^2 T^-1
      tracer_diffusivity_inside *= dimensions.diffusivity_scaling;
      tracer_diffusivity_outside *= dimensions.diffusivity_scaling;
    }
    prm.leave_subsection();
  }

  void
  ImmersedSolidGaussianParameters::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("immersed solid gaussian");
    {
      prm.declare_entry("tracer diffusivity interface",
                        "0.0",
                        Patterns::Double(),
                        "Tracer diffusivity at the immersed solid interface");
      prm.declare_entry("tracer diffusivity bulk",
                        "1.0",
                        Patterns::Double(),
                        "Tracer diffusivity in the phase bulk");
      prm.declare_entry(
        "tracer reaction constant interface",
        "0.0",
        Patterns::Double(),
        "Tracer reaction constant at the immersed solid interface");
      prm.declare_entry("tracer reaction constant bulk",
                        "0.0",
                        Patterns::Double(),
                        "Tracer reaction constant in the phase bulk");
      prm.declare_entry("thickness",
                        "1.0",
                        Patterns::Double(),
                        "Thickness to be used with the Gaussian function");
    }
    prm.leave_subsection();
  }

  void
  ImmersedSolidGaussianParameters::parse_parameters(
    ParameterHandler     &prm,
    const Dimensionality &dimensions)
  {
    prm.enter_subsection("immersed solid gaussian");
    {
      tracer_diffusivity_interface =
        prm.get_double("tracer diffusivity interface");
      tracer_diffusivity_bulk = prm.get_double("tracer diffusivity bulk");
      tracer_reaction_constant_interface =
        prm.get_double("tracer reaction constant interface");
      tracer_reaction_constant_bulk =
        prm.get_double("tracer reaction constant bulk");
      thickness = prm.get_double("thickness");

      // Diffusivity is in L^2 T^-1
      tracer_diffusivity_interface *= dimensions.diffusivity_scaling;
      tracer_diffusivity_bulk *= dimensions.diffusivity_scaling;
    }
    prm.leave_subsection();
  }

  void
  IsothermalIdealGasDensityParameters::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("isothermal_ideal_gas");
    {
      prm.declare_entry(
        "density_ref",
        "1.2", // dry air's density at normal temperature and pressure (20 °C
               // and 1 atm)
        Patterns::Double(),
        "Reference density of the gas in SI units for isothermal ideal gas equation of state in density calculation");

      prm.declare_entry(
        "R",
        "287.05", // dry air's specific gas constant as default
        Patterns::Double(),
        "Specific gas constant in SI units for isothermal ideal gas equation of state in density calculation");

      prm.declare_entry(
        "T",
        "293.15", // normal temperature (20°C) as a default
        Patterns::Double(),
        "Absolute temperature of the gas in kelvin (K) for isothermal ideal gas equation of state in density calculation");
    }
    prm.leave_subsection();
  }

  void
  IsothermalIdealGasDensityParameters::parse_parameters(
    ParameterHandler     &prm,
    const Dimensionality &dimensions)
  {
    prm.enter_subsection("isothermal_ideal_gas");
    {
      // Isothermal ideal gas equation of state parameters
      // The reference state density of the gas (rho_{ref}) is in M L^-3
      density_ref = prm.get_double("density_ref");
      density_ref *= dimensions.density_scaling;

      // The specific gas constant (R) is in L^2 T^-2 theta^-1
      R = prm.get_double("R");
      R *= dimensions.specific_gas_constant_scaling;

      // The absolute temperature (T) of the ideal gas is in theta^-1
      T = prm.get_double("T");
      T *= 1. / dimensions.temperature;
    }
    prm.leave_subsection();
  }

  void
  SurfaceTensionParameters::declare_parameters(dealii::ParameterHandler &prm)
  {
    prm.declare_entry(
      "surface tension coefficient",
      "0.0",
      Patterns::Double(),
      "Surface tension coefficient for the corresponding pair of fluids or fluid-solid pair");
    prm.declare_entry(
      "reference state temperature",
      "0.0",
      Patterns::Double(),
      "Temperature of the reference state corresponding to the surface tension coefficient");
    prm.declare_entry(
      "temperature-driven surface tension gradient",
      "0.0",
      Patterns::Double(),
      "Surface tension gradient with respect to the temperature for the corresponding pair of fluids or fluid-solid pair");
    prm.declare_entry(
      "solidus temperature",
      "0",
      Patterns::Double(),
      "Temperature of the solidus for the corresponding pair of fluids or fluid-solid pair");
    prm.declare_entry(
      "liquidus temperature",
      "1",
      Patterns::Double(),
      "Temperature of the liquidus for the corresponding pair of fluids or fluid-solid pair");
  }

  void
  SurfaceTensionParameters::parse_parameters(
    const ParameterHandler           &prm,
    const Parameters::Dimensionality &dimensions)
  {
    surface_tension_coefficient = prm.get_double("surface tension coefficient");
    surface_tension_coefficient *= dimensions.surface_tension_scaling;
    T_0 = prm.get_double("reference state temperature");
    T_0 *= 1. / dimensions.temperature;
    surface_tension_gradient =
      prm.get_double("temperature-driven surface tension gradient");
    surface_tension_gradient *= dimensions.surface_tension_gradient_scaling;
    T_solidus = prm.get_double("solidus temperature");
    T_solidus *= 1. / dimensions.temperature;
    T_liquidus = prm.get_double("liquidus temperature");
    T_liquidus *= 1. / dimensions.temperature;
    Assert(T_liquidus > T_solidus,
           PhaseChangeIntervalError(T_liquidus, T_solidus));
  }

  void
  MobilityCahnHilliardParameters::declare_parameters(
    dealii::ParameterHandler &prm)
  {
    prm.declare_entry(
      "cahn hilliard mobility constant",
      "1e-7",
      Patterns::Double(),
      "Cahn-Hilliard mobility constant for the corresponding pair of fluids");
  }

  void
  MobilityCahnHilliardParameters::parse_parameters(
    const ParameterHandler           &prm,
    const Parameters::Dimensionality &dimensions)
  {
    mobility_cahn_hilliard_constant =
      prm.get_double("cahn hilliard mobility constant");
    mobility_cahn_hilliard_constant *=
      dimensions.cahn_hilliard_mobility_scaling;
  }

  template <int dim>
  void
  ConstrainSolidDomain<dim>::declare_parameters(
    dealii::ParameterHandler &prm,
    const unsigned int        max_number_of_constraints)
  {
    prm.enter_subsection("constrain stasis");
    {
      prm.declare_entry(
        "enable",
        "false",
        Patterns::Bool(),
        "Enable/disable (true/false) the solid domain constraining feature.");

      prm.declare_entry(
        "enable domain restriction with plane",
        "false",
        Patterns::Bool(),
        "Enable/disable (true/false) the definition of a plane for geometrical\n"
        " restrictions on the domain where the solid domain constraining feature\n"
        " is applied.");
      std::string default_entry_sting = (dim == 2) ? "0., 0." : "0., 0., 0.";
      prm.declare_entry("restriction plane point",
                        default_entry_sting,
                        Patterns::List(Patterns::Double()),
                        "Domain restriction plane point coordinates.");
      prm.declare_entry(
        "restriction plane normal vector",
        default_entry_sting,
        Patterns::List(Patterns::Double()),
        "Domain restriction plane outward pointing normal vector.");

      prm.declare_entry(
        "number of constraints",
        "0",
        Patterns::Integer(),
        "Number of solid constraints (maximum of 1 per fluid).");

      // Resize vectors
      this->fluid_ids.resize(max_number_of_constraints);
      this->filtered_phase_fraction_tolerance.resize(max_number_of_constraints);
      this->temperature_min_values.resize(max_number_of_constraints);
      this->temperature_max_values.resize(max_number_of_constraints);

      // Declare default entries
      for (unsigned int c_id = 0; c_id < max_number_of_constraints; ++c_id)
        {
          prm.enter_subsection("constraint " + std::to_string(c_id));
          {
            declare_default_entries(prm);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  ConstrainSolidDomain<dim>::declare_default_entries(
    dealii::ParameterHandler &prm)
  {
    prm.declare_entry("fluid id",
                      "0",
                      Patterns::Integer(),
                      "Identifier of the fluid material that is constrained.");
    prm.declare_entry("phase fraction tolerance",
                      "1e-4",
                      Patterns::Double(),
                      "Absolute filtered phase fraction tolerance used in "
                      "conjunction with VOF simulations to select the cells "
                      "on which the constraint is applied.");
    prm.declare_entry("min temperature",
                      "-999",
                      Patterns::Double(),
                      "Minimum temperature value of the fluid for it be "
                      "considered as a solid.");
    prm.declare_entry("max temperature",
                      "0",
                      Patterns::Double(),
                      "Maximum temperature value of the fluid for it be "
                      "considered as a solid.");
  }

  template <int dim>
  void
  ConstrainSolidDomain<dim>::parse_parameters(dealii::ParameterHandler &prm)
  {
    prm.enter_subsection("constrain stasis");
    {
      this->enable = prm.get_bool("enable");

      // Restriction plane parameters
      this->enable_domain_restriction_with_plane =
        prm.get_bool("enable domain restriction with plane");
      this->restriction_plane_point =
        value_string_to_tensor<dim>(prm.get("restriction plane point"));
      this->restriction_plane_normal_vector =
        value_string_to_tensor<dim>(prm.get("restriction plane normal vector"));

      this->number_of_constraints = prm.get_integer("number of constraints");

      // Resize vectors
      this->fluid_ids.resize(number_of_constraints);
      this->filtered_phase_fraction_tolerance.resize(number_of_constraints);
      this->temperature_min_values.resize(number_of_constraints);
      this->temperature_max_values.resize(number_of_constraints);

      // Parse parameters for each constraint
      for (unsigned int c_id = 0; c_id < number_of_constraints; ++c_id)
        {
          prm.enter_subsection("constraint " + std::to_string(c_id));
          {
            parse_constraint_parameters(prm, c_id);
          }
          prm.leave_subsection();
        }
      prm.leave_subsection();
    }
  }

  template <int dim>
  void
  ConstrainSolidDomain<dim>::parse_constraint_parameters(
    const dealii::ParameterHandler &prm,
    const unsigned int              constraint_id)
  {
    this->fluid_ids[constraint_id] = prm.get_integer("fluid id");
    this->filtered_phase_fraction_tolerance[constraint_id] =
      prm.get_double("phase fraction tolerance");
    this->temperature_min_values[constraint_id] =
      prm.get_double("min temperature");
    this->temperature_max_values[constraint_id] =
      prm.get_double("max temperature");
  }

  void
  Stabilization::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("stabilization");
    {
      prm.declare_entry(
        "use default stabilization",
        "true",
        Patterns::Bool(),
        "Use the default stabilization method provided by the solver");
      prm.declare_entry(
        "stabilization",
        "pspg_supg",
        Patterns::Selection("pspg_supg|gls|grad_div"),
        "Type of stabilization used for the Navier-Stokes equations"
        "Choices are <pspg_supg|gls|grad_div>.");

      prm.declare_entry(
        "heat transfer dcdd stabilization",
        "false",
        Patterns::Bool(),
        "Apply Discontinuity-Capturing Directional Dissipation (DCDD) "
        "stabilization term on heat transfer <true|false>");

      prm.declare_entry(
        "vof dcdd stabilization",
        "true",
        Patterns::Bool(),
        "Apply Discontinuity-Capturing Directional Dissipation (DCDD) "
        "stabilization term on the VOF phase fraction <true|false>");

      prm.declare_entry(
        "pressure scaling factor",
        "1",
        Patterns::Double(),
        "This parameter can be used to change the scale of pressure in the "
        "Navier-Stokes equations. When the velocity and pressure scales are very "
        "different, using this parameter allows to reduce the condition number"
        " and reach a solution.");

      prm.declare_entry(
        "scalar limiter",
        "none",
        Patterns::Selection("none|moe"),
        "Type of scalar limiter. The limiters are only appropriate with the DG versions of the solvers and should only be used for advection-dominated problem.");
    }
    prm.leave_subsection();
  }

  void
  Stabilization::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("stabilization");
    {
      use_default_stabilization = prm.get_bool("use default stabilization");
      {
        std::string op = prm.get("stabilization");
        if (op == "pspg_supg")
          stabilization = NavierStokesStabilization::pspg_supg;
        else if (op == "gls")
          stabilization = NavierStokesStabilization::gls;
        else if (op == "grad_div")
          stabilization = NavierStokesStabilization::grad_div;
        else
          throw(std::runtime_error("Invalid stabilization strategy"));
      }
      {
        std::string op = prm.get("scalar limiter");
        if (op == "none")
          scalar_limiter = ScalarLimiters::none;
        else if (op == "moe")
          scalar_limiter = ScalarLimiters::moe;
        else
          throw(std::runtime_error("Invalid scalar limiter"));
      }

      // DCDD stabilization activation parameters
      heat_transfer_dcdd_stabilization =
        prm.get_bool("heat transfer dcdd stabilization");
      vof_dcdd_stabilization = prm.get_bool("vof dcdd stabilization");

      pressure_scaling_factor = prm.get_double("pressure scaling factor");
    }
    prm.leave_subsection();
  }

  void
  PhaseChange::parse_parameters(ParameterHandler     &prm,
                                const Dimensionality &dimensions)
  {
    prm.enter_subsection("phase change");
    {
      T_solidus = prm.get_double("solidus temperature");
      // T_solidus has units of theta
      T_solidus *= 1 / dimensions.temperature;

      T_liquidus = prm.get_double("liquidus temperature");
      // T_liquidus has units of theta
      T_liquidus *= 1 / dimensions.temperature;

      latent_enthalpy = prm.get_double("latent enthalpy");
      // latent_enthalpy  in M L^2 T^-2
      latent_enthalpy *= dimensions.enthalpy_scaling;
      cp_l = prm.get_double("specific heat liquid");
      // cp_l  in L^2 theta^-1 T^-2
      cp_l *= dimensions.specific_heat_scaling;
      cp_s = prm.get_double("specific heat solid");
      // cp_s  in L^2 theta^-1 T^-2
      cp_s *= dimensions.specific_heat_scaling;
      kinematic_viscosity_l = prm.get_double("viscosity liquid");
      // viscosity_l  in L^2 T^-1
      kinematic_viscosity_l *= dimensions.viscosity_scaling;
      kinematic_viscosity_s = prm.get_double("viscosity solid");
      // viscosity_l  in L^2 T^-1
      kinematic_viscosity_s *= dimensions.viscosity_scaling;
      thermal_conductivity_l = prm.get_double("thermal conductivity liquid");
      // thermal_conductivity_l is in M L T^-3 theta ^-1
      thermal_conductivity_l *= dimensions.thermal_conductivity_scaling;
      thermal_conductivity_s = prm.get_double("thermal conductivity solid");
      // thermal_conductivity_s is in M L T^-3 theta ^-1
      thermal_conductivity_s *= dimensions.thermal_conductivity_scaling;
      thermal_expansion_l = prm.get_double("thermal expansion liquid");
      // thermal_expansion_l is in theta^-1
      thermal_expansion_l *= dimensions.thermal_expansion_scaling;
      thermal_expansion_s = prm.get_double("thermal expansion solid");
      // thermal_expansion_l is in theta^-1
      thermal_expansion_s *= dimensions.thermal_expansion_scaling;

      // Darcy penalty terms
      penalty_l = prm.get_double("Darcy penalty liquid");
      penalty_s = prm.get_double("Darcy penalty solid");
    }

    Assert(T_liquidus > T_solidus,
           PhaseChangeIntervalError(T_liquidus, T_solidus));

    prm.leave_subsection();
  }


  void
  PhaseChange::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("phase change");
    {
      prm.declare_entry("solidus temperature",
                        "0",
                        Patterns::Double(),
                        "Temperature of the solidus");
      prm.declare_entry("liquidus temperature",
                        "1",
                        Patterns::Double(),
                        "Temperature of the liquidus");
      prm.declare_entry("latent enthalpy",
                        "1",
                        Patterns::Double(),
                        "Enthalpy of the phase change");

      prm.declare_entry("specific heat liquid",
                        "1",
                        Patterns::Double(),
                        "Specific heat of the liquid phase");

      prm.declare_entry("specific heat solid",
                        "1",
                        Patterns::Double(),
                        "Specific heat of the solid phase");

      prm.declare_entry("thermal conductivity liquid",
                        "1",
                        Patterns::Double(),
                        "Thermal conductivity of the liquid phase");

      prm.declare_entry("thermal conductivity solid",
                        "1",
                        Patterns::Double(),
                        "Thermal conductivity of the solid phase");

      prm.declare_entry("thermal expansion liquid",
                        "1",
                        Patterns::Double(),
                        "Thermal expansion coefficient of the liquid phase");

      prm.declare_entry("thermal expansion solid",
                        "0",
                        Patterns::Double(),
                        "Thermal expansion coefficient of the solid phase");

      prm.declare_entry("viscosity liquid",
                        "1",
                        Patterns::Double(),
                        "Kinematic viscosity of the liquid phase");

      prm.declare_entry("viscosity solid",
                        "1",
                        Patterns::Double(),
                        "Kinematic viscosity of the solid phase");

      prm.declare_entry("Darcy penalty liquid",
                        "0",
                        Patterns::Double(),
                        "Darcy penalty of the liquid phase");

      prm.declare_entry("Darcy penalty solid",
                        "0",
                        Patterns::Double(),
                        "Darcy penalty of the solid phase");
    }
    prm.leave_subsection();
  }



  void
  PhysicalProperties::declare_parameters(ParameterHandler &prm)
  {
    fluids.resize(max_fluids);
    number_of_fluids = 1;
    solids.resize(max_solids);
    number_of_solids = 0;
    material_interactions.resize(max_material_interactions);
    number_of_material_interactions = 0;

    prm.enter_subsection("physical properties");
    {
      prm.declare_entry("number of fluids",
                        "1",
                        Patterns::Integer(),
                        "Number of fluids");

      // Definition of the fluids in the simulation
      for (unsigned int i_fluid = 0; i_fluid < max_fluids; ++i_fluid)
        {
          fluids[i_fluid].declare_parameters(prm, "fluid", i_fluid);
        }

      prm.declare_entry("number of solids",
                        "0",
                        Patterns::Integer(),
                        "Number of solids");

      // Definition of the solids in the simulation
      for (unsigned int i_solid = 0; i_solid < max_solids; ++i_solid)
        {
          solids[i_solid].declare_parameters(prm, "solid", i_solid);
        }

      prm.declare_entry(
        "reference temperature",
        "0",
        Patterns::Double(),
        "Reference temperature used for the calculation of physical properties and thermal expansion");
    }

    // Definition of interactions between materials
    prm.declare_entry(
      "number of material interactions",
      "0",
      Patterns::Integer(),
      "Number of material interactions (either fluid-fluid or fluid-solid)");

    for (unsigned int i_material_interaction = 0;
         i_material_interaction < max_material_interactions;
         ++i_material_interaction)
      {
        material_interactions[i_material_interaction].declare_parameters(
          prm, i_material_interaction);
      }

    prm.leave_subsection();
  }

  void
  PhysicalProperties::parse_parameters(ParameterHandler     &prm,
                                       const Dimensionality &dimensions)
  {
    prm.enter_subsection("physical properties");
    {
      // Multiphase simulations parameters definition
      number_of_fluids = prm.get_integer("number of fluids");
      AssertThrow(number_of_fluids <= max_fluids,
                  NumberOfFluidsError(number_of_fluids));

      for (unsigned int i_fluid = 0; i_fluid < number_of_fluids; ++i_fluid)
        {
          fluids[i_fluid].parse_parameters(prm, "fluid", i_fluid, dimensions);
        }

      // Multiphase simulations parameters definition
      number_of_solids = prm.get_integer("number of solids");
      for (unsigned int i_solid = 0; i_solid < number_of_solids; ++i_solid)
        {
          solids[i_solid].parse_parameters(prm, "solid", i_solid, dimensions);
        }
      AssertThrow(number_of_solids <= max_solids,
                  NumberOfSolidsError(number_of_solids));

      // Definition of interactions between materials
      number_of_material_interactions =
        prm.get_integer("number of material interactions");
      AssertThrow(number_of_material_interactions <= max_material_interactions,
                  NumberOfMaterialInteractionsError(
                    number_of_material_interactions));
      material_interactions.resize(number_of_material_interactions);
      for (unsigned int i_material_interaction = 0;
           i_material_interaction < number_of_material_interactions;
           ++i_material_interaction)
        {
          material_interactions[i_material_interaction].parse_parameters(
            prm, i_material_interaction, dimensions);
          if (material_interactions[i_material_interaction]
                .material_interaction_type ==
              MaterialInteractions::MaterialInteractionsType::fluid_fluid)
            fluid_fluid_interactions_with_material_interaction_ids.insert(
              material_interactions[i_material_interaction]
                .fluid_fluid_interaction_with_material_interaction_id);
          else // fluid-solid interaction
            fluid_solid_interactions_with_material_interaction_ids.insert(
              material_interactions[i_material_interaction]
                .fluid_solid_interaction_with_material_interaction_id);
        }

      reference_temperature = prm.get_double("reference temperature");
    }
    prm.leave_subsection();
  }

  void
  Material::declare_parameters(ParameterHandler  &prm,
                               const std::string &material_prefix,
                               const unsigned int id) const
  {
    prm.enter_subsection(material_prefix + " " +
                         Utilities::int_to_string(id, 1));
    {
      prm.declare_entry("density",
                        "1",
                        Patterns::Double(),
                        "Density for the fluid corresponding to Phase = " +
                          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "kinematic viscosity",
        "1",
        Patterns::Double(),
        "Kinematic viscosity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "specific heat",
        "1",
        Patterns::Double(),
        "Specific heat for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "thermal conductivity",
        "1",
        Patterns::Double(),
        "Thermal conductivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "thermal expansion",
        "1",
        Patterns::Double(),
        "Thermal expansion coefficient for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer diffusivity model",
        "constant",
        Patterns::Selection(
          "constant|immersed solid tanh|immersed solid gaussian"),
        "Model used for the calculation of the tracer diffusivity"
        "Choices are <constant|immersed solid tanh|immersed solid gaussian>.");

      prm.declare_entry(
        "tracer diffusivity",
        "0",
        Patterns::Double(),
        "Tracer diffusivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer reaction constant model",
        "constant",
        Patterns::Selection(
          "none|constant|immersed solid tanh|immersed solid gaussian"),
        "Model used for the calculation of the tracer reaction constant"
        "Choices are <none|constant|immersed solid tanh|immersed solid gaussian>.");

      prm.declare_entry(
        "tracer reaction constant",
        "0",
        Patterns::Double(),
        "Tracer reaction constant for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer reaction order",
        "1",
        Patterns::Double(),
        "Tracer reaction order for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer reaction threshold",
        "1e-8",
        Patterns::Double(),
        "Tracer reaction threshold for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1) + ". " +
          "Lower values are more realistic but lead to stiffer linear systems when n < 1. ");

      // Declaration of the immersed solids models parameters
      immersed_solid_tanh_parameters.declare_parameters(prm);
      immersed_solid_gaussian_parameters.declare_parameters(prm);

      prm.declare_entry(
        "rheological model",
        "newtonian",
        Patterns::Selection("newtonian|power-law|carreau|phase_change"),
        "Rheological model "
        "Choices are <newtonian|power-law|carreau|phase_change>.");

      non_newtonian_parameters.declare_parameters(prm);


      prm.declare_entry("density model",
                        "constant",
                        Patterns::Selection("constant|isothermal_ideal_gas"),
                        "Model used for the calculation of the density"
                        "Choices are <constant|isothermal_ideal_gas>.");

      isothermal_ideal_gas_density_parameters.declare_parameters(prm);

      prm.declare_entry("specific heat model",
                        "constant",
                        Patterns::Selection("constant|phase_change"),
                        "Model used for the calculation of the specific heat"
                        "Choices are <constant|phase_change>.");

      phase_change_parameters.declare_parameters(prm);

      prm.declare_entry(
        "thermal conductivity model",
        "constant",
        Patterns::Selection("constant|linear|phase_change"),
        "Model used for the calculation of the thermal conductivity"
        "Choices are <constant|linear|phase_change>.");

      prm.declare_entry(
        "thermal expansion model",
        "constant",
        Patterns::Selection("constant|phase_change"),
        "Model used for the calculation of the thermal expansion coefficient"
        "Choices are <constant|phase_change>.");

      prm.declare_entry("k_A0",
                        "0",
                        Patterns::Double(),
                        "k_A0 parameter for linear conductivity model");

      prm.declare_entry("k_A1",
                        "0",
                        Patterns::Double(),
                        "k_A1 parameter for linear conductivity model");
    }
    prm.leave_subsection();
  }

  void
  Material::parse_parameters(ParameterHandler                 &prm,
                             const std::string                &material_prefix,
                             const unsigned int                id,
                             const Parameters::Dimensionality &dimensions)
  {
    prm.enter_subsection(material_prefix + " " +
                         Utilities::int_to_string(id, 1));
    {
      // String that will be used to parse the models
      std::string op;

      //---------------------------------------------------
      // Density
      //---------------------------------------------------
      op = prm.get("density model");
      if (op == "constant")
        {
          density_model = DensityModel::constant;
        }
      else if (op == "isothermal_ideal_gas")
        {
          density_model = DensityModel::isothermal_ideal_gas;
        }
      density = prm.get_double("density");
      // Density is in M L^-3, rescale
      density *= dimensions.density_scaling;
      isothermal_ideal_gas_density_parameters.parse_parameters(prm, dimensions);

      //---------------------------------------------------
      // Kinematic viscosity and Rheology
      //---------------------------------------------------
      op = prm.get("rheological model");
      if (op == "power-law")
        {
          rheological_model = RheologicalModel::powerlaw;
        }
      else if (op == "carreau")
        {
          rheological_model = RheologicalModel::carreau;
        }
      else if (op == "newtonian")
        {
          rheological_model = RheologicalModel::newtonian;
        }
      else if (op == "phase_change")
        {
          rheological_model = RheologicalModel::phase_change;
        }

      kinematic_viscosity = prm.get_double("kinematic viscosity");
      // Kinematic viscosity is in L^2 T^-1, rescale
      kinematic_viscosity *= dimensions.viscosity_scaling;
      non_newtonian_parameters.parse_parameters(prm, dimensions);

      //--------------
      // Specific heat
      //--------------
      op = prm.get("specific heat model");
      if (op == "constant")
        specific_heat_model = SpecificHeatModel::constant;
      else if (op == "phase_change")
        specific_heat_model = SpecificHeatModel::phase_change;
      specific_heat = prm.get_double("specific heat");

      // specific heat is in L^2 T^-2 theta^-1
      specific_heat *= dimensions.specific_heat_scaling;


      //----------------------
      // Thermal conductivity
      //----------------------
      op = prm.get("thermal conductivity model");
      if (op == "constant")
        thermal_conductivity_model = ThermalConductivityModel::constant;
      else if (op == "linear")
        thermal_conductivity_model = ThermalConductivityModel::linear;
      else if (op == "phase_change")
        thermal_conductivity_model = ThermalConductivityModel::phase_change;

      thermal_conductivity = prm.get_double("thermal conductivity");
      // thermal conductivity is in M L T^-3 theta^-1
      thermal_conductivity *= dimensions.thermal_conductivity_scaling;


      // Linear conductivity model parameters
      k_A0 = prm.get_double("k_A0");
      // k_A0 is in M L T^-3 theta^-1
      k_A0 *= dimensions.thermal_conductivity_scaling;

      k_A1 = prm.get_double("k_A1");
      // k_A1 is in M L T^-3 theta ^-2
      k_A1 *= dimensions.thermal_conductivity_scaling * dimensions.temperature;


      //------------------
      // Thermal expansion
      //------------------
      op = prm.get("thermal expansion model");
      if (op == "constant")
        thermal_expansion_model = ThermalExpansionModel::constant;
      else if (op == "phase_change")
        thermal_expansion_model = ThermalExpansionModel::phase_change;

      thermal_expansion = prm.get_double("thermal expansion");
      // thermal expansion is in theta^-1
      thermal_expansion *= dimensions.temperature;

      //-------------------
      // Tracer diffusivity
      //-------------------
      op = prm.get("tracer diffusivity model");
      if (op == "immersed solid tanh")
        tracer_diffusivity_model =
          TracerDiffusivityModel::immersed_boundary_tanh;
      else if (op == "immersed solid gaussian")
        tracer_diffusivity_model =
          TracerDiffusivityModel::immersed_boundary_gaussian;
      else
        tracer_diffusivity_model = TracerDiffusivityModel::constant;
      tracer_diffusivity = prm.get_double("tracer diffusivity");
      // Diffusivity is in L^2 T^-1
      tracer_diffusivity *= dimensions.diffusivity_scaling;

      //-------------------
      // Tracer reaction constant
      //-------------------
      op = prm.get("tracer reaction constant model");
      if (op == "none")
        tracer_reaction_prefactor_model = TracerReactionPrefactorModel::none;
      else if (op == "immersed solid tanh")
        tracer_reaction_prefactor_model =
          TracerReactionPrefactorModel::immersed_boundary_tanh;
      else if (op == "immersed solid gaussian")
        tracer_reaction_prefactor_model =
          TracerReactionPrefactorModel::immersed_boundary_gaussian;
      else
        tracer_reaction_prefactor_model =
          TracerReactionPrefactorModel::constant;
      tracer_reaction_constant  = prm.get_double("tracer reaction constant");
      tracer_reaction_order     = prm.get_double("tracer reaction order");
      tracer_reaction_threshold = prm.get_double("tracer reaction threshold");

      // Parsing of the immersed solids models parameters
      immersed_solid_tanh_parameters.parse_parameters(prm, dimensions);
      immersed_solid_gaussian_parameters.parse_parameters(prm, dimensions);

      //--------------------------------
      // Phase change properties
      //--------------------------------
      phase_change_parameters.parse_parameters(prm, dimensions);
    }
    prm.leave_subsection();
  }

  void
  MaterialInteractions::declare_parameters(ParameterHandler  &prm,
                                           const unsigned int id) const
  {
    prm.enter_subsection("material interaction " +
                         Utilities::int_to_string(id, 1));
    {
      prm.declare_entry(
        "type",
        "fluid-fluid",
        Patterns::Selection("fluid-fluid|fluid-solid"),
        "Type of materials interacting. The choices are <fluid-fluid|fluid-solid>");

      // Fluid-fluid interactions
      prm.enter_subsection("fluid-fluid interaction");
      {
        prm.declare_entry(
          "first fluid id",
          "0",
          Patterns::Integer(),
          "ID of the first fluid interacting with the second fluid. This value should be lower than the second fluid's id.");
        prm.declare_entry(
          "second fluid id",
          "1",
          Patterns::Integer(),
          "ID of the second fluid interacting with the first fluid. This value should be greater than the first fluid's id.");

        // Surface tension interactions
        prm.declare_entry(
          "surface tension model",
          "constant",
          Patterns::Selection("constant|linear|phase change"),
          "Model used for the calculation of the surface tension coefficient\n"
          "The choices are <constant|linear|phase change>.");
        surface_tension_parameters.declare_parameters(prm);

        // Cahn-Hilliard mobility
        prm.declare_entry(
          "cahn hilliard mobility model",
          "constant",
          Patterns::Selection("constant|quartic"),
          "Model used for the calculation of the mobility in the Cahn-Hilliard equations"
          "\n"
          "The choices are <constant|quartic>.");
        mobility_cahn_hilliard_parameters.declare_parameters(prm);
      }
      prm.leave_subsection();

      // Fluid-solid interactions
      prm.enter_subsection("fluid-solid interaction");
      {
        prm.declare_entry("fluid id",
                          "0",
                          Patterns::Integer(),
                          "ID of the fluid interacting with the solid");
        prm.declare_entry("solid id",
                          "0",
                          Patterns::Integer(),
                          "ID of the solid interacting with the fluid");

        // Surface tension interactions
        prm.declare_entry(
          "surface tension model",
          "constant",
          Patterns::Selection("constant|linear|phase change"),
          "Model used for the calculation of the surface tension coefficient\n"
          "The choices are <constant|linear|phase change>.");
        surface_tension_parameters.declare_parameters(prm);

        // Cahn-Hilliard mobility
        prm.declare_entry(
          "cahn hilliard mobility model",
          "constant",
          Patterns::Selection("constant|quartic"),
          "Model used for the calculation of the mobility in the Cahn-Hilliard equations"
          "\n"
          "The choices are <constant|quartic>.");
        mobility_cahn_hilliard_parameters.declare_parameters(prm);
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  MaterialInteractions::parse_parameters(
    ParameterHandler                 &prm,
    unsigned int                      id,
    const Parameters::Dimensionality &dimensions)
  {
    prm.enter_subsection("material interaction " +
                         Utilities::int_to_string(id, 1));
    {
      std::string op;
      op = prm.get("type");
      if (op == "fluid-fluid")
        material_interaction_type = MaterialInteractionsType::fluid_fluid;
      else if (op == "fluid-solid")
        material_interaction_type = MaterialInteractionsType::fluid_solid;
      else
        throw(std::runtime_error(
          "Invalid material interaction type. The choices are <fluid-fluid|fluid-solid>."));

      if (material_interaction_type == MaterialInteractionsType::fluid_fluid)
        {
          prm.enter_subsection("fluid-fluid interaction");
          {
            std::pair<unsigned int, unsigned int> fluid_fluid_interaction;
            fluid_fluid_interaction.first  = prm.get_integer("first fluid id");
            fluid_fluid_interaction.second = prm.get_integer("second fluid id");
            AssertThrow(fluid_fluid_interaction.first <=
                          fluid_fluid_interaction.second,
                        OrderOfFluidIDsError(fluid_fluid_interaction.first,
                                             fluid_fluid_interaction.second));
            fluid_fluid_interaction_with_material_interaction_id.first =
              fluid_fluid_interaction;
            fluid_fluid_interaction_with_material_interaction_id.second = id;

            // Surface tension
            op = prm.get("surface tension model");
            if (op == "constant")
              {
                surface_tension_model = SurfaceTensionModel::constant;
              }
            else if (op == "linear")
              {
                surface_tension_model = SurfaceTensionModel::linear;
              }
            else if (op == "phase change")
              {
                surface_tension_model = SurfaceTensionModel::phase_change;
              }
            else
              throw(std::runtime_error(
                "Invalid surface tension model. The choices are <constant|linear|phase change>."));
            surface_tension_parameters.parse_parameters(prm, dimensions);
            // Cahn-Hilliard mobility
            op = prm.get("cahn hilliard mobility model");
            if (op == "constant")
              {
                mobility_cahn_hilliard_model =
                  MobilityCahnHilliardModel::constant;
              }
            else if (op == "quartic")
              {
                mobility_cahn_hilliard_model =
                  MobilityCahnHilliardModel::quartic;
              }
            else
              throw(std::runtime_error(
                "Invalid mobility model. The choices are <constant|quartic>."));

            mobility_cahn_hilliard_parameters.parse_parameters(prm, dimensions);
          }
          prm.leave_subsection();
        }
      else // Solid-fluid interactions
        {
          prm.enter_subsection("fluid-solid interaction");
          std::pair<unsigned int, unsigned int> fluid_solid_interaction;
          fluid_solid_interaction.first  = prm.get_integer("fluid id");
          fluid_solid_interaction.second = prm.get_integer("solid id");
          fluid_solid_interaction_with_material_interaction_id.first =
            fluid_solid_interaction;
          fluid_solid_interaction_with_material_interaction_id.second = id;

          // Surface tension
          op = prm.get("surface tension model");
          if (op == "constant")
            {
              surface_tension_model = SurfaceTensionModel::constant;
              surface_tension_parameters.parse_parameters(prm, dimensions);
            }
          else if (op == "linear")
            {
              surface_tension_model = SurfaceTensionModel::linear;
              surface_tension_parameters.parse_parameters(prm, dimensions);
            }
          else if (op == "phase change")
            {
              surface_tension_model = SurfaceTensionModel::phase_change;
              surface_tension_parameters.parse_parameters(prm, dimensions);
            }
          else
            AssertThrow(
              false,
              ExcMessage(
                "Invalid surface tension model. The choices are <constant|linear|phase change>."));

          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  void
  FEM::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("FEM");
    {
      prm.declare_entry("velocity order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order velocity");
      prm.declare_entry("pressure order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order pressure");
      prm.declare_entry("void fraction order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order void fraction");
      prm.declare_entry("temperature order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order temperature");
      prm.declare_entry("tracer order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order tracer");
      prm.declare_entry("VOF order",
                        "1",
                        Patterns::Integer(),
                        "interpolation order tracer");
      prm.declare_entry(
        "phase cahn hilliard order",
        "1",
        Patterns::Integer(),
        "interpolation order phase parameter in the Cahn-Hilliard equations");
      prm.declare_entry(
        "potential cahn hilliard order",
        "1",
        Patterns::Integer(),
        "interpolation order chemical potential in the Cahn-Hilliard equations");

      prm.declare_entry(
        "tracer uses dg",
        "false",
        Patterns::Bool(),
        "Switch tracer to Discontinuous Galerkin (DG) formulation");

      prm.declare_entry(
        "VOF uses dg",
        "false",
        Patterns::Bool(),
        "Switch VOF to Discontinuous Galerkin (DG) formulation");

      prm.declare_entry("enable bubble function velocity",
                        "false",
                        Patterns::Bool(),
                        "Enable bubble enrichment function for velocity field");

      prm.declare_entry("enable bubble function pressure",
                        "false",
                        Patterns::Bool(),
                        "Enable bubble enrichment function for pressure field");
    }
    prm.leave_subsection();
  }

  void
  FEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("FEM");
    {
      velocity_order            = prm.get_integer("velocity order");
      pressure_order            = prm.get_integer("pressure order");
      void_fraction_order       = prm.get_integer("void fraction order");
      temperature_order         = prm.get_integer("temperature order");
      tracer_order              = prm.get_integer("tracer order");
      tracer_uses_dg            = prm.get_bool("tracer uses dg");
      VOF_order                 = prm.get_integer("VOF order");
      VOF_uses_dg               = prm.get_bool("VOF uses dg");
      phase_cahn_hilliard_order = prm.get_integer("phase cahn hilliard order");
      potential_cahn_hilliard_order =
        prm.get_integer("potential cahn hilliard order");
      enable_bubble_function_velocity =
        prm.get_bool("enable bubble function velocity");
      enable_bubble_function_pressure =
        prm.get_bool("enable bubble function pressure");
    }
    prm.leave_subsection();
  }

  void
  Forces::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("forces");
    {
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the non-linear solver should be printed "
        "Choices are <quiet|verbose>.");
      prm.declare_entry("calculate force",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of force");
      prm.declare_entry("calculate torque",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of torque");
      prm.declare_entry("force name",
                        "force",
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("torque name",
                        "torque",
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("output precision",
                        "10",
                        Patterns::Integer(),
                        "Calculation frequency");
      prm.declare_entry("calculation frequency",
                        "1",
                        Patterns::Integer(),
                        "Calculation frequency");
      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output frequency");
    }
    prm.leave_subsection();
  }

  void
  Forces::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("forces");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;
      calculate_force       = prm.get_bool("calculate force");
      calculate_torque      = prm.get_bool("calculate torque");
      force_output_name     = prm.get("force name");
      torque_output_name    = prm.get("torque name");
      output_precision      = prm.get_integer("output precision");
      calculation_frequency = prm.get_integer("calculation frequency");
      output_frequency      = prm.get_integer("output frequency");
    }
    prm.leave_subsection();
  }

  void
  Laser_FreeSurfaceRadiation::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("free surface radiation");
    {
      prm.declare_entry(
        "enable",
        "false",
        Patterns::Bool(),
        "Enable radiation at the free surface (air/metal interface) <true|false>");
      prm.declare_entry("Stefan-Boltzmann constant",
                        "5.6703e-8",
                        Patterns::Double(),
                        "Stefan-Boltzmann constant");
      prm.declare_entry("emissivity",
                        "0.6",
                        Patterns::Double(),
                        "Emissivity of the free surface (air/metal interface)");
      prm.declare_entry(
        "Tinf",
        "0.0",
        Patterns::Double(),
        "Temperature (Double) of environment for radiation term at the free surface (air/metal interface)");
    }
    prm.leave_subsection();
  }

  void
  Laser_FreeSurfaceRadiation::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("free surface radiation");
    {
      enable_radiation          = prm.get_bool("enable");
      Stefan_Boltzmann_constant = prm.get_double("Stefan-Boltzmann constant");
      emissivity                = prm.get_double("emissivity");
      Tinf                      = prm.get_double("Tinf");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Laser<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("laser parameters");
    {
      prm.declare_entry("enable", "false", Patterns::Bool(), "Activate laser");
      prm.declare_entry(
        "type",
        "gaussian_heat_flux_vof_interface",
        Patterns::Selection(
          "exponential_decay|gaussian_heat_flux_vof_interface|uniform_heat_flux_vof_interface"),
        "Type of laser model used."
        "Choices are <exponential_decay|gaussian_heat_flux_vof_interface|uniform_heat_flux_vof_interface>.");
      prm.declare_entry("concentration factor",
                        "2.0",
                        Patterns::Double(),
                        "Concentration factor");
      prm.declare_entry("power", "0.0", Patterns::Double(), "Laser power");
      prm.declare_entry("absorptivity",
                        "0.5",
                        Patterns::Double(),
                        "Laser absorptivity");
      prm.declare_entry("penetration depth",
                        "1.0",
                        Patterns::Double(),
                        "Penetration depth");
      prm.declare_entry("beam radius",
                        "0.0",
                        Patterns::Double(),
                        "Laser beam radius");
      radiation.declare_parameters(prm);

      prm.enter_subsection("path");
      laser_scan_path = std::make_shared<Functions::ParsedFunction<dim>>(dim);
      laser_scan_path->declare_parameters(prm, dim);
      prm.leave_subsection();

      prm.declare_entry("start time",
                        "0.0",
                        Patterns::Double(),
                        "Start time of laser");
      prm.declare_entry("end time",
                        "1.0",
                        Patterns::Double(),
                        "End time of laser");

      prm.declare_entry("beam orientation",
                        "z-",
                        Patterns::Selection("x+|x-|y+|y-|z+|z-"),
                        "Laser beam orientation "
                        "Choices are <x+|x-|y+|y-|z+|z->.");

      prm.declare_entry(
        "beam rotation angle",
        "0.0",
        Patterns::Double(),
        "Angle of rotation in rad of the beam axis with respect to the axis defined by the beam orientation parameter");

      prm.declare_entry(
        "beam rotation axis",
        "0.0, 0.0, 1.0",
        Patterns::List(Patterns::Double()),
        "Axis around which the laser beam is rotated, only used in 3D");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Laser<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("laser parameters");
    {
      activate_laser                = prm.get_bool("enable");
      const std::string type_string = prm.get("type");
      if (type_string == "exponential_decay")
        laser_type = LaserType::exponential_decay;
      else if (type_string == "gaussian_heat_flux_vof_interface")
        laser_type = LaserType::gaussian_heat_flux_vof_interface;
      else
        laser_type = LaserType::uniform_heat_flux_vof_interface;
      concentration_factor = prm.get_double("concentration factor");
      laser_power          = prm.get_double("power");
      laser_absorptivity   = prm.get_double("absorptivity");
      penetration_depth    = prm.get_double("penetration depth");

      // Check if penetration depth is a strictly positive double.
      if (activate_laser && laser_type == LaserType::exponential_decay)
        {
          AssertThrow(laser_type == LaserType::exponential_decay &&
                        penetration_depth > 0.0,
                      ParameterStrictlyGreaterThanError("penetration depth",
                                                        penetration_depth,
                                                        0.0));
        }

      beam_radius = prm.get_double("beam radius");
      radiation.parse_parameters(prm);

      prm.enter_subsection("path");
      laser_scan_path->parse_parameters(prm);
      laser_scan_path->set_time(0);
      prm.leave_subsection();

      start_time = prm.get_double("start time");
      end_time   = prm.get_double("end time");

      std::string op;
      op = prm.get("beam orientation");
      if (op == "x+")
        {
          beam_direction                     = true;
          beam_orientation                   = BeamOrientation::x_plus;
          beam_orientation_coordinate        = 0;
          perpendicular_plane_coordinate_one = 1;
          beam_axis[0]                       = 1;
          beam_axis[1]                       = 0;
          if constexpr (dim == 3)
            {
              perpendicular_plane_coordinate_two = 2;
              beam_axis[2]                       = 0;
            }
        }
      else if (op == "x-")
        {
          beam_direction                     = false;
          beam_orientation                   = BeamOrientation::x_minus;
          beam_orientation_coordinate        = 0;
          perpendicular_plane_coordinate_one = 1;
          beam_axis[0]                       = -1;
          beam_axis[1]                       = 0;
          if constexpr (dim == 3)
            {
              perpendicular_plane_coordinate_two = 2;
              beam_axis[2]                       = 0;
            }
        }
      else if (op == "y+")
        {
          beam_direction                     = true;
          beam_orientation                   = BeamOrientation::y_plus;
          perpendicular_plane_coordinate_one = 0;
          beam_orientation_coordinate        = 1;
          beam_axis[0]                       = 0;
          beam_axis[1]                       = 1;
          if constexpr (dim == 3)
            {
              perpendicular_plane_coordinate_two = 2;
              beam_axis[2]                       = 0;
            }
        }
      else if (op == "y-")
        {
          beam_direction                     = false;
          beam_orientation                   = BeamOrientation::y_minus;
          perpendicular_plane_coordinate_one = 0;
          beam_orientation_coordinate        = 1;
          beam_axis[0]                       = 0;
          beam_axis[1]                       = -1;
          if constexpr (dim == 3)
            {
              perpendicular_plane_coordinate_two = 2;
              beam_axis[2]                       = 0;
            }
        }
      else if (op == "z+")
        {
          if constexpr (dim == 3)
            {
              beam_direction                     = true;
              beam_orientation                   = BeamOrientation::z_plus;
              perpendicular_plane_coordinate_one = 0;
              perpendicular_plane_coordinate_two = 1;
              beam_orientation_coordinate        = 2;
              beam_axis[0]                       = 0;
              beam_axis[1]                       = 0;
              beam_axis[2]                       = 1;
            }
          else if constexpr (dim == 2)
            Assert(dim == 2, TwoDimensionalLaserError(dim));
        }
      else if (op == "z-")
        {
          if constexpr (dim == 3)
            {
              beam_direction                     = false;
              beam_orientation                   = BeamOrientation::z_minus;
              perpendicular_plane_coordinate_one = 0;
              perpendicular_plane_coordinate_two = 1;
              beam_orientation_coordinate        = 2;
              beam_axis[0]                       = 0;
              beam_axis[1]                       = 0;
              beam_axis[2]                       = -1;
            }
          else if constexpr (dim == 2)
            Assert(dim == 2, TwoDimensionalLaserError(dim));
        }
      // Initial rotation axis and angle
      if constexpr (dim == 3)
        rotation_axis =
          value_string_to_tensor<dim>(prm.get("beam rotation axis"));

      rotation_angle = prm.get_double("beam rotation angle");

      if constexpr (dim == 2)
        rotation_matrix =
          Physics::Transformations::Rotations::rotation_matrix_2d(
            rotation_angle);
      if constexpr (dim == 3)
        rotation_matrix =
          Physics::Transformations::Rotations::rotation_matrix_3d(
            rotation_axis, rotation_angle);

      beam_axis = rotation_matrix * beam_axis;
    }
    prm.leave_subsection();
  }

  void
  PostProcessing::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("post-processing");
    {
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the post-processing values should be printed "
        "Choices are <quiet|verbose>.");

      prm.declare_entry(
        "calculate kinetic energy",
        "false",
        Patterns::Bool(),
        "Enable calculation of total kinetic energy. The total kinetic "
        "energy is calculated from the volumetric integral of the kinetic energy over the domain.");

      prm.declare_entry(
        "calculate enstrophy",
        "false",
        Patterns::Bool(),
        "Enable calculation of total enstrophy. The total enstrophy "
        "is calculated from the volumetric integral of the enstrophy over the domain.");

      prm.declare_entry(
        "calculate pressure power",
        "false",
        Patterns::Bool(),
        "Enable calculation of the pressure power. The pressure power "
        "is calculated from the volumetric integral of u.grad(p) over the domain.");

      prm.declare_entry(
        "calculate viscous dissipation",
        "false",
        Patterns::Bool(),
        "Enable calculation of the viscous dissipation. The viscous dissipation "
        "is calculated from the volumetric integral of grad(u).tau over the domain.");

      prm.declare_entry("calculate apparent viscosity",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of apparent viscosity");

      prm.declare_entry("calculate average velocities",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of average velocities.");

      prm.declare_entry(
        "calculate average temperature and heat flux",
        "false",
        Patterns::Bool(),
        "Enable calculation of time average temperature and time average heat flux");

      prm.declare_entry(
        "calculate pressure drop",
        "false",
        Patterns::Bool(),
        "Enable calculation of pressure drop between two boundaries.");

      prm.declare_entry("inlet boundary id",
                        "0",
                        Patterns::Integer(),
                        "Inlet boundary ID for pressure drop calculation");

      prm.declare_entry("outlet boundary id",
                        "1",
                        Patterns::Integer(),
                        "Outlet boundary ID for pressure drop calculation");

      prm.declare_entry("calculate flow rate",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of flow rate at boundaries.");

      prm.declare_entry(
        "calculate tracer flow rate",
        "false",
        Patterns::Bool(),
        "Enable calculation of tracer flow rate at boundaries.");

      prm.declare_entry(
        "initial time for average velocity",
        "0.0",
        Patterns::Double(),
        "Initial time to start calculations for average velocities");

      prm.declare_entry(
        "initial time for average temperature and heat flux",
        "0.0",
        Patterns::Double(),
        "Initial time to start calculations for average temperature");

      prm.declare_entry("kinetic energy name",
                        "kinetic_energy",
                        Patterns::FileName(),
                        "File output kinetic energy");

      prm.declare_entry("pressure drop name",
                        "pressure_drop",
                        Patterns::FileName(),
                        "File output pressure drop");

      prm.declare_entry("flow rate name",
                        "flow_rate",
                        Patterns::FileName(),
                        "File output volumetric flux");

      prm.declare_entry("tracer flow rate name",
                        "tracer_flow_rate",
                        Patterns::FileName(),
                        "Output file name for tracer flow rate");

      prm.declare_entry("enstrophy name",
                        "enstrophy",
                        Patterns::FileName(),
                        "File output enstrophy");

      prm.declare_entry("pressure power name",
                        "pressure_power",
                        Patterns::FileName(),
                        "File output pressure power");

      prm.declare_entry("viscous dissipation name",
                        "viscous_dissipation",
                        Patterns::FileName(),
                        "File output viscous dissipation");

      prm.declare_entry("apparent viscosity name",
                        "apparent_viscosity",
                        Patterns::FileName(),
                        "File output apparent viscosity");

      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output frequency");

      prm.declare_entry("calculate tracer statistics",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of tracer statistics.");

      prm.declare_entry("tracer statistics name",
                        "tracer_statistics",
                        Patterns::FileName(),
                        "File name output tracer statistics");

      prm.declare_entry(
        "calculate phase statistics",
        "false",
        Patterns::Bool(),
        "Enable calculation of phase statistics: maximum, minimum, average and integral over the domain (Cahn-Hilliard).");

      prm.declare_entry("phase statistics name",
                        "phase_statistics",
                        Patterns::FileName(),
                        "File name output phase statistics (Cahn-Hilliard)");

      prm.declare_entry("calculate temperature statistics",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of temperature statistics.");

      prm.declare_entry("temperature statistics name",
                        "temperature_statistics",
                        Patterns::FileName(),
                        "File name output temperature statistics");

      prm.declare_entry(
        "calculate liquid fraction",
        "false",
        Patterns::Bool(),
        "Enable calculation of the liquid fraction. The liquid fraction "
        "is calculated from the volume integral of the liquid fraction divided by the volume of the domain.");

      prm.declare_entry("liquid fraction name",
                        "liquid_fraction",
                        Patterns::FileName(),
                        "File name output liquid fraction");

      prm.declare_entry("calculate heat flux",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of heat flux.");

      prm.declare_entry("heat flux name",
                        "heat_flux",
                        Patterns::FileName(),
                        "File name output for the heat flux");

      prm.declare_entry("convective flux name",
                        "convective_flux",
                        Patterns::FileName(),
                        "File name output for the convective flux");

      prm.declare_entry("nitsche flux name",
                        "nitsche_heat_flux",
                        Patterns::FileName(),
                        "File name output for the convective flux");

      prm.declare_entry("postprocessed fluid",
                        "both",
                        Patterns::Selection("fluid 0|fluid 1|both"),
                        "Fluid domain used for thermal postprocesses "
                        "in the heat equation <fluid 0|fluid 1|both>");

      prm.declare_entry(
        "calculate barycenter",
        "false",
        Patterns::Bool(),
        "Enable calculation of the barycenter location and velocity of fluid 1 in VOF and Cahn-Hilliard simulations.");

      prm.declare_entry(
        "barycenter name",
        "barycenter_information",
        Patterns::FileName(),
        "Name of barycenter information output file in VOF or Cahn-Hilliard simulations");

      prm.declare_entry(
        "calculate mass conservation",
        "true",
        Patterns::Bool(),
        "Enable calculation of the mass and momentum of both fluids in VOF simulations.");

      prm.declare_entry(
        "mass conservation name",
        "mass_conservation_information",
        Patterns::FileName(),
        "Name of mass conservation output file in VOF simulations");

      prm.declare_entry(
        "calculate phase energy",
        "false",
        Patterns::Bool(),
        "Enable calculation of phase energies, including: total energy, bulk energy, and interface energy");

      prm.declare_entry(
        "phase energy name",
        "phase_energy",
        Patterns::FileName(),
        "Name of energy output file in Cahn-Hilliard simulations. The file is stored in the output folder specified in the simulation control subsection");

      prm.declare_entry(
        "calculate phase volumes",
        "false",
        Patterns::Bool(),
        "Enable calculation of total volume each phases in cfd-dem simulation, including: total volume of fluid, and total volume of particles");

      prm.declare_entry(
        "phase volumes name",
        "phase_volumes",
        Patterns::FileName(),
        "Name of phases volume output file in cfd-dem simulations. The file is stored in the output folder specified in the simulation control subsection");
    }
    prm.leave_subsection();
  }

  void
  PostProcessing::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("post-processing");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;

      calculate_kinetic_energy = prm.get_bool("calculate kinetic energy");
      calculate_enstrophy      = prm.get_bool("calculate enstrophy");
      calculate_pressure_power = prm.get_bool("calculate pressure power");
      calculate_viscous_dissipation =
        prm.get_bool("calculate viscous dissipation");
      calculate_apparent_viscosity =
        prm.get_bool("calculate apparent viscosity");
      calculate_average_velocities =
        prm.get_bool("calculate average velocities");
      calculate_average_temp_and_hf =
        prm.get_bool("calculate average temperature and heat flux");
      calculate_pressure_drop    = prm.get_bool("calculate pressure drop");
      inlet_boundary_id          = prm.get_integer("inlet boundary id");
      outlet_boundary_id         = prm.get_integer("outlet boundary id");
      calculate_flow_rate        = prm.get_bool("calculate flow rate");
      calculate_tracer_flow_rate = prm.get_bool("calculate tracer flow rate");
      initial_time_for_average_velocities =
        prm.get_double("initial time for average velocity");
      initial_time_for_average_temp_and_hf =
        prm.get_double("initial time for average temperature and heat flux");
      kinetic_energy_output_name      = prm.get("kinetic energy name");
      pressure_drop_output_name       = prm.get("pressure drop name");
      flow_rate_output_name           = prm.get("flow rate name");
      tracer_flow_rate_output_name    = prm.get("tracer flow rate name");
      enstrophy_output_name           = prm.get("enstrophy name");
      pressure_power_output_name      = prm.get("pressure power name");
      viscous_dissipation_output_name = prm.get("viscous dissipation name");
      apparent_viscosity_output_name  = prm.get("apparent viscosity name");
      output_frequency                = prm.get_integer("output frequency");
      calculate_tracer_statistics = prm.get_bool("calculate tracer statistics");
      tracer_output_name          = prm.get("tracer statistics name");
      calculate_phase_statistics  = prm.get_bool("calculate phase statistics");
      phase_output_name           = prm.get("phase statistics name");
      calculate_temperature_statistics =
        prm.get_bool("calculate temperature statistics");
      calculate_liquid_fraction   = prm.get_bool("calculate liquid fraction");
      liquid_fraction_output_name = prm.get("liquid fraction name");
      temperature_output_name     = prm.get("temperature statistics name");
      calculate_heat_flux         = prm.get_bool("calculate heat flux");
      heat_flux_output_name       = prm.get("heat flux name");
      calculate_barycenter        = prm.get_bool("calculate barycenter");
      barycenter_output_name      = prm.get("barycenter name");
      calculate_mass_conservation = prm.get_bool("calculate mass conservation");
      mass_conservation_output_name = prm.get("mass conservation name");
      calculate_phase_energy        = prm.get_bool("calculate phase energy");
      phase_energy_output_name      = prm.get("phase energy name");
      calculate_phase_volumes       = prm.get_bool("calculate phase volumes");
      phase_volumes_output_name     = prm.get("phase volumes name");

      // Viscous dissipative fluid
      const std::string op_fluid = prm.get("postprocessed fluid");
      if (op_fluid == "fluid 1")
        postprocessed_fluid = Parameters::FluidIndicator::fluid1;
      else if (op_fluid == "fluid 0")
        postprocessed_fluid = Parameters::FluidIndicator::fluid0;
      else if (op_fluid == "both")
        postprocessed_fluid = Parameters::FluidIndicator::both;
      else
        throw(
          std::runtime_error("Invalid postprocessed fluid. "
                             "Options are 'fluid 0', 'fluid 1' or 'both'."));
    }
    prm.leave_subsection();
  }

  void
  NonLinearSolver::declare_parameters(ParameterHandler  &prm,
                                      const std::string &physics_name)
  {
    prm.enter_subsection("non-linear solver");
    {
      prm.enter_subsection(physics_name);
      {
        prm.declare_entry(
          "verbosity",
          "verbose",
          Patterns::Selection("quiet|verbose"),
          "State whether the outputs from the non-linear solver should be printed. "
          "Choices are <quiet|verbose>.");

        prm.declare_entry(
          "solver",
          "newton",
          Patterns::Selection("newton|kinsol_newton|inexact_newton"),
          "Non-linear solver that will be used "
          "Choices are <newton|kinsol_newton|inexact_newton>."
          " The newton solver is a traditional newton solver with"
          "an analytical jacobian formulation. The jacobian matrix and the preconditioner"
          "are assembled every iteration. In the kinsol_newton method, the nonlinear solver"
          "Kinsol from the SUNDIALS library is used. This solver has an internal algorithm"
          "that decides whether to reassemble the Jacobian matrix or not.");

        prm.declare_entry(
          "kinsol strategy",
          "line_search",
          Patterns::Selection("normal_newton|line_search|fixed_point|picard"),
          "Strategy that will be used by the kinsol newton solver");

        prm.declare_entry("tolerance",
                          "1e-6",
                          Patterns::Double(),
                          "Newton solver tolerance");
        prm.declare_entry("max iterations",
                          "10",
                          Patterns::Integer(),
                          "Maximum number of Newton Iterations");

        prm.declare_entry(
          "step tolerance",
          "0.9",
          Patterns::Double(),
          "Newton solver relative tolerance between steps."
          " If a newton iteration leads to a residual > step tolerance"
          " * previous residual then the theta relaxation"
          " is applied until this criteria is satisfied");

        prm.declare_entry(
          "matrix tolerance",
          "0.1",
          Patterns::Double(),
          "This parameter controls the frequency at which the matrix is refreshed in the inexact Newton solvers"
          "If the residual after a newton step < previous residual * matrix tolerance, the matrix is not re-assembled");

        prm.declare_entry(
          "force rhs calculation",
          "false",
          Patterns::Bool(),
          "This is required if there is a fixed point component to the non-linear"
          "solver that is changed at the beginning of every newton iteration."
          "This is notably the case of the sharp edge method."
          "The default value of this parameter is false.");


        prm.declare_entry("residual precision",
                          "4",
                          Patterns::Integer(),
                          "Number of digits displayed when showing residuals");
        prm.declare_entry(
          "reuse matrix",
          "false",
          Patterns::Bool(),
          "Reuse the last jacobian matrix for the next non-linear problem solution");

        prm.declare_entry(
          "reuse preconditioner",
          "false",
          Patterns::Bool(),
          "Reuse the last preconditioner for the next non-linear problem solution");

        prm.declare_entry(
          "abort at convergence failure",
          "false",
          Patterns::Bool(),
          "Aborts Lethe by throwing an exception if non-linear solver convergence has failed");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  NonLinearSolver::parse_parameters(ParameterHandler  &prm,
                                    const std::string &physics_name)
  {
    prm.enter_subsection("non-linear solver");
    {
      prm.enter_subsection(physics_name);
      {
        const std::string op = prm.get("verbosity");
        if (op == "verbose")
          verbosity = Parameters::Verbosity::verbose;
        else if (op == "quiet")
          verbosity = Parameters::Verbosity::quiet;
        else
          throw(std::runtime_error("Invalid verbosity level"));

        const std::string str_solver = prm.get("solver");
        if (str_solver == "newton")
          solver = SolverType::newton;
        else if (str_solver == "kinsol_newton")
          solver = SolverType::kinsol_newton;
        else if (str_solver == "inexact_newton")
          solver = SolverType::inexact_newton;
        else
          throw(std::runtime_error("Invalid non-linear solver "));

        const std::string str_kinsol_strategy = prm.get("kinsol strategy");
        if (str_kinsol_strategy == "normal_newton")
          kinsol_strategy = KinsolStrategy::normal_newton;
        else if (str_kinsol_strategy == "line_search")
          kinsol_strategy = KinsolStrategy::line_search;
        else if (str_kinsol_strategy == "picard")
          kinsol_strategy = KinsolStrategy::picard;
        else
          throw(std::runtime_error(
            "Invalid strategy for kinsol non-linear solver "));

        tolerance             = prm.get_double("tolerance");
        step_tolerance        = prm.get_double("step tolerance");
        matrix_tolerance      = prm.get_double("matrix tolerance");
        max_iterations        = prm.get_integer("max iterations");
        display_precision     = prm.get_integer("residual precision");
        force_rhs_calculation = prm.get_bool("force rhs calculation");
        reuse_matrix          = prm.get_bool("reuse matrix");
        reuse_preconditioner  = prm.get_bool("reuse preconditioner");
        abort_at_convergence_failure =
          prm.get_bool("abort at convergence failure");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
  void
  Mesh::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh");
    {
      prm.declare_entry("type",
                        "dealii",
                        Patterns::Selection(
                          "gmsh|dealii|periodic_hills|cylinder"),
                        "Type of mesh "
                        "Choices are <gmsh|dealii|periodic_hills|cylinder>.");

      prm.declare_entry("file name",
                        "none",
                        Patterns::FileName(),
                        "GMSH file name");

      prm.declare_entry("initial refinement",
                        "0",
                        Patterns::Integer(),
                        "Initial refinement of the mesh");

      prm.declare_entry(
        "initial boundary refinement",
        "0",
        Patterns::Integer(),
        "Initial refinement of the mesh at the boundaries specified by the user");

      prm.declare_entry(
        "boundaries refined",
        "",
        Patterns::List(Patterns::Integer()),
        "Boundary ids of the boundaries to be initially refined");

      if (prm.get("type") == "periodic_hills")
        {
          prm.declare_entry("grid arguments", "1 ; 1 ; 1 ; 1 ; 1");
        }
      else
        {
          prm.declare_entry("grid type", "hyper_cube");
          prm.declare_entry("grid arguments", "-1 : 1 : false");
        }
      prm.declare_entry(
        "enable target size",
        "false",
        Patterns::Bool(),
        "Enable initial refinement until target size is reached.");

      prm.declare_entry(
        "simplex",
        "false",
        Patterns::Bool(),
        "Indicates that the mesh used is a mesh made of only simplex elements.");

      prm.declare_entry(
        "check diamond cells",
        "false",
        Patterns::Bool(),
        "Enables checking the input grid for diamond-shaped cells.");

      prm.declare_entry(
        "expand particle-wall contact search",
        "false",
        Patterns::Bool(),
        "Enables adding the boundary neighbor cells of boundary cells to the"
        "particle-wall contact search list. This feature should only be "
        "activated in geometries with concave boundaries. (For example, for "
        "particles flow inside a cylinder or sphere). In geometries with "
        "convex boundaries, this feature MUST NOT be activated");

      prm.declare_entry("target size",
                        "1",
                        Patterns::Double(),
                        "Target size of the initial refinement");


      prm.declare_entry("grid type", "hyper_cube");
      prm.declare_entry("grid arguments", "-1 : 1 : false");

      prm.declare_entry(
        "initial translation",
        "0, 0, 0",
        Patterns::List(Patterns::Double()),
        "Component of the desired translation of the mesh at initialization");

      prm.declare_entry(
        "initial rotation axis",
        "1, 0, 0",
        Patterns::List(Patterns::Double()),
        "Component of the desired rotation of the mesh at initialization");

      prm.declare_entry(
        "initial rotation angle",
        "0",
        Patterns::Double(),
        "Angle of rotation of the mesh at initialization around the axis in radian");

      prm.declare_entry("scale",
                        "1",
                        Patterns::Double(),
                        "Scaling factor used for the mesh.");
    }
    prm.leave_subsection();
  }

  void
  Mesh::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh");
    {
      {
        const std::string op = prm.get("type");
        if (op == "gmsh")
          type = Type::gmsh;
        else if (op == "dealii")
          type = Type::dealii;
        else if (op == "periodic_hills")
          type = Type::periodic_hills;
        else if (op == "cylinder")
          type = Type::cylinder;
        else
          throw std::logic_error(
            "Error, invalid mesh type. Choices are gmsh and dealii");
      }

      file_name = prm.get("file name");

      initial_refinement = prm.get_integer("initial refinement");
      initial_refinement_at_boundaries =
        prm.get_integer("initial boundary refinement");

      boundaries_to_refine =
        convert_string_to_vector<int>(prm, "boundaries refined");

      grid_type      = prm.get("grid type");
      grid_arguments = prm.get("grid arguments");

      refine_until_target_size = prm.get_bool("enable target size");
      simplex                  = prm.get_bool("simplex");
      check_for_diamond_cells  = prm.get_bool("check diamond cells");
      expand_particle_wall_contact_search =
        prm.get_bool("expand particle-wall contact search");
      target_size = prm.get_double("target size");

      // Initial translation
      translation = value_string_to_tensor<3>(prm.get("initial translation"));

      // Initial rotation axis and angle
      rotation_axis =
        value_string_to_tensor<3>(prm.get("initial rotation axis"));
      rotation_angle = prm.get_double("initial rotation angle");

      // Scaling factor
      scale = prm.get_double("scale");
    }
    prm.leave_subsection();
  }

  void
  MeshBoxRefinement::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("box refinement");
    {
      box_mesh = std::make_shared<Mesh>();
      box_mesh->declare_parameters(prm);

      prm.declare_entry("initial refinement",
                        "0",
                        Patterns::Integer(),
                        "Initial refinement of the principal mesh");
    }
    prm.leave_subsection();
  }

  void
  MeshBoxRefinement::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("box refinement");
    {
      box_mesh->parse_parameters(prm);

      initial_refinement = prm.get_integer("initial refinement");
    }
    prm.leave_subsection();
  }

  void
  LinearSolver::declare_parameters(ParameterHandler  &prm,
                                   const std::string &physics_name)
  {
    prm.enter_subsection("linear solver");
    {
      prm.enter_subsection(physics_name);
      {
        prm.declare_entry(
          "verbosity",
          "verbose",
          Patterns::Selection("quiet|verbose|extra verbose"),
          "State whether output from solver runs should be printed. "
          "Choices are <quiet|verbose|extra verbose>.");
        prm.declare_entry(
          "method",
          "gmres",
          Patterns::Selection("gmres|bicgstab|direct"),
          "The iterative solver for the linear system of equations. "
          "Choices are <gmres|bicgstab|direct>.");

        prm.declare_entry(
          "rescale residual",
          "false",
          Patterns::Bool(),
          "Rescale the residual by the square root of the volume of the triangulation");
        prm.declare_entry("relative residual",
                          "1e-3",
                          Patterns::Double(),
                          "Linear solver residual");
        prm.declare_entry("minimum residual",
                          "1e-12",
                          Patterns::Double(),
                          "Linear solver minimum residual");
        prm.declare_entry("max iters",
                          "1000",
                          Patterns::Integer(),
                          "Maximum solver iterations");

        prm.declare_entry("max krylov vectors",
                          "100",
                          Patterns::Integer(),
                          "Maximum number of krylov vectors for GMRES");

        prm.declare_entry(
          "enable hessians in jacobian",
          "true",
          Patterns::Bool(),
          "Turns off the terms involving the hessian in the Jacobian");

        prm.declare_entry(
          "enable hessians in residual",
          "true",
          Patterns::Bool(),
          "Turns off the terms involving the hessian in the rhs");

        prm.declare_entry("preconditioner",
                          "ilu",
                          Patterns::Selection("amg|ilu|lsmg|gcmg"),
                          "The preconditioner for the linear solver."
                          "Choices are <amg|ilu|lsmg|gcmg>.");


        prm.declare_entry("ilu preconditioner fill",
                          "0",
                          Patterns::Double(),
                          "Ilu preconditioner fill");

        prm.declare_entry("ilu preconditioner absolute tolerance",
                          "1e-12",
                          Patterns::Double(),
                          "Ilu preconditioner tolerance");

        prm.declare_entry("ilu preconditioner relative tolerance",
                          "1.00",
                          Patterns::Double(),
                          "Ilu relative tolerance");

        prm.declare_entry("amg preconditioner ilu fill",
                          "0",
                          Patterns::Double(),
                          "amg preconditioner ilu smoother fill");

        prm.declare_entry("amg preconditioner ilu absolute tolerance",
                          "1e-12",
                          Patterns::Double(),
                          "amg preconditioner ilu smoother absolute tolerance");

        prm.declare_entry("amg preconditioner ilu relative tolerance",
                          "1.00",
                          Patterns::Double(),
                          "amg preconditioner ilu smoother relative tolerance");

        prm.declare_entry("amg aggregation threshold",
                          "1e-14",
                          Patterns::Double(),
                          "amg aggregation threshold");
        prm.declare_entry("amg n cycles",
                          "1",
                          Patterns::Integer(),
                          "amg number of cycles");
        prm.declare_entry("amg w cycles",
                          "false",
                          Patterns::Bool(),
                          "amg w cycling. If this is set to true, W cycling is "
                          "used. Otherwise, V cycling is used.");
        prm.declare_entry("amg smoother sweeps",
                          "2",
                          Patterns::Integer(),
                          "amg smoother sweeps");
        prm.declare_entry("amg smoother overlap",
                          "1",
                          Patterns::Integer(),
                          "amg smoother overlap");
        prm.declare_entry(
          "force linear solver continuation",
          "false",
          Patterns::Bool(),
          "A boolean that will force the linear solver to continue even if it fails");

        prm.declare_entry("mg min level",
                          "-1",
                          Patterns::Integer(),
                          "mg min level");

        prm.declare_entry("mg level min cells",
                          "-1",
                          Patterns::Integer(),
                          "mg minimum number of cells for coarse level");

        prm.declare_entry("mg int level",
                          "-1",
                          Patterns::Integer(),
                          "mg int level");

        prm.declare_entry(
          "mg enable hessians in jacobian",
          "true",
          Patterns::Bool(),
          "Turns off the terms involving the hessian in the Jacobian of mg operators");

        prm.declare_entry("mg smoother iterations",
                          "10",
                          Patterns::Integer(),
                          "mg smoother iterations for lsmg or gcmg");

        prm.declare_entry("mg smoother relaxation",
                          "0.5",
                          Patterns::Double(),
                          "mg smoother relaxation for lsmg or gcmg");

        prm.declare_entry("mg smoother preconditioner type",
                          "inverse diagonal",
                          Patterns::Selection(
                            "inverse diagonal|additive schwarz method"),
                          "Preconditioner of smoother");

        prm.declare_entry("mg smoother eig estimation",
                          "true",
                          Patterns::Bool(),
                          "estimate eigenvalues for relaxation parameter");

        prm.declare_entry("eig estimation smoothing range",
                          "10",
                          Patterns::Double(),
                          "sets range between largest and smallest eig");

        prm.declare_entry("eig estimation cg n iterations",
                          "10",
                          Patterns::Integer(),
                          "cg iterations performed to find eigenvalue");

        prm.declare_entry("eig estimation verbosity",
                          "verbose",
                          Patterns::Selection("quiet|verbose"),
                          "State whether MG should print max and min eigenvalue"
                          "Choices are <quiet|verbose>.");

        prm.declare_entry("mg coarse grid solver",
                          "direct",
                          Patterns::Selection("gmres|amg|ilu|direct"),
                          "The coarse grid solver for lsmg or gcmg"
                          "Choices are <gmres|amg|ilu|direct>.");

        prm.declare_entry(
          "mg coarse grid use fe q iso q1",
          "false",
          Patterns::Bool(),
          "use elements with linear interpolation for coarse grid");

        prm.declare_entry("mg coarsening type",
                          "h",
                          Patterns::Selection("h|p|hp|ph"),
                          "mg coarsening type for gcmg");

        prm.declare_entry("mg p coarsening type",
                          "decrease by one",
                          Patterns::Selection(
                            "decrease by one|bisect|go to one"),
                          "mg p coarsening type for gcmg");

        prm.declare_entry("mg gmres max iterations",
                          "2000",
                          Patterns::Integer(),
                          "mg gmres iterations for lsmg or gcmg");

        prm.declare_entry("mg gmres tolerance",
                          "1e-14",
                          Patterns::Double(),
                          "mg gmres tolerance n for lsmg or gcmg");

        prm.declare_entry("mg gmres reduce",
                          "1e-4",
                          Patterns::Double(),
                          "mg gmres reduce for lsmg or gcmg");

        prm.declare_entry("mg gmres max krylov vectors",
                          "30",
                          Patterns::Integer(),
                          "mg gmres max krylov vectors for lsmg or gcmg");

        prm.declare_entry("mg gmres preconditioner",
                          "amg",
                          Patterns::Selection("amg|ilu"),
                          "The preconditioner for the mg gmres solver"
                          "Choices are <amg|ilu>.");

        prm.declare_entry("mg amg use default parameters",
                          "false",
                          Patterns::Bool(),
                          "Use default parameters for Trilinos AMG");

        prm.declare_entry(
          "mg verbosity",
          "verbose",
          Patterns::Selection("quiet|verbose|extra verbose"),
          "State whether LSMG or GCMG should print information about levels "
          "Choices are <quiet|verbose|extra verbose>.");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
  void
  LinearSolver::parse_parameters(ParameterHandler  &prm,
                                 const std::string &physics_name)
  {
    prm.enter_subsection("linear solver");
    {
      prm.enter_subsection(physics_name);
      {
        const std::string sv = prm.get("method");
        if (sv == "gmres")
          solver = SolverType::gmres;
        else if (sv == "bicgstab")
          solver = SolverType::bicgstab;
        else if (sv == "direct")
          solver = SolverType::direct;
        else
          throw std::logic_error(
            "Error, invalid iterative solver type. Choices are amg, gmres, bicgstab or direct");

        const std::string op = prm.get("verbosity");
        if (op == "verbose")
          verbosity = Parameters::Verbosity::verbose;
        else if (op == "quiet")
          verbosity = Parameters::Verbosity::quiet;
        else if (op == "extra verbose")
          verbosity = Parameters::Verbosity::extra_verbose;
        else
          throw(
            std::runtime_error("Unknown verbosity mode for the linear solver"));

        rescale_residual_by_volume = prm.get_bool("rescale residual");
        relative_residual          = prm.get_double("relative residual");
        minimum_residual           = prm.get_double("minimum residual");
        max_iterations             = prm.get_integer("max iters");
        max_krylov_vectors         = prm.get_integer("max krylov vectors");
        enable_hessians_jacobian = prm.get_bool("enable hessians in jacobian");
        enable_hessians_residual = prm.get_bool("enable hessians in residual");

        Assert(enable_hessians_residual || !enable_hessians_jacobian,
               ExcNotImplemented());

        const std::string precond = prm.get("preconditioner");
        if (precond == "amg")
          preconditioner = PreconditionerType::amg;
        else if (precond == "ilu")
          preconditioner = PreconditionerType::ilu;
        else if (precond == "lsmg")
          preconditioner = PreconditionerType::lsmg;
        else if (precond == "gcmg")
          preconditioner = PreconditionerType::gcmg;
        else
          throw std::logic_error(
            "Error, invalid preconditioner type. Choices are amg, ilu, lsmg or gcmg.");


        ilu_precond_fill = prm.get_integer("ilu preconditioner fill");
        ilu_precond_atol =
          prm.get_double("ilu preconditioner absolute tolerance");
        ilu_precond_rtol =
          prm.get_double("ilu preconditioner relative tolerance");

        amg_precond_ilu_fill = prm.get_integer("amg preconditioner ilu fill");
        amg_precond_ilu_atol =
          prm.get_double("amg preconditioner ilu absolute tolerance");
        amg_precond_ilu_rtol =
          prm.get_double("amg preconditioner ilu relative tolerance");
        amg_aggregation_threshold = prm.get_double("amg aggregation threshold");
        amg_n_cycles              = prm.get_integer("amg n cycles");
        amg_w_cycles              = prm.get_bool("amg w cycles");
        amg_smoother_sweeps       = prm.get_integer("amg smoother sweeps");
        amg_smoother_overlap      = prm.get_integer("amg smoother overlap");

        force_linear_solver_continuation =
          prm.get_bool("force linear solver continuation");

        mg_min_level       = prm.get_integer("mg min level");
        mg_level_min_cells = prm.get_integer("mg level min cells");
        mg_int_level       = prm.get_integer("mg int level");
        mg_enable_hessians_jacobian =
          prm.get_bool("mg enable hessians in jacobian");
        Assert(enable_hessians_jacobian || !mg_enable_hessians_jacobian,
               ExcNotImplemented());

        mg_smoother_iterations = prm.get_integer("mg smoother iterations");
        mg_smoother_relaxation = prm.get_double("mg smoother relaxation");

        const std::string mg_smoother_preconditioner_type_str =
          prm.get("mg smoother preconditioner type");
        if (mg_smoother_preconditioner_type_str == "inverse diagonal")
          this->mg_smoother_preconditioner_type =
            MultigridSmootherPreconditionerType::InverseDiagonal;
        else if (mg_smoother_preconditioner_type_str ==
                 "additive schwarz method")
          this->mg_smoother_preconditioner_type =
            MultigridSmootherPreconditionerType::AdditiveSchwarzMethod;
        else
          AssertThrow(false, ExcNotImplemented());

        mg_smoother_eig_estimation = prm.get_bool("mg smoother eig estimation");
        eig_estimation_smoothing_range =
          prm.get_double("eig estimation smoothing range");
        eig_estimation_cg_n_iterations =
          prm.get_integer("eig estimation cg n iterations");

        const std::string eig_estimation_v =
          prm.get("eig estimation verbosity");
        if (eig_estimation_v == "verbose")
          eig_estimation_verbose = Parameters::Verbosity::verbose;
        else if (eig_estimation_v == "quiet")
          eig_estimation_verbose = Parameters::Verbosity::quiet;
        else
          throw(std::runtime_error(
            "Unknown verbosity mode for the eigenvalue estimation"));

        const std::string cg_solver = prm.get("mg coarse grid solver");
        if (cg_solver == "gmres")
          mg_coarse_grid_solver = CoarseGridSolverType::gmres;
        else if (cg_solver == "amg")
          mg_coarse_grid_solver = CoarseGridSolverType::amg;
        else if (cg_solver == "ilu")
          mg_coarse_grid_solver = CoarseGridSolverType::ilu;
        else if (cg_solver == "direct")
          mg_coarse_grid_solver = CoarseGridSolverType::direct;
        else
          throw std::logic_error(
            "Error, invalid coarse grid solver type. Choices are gmres, amg, ilu or direct.");

        mg_use_fe_q_iso_q1 = prm.get_bool("mg coarse grid use fe q iso q1");

        const std::string mg_coarsening_type_str =
          prm.get("mg coarsening type");
        if (mg_coarsening_type_str == "h")
          this->mg_coarsening_type = MultigridCoarseningSequenceType::h;
        else if (mg_coarsening_type_str == "p")
          this->mg_coarsening_type = MultigridCoarseningSequenceType::p;
        else if (mg_coarsening_type_str == "hp")
          this->mg_coarsening_type = MultigridCoarseningSequenceType::hp;
        else if (mg_coarsening_type_str == "ph")
          this->mg_coarsening_type = MultigridCoarseningSequenceType::ph;
        else
          AssertThrow(false, ExcNotImplemented());

        const std::string mg_p_coarsening_type_str =
          prm.get("mg p coarsening type");
        if (mg_p_coarsening_type_str == "decrease by one")
          this->mg_p_coarsening_type = MGTransferGlobalCoarseningTools::
            PolynomialCoarseningSequenceType::decrease_by_one;
        else if (mg_p_coarsening_type_str == "bisect")
          this->mg_p_coarsening_type = MGTransferGlobalCoarseningTools::
            PolynomialCoarseningSequenceType::bisect;
        else if (mg_p_coarsening_type_str == "go to one")
          this->mg_p_coarsening_type = MGTransferGlobalCoarseningTools::
            PolynomialCoarseningSequenceType::go_to_one;
        else
          AssertThrow(false, ExcNotImplemented());

        AssertThrow((!mg_use_fe_q_iso_q1) ||
                      (this->mg_coarsening_type ==
                       MultigridCoarseningSequenceType::h),
                    ExcNotImplemented());

        AssertThrow(
          (preconditioner != PreconditionerType::lsmg) ||
            (this->mg_coarsening_type == MultigridCoarseningSequenceType::h ||
             this->mg_coarsening_type == MultigridCoarseningSequenceType::hp),
          ExcNotImplemented());

        mg_gmres_max_iterations = prm.get_integer("mg gmres max iterations");
        mg_gmres_tolerance      = prm.get_double("mg gmres tolerance");
        mg_gmres_reduce         = prm.get_double("mg gmres reduce");
        mg_gmres_max_krylov_vectors =
          prm.get_integer("mg gmres max krylov vectors");

        const std::string cg_precond = prm.get("mg gmres preconditioner");
        if (cg_precond == "amg")
          mg_gmres_preconditioner = PreconditionerType::amg;
        else if (cg_precond == "ilu")
          mg_gmres_preconditioner = PreconditionerType::ilu;
        else
          throw std::logic_error(
            "Error, invalid preconditioner type for mg gmres solver. Choices are amg or ilu.");

        mg_amg_use_default_parameters =
          prm.get_bool("mg amg use default parameters");

        const std::string mg_op = prm.get("mg verbosity");
        if (mg_op == "verbose")
          mg_verbosity = Parameters::Verbosity::verbose;
        else if (mg_op == "extra verbose")
          mg_verbosity = Parameters::Verbosity::extra_verbose;
        else if (mg_op == "quiet")
          mg_verbosity = Parameters::Verbosity::quiet;
        else
          throw(std::runtime_error(
            "Unknown verbosity mode for the LSMG or GCMG preconditioners"));
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  void
  MeshAdaptation::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh adaptation");
    {
      prm.declare_entry("initial refinement steps",
                        "0",
                        Patterns::Integer(),
                        "Number of pre-solve adaptive mesh refinement steps");

      prm.declare_entry("type",
                        "none",
                        Patterns::Selection("none|uniform|kelly"),
                        "Type of mesh adaptation"
                        "Choices are <none|uniform|kelly>.");

      prm.declare_entry(
        "fraction refinement",
        "0.1",
        Patterns::List(Patterns::Double()),
        "Fraction of refined elements"
        "For multi-variables refinement, separate the different fractions with a comma "
        "(ex/ 'set fraction refinement = 0.1,0.1')");

      prm.declare_entry(
        "fraction coarsening",
        "0.05",
        Patterns::List(Patterns::Double()),
        "Fraction of coarsened elements"
        "For multi-variables refinement, separate the different fractions with a comma "
        "(ex/ 'set fraction coarsening = 0.05,0.05')");

      prm.declare_entry(
        "variable",
        "velocity",
        Patterns::List(Patterns::Selection(
          "velocity|pressure|phase|temperature|phase_cahn_hilliard|chemical_potential_cahn_hilliard|tracer")),
        "Variable(s) for kelly estimation"
        "Choices are <velocity|pressure|phase|temperature|phase_cahn_hilliard|chemical_potential_cahn_hilliard|tracer>."
        "For multi-variables refinement, separate the different variables with a comma "
        "(ex/ 'set variables = velocity,temperature')");

      prm.declare_entry(
        "fraction type",
        "number",
        Patterns::Selection("number|fraction"),
        "How the fraction of refinement/coarsening are interpreted"
        "Choices are <number|fraction>.");
      prm.declare_entry("max number elements",
                        "100000000",
                        Patterns::Integer(),
                        "Maximum number of elements");
      prm.declare_entry("max refinement level",
                        "10",
                        Patterns::Integer(),
                        "Maximum refinement level");
      prm.declare_entry("min refinement level",
                        "0",
                        Patterns::Integer(),
                        "Minimum refinement level");
      prm.declare_entry("frequency",
                        "1",
                        Patterns::Integer(),
                        "Frequency of the mesh refinement");
      prm.declare_entry(
        "mesh refinement controller",
        "false",
        Patterns::Bool(),
        "Fraction of refined elements"
        "Enable a controller that will target a specific number of elements in the mesh equal to the maximum number of elements");
      prm.declare_entry("fix boundary refinement",
                        "false",
                        Patterns::Bool(),
                        "Enable fix boundary refinement");
      prm.declare_entry("boundaries fixed",
                        "",
                        Patterns::List(Patterns::Integer()),
                        "Boundary ids of the boundaries to be fixed");
    }
    prm.leave_subsection();
  }

  void
  MeshAdaptation::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mesh adaptation");
    {
      initial_refinement = prm.get_integer("initial refinement steps");

      const std::string op = prm.get("type");
      if (op == "none")
        type = Type::none;
      if (op == "uniform")
        type = Type::uniform;
      if (op == "kelly")
        type = Type::kelly;

      // Getting multivariables refinement parameters
      const std::string        var_op   = prm.get("variable");
      std::vector<std::string> var_vec  = Utilities::split_string_list(var_op);
      const std::string        coars_op = prm.get("fraction coarsening");
      std::vector<std::string> coars_vec =
        Utilities::split_string_list(coars_op);
      const std::string        refin_op = prm.get("fraction refinement");
      std::vector<std::string> refin_vec =
        Utilities::split_string_list(refin_op);

      // Checking that the sizes are coherent
      Assert(coars_vec.size() == var_vec.size(),
             MultipleAdaptationSizeError("fraction coarsening",
                                         coars_vec.size(),
                                         var_vec.size()));
      Assert(refin_vec.size() == var_vec.size(),
             MultipleAdaptationSizeError("fraction refinement",
                                         refin_vec.size(),
                                         var_vec.size()));

      // Create map of refinement variables
      for (std::vector<int>::size_type i = 0; i != var_vec.size(); ++i)
        {
          // Parsing variable for this index
          if (var_vec[i] == "velocity")
            vars = Variable::velocity;
          else if (var_vec[i] == "pressure")
            vars = Variable::pressure;
          else if (var_vec[i] == "phase")
            vars = Variable::phase;
          else if (var_vec[i] == "temperature")
            vars = Variable::temperature;
          else if (var_vec[i] == "phase_cahn_hilliard")
            vars = Variable::phase_cahn_hilliard;
          else if (var_vec[i] == "chemical_potential_cahn_hilliard")
            vars = Variable::chemical_potential_cahn_hilliard;
          else if (var_vec[i] == "tracer")
            vars = Variable::tracer;
          else
            throw std::logic_error(
              "Error, invalid mesh adaptation variable. Choices are velocity, pressure, phase, temperature, phase_cahn_hilliard, chemical_potential_cahn_hilliard or tracer");

          var_adaptation_param.coarsening_fraction = std::stod(coars_vec[i]);
          var_adaptation_param.refinement_fraction = std::stod(refin_vec[i]);

          // defining adaptation map for this variable
          variables[vars] = var_adaptation_param;
        }

      const std::string fop = prm.get("fraction type");
      if (fop == "number")
        fractionType = FractionType::number;
      if (fop == "fraction")
        fractionType = FractionType::fraction;
      maximum_number_elements      = prm.get_integer("max number elements");
      maximum_refinement_level     = prm.get_integer("max refinement level");
      minimum_refinement_level     = prm.get_integer("min refinement level");
      frequency                    = prm.get_integer("frequency");
      refinement_at_frequency      = frequency != 0;
      mesh_controller_is_enabled   = prm.get_bool("mesh refinement controller");
      is_boundary_refinement_fixed = prm.get_bool("fix boundary refinement");
      boundaries_to_fix =
        convert_string_to_vector<int>(prm, "boundaries fixed");
    }
    prm.leave_subsection();
  }

  void
  Testing::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("test");
    {
      prm.declare_entry(
        "enable",
        "false",
        Patterns::Bool(),
        "Enable testing mode of a solver. Some solvers have a specific"
        "testing mode which enables the output of debug variables. This"
        "testing mode is generally used only for the automatic testing bench using ctest.");
      prm.declare_entry(
        "type",
        "particles",
        Patterns::Selection("particles|mobility_status|subdomain"),
        "Output type for testing mode. Currently, particles type will output "
        "each particle with some information and mobility_status or subdomain output results "
        "in deal.II format.");
    }
    prm.leave_subsection();
  }

  void
  Testing::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("test");
    {
      enabled = prm.get_bool("enable");
      if (enabled)
        {
          const std::string op = prm.get("type");
          if (op == "particles")
            test_type = TestType::particles;
          else if (op == "mobility_status")
            test_type = TestType::mobility_status;
          else if (op == "subdomain")
            test_type = TestType::subdomain;
          else
            throw std::logic_error(
              "Error, invalid testing type. Current choices are particles, "
              "mobility_status or subdomain");
        }
    }
    prm.leave_subsection();
  }

  void
  Restart::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("restart");
    {
      prm.declare_entry("filename",
                        "restart",
                        Patterns::FileName(),
                        "Prefix for the filename of checkpoints");
      prm.declare_entry("restart",
                        "false",
                        Patterns::Bool(),
                        "Frequency for checkpointing");
      prm.declare_entry(
        "checkpoint",
        "false",
        Patterns::Bool(),
        "Enable checkpointing. Checkpointing creates a restart"
        "point from which the simulation can be restarted from.");

      prm.declare_entry("frequency",
                        "1",
                        Patterns::Integer(),
                        "Frequency for checkpointing");
    }
    prm.leave_subsection();
  }

  void
  Restart::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("restart");
    {
      filename   = prm.get("filename");
      checkpoint = prm.get_bool("checkpoint");
      restart    = prm.get_bool("restart");
      frequency  = prm.get_integer("frequency");
    }
    prm.leave_subsection();
  }

  void
  VelocitySource::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity source");
    {
      prm.declare_entry(
        "rotating frame type",
        "none",
        Patterns::Selection("none|srf"),
        "Rotating frame velocity-dependent source terms"
        "Choices are <none|srf>. The srf stands"
        "for single rotating frame and adds"
        "the coriolis and the centrifugal force to the Navier-Stokes equations");

      prm.declare_entry(
        "omega_x",
        "0 ",
        Patterns::Double(),
        "X component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "omega_y",
        "0 ",
        Patterns::Double(),
        "Y component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "omega_z",
        "0 ",
        Patterns::Double(),
        "Z component of the angular velocity vector of the frame of reference");

      prm.declare_entry("Darcy type",
                        "none",
                        Patterns::Selection("none|phase_change"),
                        "Darcy-like permeability term"
                        "Choices are <none|phase_change>.");
    }
    prm.leave_subsection();
  }

  void
  VelocitySource::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity source");
    {
      const std::string op = prm.get("rotating frame type");
      if (op == "none")
        rotating_frame_type = RotatingFrameType::none;
      else if (op == "srf")
        rotating_frame_type = RotatingFrameType::srf;
      else
        throw std::logic_error("Error, invalid velocity source type");

      const std::string darcy = prm.get("Darcy type");
      if (darcy == "none")
        darcy_type = DarcySourceType::none;
      else if (darcy == "phase_change")
        darcy_type = DarcySourceType::phase_change;
      else
        throw std::logic_error("Error, invalid Darcy source type");

      omega_x = prm.get_double("omega_x");
      omega_y = prm.get_double("omega_y");
      omega_z = prm.get_double("omega_z");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  IBParticles<dim>::declare_default_entry(ParameterHandler &prm,
                                          unsigned int      index)
  {
    prm.declare_entry(
      "integrate motion",
      "false",
      Patterns::Bool(),
      "Bool to define if the particle trajectory is integrated meaning its velocity and position will be updated at each time step according to the hydrodynamic force applied to it");
    prm.declare_entry(
      "mesh-based precalculations",
      "true",
      Patterns::Bool(),
      "Bool to define if precalculations should be performed between refinements. Precalculations can introduce shape deformation when the type is RBF and some nodes are located outside the background mesh.");

    prm.enter_subsection("position");
    particles[index].f_position =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    particles[index].f_position->declare_parameters(prm, dim);
    prm.leave_subsection();

    prm.enter_subsection("orientation");
    particles[index].f_orientation =
      std::make_shared<Functions::ParsedFunction<dim>>(3);
    particles[index].f_orientation->declare_parameters(prm, 3);
    prm.leave_subsection();

    prm.enter_subsection("velocity");
    particles[index].f_velocity =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    particles[index].f_velocity->declare_parameters(prm, dim);
    prm.leave_subsection();

    prm.enter_subsection("omega");
    particles[index].f_omega =
      std::make_shared<Functions::ParsedFunction<dim>>(3);
    particles[index].f_omega->declare_parameters(prm, 3);
    prm.leave_subsection();

    prm.declare_entry(
      "type",
      "sphere",
      Patterns::Selection(
        "sphere|hyper rectangle|ellipsoid|torus|cone|cylinder|cylindrical tube|cylindrical helix|cut hollow sphere|death star|superquadric|rbf|opencascade|plane|composite"),
      "The type of shape considered."
      "Choices are <sphere|hyper rectangle|ellipsoid|torus|cone|cylinder|cylindrical tube|cylindrical helix|cut hollow sphere|death star|superquadric|rbf|opencascade|composite>."
      "The parameter for a sphere is: radius. "
      "The parameters for a hyper rectangle are, in order: x half length,"
      "y half length, z half length."
      "The parameters for an ellipsoid are, in order: x radius,"
      "y radius, z radius. "
      "The parameters for a torus are, in order: torus radius,"
      "torus thickness radius. "
      "The parameters for a cone are, in order: tan(base angle),"
      " height. "
      "The parameters for a cylinder are, in order: radius, half-length."
      "It is aligned to the z axis by default."
      "The parameters for a cylindrical tube are, in order: inside radius,"
      "outside radius, half-length. It is aligned to the z axis by default."
      "The parameters for a cylindrical helix are, in order: helix radius,"
      "tube radius, helix total height, pitch (height between two consecutive "
      "loops)."
      "The parameters for a cut hollow sphere are, in order: sphere radius,"
      "cut thickness, wall thickness. "
      "The parameters for a death star are, in order: sphere radius,"
      "smaller sphere radius, distance between centers."
      "The parameters for a superquadric are, in order: "
      "a, b, c, r, s, t, epsilon. "
      "The first three are half-lengths in x, y, and z."
      "The next three are the blockiness in the x, y, and z directions."
      "The last is the tolerance of the found surface."
      "The parameter for an rbf is the file name."
      "The parameter for a composite is the file name.");

    prm.declare_entry("shape arguments",
                      "1",
                      Patterns::Anything(),
                      "Arguments defining the geometry");

    prm.declare_entry(
      "layer thickening",
      "0",
      Patterns::Double(),
      "Thickness (positive or negative) of uniform additional layer of solid on particle."
      "A negative value will decrease the particle's thickness by subtracting a layer of specified width.");

    prm.declare_entry(
      "pressure location",
      "0; 0; 0",
      Patterns::Anything(),
      "position relative to the center of the particle for the location of the point where the pressure is imposed inside the particle");



    prm.enter_subsection("physical properties");
    {
      prm.declare_entry("density",
                        "1",
                        Patterns::Double(),
                        "density of the particle ");
      prm.declare_entry(
        "volume",
        "0",
        Patterns::Double(),
        "The volume occupied by the particle. If it is left empty, the volume is automatically calculated if possible otherwise the volume of a sphere is used instead");
      prm.declare_entry(
        "inertia",
        "1 ;0 ;0 ;0 ;1 ;0 ;0 ;0 ;1",
        Patterns::Anything(),
        "Moments of inertia of the particle in the reference frame of the fluid. The entry sequence corresponds to : I_xx ;I_xy ;I_xz ;I_yx ;I_yy ;I_yz ;I_zx ;I_zy ;I_zz");
      prm.declare_entry("youngs modulus",
                        "100000000",
                        Patterns::Double(),
                        "The Young's modulus of particle in case of contact");
      prm.declare_entry("poisson ratio",
                        "0.3",
                        Patterns::Double(),
                        "The poisson ratio of particle in case of contact");
      prm.declare_entry(
        "restitution coefficient",
        "1",
        Patterns::Double(),
        "The restitution coefficient of particle in case of contact");
      prm.declare_entry(
        "friction coefficient",
        "0",
        Patterns::Double(),
        "The friction coefficient of particle in case of contact");
      prm.declare_entry(
        "rolling friction coefficient",
        "0",
        Patterns::Double(),
        "The rolling friction coefficient of particle in case of contact");
      prm.leave_subsection();
    }
  }

  template <int dim>
  void
  IBParticles<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("particles");
    {
      prm.declare_entry(
        "number of particles",
        "1",
        Patterns::Integer(),
        "The number of particles represented by IB. The maximal number of particles is equal to 10 when defined individually. If particles are loaded from a file, this parameter is overridden, and there is no limit to the number of particles.");

      prm.declare_entry(
        "assemble Navier-Stokes inside particles",
        "false",
        Patterns::Bool(),
        "Bool to define if you assemble the inside of particles with the NS equation.");

      prm.enter_subsection("extrapolation function");
      {
        prm.declare_entry(
          "stencil order",
          "2",
          Patterns::Integer(),
          "The polynomial order used in the extrapolation function");
        prm.declare_entry(
          "length ratio",
          "4",
          Patterns::Double(),
          "The length ratio used to define the points for the IB stencil. See definition of epsilon_n in the paper on sharp IB.");
        prm.declare_entry(
          "enable extrapolation",
          "true",
          Patterns::Bool(),
          "Bool to define if extrapolation should be enabled (default). If disabled, all velocity degrees of freedom "
          "in a cell will be set to the particle velocity if that cell is cut. Setting to false is intended for "
          "debugging purposes.");

        prm.leave_subsection();
      }

      prm.enter_subsection("input file");
      {
        prm.declare_entry(
          "load particles from file",
          "false",
          Patterns::Bool(),
          "Bool to define if particles are loaded from an external file");
        prm.declare_entry("particles file",
                          "particles",
                          Patterns::FileName(),
                          "The file name from which we load the particles");
        prm.leave_subsection();
      }

      prm.enter_subsection("local mesh refinement");
      {
        prm.declare_entry(
          "initial refinement",
          "0",
          Patterns::Integer(),
          "Number of refinements around the particles before the start of the simulation ");
        prm.declare_entry(
          "refine mesh inside radius factor",
          "0.5",
          Patterns::Double(),
          "The factor that multiplies the radius to define the inside bound for the refinement of the mesh");
        prm.declare_entry(
          "refine mesh outside radius factor",
          "1.5",
          Patterns::Double(),
          "The factor that multiplies the radius to define the outside bound for the refinement of the mesh");
        prm.declare_entry(
          "refinement zone extrapolation",
          "false",
          Patterns::Bool(),
          "This parameter enables the extrapolation in time of the refinement zone. This means that it will try to refine where the particle will be at the end of the time step instead of the initial position.");
        prm.leave_subsection();
      }

      prm.enter_subsection("output");
      {
        prm.declare_entry(
          "calculate force",
          "true",
          Patterns::Bool(),
          "Bool to define if the force is evaluated on each particle ");
        prm.declare_entry(
          "ib force output file",
          "ib_force",
          Patterns::FileName(),
          "The name of the file where the data on the force of each particle is stored");
        prm.declare_entry(
          "ib particles pvd file",
          "ib_particles_data",
          Patterns::FileName(),
          "The output files of the pvd data for the ib particles");
        prm.declare_entry(
          "print DEM",
          "true",
          Patterns::Bool(),
          "Bool to define if particles' information are printed on the terminal when particles' time-step is finished");
        prm.declare_entry(
          "enable extra sharp interface vtu output field",
          "false",
          Patterns::Bool(),
          "This parameter enables the output of more information related to the particles in the vtu file.");
        prm.leave_subsection();
      }

      prm.enter_subsection("DEM");
      {
        prm.declare_entry(
          "contact search radius factor",
          "3",
          Patterns::Double(),
          "The factor that multiplies the radius to define the region of contact search around the particle");
        prm.declare_entry(
          "contact search frequency",
          "1",
          Patterns::Integer(),
          "The frequency of update in the contact candidates list");
        prm.declare_entry(
          "particle nonlinear tolerance",
          "1e-6",
          Patterns::Double(),
          "The nonlinear tolerance for the coupling of the particle dynamics and the fluid");
        prm.declare_entry("DEM coupling frequency",
                          "1000",
                          Patterns::Integer(),
                          "The number of DEM time steps per CFD time step");
        prm.declare_entry("alpha",
                          "1",
                          Patterns::Double(),
                          "relaxation parameter");

        prm.declare_entry("enable lubrication force",
                          "true",
                          Patterns::Bool(),
                          "Bool to enable or disable the lubrication force");
        prm.declare_entry(
          "lubrication range max",
          "2",
          Patterns::Double(),
          "Gap require to consider the lubrication force. This value is multiplied the smallest cell size");
        prm.declare_entry(
          "lubrication range min",
          "0.1",
          Patterns::Double(),
          "Smallest gap considered for the lubrification force calculation. This value is multiplied by the smallest cell size");

        prm.declare_entry(
          "explicit contact impulsion",
          "false",
          Patterns::Bool(),
          "Bool to enable or disable the use of explicit contact impulsion evaluation in the resolution of the coupling of the particle. When it is set to true, this parameter results in the code only performing the DEM calculation once per CFD time step and using the resulting contact impulsion to evaluate all the other Newton's iterations. This reduces the number of times the DEM calculation is made.");

        prm.declare_entry(
          "explicit position integration",
          "false",
          Patterns::Bool(),
          "Bool to enable or disable the explicit position integration. This means that the particle position is obtained directly by the integration of the previous velocities only. This avoids multiple cut cell mapping for each newton iteration. Note that this limits the order of convergence in time to one.");

        prm.declare_entry(
          "approximate radius for contact",
          "false",
          Patterns::Bool(),
          "Bool to turn on or off using the approximate radius of the particles during contact. If activated, the radius used in the contact calculation is constant and fixed to the effective radius of the shape. If not, the radius of curvature of the shape at the contact point is evaluated. For some shapes, this can be numerically expensive to evaluate.");



        prm.enter_subsection("wall physical properties");
        {
          prm.declare_entry(
            "wall youngs modulus",
            "100000000",
            Patterns::Double(),
            "The wall Young's modulus if IB particles are in contact with it");

          prm.declare_entry(
            "wall poisson ratio",
            "0.3",
            Patterns::Double(),
            "The wall poisson ratio if IB particles are in contact with it");

          prm.declare_entry(
            "wall rolling friction coefficient",
            "0",
            Patterns::Double(),
            "The wall rolling friction coefficient if IB particles are in contact with it");

          prm.declare_entry(
            "wall friction coefficient",
            "0",
            Patterns::Double(),
            "The wall friction coefficient if IB particles are in contact with it");

          prm.declare_entry(
            "wall restitution coefficient",
            "1",
            Patterns::Double(),
            "The wall restitution coefficient if IB particles are in contact with it");
          prm.leave_subsection();
        }

        prm.enter_subsection("gravity");
        f_gravity->declare_parameters(prm, dim);
        prm.leave_subsection();

        prm.leave_subsection();
      }

      unsigned int max_ib_particles = 10;
      particles.resize(max_ib_particles);
      for (unsigned int i = 0; i < max_ib_particles; ++i)
        {
          std::string section = "particle info " + std::to_string(i);
          prm.enter_subsection(section);
          {
            declare_default_entry(prm, i);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  IBParticles<dim>::parse_parameters(ParameterHandler &prm)
  {
    using numbers::PI;
    prm.enter_subsection("particles");
    {
      prm.enter_subsection("extrapolation function");
      {
        order                = prm.get_integer("stencil order");
        length_ratio         = prm.get_double("length ratio");
        enable_extrapolation = prm.get_bool("enable extrapolation");
        prm.leave_subsection();
      }

      prm.enter_subsection("input file");
      {
        load_particles_from_file = prm.get_bool("load particles from file");
        particles_file           = prm.get("particles file");
        prm.leave_subsection();
      }

      prm.enter_subsection("local mesh refinement");
      {
        initial_refinement = prm.get_integer("initial refinement");
        inside_radius      = prm.get_double("refine mesh inside radius factor");
        outside_radius = prm.get_double("refine mesh outside radius factor");
        time_extrapolation_of_refinement_zone =
          prm.get_bool("refinement zone extrapolation");
        prm.leave_subsection();
      }

      prm.enter_subsection("output");
      {
        calculate_force_ib   = prm.get_bool("calculate force");
        ib_force_output_file = prm.get("ib force output file");
        print_dem            = prm.get_bool("print DEM");
        enable_extra_sharp_interface_vtu_output_field =
          prm.get_bool("enable extra sharp interface vtu output field");
        ib_particles_pvd_file = prm.get("ib particles pvd file");
        prm.leave_subsection();
      }
      prm.enter_subsection("DEM");
      {
        alpha = prm.get_double("alpha");

        contact_search_radius_factor =
          prm.get_double("contact search radius factor");
        if (contact_search_radius_factor < 1.)
          throw(std::logic_error(
            "Error, the parameter 'contact search radius factor' cannot be < 1."));
        contact_search_frequency = prm.get_integer("contact search frequency");
        particle_nonlinear_tolerance =
          prm.get_double("particle nonlinear tolerance");
        coupling_frequency       = prm.get_integer("DEM coupling frequency");
        enable_lubrication_force = prm.get_bool("enable lubrication force");
        lubrication_range_max    = prm.get_double("lubrication range max");
        lubrication_range_min    = prm.get_double("lubrication range min");
        explicit_contact_impulsion_calculation =
          prm.get_bool("explicit contact impulsion");
        explicit_position_integration_calculation =
          prm.get_bool("explicit position integration");
        approximate_radius_for_contact =
          prm.get_bool("approximate radius for contact");

        prm.enter_subsection("wall physical properties");
        {
          wall_youngs_modulus = prm.get_double("wall youngs modulus");
          wall_poisson_ratio  = prm.get_double("wall poisson ratio");
          wall_rolling_friction_coefficient =
            prm.get_double("wall rolling friction coefficient");
          wall_friction_coefficient =
            prm.get_double("wall friction coefficient");
          wall_restitution_coefficient =
            prm.get_double("wall restitution coefficient");
          prm.leave_subsection();
        }
        f_gravity = std::make_shared<Functions::ParsedFunction<dim>>(dim);
        prm.enter_subsection("gravity");
        f_gravity->parse_parameters(prm);
        prm.leave_subsection();
        prm.leave_subsection();
      }

      nb_particles = prm.get_integer("number of particles");

      assemble_navier_stokes_inside =
        prm.get_bool("assemble Navier-Stokes inside particles");

      particles.resize(nb_particles);
      for (unsigned int i = 0; i < nb_particles; ++i)
        {
          particles[i].initialize_all();
          std::string section = "particle info " + std::to_string(i);
          prm.enter_subsection(section);

          particles[i].integrate_motion = prm.get_bool("integrate motion");
          particles[i].mesh_based_precalculations =
            prm.get_bool("mesh-based precalculations");

          prm.enter_subsection("position");
          particles[i].f_position->parse_parameters(prm);
          particles[i].f_position->set_time(0);
          prm.leave_subsection();
          prm.enter_subsection("orientation");
          particles[i].f_orientation->parse_parameters(prm);
          particles[i].f_orientation->set_time(0);
          prm.leave_subsection();

          prm.enter_subsection("velocity");
          particles[i].f_velocity->parse_parameters(prm);
          particles[i].f_velocity->set_time(0);
          prm.leave_subsection();
          prm.enter_subsection("omega");
          particles[i].f_omega->parse_parameters(prm);
          particles[i].f_omega->set_time(0);
          prm.leave_subsection();
          particles[i].position[0] =
            particles[i].f_position->value(particles[i].position, 0);
          particles[i].position[1] =
            particles[i].f_position->value(particles[i].position, 1);
          particles[i].orientation[0] =
            particles[i].f_orientation->value(particles[i].position, 0);
          particles[i].orientation[1] =
            particles[i].f_orientation->value(particles[i].position, 1);
          particles[i].orientation[2] =
            particles[i].f_orientation->value(particles[i].position, 2);
          particles[i].velocity[0] =
            particles[i].f_velocity->value(particles[i].position, 0);
          particles[i].velocity[1] =
            particles[i].f_velocity->value(particles[i].position, 1);
          particles[i].omega[0] =
            particles[i].f_omega->value(particles[i].position, 0);
          particles[i].omega[1] =
            particles[i].f_omega->value(particles[i].position, 1);
          particles[i].omega[2] =
            particles[i].f_omega->value(particles[i].position, 2);

          particles[i].particle_id = i;

          std::string pressure_location_str = prm.get("pressure location");
          std::vector<std::string> pressure_location_str_list(
            Utilities::split_string_list(pressure_location_str, ";"));
          std::vector<double> pressure_list =
            Utilities::string_to_double(pressure_location_str_list);
          particles[i].pressure_location[0] = pressure_list[0];
          particles[i].pressure_location[1] = pressure_list[1];


          if (dim == 3)
            {
              particles[i].position[2] =
                particles[i].f_position->value(particles[i].position, 2);
              particles[i].velocity[2] =
                particles[i].f_velocity->value(particles[i].position, 2);
              particles[i].pressure_location[2] = pressure_list[2];
            }

          std::string shape_type          = prm.get("type");
          std::string shape_arguments_str = prm.get("shape arguments");
          particles[i].initialize_shape(shape_type, shape_arguments_str);

          particles[i].set_layer_thickening(prm.get_double("layer thickening"));

          particles[i].radius = particles[i].shape->effective_radius;
          prm.enter_subsection("physical properties");
          {
            std::string              inertia_str = prm.get("inertia");
            std::vector<std::string> inertia_str_list(
              Utilities::split_string_list(inertia_str, ";"));
            const std::vector<double> inertia_list =
              Utilities::string_to_double(inertia_str_list);
            if (inertia_str_list.size() == 9)
              {
                particles[i].inertia[0][0] = inertia_list[0];
                particles[i].inertia[0][1] = inertia_list[1];
                particles[i].inertia[0][2] = inertia_list[2];
                particles[i].inertia[1][0] = inertia_list[3];
                particles[i].inertia[1][1] = inertia_list[4];
                particles[i].inertia[1][2] = inertia_list[5];
                particles[i].inertia[2][0] = inertia_list[6];
                particles[i].inertia[2][1] = inertia_list[7];
                particles[i].inertia[2][2] = inertia_list[8];
              }
            else if (inertia_str_list.size() == 1)
              {
                // If only one inertia value is given, we assume that the
                // inertia is uniform in all axes.
                particles[i].inertia[0][0] = inertia_list[0];
                particles[i].inertia[0][1] = 0;
                particles[i].inertia[0][2] = 0;
                particles[i].inertia[1][0] = 0;
                particles[i].inertia[1][1] = inertia_list[0];
                particles[i].inertia[1][2] = 0;
                particles[i].inertia[2][0] = 0;
                particles[i].inertia[2][1] = 0;
                particles[i].inertia[2][2] = inertia_list[0];
              }
            else
              {
                throw(std::runtime_error(
                  " Invalid inertia matrix. The inertia is given as a 3 by 3 matrix or a single value if the inertia is uniform around each axis."));
              }

            particles[i].youngs_modulus = prm.get_double("youngs modulus");
            particles[i].restitution_coefficient =
              prm.get_double("restitution coefficient");
            particles[i].friction_coefficient =
              prm.get_double("friction coefficient");
            particles[i].poisson_ratio = prm.get_double("poisson ratio");
            particles[i].rolling_friction_coefficient =
              prm.get_double("rolling friction coefficient");

            double volume = prm.get_double("volume");
            if (volume == 0)
              {
                // value is automatically defined.
                volume = particles[i].shape->displaced_volume();
                if (volume == 0)
                  {
                    if (dim == 2)
                      {
                        volume = PI * particles[i].radius * particles[i].radius;
                      }
                    else if (dim == 3)
                      {
                        volume = 4.0 / 3.0 * PI * particles[i].radius *
                                 particles[i].radius * particles[i].radius;
                      }
                  }
              }
            particles[i].volume = volume;
            particles[i].mass = particles[i].volume * prm.get_double("density");

            particles[i].initialize_previous_solution();
            particles[i].set_position(particles[i].position);
            particles[i].set_orientation(particles[i].orientation);

            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
      prm.leave_subsection();
    }
  }

  void
  DynamicFlowControl::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("flow control");
    {
      prm.declare_entry("enable",
                        "false",
                        Patterns::Bool(),
                        "Enable flow rate control");
      prm.declare_entry("average velocity",
                        "0",
                        Patterns::Double(),
                        "The target average velocity");
      prm.declare_entry("inlet boundary id",
                        "0",
                        Patterns::Integer(),
                        "Boundary id of the inlet flow");
      prm.declare_entry("flow direction",
                        "0",
                        Patterns::Integer(),
                        "Flow direction at flow inlet");
      prm.declare_entry("initial beta",
                        "0",
                        Patterns::Double(),
                        "Beta coefficient value for the first step time");
      prm.declare_entry("alpha",
                        "1",
                        Patterns::Double(),
                        "Relaxation coefficient for flow controller");
      prm.declare_entry("enable beta particle",
                        "false",
                        Patterns::Bool(),
                        "Enable beta force for particles");
      prm.declare_entry("beta threshold",
                        "0.0",
                        Patterns::Double(),
                        "Enable beta force for particles");
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the flow control information should be printed "
        "Choices are <quiet|verbose>.");
    }
    prm.leave_subsection();
  }

  void
  DynamicFlowControl::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("flow control");
    {
      // Enable flow control
      enable_flow_control = prm.get_bool("enable");

      // Set the target value for the flow control and flow direction
      average_velocity_0 = prm.get_double("average velocity");
      flow_direction     = prm.get_integer("flow direction");
      boundary_flow_id   = prm.get_integer("inlet boundary id");

      // Tuning parameters for the flow controller
      beta_0         = prm.get_double("initial beta");
      alpha          = prm.get_double("alpha");
      beta_threshold = prm.get_double("beta threshold");

      // Enable printing of flow control information
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;

      // Enable beta force for particles (CFD-DEM)
      enable_beta_particle = prm.get_bool("enable beta particle");
    }

    prm.leave_subsection();
  }

  void
  Evaporation::declare_parameters(dealii::ParameterHandler &prm)
  {
    prm.enter_subsection("evaporation");
    {
      prm.declare_entry(
        "evaporation mass flux model",
        "constant",
        Patterns::Selection("constant|temperature_dependent"),
        "Model used for the calculation of the evaporative mass flux"
        "Choices are <constant|temperature_dependent>.");
      prm.declare_entry(
        "enable evaporative cooling",
        "false",
        Patterns::Bool(),
        "Enable the evaporative cooling at the free surface (gas/liquid interface) in the energy equation <true|false>");
      prm.declare_entry(
        "enable recoil pressure",
        "false",
        Patterns::Bool(),
        "Enable the recoil pressure due to evaporation at the free surface (gas/liquid interface) in the momentum equation <true|false>");
      prm.declare_entry(
        "evaporation mass flux",
        "0.0",
        Patterns::Double(),
        "Evaporation mass flux used if the constant evaporation model is selected in M*L^-2*T^-1");
      prm.declare_entry(
        "evaporation coefficient",
        "0.82",
        Patterns::Double(),
        "Evaporation coefficient corresponding to the ratio between the net mass flux (evaporation-condensation) and the mass flux of evaporation");
      prm.declare_entry(
        "recoil pressure coefficient",
        "0.56",
        Patterns::Double(),
        "Recoil pressure coefficient corresponding to the factor applied to the saturation pressure to compute the recoil pressure in an out of equilibrium evaporation");
      prm.declare_entry("molar mass",
                        "1.0",
                        Patterns::Double(),
                        "Molar mass of the material in M*N^-1");
      prm.declare_entry("boiling temperature",
                        "1.0",
                        Patterns::Double(),
                        "Boiling temperature in Theta");
      prm.declare_entry("evaporation latent heat",
                        "0.0",
                        Patterns::Double(),
                        "Latent heat of evaporation in L^2*T^-2");
      prm.declare_entry("ambient pressure",
                        "101325",
                        Patterns::Double(),
                        "Pressure of the ambient gas in M*L^-1*T^-2");
      prm.declare_entry("ambient gas density",
                        "1.0",
                        Patterns::Double(),
                        "Ambient gas density in M*L^-3");
      prm.declare_entry("liquid density",
                        "10.0",
                        Patterns::Double(),
                        "Liquid density in M*L^-3");
      prm.declare_entry("universal gas constant",
                        "8.3145",
                        Patterns::Double(),
                        "Universal gas constant in M*L^2*T^-2*Theta^-1*N^-1");
    }
    prm.leave_subsection();
  }

  void
  Evaporation::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("evaporation");
    {
      std::string op;
      op = prm.get("evaporation mass flux model");
      if (op == "constant")
        {
          evaporative_mass_flux_model_type =
            EvaporativeMassFluxModelType::constant;
        }
      else if (op == "temperature_dependent")
        {
          evaporative_mass_flux_model_type =
            EvaporativeMassFluxModelType::temperature_dependent;
        }
      else
        throw(std::runtime_error(
          "Invalid evaporative mass flux model. The choices are <constant|temperature_dependent>."));

      enable_evaporation_cooling = prm.get_bool("enable evaporative cooling");
      enable_recoil_pressure     = prm.get_bool("enable recoil pressure");

      evaporation_mass_flux   = prm.get_double("evaporation mass flux");
      evaporation_coefficient = prm.get_double("evaporation coefficient");
      recoil_pressure_coefficient =
        prm.get_double("recoil pressure coefficient");
      molar_mass              = prm.get_double("molar mass");
      boiling_temperature     = prm.get_double("boiling temperature");
      latent_heat_evaporation = prm.get_double("evaporation latent heat");
      ambient_pressure        = prm.get_double("ambient pressure");
      ambient_gas_density     = prm.get_double("ambient gas density");
      liquid_density          = prm.get_double("liquid density");
      universal_gas_constant  = prm.get_double("universal gas constant");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Mortar<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mortar");
    {
      prm.declare_entry("enable",
                        "false",
                        Patterns::Bool(),
                        "Enable mortar interface <true|false>");
      rotor_mesh = std::make_shared<Mesh>();
      rotor_mesh->declare_parameters(prm);
      prm.declare_entry("rotor boundary id",
                        "1",
                        Patterns::Integer(),
                        "Rotor boundary ID # of the mortar matching interface");
      prm.declare_entry(
        "stator boundary id",
        "2",
        Patterns::Integer(),
        "Stator boundary ID # of the mortar matching interface");
      std::string default_entry_point = (dim == 2) ? "0., 0." : "0., 0., 0.";
      prm.declare_entry("center of rotation",
                        default_entry_point,
                        Patterns::List(Patterns::Double()),
                        "Center of rotation coordinates of rotor domain");

      prm.declare_entry("rotation axis",
                        "0.0, 0.0, 1.0",
                        Patterns::List(Patterns::Double()),
                        "Unit vector representing the rotor rotation axis."
                        "Choices are <1, 0, 0 | 0, 1, 0 | 0, 0, 1>.");

      prm.enter_subsection("rotor rotation angle");
      rotor_rotation_angle = std::make_shared<Functions::ParsedFunction<dim>>();
      rotor_rotation_angle->declare_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("rotor angular velocity");
      rotor_angular_velocity =
        std::make_shared<Functions::ParsedFunction<dim>>();
      rotor_angular_velocity->declare_parameters(prm);
      prm.leave_subsection();

      prm.declare_entry("penalty factor",
                        "1.",
                        Patterns::Double(),
                        "Penalty factor for mortar elements");
      prm.declare_entry("oversampling factor",
                        "1",
                        Patterns::Integer(),
                        "Oversampling factor for quadrature points");
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the mortar information should be printed "
        "Choices are <quiet|verbose>.");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Mortar<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("mortar");
    {
      enable = prm.get_bool("enable");
      rotor_mesh->parse_parameters(prm);
      rotor_boundary_id  = prm.get_integer("rotor boundary id");
      stator_boundary_id = prm.get_integer("stator boundary id");
      center_of_rotation =
        value_string_to_tensor<dim>(prm.get("center of rotation"));
      rotation_axis = value_string_to_tensor<dim>(prm.get("rotation axis"));
      prm.enter_subsection("rotor rotation angle");
      rotor_rotation_angle->parse_parameters(prm);
      rotor_rotation_angle->set_time(0);
      prm.leave_subsection();

      prm.enter_subsection("rotor angular velocity");
      rotor_angular_velocity->parse_parameters(prm);
      rotor_angular_velocity->set_time(0);
      prm.leave_subsection();

      sip_factor          = prm.get_double("penalty factor");
      oversampling_factor = prm.get_integer("oversampling factor");

      // Enable printing of mortar information
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;
    }
    prm.leave_subsection();
  }

  // Explicitly instantiate template classes and structs
  template class Laser<2>;
  template class Laser<3>;
  template class IBParticles<2>;
  template class IBParticles<3>;
  template struct ConstrainSolidDomain<2>;
  template struct ConstrainSolidDomain<3>;
  template struct Mortar<2>;
  template struct Mortar<3>;

} // namespace Parameters
