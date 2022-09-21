#include <core/parameters.h>

#include <deal.II/base/exceptions.h>

DeclException2(
  PhaseChangeIntervalError,
  double,
  double,
  << "Liquidus temperature : " << arg1
  << " is not strictly superior to Solidus temperature: " << arg2
  << " The liquidus temperature specific is below or equal to the solidus temperature."
  << " The phase change specific heat model requires that T_liquidus>T_solidus.");

DeclException1(NumberOfFluidsError,
               int,
               << "Number of fluids: " << arg1
               << " is not 1 (single phase simulation) or 2 (VOF simulation)");

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

namespace Parameters
{
  void
  SimulationControl::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("simulation control");
    {
      prm.declare_entry(
        "method",
        "steady",
        Patterns::Selection("steady|steady_bdf|bdf1|bdf2|bdf3|sdirk2|sdirk3"),
        "The kind of solver for the linear system. "
        "Choices are <steady|steady_bdf|bdf1|bdf2|bdf3|sdirk2|sdirk3>.");

      prm.declare_entry(
        "bdf startup method",
        "multiple step bdf",
        Patterns::Selection("multiple step bdf|sdirk step|initial solution"),
        "The kind of method used to startup high order bdf methods "
        "Choices are <multiple step bdf|sdirk step|initial solution>.");

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
      prm.declare_entry("number mesh adapt",
                        "0",
                        Patterns::Integer(),
                        "Number of mesh adaptation (for steady simulations)");
      prm.declare_entry("max cfl",
                        "1",
                        Patterns::Double(),
                        "Maximum CFL value");
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
                        "Output frequency");

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


      prm.declare_entry("output time", "1", Patterns::Double(), "Output time");

      prm.declare_entry(
        "output control",
        "iteration",
        Patterns::Selection("iteration|time"),
        "The control for the output of the simulation results"
        "Results can be either outputted at constant iteration frequency or at constant time");

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
      else if (sv == "sdirk2")
        method = TimeSteppingMethod::sdirk22;
      else if (sv == "sdirk3")
        method = TimeSteppingMethod::sdirk33;
      else
        {
          std::runtime_error("Invalid time stepping scheme");
        }
      const std::string bdf_startup_string = prm.get("bdf startup method");
      if (bdf_startup_string == "multiple step bdf")
        bdf_startup_method = BDFStartupMethods::multiple_step_bdf;
      else if (bdf_startup_string == "sdirk step")
        bdf_startup_method = BDFStartupMethods::sdirk_step;
      else if (bdf_startup_string == "initial solution")
        bdf_startup_method = BDFStartupMethods::initial_solution;
      else
        {
          std::runtime_error("Invalid bdf startup scheme");
        }

      const std::string osv = prm.get("output control");
      if (osv == "iteration")
        output_control = OutputControl::iteration;
      else if (osv == "time")
        output_control = OutputControl::time;
      else
        {
          std::runtime_error("Invalid output control scheme");
        }
      dt             = prm.get_double("time step");
      timeEnd        = prm.get_double("time end");
      adapt          = prm.get_bool("adapt");
      maxCFL         = prm.get_double("max cfl");
      stop_tolerance = prm.get_double("stop tolerance");
      adaptative_time_step_scaling =
        prm.get_double("adaptative time step scaling");
      startup_timestep_scaling = prm.get_double("startup time scaling");
      number_mesh_adaptation   = prm.get_integer("number mesh adapt");

      output_folder     = prm.get("output path");
      output_name       = prm.get("output name");
      output_frequency  = prm.get_integer("output frequency");
      output_time       = prm.get_double("output time");
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
  PowerLawParameters::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("power-law");
    {
      K              = prm.get_double("K");
      n              = prm.get_double("n");
      shear_rate_min = prm.get_double("shear rate min");
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
                        "Viscosity at rest");
      prm.declare_entry("viscosity_inf",
                        "1.0",
                        Patterns::Double(),
                        "Viscosity for an infinite shear rate");
      prm.declare_entry("lambda", "1.0", Patterns::Double(), "Relaxation time");
      prm.declare_entry("a", "2.0", Patterns::Double(), "Carreau parameter");
      prm.declare_entry("n", "0.5", Patterns::Double(), "Power parameter");
    }
    prm.leave_subsection();
  }

  void
  CarreauParameters::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("carreau");
    {
      viscosity_0   = prm.get_double("viscosity_0");
      viscosity_inf = prm.get_double("viscosity_inf");
      lambda        = prm.get_double("lambda");
      a             = prm.get_double("a");
      n             = prm.get_double("n");
    }
    prm.leave_subsection();
  }

  void
  NonNewtonian::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("non newtonian");
    {
      powerlaw_parameters.declare_parameters(prm);
      carreau_parameters.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  void
  NonNewtonian::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("non newtonian");
    {
      powerlaw_parameters.parse_parameters(prm);
      carreau_parameters.parse_parameters(prm);
    }
    prm.leave_subsection();
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
    }
    prm.leave_subsection();
  }

  void
  Stabilization::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("stabilization");
    {
      use_default_stabilization = prm.get_bool("use default stabilization");
      std::string op            = prm.get("stabilization");
      if (op == "pspg_supg")
        stabilization = NavierStokesStabilization::pspg_supg;
      else if (op == "gls")
        stabilization = NavierStokesStabilization::gls;
      else if (op == "grad_div")
        stabilization = NavierStokesStabilization::grad_div;
      else
        throw(std::runtime_error("Invalid stabilization strategy"));
    }
    prm.leave_subsection();
  }

  void
  PhaseChange::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("phase change");
    {
      T_solidus              = prm.get_double("solidus temperature");
      T_liquidus             = prm.get_double("liquidus temperature");
      latent_enthalpy        = prm.get_double("latent enthalpy");
      cp_l                   = prm.get_double("specific heat liquid");
      cp_s                   = prm.get_double("specific heat solid");
      viscosity_l            = prm.get_double("viscosity liquid");
      viscosity_s            = prm.get_double("viscosity solid");
      thermal_conductivity_l = prm.get_double("thermal conductivity liquid");
      thermal_conductivity_s = prm.get_double("thermal conductivity solid");
      thermal_expansion_l    = prm.get_double("thermal expansion liquid");
      thermal_expansion_s    = prm.get_double("thermal expansion solid");
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
                        "1",
                        Patterns::Double(),
                        "Thermal expansion coefficient of the solid phase");

      prm.declare_entry("viscosity liquid",
                        "1",
                        Patterns::Double(),
                        "Viscosity of the liquid phase");

      prm.declare_entry("viscosity solid",
                        "1",
                        Patterns::Double(),
                        "Viscosity of the solid phase");
    }
    prm.leave_subsection();
  }

  void
  PhysicalProperties::declare_parameters(ParameterHandler &prm)
  {
    fluids.resize(max_fluids);
    number_of_fluids = 1;

    prm.enter_subsection("physical properties");
    {
      prm.declare_entry("number of fluids",
                        "1",
                        Patterns::Integer(),
                        "Number of fluids");

      // Multiphasic simulations parameters definition
      for (unsigned int i_fluid = 0; i_fluid < max_fluids; ++i_fluid)
        {
          fluids[i_fluid].declare_parameters(prm, i_fluid);
        }
    }
    prm.leave_subsection();
  }

  void
  PhysicalProperties::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("physical properties");
    {
      // Multiphasic simulations parameters definition
      number_of_fluids = prm.get_integer("number of fluids");
      Assert(number_of_fluids == 1 || number_of_fluids == 2,
             NumberOfFluidsError(number_of_fluids));

      for (unsigned int i_fluid = 0; i_fluid < number_of_fluids; ++i_fluid)
        {
          fluids[i_fluid].parse_parameters(prm, i_fluid);
        }
    }
    prm.leave_subsection();
  }

  void
  Fluid::declare_parameters(ParameterHandler &prm, unsigned int id)
  {
    prm.enter_subsection("fluid " + Utilities::int_to_string(id, 1));
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
        "tracer diffusivity",
        "0",
        Patterns::Double(),
        "Tracer diffusivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "rheological model",
        "newtonian",
        Patterns::Selection("newtonian|power-law|carreau|phase_change"),
        "Rheological model "
        "Choices are <newtonian|power-law|carreau|phase_change>.");

      non_newtonian_parameters.declare_parameters(prm);


      prm.declare_entry("density model",
                        "constant",
                        Patterns::Selection("constant"),
                        "Model used for the calculation of the density"
                        "Choices are <constant>.");

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
  Fluid::parse_parameters(ParameterHandler &prm, unsigned int id)
  {
    prm.enter_subsection("fluid " + Utilities::int_to_string(id, 1));
    {
      density              = prm.get_double("density");
      viscosity            = prm.get_double("kinematic viscosity");
      specific_heat        = prm.get_double("specific heat");
      thermal_conductivity = prm.get_double("thermal conductivity");
      thermal_expansion    = prm.get_double("thermal expansion");
      tracer_diffusivity   = prm.get_double("tracer diffusivity");


      // Parse models for the physical properties
      std::string op;

      // Density
      op = prm.get("density model");
      if (op == "constant")
        density_model = DensityModel::constant;

      // Thermal conductivity
      op = prm.get("thermal conductivity model");
      if (op == "constant")
        thermal_conductivity_model = ThermalConductivityModel::constant;
      else if (op == "linear")
        thermal_conductivity_model = ThermalConductivityModel::linear;
      else if (op == "phase_change")
        thermal_conductivity_model = ThermalConductivityModel::phase_change;

      // Linear conductivity model parameters
      k_A0 = prm.get_double("k_A0");
      k_A1 = prm.get_double("k_A1");

      // Thermal expansion
      op = prm.get("thermal expansion model");
      if (op == "constant")
        thermal_expansion_model = ThermalExpansionModel::constant;
      else if (op == "phase_change")
        thermal_expansion_model = ThermalExpansionModel::phase_change;

      // Specific heat
      op = prm.get("specific heat model");
      if (op == "constant")
        specific_heat_model = SpecificHeatModel::constant;
      else if (op == "phase_change")
        specific_heat_model = SpecificHeatModel::phase_change;

      phase_change_parameters.parse_parameters(prm);

      // Rheology
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

      non_newtonian_parameters.parse_parameters(prm);
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
      prm.declare_entry("qmapping all",
                        "false",
                        Patterns::Bool(),
                        "Apply high order mapping everywhere");
    }
    prm.leave_subsection();
  }

  void
  FEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("FEM");
    {
      velocity_order      = prm.get_integer("velocity order");
      pressure_order      = prm.get_integer("pressure order");
      void_fraction_order = prm.get_integer("void fraction order");
      temperature_order   = prm.get_integer("temperature order");
      tracer_order        = prm.get_integer("tracer order");
      VOF_order           = prm.get_integer("VOF order");
      qmapping_all        = prm.get_bool("qmapping all");
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

  template <int dim>
  void
  Laser<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("laser parameters");
    {
      prm.declare_entry("enable", "false", Patterns::Bool(), "Activate laser");
      prm.declare_entry("concentration factor",
                        "2.0",
                        Patterns::Double(),
                        "Concentration factor");
      prm.declare_entry("power", "100.0", Patterns::Double(), "Laser power");
      prm.declare_entry("absorptivity",
                        "0.5",
                        Patterns::Double(),
                        "Laser absorptivity");
      prm.declare_entry("penetration depth",
                        "0.0",
                        Patterns::Double(),
                        "Penetration depth");
      prm.declare_entry("beam radius",
                        "0.0",
                        Patterns::Double(),
                        "Laser beam radius");


      prm.enter_subsection("path");
      laser_scan_path = std::make_shared<Functions::ParsedFunction<dim>>(dim);
      laser_scan_path->declare_parameters(prm, dim);
      if (dim == 2)
        prm.set("Function expression", "0; 0");
      if (dim == 3)
        prm.set("Function expression", "0; 0; 0");
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
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Laser<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("laser parameters");
    {
      activate_laser       = prm.get_bool("enable");
      concentration_factor = prm.get_double("concentration factor");
      laser_power          = prm.get_double("power");
      laser_absorptivity   = prm.get_double("absorptivity");
      penetration_depth    = prm.get_double("penetration depth");
      beam_radius          = prm.get_double("beam radius");

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
          beam_direction                     = 1;
          beam_orientation                   = BeamOrientation::x_plus;
          beam_orientation_coordinate        = 0;
          perpendicular_plane_coordinate_one = 1;
          if constexpr (dim == 3)
            perpendicular_plane_coordinate_two = 2;
        }
      else if (op == "x-")
        {
          beam_direction                     = 0;
          beam_orientation                   = BeamOrientation::x_minus;
          beam_orientation_coordinate        = 0;
          perpendicular_plane_coordinate_one = 1;
          if constexpr (dim == 3)
            perpendicular_plane_coordinate_two = 2;
        }
      else if (op == "y+")
        {
          beam_direction                     = 1;
          beam_orientation                   = BeamOrientation::y_plus;
          perpendicular_plane_coordinate_one = 0;
          beam_orientation_coordinate        = 1;
          if constexpr (dim == 3)
            perpendicular_plane_coordinate_two = 2;
        }
      else if (op == "y-")
        {
          beam_direction                     = 0;
          beam_orientation                   = BeamOrientation::y_minus;
          perpendicular_plane_coordinate_one = 0;
          beam_orientation_coordinate        = 1;
          if constexpr (dim == 3)
            perpendicular_plane_coordinate_two = 2;
        }
      else if (op == "z+")
        {
          if constexpr (dim == 3)
            {
              beam_direction                     = 1;
              beam_orientation                   = BeamOrientation::z_plus;
              perpendicular_plane_coordinate_one = 0;
              perpendicular_plane_coordinate_two = 1;
              beam_orientation_coordinate        = 2;
            }
          else if constexpr (dim == 2)
            Assert(dim == 2, TwoDimensionalLaserError(dim));
        }
      else if (op == "z-")
        {
          if constexpr (dim == 3)
            {
              beam_direction                     = 0;
              beam_orientation                   = BeamOrientation::z_minus;
              perpendicular_plane_coordinate_one = 0;
              perpendicular_plane_coordinate_two = 1;
              beam_orientation_coordinate        = 2;
            }
          else if constexpr (dim == 2)
            Assert(dim == 2, TwoDimensionalLaserError(dim));
        }
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

      prm.declare_entry("calculate apparent viscosity",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of apparent viscosity");

      prm.declare_entry("calculate average velocities",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of average velocities.");

      prm.declare_entry("calculate pressure drop",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of between two boundaries.");

      prm.declare_entry(
        "calculate temperature range",
        "false",
        Patterns::Bool(),
        "Enable calculation of minimum and maximum temperature.");

      prm.declare_entry("inlet boundary id",
                        "0",
                        Patterns::Integer(),
                        "Inlet boundary ID for pressure drop calculation");

      prm.declare_entry("outlet boundary id",
                        "1",
                        Patterns::Integer(),
                        "Outlet boundary ID for pressure drop calculation");

      prm.declare_entry(
        "initial time",
        "0.0",
        Patterns::Double(),
        "Initial time to start calculations for average velocities");

      prm.declare_entry("kinetic energy name",
                        "kinetic_energy",
                        Patterns::FileName(),
                        "File output kinetic energy");

      prm.declare_entry("pressure drop name",
                        "pressure_drop",
                        Patterns::FileName(),
                        "File output pressure drop");

      prm.declare_entry("enstrophy name",
                        "enstrophy",
                        Patterns::FileName(),
                        "File output enstrophy");

      prm.declare_entry("apparent viscosity name",
                        "apparent_viscosity",
                        Patterns::FileName(),
                        "File output apparent viscosity");

      prm.declare_entry("calculation frequency",
                        "1",
                        Patterns::Integer(),
                        "Calculation frequency");

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
      calculate_apparent_viscosity =
        prm.get_bool("calculate apparent viscosity");
      calculate_average_velocities =
        prm.get_bool("calculate average velocities");
      calculate_pressure_drop = prm.get_bool("calculate pressure drop");
      calculate_min_max_temperature =
        prm.get_bool("calculate temperature range");
      inlet_boundary_id              = prm.get_integer("inlet boundary id");
      outlet_boundary_id             = prm.get_integer("outlet boundary id");
      initial_time                   = prm.get_double("initial time");
      kinetic_energy_output_name     = prm.get("kinetic energy name");
      pressure_drop_output_name      = prm.get("pressure drop name");
      enstrophy_output_name          = prm.get("enstrophy name");
      apparent_viscosity_output_name = prm.get("apparent viscosity name");
      calculation_frequency          = prm.get_integer("calculation frequency");
      output_frequency               = prm.get_integer("output frequency");
      calculate_tracer_statistics = prm.get_bool("calculate tracer statistics");
      tracer_output_name          = prm.get("tracer statistics name");
    }
    prm.leave_subsection();
  }

  void
  NonLinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("non-linear solver");
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
        "abort at convergence failure",
        "false",
        Patterns::Bool(),
        "Aborts Lethe by throwing an exception if non-linear solver convergence has failed");
    }
    prm.leave_subsection();
  }

  void
  NonLinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("non-linear solver");
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
        throw(
          std::runtime_error("Invalid strategy for kinsol non-linear solver "));

      tolerance             = prm.get_double("tolerance");
      step_tolerance        = prm.get_double("step tolerance");
      matrix_tolerance      = prm.get_double("matrix tolerance");
      max_iterations        = prm.get_integer("max iterations");
      display_precision     = prm.get_integer("residual precision");
      force_rhs_calculation = prm.get_bool("force rhs calculation");
      reuse_matrix          = prm.get_bool("reuse matrix");
      abort_at_convergence_failure =
        prm.get_bool("abort at convergence failure");
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
                        Patterns::Selection("gmsh|dealii|periodic_hills"),
                        "Type of mesh "
                        "Choices are <gmsh|dealii|periodic_hills>.");

      prm.declare_entry("file name",
                        "none",
                        Patterns::FileName(),
                        "GMSH file name");

      prm.declare_entry("initial refinement",
                        "0",
                        Patterns::Integer(),
                        "Initial refinement of the mesh");

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
        "translate",
        "false",
        Patterns::Bool(),
        "Indicates that the mesh should be translated. To be "
        "used to reposition solid grids relative to the fluid grid.");

      prm.declare_entry("delta_x",
                        "0 ",
                        Patterns::Double(),
                        "X component of the desired translation of the mesh");

      prm.declare_entry("delta_y",
                        "0 ",
                        Patterns::Double(),
                        "Y component of the desired translation of the mesh");

      prm.declare_entry("delta_z",
                        "0 ",
                        Patterns::Double(),
                        "Z component of the desired translation of the mesh");

      // Grid rotation at initiation
      prm.declare_entry(
        "rotate",
        "false",
        Patterns::Bool(),
        "Indicates that the mesh should be rotated. To be "
        "used to reposition solid grids relative to the fluid grid.");

      prm.declare_entry(
        "axis",
        "0",
        Patterns::Integer(),
        "Axis around which the rotation should occur (0 for x, 1 for y, 2 for z)");

      prm.declare_entry("angle",
                        "0 ",
                        Patterns::Double(),
                        "Angle of rotation around the axis in radian");
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
        else
          throw std::logic_error(
            "Error, invalid mesh type. Choices are gmsh and dealii");
      }

      file_name = prm.get("file name");

      initial_refinement = prm.get_integer("initial refinement");

      grid_type      = prm.get("grid type");
      grid_arguments = prm.get("grid arguments");

      refine_until_target_size = prm.get_bool("enable target size");
      simplex                  = prm.get_bool("simplex");
      check_for_diamond_cells  = prm.get_bool("check diamond cells");
      expand_particle_wall_contact_search =
        prm.get_bool("expand particle-wall contact search");
      target_size = prm.get_double("target size");


      translate = prm.get_bool("translate");
      delta_x   = prm.get_double("delta_x");
      delta_y   = prm.get_double("delta_y");
      delta_z   = prm.get_double("delta_z");

      rotate = prm.get_bool("rotate");
      axis   = prm.get_integer("axis");
      angle  = prm.get_double("angle");
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
  LinearSolver::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("linear solver");
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
        Patterns::Selection("gmres|bicgstab|amg|direct"),
        "The iterative solver for the linear system of equations. "
        "Choices are <gmres|bicgstab|amg|tfqmr|direct>. gmres is a GMRES iterative "
        "solver "
        "with ILU preconditioning. bicgstab is a BICGSTAB iterative solver "
        "with ILU preconditioning. "
        "amg is GMRES + AMG preconditioning with an ILU coarsener and "
        "smoother. On coarse meshes, the gmres/bicgstab solver with ILU "
        "preconditioning is more efficient. "
        "As the number of mesh elements increase, the amg solver is the most "
        "efficient. Generally, at 1M elements, the amg solver always "
        "outperforms the gmres or bicgstab");
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
                        "amg preconditioner ilu smoother/coarsener fill");

      prm.declare_entry(
        "amg preconditioner ilu absolute tolerance",
        "1e-12",
        Patterns::Double(),
        "amg preconditioner ilu smoother/coarsener absolute tolerance");

      prm.declare_entry(
        "amg preconditioner ilu relative tolerance",
        "1.00",
        Patterns::Double(),
        "amg preconditioner ilu smoother/coarsener relative tolerance");

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
    }
    prm.leave_subsection();
  }
  void
  LinearSolver::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("linear solver");
    {
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

      const std::string sv = prm.get("method");
      if (sv == "amg")
        solver = SolverType::amg;
      else if (sv == "gmres")
        solver = SolverType::gmres;
      else if (sv == "bicgstab")
        solver = SolverType::bicgstab;
      else if (sv == "direct")
        solver = SolverType::direct;
      else
        throw std::logic_error(
          "Error, invalid iterative solver type. Choices are amg, gmres, bicgstab or direct");

      relative_residual  = prm.get_double("relative residual");
      minimum_residual   = prm.get_double("minimum residual");
      max_iterations     = prm.get_integer("max iters");
      max_krylov_vectors = prm.get_integer("max krylov vectors");

      ilu_precond_fill = prm.get_double("ilu preconditioner fill");
      ilu_precond_atol =
        prm.get_double("ilu preconditioner absolute tolerance");
      ilu_precond_rtol =
        prm.get_double("ilu preconditioner relative tolerance");
      amg_precond_ilu_fill = prm.get_double("amg preconditioner ilu fill");
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
        Patterns::List(
          Patterns::Selection("velocity|pressure|phase|temperature")),
        "Variable(s) for kelly estimation"
        "Choices are <velocity|pressure|phase|temperature>."
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
          else
            throw std::logic_error(
              "Error, invalid mesh adaptation variable. Choices are velocity, pressure, phase or temperature");

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
      maximum_number_elements  = prm.get_integer("max number elements");
      maximum_refinement_level = prm.get_integer("max refinement level");
      minimum_refinement_level = prm.get_integer("min refinement level");
      frequency                = prm.get_integer("frequency");
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
    }
    prm.leave_subsection();
  }

  void
  Testing::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("test");
    {
      enabled = prm.get_bool("enable");
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
        "type",
        "none",
        Patterns::Selection("none|srf"),
        "Velocity-dependent source terms"
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
    }
    prm.leave_subsection();
  }

  void
  VelocitySource::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("velocity source");
    {
      const std::string op = prm.get("type");
      if (op == "none")
        type = VelocitySourceType::none;
      else if (op == "srf")
        type = VelocitySourceType::srf;
      else
        throw std::logic_error("Error, invalid velocity source type");

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
    prm.enter_subsection("position");
    particles[index].f_position =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    particles[index].f_position->declare_parameters(prm, dim);
    if (dim == 2)
      prm.set("Function expression", "0; 0");
    if (dim == 3)
      prm.set("Function expression", "0; 0; 0");
    prm.leave_subsection();

    prm.enter_subsection("orientation");
    particles[index].f_orientation =
      std::make_shared<Functions::ParsedFunction<dim>>(3);
    particles[index].f_orientation->declare_parameters(prm, dim);
    prm.set("Function expression", "0; 0; 0");
    prm.leave_subsection();

    prm.enter_subsection("velocity");
    particles[index].f_velocity =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    particles[index].f_velocity->declare_parameters(prm, dim);
    if (dim == 2)
      prm.set("Function expression", "0; 0");
    if (dim == 3)
      prm.set("Function expression", "0; 0; 0");
    prm.leave_subsection();

    prm.enter_subsection("omega");
    particles[index].f_omega =
      std::make_shared<Functions::ParsedFunction<dim>>(3);
    particles[index].f_omega->declare_parameters(prm, dim);
    prm.set("Function expression", "0; 0; 0");
    prm.leave_subsection();

    prm.declare_entry(
      "type",
      "sphere",
      Patterns::Selection(
        "sphere|rectangle|ellipsoid|torus|cone|cut hollow sphere|death star|rbf"),
      "The type of shape considered."
      "Choices are <sphere|rectangle|ellipsoid|torus|cone|cut hollow sphere|death star|rbf>."
      "The parameter for a sphere is: radius. "
      "The parameters for a rectangle are, in order: x half length,"
      "y half length, z half length."
      "The parameters for an ellipsoid are, in order: x radius,"
      "y radius, z radius. "
      "The parameters for a torus are, in order: torus radius,"
      "torus thickness radius. "
      "The parameters for a cone are, in order: tan(base angle),"
      " height. "
      "The parameters for a cut hollow sphere are, in order: sphere radius,"
      "cut thickness, wall thickness. "
      "The parameters for a death star are, in order: sphere radius,"
      "smaller sphere radius, distance between centers."
      "The parameter for an rbf is the file name.");
    prm.declare_entry("shape arguments",
                      "1",
                      Patterns::Anything(),
                      "Arguments defining the geometry");

    prm.declare_entry(
      "pressure location",
      "0; 0; 0",
      Patterns::Anything(),
      "position relative to the center of the particle for the location of the point where the pressure is imposed inside the particle");
    prm.declare_entry("density",
                      "1",
                      Patterns::Double(),
                      "density of the particle ");
    prm.declare_entry("inertia",
                      "1",
                      Patterns::Double(),
                      "uniform rotational moment of inertia");
    prm.declare_entry("youngs modulus",
                      "100000000",
                      Patterns::Double(),
                      "The youngs modulus of particle in case of contact");
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
        "Number of particles reprensented by IB max number of particles = 10 ");
      prm.declare_entry(
        "initial refinement",
        "0",
        Patterns::Integer(),
        "number of refinement around the particles before the start of the simulation ");
      prm.declare_entry(
        "stencil order",
        "2",
        Patterns::Integer(),
        "Number of particles reprensented by IB max number of particles = 10 ");
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
        "load particles from file",
        "false",
        Patterns::Bool(),
        "Bool to define if particles are loaded from an external file");
      prm.declare_entry("particles file",
                        "particles",
                        Patterns::FileName(),
                        "The file name from which we load the particles");
      prm.declare_entry(
        "ib particles pvd file",
        "ib_particles_data",
        Patterns::FileName(),
        "The output files of the pvd data for the ib particles");
      prm.declare_entry(
        "integrate motion",
        "false",
        Patterns::Bool(),
        "Bool to define if the particle trajectory is integrated meaning it's velocity and position will be updated at each time step according to the hydrodynamic force applied to it");
      prm.declare_entry(
        "print DEM",
        "true",
        Patterns::Bool(),
        "Bool to define if particles' informations are printed on the terminal when particles' time-step is finished");
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
        "assemble Navier-Stokes inside particles",
        "false",
        Patterns::Bool(),
        "Bool to define if you assemble the inside of particles with the NS equation.");
      prm.declare_entry(
        "length ratio",
        "4",
        Patterns::Double(),
        "The length ratio used to define the points for the IB stencil. See definition of epsilon_n in the paper on sharp IB.");
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
        "wall youngs modulus",
        "100000000",
        Patterns::Double(),
        "The wall youngs modulus if IB particles are in contact with it");

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


      prm.enter_subsection("gravity");
      f_gravity->declare_parameters(prm, dim);
      if (dim == 2)
        prm.set("Function expression", "0; 0");
      if (dim == 3)
        prm.set("Function expression", "0; 0; 0");
      prm.leave_subsection();

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
      nb                 = prm.get_integer("number of particles");
      order              = prm.get_integer("stencil order");
      initial_refinement = prm.get_integer("initial refinement");
      inside_radius      = prm.get_double("refine mesh inside radius factor");
      outside_radius     = prm.get_double("refine mesh outside radius factor");
      calculate_force_ib = prm.get_bool("calculate force");
      ib_force_output_file = prm.get("ib force output file");
      integrate_motion     = prm.get_bool("integrate motion");
      print_dem            = prm.get_bool("print DEM");
      alpha                = prm.get_double("alpha");
      contact_search_radius_factor =
        prm.get_double("contact search radius factor");
      if (contact_search_radius_factor < 1.)
        throw(std::logic_error(
          "Error, the parameter 'contact search radius factor' cannot be < 1."));

      contact_search_frequency = prm.get_integer("contact search frequency");

      load_particles_from_file = prm.get_bool("load particles from file");
      particles_file           = prm.get("particles file");

      assemble_navier_stokes_inside =
        prm.get_bool("assemble Navier-Stokes inside particles");
      length_ratio          = prm.get_double("length ratio");
      ib_particles_pvd_file = prm.get("ib particles pvd file");
      particle_nonlinear_tolerance =
        prm.get_double("particle nonlinear tolerance");

      f_gravity = std::make_shared<Functions::ParsedFunction<dim>>(dim);
      prm.enter_subsection("gravity");
      f_gravity->parse_parameters(prm);
      prm.leave_subsection();

      wall_youngs_modulus = prm.get_double("wall youngs modulus");
      wall_poisson_ratio  = prm.get_double("wall poisson ratio");
      wall_rolling_friction_coefficient =
        prm.get_double("wall rolling friction coefficient");
      wall_friction_coefficient = prm.get_double("wall friction coefficient");
      wall_restitution_coefficient =
        prm.get_double("wall restitution coefficient");
      coupling_frequency       = prm.get_double("DEM coupling frequency");
      enable_lubrication_force = prm.get_bool("enable lubrication force");
      lubrication_range_max    = prm.get_double("lubrication range max");
      lubrication_range_min    = prm.get_double("lubrication range min");

      particles.resize(nb);
      for (unsigned int i = 0; i < nb; ++i)
        {
          particles[i].initialize_all();
          std::string section = "particle info " + std::to_string(i);
          prm.enter_subsection(section);

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

          particles[i].particle_id   = i;
          particles[i].inertia[0][0] = prm.get_double("inertia");
          particles[i].inertia[1][1] = prm.get_double("inertia");
          particles[i].inertia[2][2] = prm.get_double("inertia");

          std::string pressure_location_str = prm.get("pressure location");
          std::vector<std::string> pressure_location_str_list(
            Utilities::split_string_list(pressure_location_str, ";"));
          std::vector<double> pressure_list =
            Utilities::string_to_double(pressure_location_str_list);
          particles[i].pressure_location[0] = pressure_list[0];
          particles[i].pressure_location[1] = pressure_list[1];


          particles[i].youngs_modulus = prm.get_double("youngs modulus");
          particles[i].restitution_coefficient =
            prm.get_double("restitution coefficient");
          particles[i].friction_coefficient =
            prm.get_double("friction coefficient");
          particles[i].poisson_ratio = prm.get_double("poisson ratio");
          particles[i].rolling_friction_coefficient =
            prm.get_double("rolling friction coefficient");

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
          std::vector<std::string> shape_arguments_str_list(
            Utilities::split_string_list(shape_arguments_str, ";"));
          std::vector<double> shape_arguments;
          if (shape_type == "rbf")
            {
              if constexpr (dim == 3)
                {
                  std::string shape_name = shape_arguments_str_list[0];
                  std::map<std::string, std::vector<double>> rbf_data;
                  fill_vectors_from_file(rbf_data, shape_name, ",");
                  size_t number_of_nodes = rbf_data["weight"].size();

                  std::vector<double> weights;
                  std::vector<double> basis_functions;
                  std::vector<double> support_radii;
                  std::vector<double> nodes_x;
                  std::vector<double> nodes_y;
                  std::vector<double> nodes_z;

                  weights.resize(number_of_nodes);
                  basis_functions.resize(number_of_nodes);
                  support_radii.resize(number_of_nodes);
                  nodes_x.resize(number_of_nodes);
                  nodes_y.resize(number_of_nodes);
                  nodes_z.resize(number_of_nodes);

                  for (size_t n_i = 0; n_i < number_of_nodes; n_i++)
                    {
                      weights[n_i]         = rbf_data["weight"][n_i];
                      support_radii[n_i]   = rbf_data["support_radius"][n_i];
                      basis_functions[n_i] = rbf_data["basis_function"][n_i];
                      nodes_x[n_i]         = rbf_data["xnode"][n_i];
                      nodes_y[n_i]         = rbf_data["ynode"][n_i];
                      nodes_z[n_i]         = rbf_data["znode"][n_i];
                    }
                  shape_arguments.reserve((dim + 3) * number_of_nodes);
                  shape_arguments.insert(shape_arguments.end(),
                                         weights.begin(),
                                         weights.end());
                  shape_arguments.insert(shape_arguments.end(),
                                         support_radii.begin(),
                                         support_radii.end());
                  shape_arguments.insert(shape_arguments.end(),
                                         basis_functions.begin(),
                                         basis_functions.end());
                  shape_arguments.insert(shape_arguments.end(),
                                         nodes_x.begin(),
                                         nodes_x.end());
                  shape_arguments.insert(shape_arguments.end(),
                                         nodes_y.begin(),
                                         nodes_y.end());
                  shape_arguments.insert(shape_arguments.end(),
                                         nodes_z.begin(),
                                         nodes_z.end());
                }
            }
          else
            {
              shape_arguments =
                Utilities::string_to_double(shape_arguments_str_list);
            }
          particles[i].initialize_shape(shape_type, shape_arguments);

          particles[i].radius = particles[i].shape->effective_radius;

          if (dim == 2)
            {
              particles[i].mass = PI * particles[i].radius *
                                  particles[i].radius *
                                  prm.get_double("density");
            }
          else if (dim == 3)
            {
              particles[i].mass = 4.0 / 3.0 * PI * particles[i].radius *
                                  particles[i].radius * particles[i].radius *
                                  prm.get_double("density");
            }
          particles[i].initialize_previous_solution();
          prm.leave_subsection();
        }
      prm.leave_subsection();
    }
  }

  template <int dim>
  void
  IBParticles<dim>::initialize_shape(const unsigned int        i,
                                     const std::string         type,
                                     const std::vector<double> shape_arguments)
  {
    if (type == "sphere")
      particles[i].shape =
        std::make_shared<Sphere<dim>>(shape_arguments[0],
                                      particles[i].position,
                                      particles[i].orientation);
    else if (type == "rectangle")
      {
        Tensor<1, dim> half_lengths;
        for (unsigned int i = 0; i < dim; ++i)
          {
            half_lengths[i] = shape_arguments[i];
          }
        particles[i].shape =
          std::make_shared<Rectangle<dim>>(half_lengths,
                                           particles[i].position,
                                           particles[i].orientation);
      }
    else if (type == "ellipsoid")
      {
        Tensor<1, dim> radii;
        for (unsigned int i = 0; i < dim; ++i)
          {
            radii[i] = shape_arguments[i];
          }
        particles[i].shape =
          std::make_shared<Ellipsoid<dim>>(radii,
                                           particles[i].position,
                                           particles[i].orientation);
      }
    else if (type == "torus")
      {
        if constexpr (dim == 3)

          particles[i].shape =
            std::make_shared<Torus<dim>>(shape_arguments[0],
                                         shape_arguments[1],
                                         particles[i].position,
                                         particles[i].orientation);
      }
    else if (type == "cone")
      {
        if constexpr (dim == 3)
          particles[i].shape =
            std::make_shared<Cone<dim>>(shape_arguments[0],
                                        shape_arguments[1],
                                        particles[i].position,
                                        particles[i].orientation);
      }
    else if (type == "cut hollow sphere")
      {
        if constexpr (dim == 3)
          particles[i].shape =
            std::make_shared<CutHollowSphere<dim>>(shape_arguments[0],
                                                   shape_arguments[1],
                                                   shape_arguments[2],
                                                   particles[i].position,
                                                   particles[i].orientation);
      }
    else if (type == "death star")
      {
        if constexpr (dim == 3)
          particles[i].shape =
            std::make_shared<DeathStar<dim>>(shape_arguments[0],
                                             shape_arguments[1],
                                             shape_arguments[2],
                                             particles[i].position,
                                             particles[i].orientation);
      }
    else if (type == "rbf")
      {
        if constexpr (dim == 3)
          {
            constexpr unsigned int numbers_per_nodes = dim + 3;
            size_t number_of_nodes = shape_arguments.size() / numbers_per_nodes;
            std::vector<double>                               support_radii;
            std::vector<enum RBFShape<dim>::RBFBasisFunction> basis_functions;
            std::vector<double>                               weights;
            std::vector<Point<dim>>                           nodes;
            support_radii.resize(number_of_nodes);
            basis_functions.resize(number_of_nodes);
            weights.resize(number_of_nodes);
            nodes.resize(number_of_nodes);
            for (size_t n_i = 0; n_i < number_of_nodes; n_i++)
              {
                weights[n_i]       = shape_arguments[0 * number_of_nodes + n_i];
                support_radii[n_i] = shape_arguments[1 * number_of_nodes + n_i];
                basis_functions[n_i] =
                  static_cast<enum RBFShape<dim>::RBFBasisFunction>(
                    round(shape_arguments[2 * number_of_nodes + n_i]));
                nodes[n_i][0] = shape_arguments[3 * number_of_nodes + n_i];
                nodes[n_i][1] = shape_arguments[4 * number_of_nodes + n_i];
                nodes[n_i][2] = shape_arguments[5 * number_of_nodes + n_i];
              }
            particles[i].shape =
              std::make_shared<RBFShape<dim>>(support_radii,
                                              basis_functions,
                                              weights,
                                              nodes,
                                              particles[i].position,
                                              particles[i].orientation);
          }
      }
    else
      StandardExceptions::ExcNotImplemented();
  }

  void
  DynamicFlowControl::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("flow control");
    {
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether from the flow control information should be printed "
        "Choices are <quiet|verbose>.");
      prm.declare_entry("enable",
                        "false",
                        Patterns::Bool(),
                        "Enable flow rate control");
      prm.declare_entry("volumetric flow rate",
                        "0",
                        Patterns::Double(),
                        "Volumetric flow rate");
      prm.declare_entry("boundary id",
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
    }
    prm.leave_subsection();
  }

  void
  DynamicFlowControl::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("flow control");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;
      enable_flow_control = prm.get_bool("enable");
      flow_rate_0         = prm.get_double("volumetric flow rate");
      boundary_flow_id    = prm.get_integer("boundary id");
      flow_direction      = prm.get_integer("flow direction");
      beta_0              = prm.get_double("initial beta");
    }
    prm.leave_subsection();
  }

  template class Laser<2>;
  template class Laser<3>;
  template class IBParticles<2>;
  template class IBParticles<3>;
} // namespace Parameters
