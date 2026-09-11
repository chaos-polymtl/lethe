// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
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
  << " is not 1 (single phase simulation) or 2 (CLS/Cahn-Hilliard simulations). This is currently not supported.");

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

DeclException4(
  ListSizeIncoherentWithDeclaredNumber,
  std::string,
  unsigned int,
  std::string,
  unsigned int,
  << "'" << arg1 << "' was set to " << arg2 << " and the list '" << arg3
  << "' has a size of " << arg4
  << ". However, the list size should correspond to the number declared.");

DeclException4(ListsSizeMismatch,
               std::string,
               std::string,
               unsigned int,
               unsigned int,
               << "There is a size mismatch between 2 lists. The list '" << arg1
               << "' and the list '" << arg2 << "' have respectively sizes of "
               << arg3 << " and " << arg4
               << ". However, they should be of the same size.");

// DeclException

namespace Parameters
{
  namespace
  {
    // Reverse mappings for the shared Verbosity/FluidIndicator enums, used to
    // derive declare_entry's default-value string from the corresponding
    // struct member's own in-class default. Kept in sync by hand with the
    // string->enum chains in the various parse_parameters() below.
    //
    // NOTE: all per-struct `to_string(EnumType)` overloads added throughout
    // this file must also live directly inside `namespace Parameters` (in
    // their own `namespace { ... }` block is fine) rather than in the global
    // namespace: unqualified lookup stops at the first enclosing scope where
    // a `to_string` is found, so a struct-local overload declared inside
    // Parameters would otherwise hide these instead of overloading with them.
    std::string
    to_string(const Verbosity verbosity)
    {
      switch (verbosity)
        {
          case Verbosity::quiet:
            return "quiet";
          case Verbosity::verbose:
            return "verbose";
          case Verbosity::extra_verbose:
            return "extra verbose";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const FluidIndicator indicator)
    {
      switch (indicator)
        {
          case FluidIndicator::fluid0:
            return "fluid 0";
          case FluidIndicator::fluid1:
            return "fluid 1";
          case FluidIndicator::both:
            return "both";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  SizeOfSubsections
  get_size_of_subsections(const std::string &file_name,
                          const bool         require_subsection_size)
  {
    SizeOfSubsections sizes;
    sizes.boundary_conditions =
      get_max_subsection_size(file_name, require_subsection_size);
    sizes.manifolds =
      get_max_subsection_size(file_name, require_subsection_size);
    return sizes;
  }

  namespace
  {
    std::string
    to_string(const SimulationControl::TimeSteppingMethod method)
    {
      switch (method)
        {
          case SimulationControl::TimeSteppingMethod::steady:
            return "steady";
          case SimulationControl::TimeSteppingMethod::steady_bdf:
            return "steady_bdf";
          case SimulationControl::TimeSteppingMethod::bdf1:
            return "bdf1";
          case SimulationControl::TimeSteppingMethod::bdf2:
            return "bdf2";
          case SimulationControl::TimeSteppingMethod::bdf3:
            return "bdf3";
          case SimulationControl::TimeSteppingMethod::sdirk22:
            return "sdirk22";
          case SimulationControl::TimeSteppingMethod::sdirk33:
            return "sdirk33";
          case SimulationControl::TimeSteppingMethod::sdirk43:
            return "sdirk43";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const SimulationControl::BDFStartupMethods method)
    {
      switch (method)
        {
          case SimulationControl::BDFStartupMethods::multiple_step_bdf:
            return "multiple step bdf";
          case SimulationControl::BDFStartupMethods::initial_solution:
            return "initial solution";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const SimulationControl::EndControl control)
    {
      switch (control)
        {
          case SimulationControl::EndControl::iteration:
            return "iteration";
          case SimulationControl::EndControl::time:
            return "time";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const SimulationControl::OutputControl control)
    {
      switch (control)
        {
          case SimulationControl::OutputControl::iteration:
            return "iteration";
          case SimulationControl::OutputControl::time:
            return "time";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const std::vector<double> &values)
    {
      std::string result;
      for (unsigned int i = 0; i < values.size(); ++i)
        {
          if (i != 0)
            result += ", ";
          result += Patterns::Tools::Convert<double>::to_string(values[i]);
        }
      return result;
    }
  } // namespace

  void
  SimulationControl::declare_parameters(ParameterHandler &prm)
  {
    const SimulationControl defaults;
    prm.enter_subsection("simulation control");
    {
      prm.declare_entry(
        "method",
        to_string(defaults.method),
        Patterns::Selection(
          "steady|steady_bdf|bdf1|bdf2|bdf3|sdirk22|sdirk33|sdirk43"),
        "The time integration scheme. "
        "Choices are <steady|steady_bdf|bdf1|bdf2|bdf3|sdirk22|sdirk33|sdirk43>.");

      prm.declare_entry(
        "bdf startup method",
        to_string(defaults.bdf_startup_method),
        Patterns::Selection("multiple step bdf|initial solution"),
        "The kind of method used to startup high order bdf methods "
        "Choices are <multiple step bdf|initial solution>.");

      prm.declare_entry("time step",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.dt),
                        Patterns::Double(),
                        "Time step value");
      prm.declare_entry(
        "time end",
        Patterns::Tools::Convert<double>::to_string(defaults.time_end),
        Patterns::Double(),
        "Time value at which a transient simulation ends. Only used when the "
        "end control is set to time");
      prm.declare_entry(
        "iteration end",
        Patterns::Tools::Convert<unsigned int>::to_string(
          defaults.iteration_end),
        Patterns::Integer(0),
        "Transient iteration number at which a transient simulation ends. "
        "Only used when the end control is set to iteration");
      prm.declare_entry("startup time scaling",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.startup_timestep_scaling),
                        Patterns::Double(),
                        "Scaling factor used in the iterations necessary to "
                        "start-up the BDF schemes.");
      prm.declare_entry(
        "adapt time step to respect CFL",
        Patterns::Tools::Convert<bool>::to_string(defaults.adapt_with_cfl),
        Patterns::Bool(),
        "Adapt the time step to respect the maximum CFL condition. When multiple conditions are applied to the time step, this ensures that the CFL condition is also respected (Δt ≤ Δt_{CFL}). <true|false>");
      prm.declare_alias("adapt time step to respect CFL", "adapt", true);
      prm.declare_entry(
        "override time step on restart",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.override_time_step_on_restart),
        Patterns::Bool(),
        "Override checkpointed time step upon restart <true|false>");
      prm.declare_entry(
        "time step independent of end time",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.time_step_independent_of_end_time),
        Patterns::Bool(),
        "Ensures that the correct time step is kept when using adaptive time step simulations");
      prm.declare_entry("number mesh adapt",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.number_mesh_adaptation),
                        Patterns::Integer(),
                        "Number of mesh adaptation (for steady simulations)");
      prm.declare_entry("max cfl",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.maxCFL),
                        Patterns::Double(),
                        "Maximum CFL value");
      prm.declare_entry("max time step",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.max_dt),
                        Patterns::Double(),
                        "Maximum time step value");
      prm.declare_entry(
        "adapt time step to respect CTR",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.adapt_with_capillary_time_step_ratio),
        Patterns::Bool(),
        "By setting to 'true', it ensures that the imposed maximum capillary time-step ratio (CTR) is respected throughout the simulation (Δt ≤ Δt_{CTR}). <true|false>");
      prm.declare_entry(
        "max capillary time-step ratio",
        Patterns::Tools::Convert<double>::to_string(
          defaults.max_capillary_time_step_ratio),
        Patterns::Double(0),
        "The capillary time-step ratio (CTR) corresponds to the ratio of the time step over capillary time-step constraint (Δt/Δt_σ)");
      prm.declare_entry("stop tolerance",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.stop_tolerance),
                        Patterns::Double(),
                        "Tolerance at which the simulation is stopped");

      prm.declare_entry("adaptative time step scaling",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.adaptative_time_step_scaling),
                        Patterns::Double(),
                        "Adaptative time step scaling");
      prm.declare_entry("output path",
                        defaults.output_folder,
                        Patterns::FileName(),
                        "File output prefix");

      prm.declare_entry("output name",
                        defaults.output_name,
                        Patterns::FileName(),
                        "File output prefix");


      prm.declare_entry("output frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.output_iteration_frequency),
                        Patterns::Integer(),
                        "Output iteration frequency");

      prm.declare_entry("output time frequency",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.output_time_frequency),
                        Patterns::Double(),
                        "Output time frequency");

      prm.declare_entry(
        "output boundaries",
        Patterns::Tools::Convert<bool>::to_string(defaults.output_boundaries),
        Patterns::Bool(),
        "Output the boundaries of the domain along with their ID");

      prm.declare_entry("log frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.log_frequency),
                        Patterns::Integer(),
                        "log frequency");

      prm.declare_entry("log precision",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.log_precision),
                        Patterns::Integer(),
                        "Display precision when writing to log",
                        "This setting percolates to all output to the log");

      prm.declare_entry("output times",
                        to_string(defaults.output_times_vector),
                        Patterns::List(Patterns::Double()),
                        "List of specific output times separated with a comma");

      prm.declare_entry(
        "end control",
        to_string(defaults.end_control),
        Patterns::Selection("iteration|time"),
        "The control for the end of a transient simulation. The end "
        "condition is either a maximum time value (time end) or a maximum "
        "number of transient iterations (iteration end)");

      prm.declare_entry(
        "output control",
        to_string(defaults.output_control),
        Patterns::Selection("iteration|time"),
        "The control for the output of the simulation results"
        "Results can be either outputted at constant iteration frequency or at constant time");

      prm.declare_entry("output time interval",
                        to_string(defaults.output_time_interval),
                        Patterns::List(Patterns::Double()),
                        "Output files for a desired time interval");

      prm.declare_entry("subdivision",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.subdivision),
                        Patterns::Integer(),
                        "Subdivision of mesh cell in postprocessing");

      prm.declare_entry("group files",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.group_files),
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

      const std::string end_control_string = prm.get("end control");
      if (end_control_string == "iteration")
        end_control = EndControl::iteration;
      else if (end_control_string == "time")
        end_control = EndControl::time;
      else
        {
          AssertThrow(false, ExcMessage("Invalid end control scheme"));
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
      dt             = prm.get_double("time step");
      time_end       = prm.get_double("time end");
      iteration_end  = prm.get_integer("iteration end");
      adapt_with_cfl = prm.get_bool("adapt time step to respect CFL");
      time_step_independent_of_end_time =
        prm.get_bool("time step independent of end time");
      maxCFL = prm.get_double("max cfl");
      max_dt = prm.get_double("max time step");
      adapt_with_capillary_time_step_ratio =
        prm.get_bool("adapt time step to respect CTR");
      max_capillary_time_step_ratio =
        prm.get_double("max capillary time-step ratio");
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
      time_step_adaptation_required =
        adapt_with_cfl || adapt_with_capillary_time_step_ratio;
    }
    prm.leave_subsection();
  } // namespace Parameters

  namespace
  {
    std::string
    to_string(const Timer::Type type)
    {
      switch (type)
        {
          case Timer::Type::none:
            return "none";
          case Timer::Type::iteration:
            return "iteration";
          case Timer::Type::end:
            return "end";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  Timer::declare_parameters(ParameterHandler &prm)
  {
    const Timer defaults;
    prm.enter_subsection("timer");
    {
      prm.declare_entry("type",
                        to_string(defaults.type),
                        Patterns::Selection("none|iteration|end"),
                        "Clock monitoring methods "
                        "Choices are <none|iteration|end>.");
      prm.declare_entry(
        "write time in error table",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.write_time_in_error_table),
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
    const PowerLawParameters defaults;
    prm.enter_subsection("power-law");
    {
      prm.declare_entry("K",
                        Patterns::Tools::Convert<double>::to_string(defaults.K),
                        Patterns::Double(),
                        "Fluid consistency index");
      prm.declare_entry("n",
                        Patterns::Tools::Convert<double>::to_string(defaults.n),
                        Patterns::Double(),
                        "Flow behavior index");
      prm.declare_entry("shear rate min",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.shear_rate_min),
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
    const CarreauParameters defaults;
    prm.enter_subsection("carreau");
    {
      prm.declare_entry("viscosity_0",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.kinematic_viscosity_0),
                        Patterns::Double(),
                        "Kinematic viscosity at rest");
      prm.declare_entry("viscosity_inf",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.kinematic_viscosity_inf),
                        Patterns::Double(),
                        "Kinematic viscosity for an infinite shear rate");
      prm.declare_entry("lambda",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.lambda),
                        Patterns::Double(),
                        "Relaxation time");
      prm.declare_entry("a",
                        Patterns::Tools::Convert<double>::to_string(defaults.a),
                        Patterns::Double(),
                        "Carreau parameter");
      prm.declare_entry("n",
                        Patterns::Tools::Convert<double>::to_string(defaults.n),
                        Patterns::Double(),
                        "Power parameter");
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
    const ImmersedSolidTanhParameters defaults;
    prm.enter_subsection("immersed solid tanh");
    {
      prm.declare_entry("tracer diffusivity inside",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_diffusivity_inside),
                        Patterns::Double(),
                        "Tracer diffusivity inside the immersed solid");
      prm.declare_entry("tracer diffusivity outside",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_diffusivity_outside),
                        Patterns::Double(),
                        "Tracer diffusivity outside the immersed solid");
      prm.declare_entry("tracer reaction constant inside",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_reaction_constant_inside),
                        Patterns::Double(),
                        "Tracer reaction constant inside the immersed solid");
      prm.declare_entry("tracer reaction constant outside",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_reaction_constant_outside),
                        Patterns::Double(),
                        "Tracer reaction constant outside the immersed solid");
      prm.declare_entry("thickness",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.thickness),
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
    const ImmersedSolidGaussianParameters defaults;
    prm.enter_subsection("immersed solid gaussian");
    {
      prm.declare_entry("tracer diffusivity interface",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_diffusivity_interface),
                        Patterns::Double(),
                        "Tracer diffusivity at the immersed solid interface");
      prm.declare_entry("tracer diffusivity bulk",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_diffusivity_bulk),
                        Patterns::Double(),
                        "Tracer diffusivity in the phase bulk");
      prm.declare_entry(
        "tracer reaction constant interface",
        Patterns::Tools::Convert<double>::to_string(
          defaults.tracer_reaction_constant_interface),
        Patterns::Double(),
        "Tracer reaction constant at the immersed solid interface");
      prm.declare_entry("tracer reaction constant bulk",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.tracer_reaction_constant_bulk),
                        Patterns::Double(),
                        "Tracer reaction constant in the phase bulk");
      prm.declare_entry("thickness",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.thickness),
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
    // dry air's density, specific gas constant and normal temperature
    // (20 °C, 1 atm) as defaults
    const IsothermalIdealGasDensityParameters defaults;
    prm.enter_subsection("isothermal_ideal_gas");
    {
      prm.declare_entry(
        "density_ref",
        Patterns::Tools::Convert<double>::to_string(defaults.density_ref),
        Patterns::Double(),
        "Reference density of the gas in SI units for isothermal ideal gas equation of state in density calculation");

      prm.declare_entry(
        "R",
        Patterns::Tools::Convert<double>::to_string(defaults.R),
        Patterns::Double(),
        "Specific gas constant in SI units for isothermal ideal gas equation of state in density calculation");

      prm.declare_entry(
        "T",
        Patterns::Tools::Convert<double>::to_string(defaults.T),
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
    const SurfaceTensionParameters defaults;
    prm.declare_entry(
      "surface tension coefficient",
      Patterns::Tools::Convert<double>::to_string(
        defaults.surface_tension_coefficient),
      Patterns::Double(),
      "Surface tension coefficient for the corresponding pair of fluids or fluid-solid pair");
    prm.declare_entry(
      "reference state temperature",
      Patterns::Tools::Convert<double>::to_string(defaults.T_0),
      Patterns::Double(),
      "Temperature of the reference state corresponding to the surface tension coefficient");
    prm.declare_entry(
      "temperature-driven surface tension gradient",
      Patterns::Tools::Convert<double>::to_string(
        defaults.surface_tension_gradient),
      Patterns::Double(),
      "Surface tension gradient with respect to the temperature for the corresponding pair of fluids or fluid-solid pair");
    prm.declare_entry(
      "solidus temperature",
      Patterns::Tools::Convert<double>::to_string(defaults.T_solidus),
      Patterns::Double(),
      "Temperature of the solidus for the corresponding pair of fluids or fluid-solid pair");
    prm.declare_entry(
      "liquidus temperature",
      Patterns::Tools::Convert<double>::to_string(defaults.T_liquidus),
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
    const MobilityCahnHilliardParameters defaults;
    prm.declare_entry(
      "cahn hilliard mobility constant",
      Patterns::Tools::Convert<double>::to_string(
        defaults.mobility_cahn_hilliard_constant),
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
        Patterns::Tools::Convert<bool>::to_string(this->enable),
        Patterns::Bool(),
        "Enable/disable (true/false) the solid domain constraining feature.");

      prm.declare_entry(
        "enable domain restriction with plane",
        Patterns::Tools::Convert<bool>::to_string(
          this->enable_domain_restriction_with_plane),
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
        Patterns::Tools::Convert<unsigned int>::to_string(
          this->number_of_constraints),
        Patterns::Integer(),
        "Number of solid constraints (maximum of 1 per fluid).");

      // Resize vectors
      this->fluid_ids.resize(max_number_of_constraints);
      this->filtered_phase_indicator_tolerance.resize(
        max_number_of_constraints);
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
    prm.declare_entry("phase indicator tolerance",
                      "1e-4",
                      Patterns::Double(),
                      "Absolute filtered phase indicator tolerance used in "
                      "conjunction with CLS simulations to select the cells "
                      "on which the constraint is applied.");
    prm.declare_alias("phase indicator tolerance",
                      "phase fraction tolerance",
                      true);
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
      this->filtered_phase_indicator_tolerance.resize(number_of_constraints);
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
    this->filtered_phase_indicator_tolerance[constraint_id] =
      prm.get_double("phase indicator tolerance");
    this->temperature_min_values[constraint_id] =
      prm.get_double("min temperature");
    this->temperature_max_values[constraint_id] =
      prm.get_double("max temperature");
  }

  namespace
  {
    std::string
    to_string(const Stabilization::NavierStokesStabilization type)
    {
      switch (type)
        {
          case Stabilization::NavierStokesStabilization::pspg_supg:
            return "pspg_supg";
          case Stabilization::NavierStokesStabilization::gls:
            return "gls";
          case Stabilization::NavierStokesStabilization::grad_div:
            return "grad_div";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Stabilization::ScalarLimiters type)
    {
      switch (type)
        {
          case Stabilization::ScalarLimiters::none:
            return "none";
          case Stabilization::ScalarLimiters::moe:
            return "moe";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  Stabilization::declare_parameters(ParameterHandler &prm)
  {
    const Stabilization defaults;
    prm.enter_subsection("stabilization");
    {
      prm.declare_entry(
        "use default stabilization",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.use_default_stabilization),
        Patterns::Bool(),
        "Use the default stabilization method provided by the solver");
      prm.declare_entry(
        "stabilization",
        to_string(defaults.stabilization),
        Patterns::Selection("pspg_supg|gls|grad_div"),
        "Type of stabilization used for the Navier-Stokes equations"
        "Choices are <pspg_supg|gls|grad_div>.");

      prm.declare_entry(
        "heat transfer dcdd stabilization",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.heat_transfer_dcdd_stabilization),
        Patterns::Bool(),
        "Apply Discontinuity-Capturing Directional Dissipation (DCDD) "
        "stabilization term on heat transfer <true|false>");

      prm.declare_entry(
        "cls dcdd stabilization",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.cls_dcdd_stabilization),
        Patterns::Bool(),
        "Apply Discontinuity-Capturing Directional Dissipation (DCDD) "
        "stabilization term on the CLS phase indicator <true|false>");
      prm.declare_alias("cls dcdd stabilization",
                        "vof dcdd stabilization",
                        true);

      prm.declare_entry(
        "cls dcdd diffusion factor",
        Patterns::Tools::Convert<double>::to_string(
          defaults.dcdd_diffusion_coeff),
        Patterns::Double(),
        "Diffusion factor scaling the DCDD stabilization term in the CLS "
        "equation");
      prm.declare_alias("cls dcdd diffusion factor",
                        "vof dcdd diffusion factor",
                        true);

      prm.declare_entry(
        "pressure scaling factor",
        Patterns::Tools::Convert<double>::to_string(
          defaults.pressure_scaling_factor),
        Patterns::Double(),
        "This parameter can be used to change the scale of pressure in the "
        "Navier-Stokes equations. When the velocity and pressure scales are very "
        "different, using this parameter allows to reduce the condition number"
        " and reach a solution.");

      prm.declare_entry(
        "scalar limiter",
        to_string(defaults.scalar_limiter),
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
      cls_dcdd_stabilization = prm.get_bool("cls dcdd stabilization");
      dcdd_diffusion_coeff   = prm.get_double("cls dcdd diffusion factor");

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
    const PhaseChange defaults;
    prm.enter_subsection("phase change");
    {
      prm.declare_entry("solidus temperature",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.T_solidus),
                        Patterns::Double(),
                        "Temperature of the solidus");
      prm.declare_entry("liquidus temperature",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.T_liquidus),
                        Patterns::Double(),
                        "Temperature of the liquidus");
      prm.declare_entry("latent enthalpy",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.latent_enthalpy),
                        Patterns::Double(),
                        "Enthalpy of the phase change");

      prm.declare_entry("specific heat liquid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.cp_l),
                        Patterns::Double(),
                        "Specific heat of the liquid phase");

      prm.declare_entry("specific heat solid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.cp_s),
                        Patterns::Double(),
                        "Specific heat of the solid phase");

      prm.declare_entry("thermal conductivity liquid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.thermal_conductivity_l),
                        Patterns::Double(),
                        "Thermal conductivity of the liquid phase");

      prm.declare_entry("thermal conductivity solid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.thermal_conductivity_s),
                        Patterns::Double(),
                        "Thermal conductivity of the solid phase");

      prm.declare_entry("thermal expansion liquid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.thermal_expansion_l),
                        Patterns::Double(),
                        "Thermal expansion coefficient of the liquid phase");

      prm.declare_entry("thermal expansion solid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.thermal_expansion_s),
                        Patterns::Double(),
                        "Thermal expansion coefficient of the solid phase");

      prm.declare_entry("viscosity liquid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.kinematic_viscosity_l),
                        Patterns::Double(),
                        "Kinematic viscosity of the liquid phase");

      prm.declare_entry("viscosity solid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.kinematic_viscosity_s),
                        Patterns::Double(),
                        "Kinematic viscosity of the solid phase");

      prm.declare_entry("Darcy penalty liquid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.penalty_l),
                        Patterns::Double(),
                        "Darcy penalty of the liquid phase");

      prm.declare_entry("Darcy penalty solid",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.penalty_s),
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
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          number_of_fluids),
                        Patterns::Integer(),
                        "Number of fluids");

      // Definition of the fluids in the simulation
      for (unsigned int i_fluid = 0; i_fluid < max_fluids; ++i_fluid)
        {
          fluids[i_fluid].declare_parameters(prm, "fluid", i_fluid);
        }

      prm.declare_entry("number of solids",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          number_of_solids),
                        Patterns::Integer(),
                        "Number of solids");

      // Definition of the solids in the simulation
      for (unsigned int i_solid = 0; i_solid < max_solids; ++i_solid)
        {
          solids[i_solid].declare_parameters(prm, "solid", i_solid);
        }

      prm.declare_entry(
        "reference temperature",
        Patterns::Tools::Convert<double>::to_string(reference_temperature),
        Patterns::Double(),
        "Reference temperature used for the calculation of physical properties and thermal expansion");
    }

    // Definition of interactions between materials
    prm.declare_entry(
      "number of material interactions",
      Patterns::Tools::Convert<unsigned int>::to_string(
        number_of_material_interactions),
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

  namespace
  {
    std::string
    to_string(const Material::RheologicalModel model)
    {
      switch (model)
        {
          case Material::RheologicalModel::powerlaw:
            return "power-law";
          case Material::RheologicalModel::carreau:
            return "carreau";
          case Material::RheologicalModel::newtonian:
            return "newtonian";
          case Material::RheologicalModel::phase_change:
            return "phase_change";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::DensityModel model)
    {
      switch (model)
        {
          case Material::DensityModel::constant:
            return "constant";
          case Material::DensityModel::isothermal_ideal_gas:
            return "isothermal_ideal_gas";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::SpecificHeatModel model)
    {
      switch (model)
        {
          case Material::SpecificHeatModel::constant:
            return "constant";
          case Material::SpecificHeatModel::phase_change:
            return "phase_change";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::ThermalConductivityModel model)
    {
      switch (model)
        {
          case Material::ThermalConductivityModel::constant:
            return "constant";
          case Material::ThermalConductivityModel::linear:
            return "linear";
          case Material::ThermalConductivityModel::phase_change:
            return "phase_change";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::ThermalExpansionModel model)
    {
      switch (model)
        {
          case Material::ThermalExpansionModel::constant:
            return "constant";
          case Material::ThermalExpansionModel::phase_change:
            return "phase_change";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::TracerDiffusivityModel model)
    {
      switch (model)
        {
          case Material::TracerDiffusivityModel::constant:
            return "constant";
          case Material::TracerDiffusivityModel::immersed_boundary_tanh:
            return "immersed solid tanh";
          case Material::TracerDiffusivityModel::immersed_boundary_gaussian:
            return "immersed solid gaussian";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::TracerReactionPrefactorModel model)
    {
      switch (model)
        {
          case Material::TracerReactionPrefactorModel::none:
            return "none";
          case Material::TracerReactionPrefactorModel::constant:
            return "constant";
          case Material::TracerReactionPrefactorModel::immersed_boundary_tanh:
            return "immersed solid tanh";
          case Material::TracerReactionPrefactorModel::
            immersed_boundary_gaussian:
            return "immersed solid gaussian";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::ElectricConductivityModel model)
    {
      switch (model)
        {
          case Material::ElectricConductivityModel::constant:
            return "constant";
          case Material::ElectricConductivityModel::polynomial:
            return "polynomial";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::ElectricPermittivityModel model)
    {
      switch (model)
        {
          case Material::ElectricPermittivityModel::constant:
            return "constant";
          case Material::ElectricPermittivityModel::polynomial:
            return "polynomial";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Material::MagneticPermeabilityModel model)
    {
      switch (model)
        {
          case Material::MagneticPermeabilityModel::constant:
            return "constant";
          case Material::MagneticPermeabilityModel::polynomial:
            return "polynomial";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  Material::declare_parameters(ParameterHandler  &prm,
                               const std::string &material_prefix,
                               const unsigned int id) const
  {
    prm.enter_subsection(material_prefix + " " +
                         Utilities::int_to_string(id, 1));
    {
      prm.declare_entry("density",
                        Patterns::Tools::Convert<double>::to_string(density),
                        Patterns::Double(),
                        "Density for the fluid corresponding to Phase = " +
                          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "kinematic viscosity",
        Patterns::Tools::Convert<double>::to_string(kinematic_viscosity),
        Patterns::Double(),
        "Kinematic viscosity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "specific heat",
        Patterns::Tools::Convert<double>::to_string(specific_heat),
        Patterns::Double(),
        "Specific heat for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "thermal conductivity",
        Patterns::Tools::Convert<double>::to_string(thermal_conductivity),
        Patterns::Double(),
        "Thermal conductivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "thermal expansion",
        Patterns::Tools::Convert<double>::to_string(thermal_expansion),
        Patterns::Double(),
        "Thermal expansion coefficient for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer diffusivity model",
        to_string(tracer_diffusivity_model),
        Patterns::Selection(
          "constant|immersed solid tanh|immersed solid gaussian"),
        "Model used for the calculation of the tracer diffusivity"
        "Choices are <constant|immersed solid tanh|immersed solid gaussian>.");

      prm.declare_entry(
        "tracer diffusivity",
        Patterns::Tools::Convert<double>::to_string(tracer_diffusivity),
        Patterns::Double(),
        "Tracer diffusivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer reaction constant model",
        to_string(tracer_reaction_prefactor_model),
        Patterns::Selection(
          "none|constant|immersed solid tanh|immersed solid gaussian"),
        "Model used for the calculation of the tracer reaction constant"
        "Choices are <none|constant|immersed solid tanh|immersed solid gaussian>.");

      prm.declare_entry(
        "tracer reaction constant",
        Patterns::Tools::Convert<double>::to_string(tracer_reaction_constant),
        Patterns::Double(),
        "Tracer reaction constant for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer reaction order",
        Patterns::Tools::Convert<double>::to_string(tracer_reaction_order),
        Patterns::Double(),
        "Tracer reaction order for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "tracer reaction threshold",
        Patterns::Tools::Convert<double>::to_string(tracer_reaction_threshold),
        Patterns::Double(),
        "Tracer reaction threshold for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1) + ". " +
          "Lower values are more realistic but lead to stiffer linear systems when n < 1. ");

      // Declaration of the immersed solids models parameters
      immersed_solid_tanh_parameters.declare_parameters(prm);
      immersed_solid_gaussian_parameters.declare_parameters(prm);

      prm.declare_entry(
        "rheological model",
        to_string(rheological_model),
        Patterns::Selection("newtonian|power-law|carreau|phase_change"),
        "Rheological model "
        "Choices are <newtonian|power-law|carreau|phase_change>.");

      non_newtonian_parameters.declare_parameters(prm);


      prm.declare_entry("density model",
                        to_string(density_model),
                        Patterns::Selection("constant|isothermal_ideal_gas"),
                        "Model used for the calculation of the density"
                        "Choices are <constant|isothermal_ideal_gas>.");

      isothermal_ideal_gas_density_parameters.declare_parameters(prm);

      prm.declare_entry("specific heat model",
                        to_string(specific_heat_model),
                        Patterns::Selection("constant|phase_change"),
                        "Model used for the calculation of the specific heat"
                        "Choices are <constant|phase_change>.");

      phase_change_parameters.declare_parameters(prm);

      prm.declare_entry(
        "thermal conductivity model",
        to_string(thermal_conductivity_model),
        Patterns::Selection("constant|linear|phase_change"),
        "Model used for the calculation of the thermal conductivity"
        "Choices are <constant|linear|phase_change>.");

      prm.declare_entry(
        "thermal expansion model",
        to_string(thermal_expansion_model),
        Patterns::Selection("constant|phase_change"),
        "Model used for the calculation of the thermal expansion coefficient"
        "Choices are <constant|phase_change>.");

      prm.declare_entry("k_A0",
                        Patterns::Tools::Convert<double>::to_string(k_A0),
                        Patterns::Double(),
                        "k_A0 parameter for linear conductivity model");

      prm.declare_entry("k_A1",
                        Patterns::Tools::Convert<double>::to_string(k_A1),
                        Patterns::Double(),
                        "k_A1 parameter for linear conductivity model");

      // ----------------------------------
      // Electromagnetic properties
      // ----------------------------------
      prm.declare_entry(
        "electric conductivity model",
        to_string(electric_conductivity_model),
        Patterns::Selection("constant|polynomial"),
        "Model used for the calculation of the electric conductivity"
        "Choices are <constant|polynomial>.");
      prm.declare_entry(
        "electric conductivity",
        Patterns::Tools::Convert<double>::to_string(electric_conductivity),
        Patterns::Double(),
        "Electric conductivity for the material corresponding to: " +
          material_prefix + " " + Utilities::int_to_string(id, 1));

      prm.declare_entry(
        "electric conductivity polynomial coefficients",
        "0",
        Patterns::List(Patterns::Double(), 1),
        "Coefficients of the polynomial model for the electric conductivity of the material corresponding to: " +
          material_prefix + " " + Utilities::int_to_string(id, 1) +
          ". The coefficients are given in decreasing order of the polynomial degree and need to be separated by commas.");

      prm.declare_entry(
        "electric permittivity model",
        to_string(electric_permittivity_model),
        Patterns::Selection("constant|polynomial"),
        "Model used for the calculation of the electric permittivity"
        "Choices are <constant|polynomial>.");
      prm.declare_entry(
        "electric permittivity real part",
        Patterns::Tools::Convert<double>::to_string(electric_permittivity_real),
        Patterns::Double(),
        "Real part of the electric permittivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "electric permittivity imag part",
        Patterns::Tools::Convert<double>::to_string(electric_permittivity_imag),
        Patterns::Double(),
        "Imaginary part of the electric permittivity for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "electric permittivity real part polynomial coefficients",
        "1",
        Patterns::List(Patterns::Double(), 1),
        "Coefficients of the polynomial model for the electric permittivity of the material corresponding to: " +
          material_prefix + " " + Utilities::int_to_string(id, 1) +
          ". The coefficients are given in decreasing order of the polynomial degree and need to be separated by commas.");

      prm.declare_entry(
        "electric permittivity imag part polynomial coefficients",
        "0",
        Patterns::List(Patterns::Double(), 1),
        "Coefficients of the polynomial model for the electric permittivity of the material corresponding to: " +
          material_prefix + " " + Utilities::int_to_string(id, 1) +
          ". The coefficients are given in decreasing order of the polynomial degree and need to be separated by commas.");

      prm.declare_entry(
        "magnetic permeability model",
        to_string(magnetic_permeability_model),
        Patterns::Selection("constant|polynomial"),
        "Model used for the calculation of the magnetic permeability"
        "Choices are <constant|polynomial>.");
      prm.declare_entry(
        "magnetic permeability real part",
        Patterns::Tools::Convert<double>::to_string(magnetic_permeability_real),
        Patterns::Double(),
        "Real part of the magnetic permeability for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "magnetic permeability imag part",
        Patterns::Tools::Convert<double>::to_string(magnetic_permeability_imag),
        Patterns::Double(),
        "Imaginary part of the magnetic permeability for the fluid corresponding to Phase = " +
          Utilities::int_to_string(id, 1));
      prm.declare_entry(
        "magnetic permeability real part polynomial coefficients",
        "1",
        Patterns::List(Patterns::Double(), 1),
        "Coefficients of the polynomial model for the magnetic permeability of the material corresponding to: " +
          material_prefix + " " + Utilities::int_to_string(id, 1) +
          ". The coefficients are given in decreasing order of the polynomial degree and need to be separated by commas.");
      prm.declare_entry(
        "magnetic permeability imag part polynomial coefficients",
        "0",
        Patterns::List(Patterns::Double(), 1),
        "Coefficients of the polynomial model for the magnetic permeability of the material corresponding to: " +
          material_prefix + " " + Utilities::int_to_string(id, 1) +
          ". The coefficients are given in decreasing order of the polynomial degree and need to be separated by commas.");
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

      //--------------------------------
      // Electromagnetic properties
      //--------------------------------
      op = prm.get("electric conductivity model");
      if (op == "constant")
        {
          electric_conductivity_model = ElectricConductivityModel::constant;
          electric_conductivity       = prm.get_double("electric conductivity");
        }
      else if (op == "polynomial")
        {
          electric_conductivity_model = ElectricConductivityModel::polynomial;

          const std::string electric_conductivity_coefficients_string =
            prm.get("electric conductivity polynomial coefficients");
          const std::vector<std::string>
            electric_conductivity_coefficients_vec =
              Utilities::split_string_list(
                electric_conductivity_coefficients_string);

          for (const std::string &coefficient :
               electric_conductivity_coefficients_vec)
            {
              electric_conductivity_polynomial_coefficients.push_back(
                std::stod(coefficient));
            }
        }

      op = prm.get("electric permittivity model");
      if (op == "constant")
        {
          electric_permittivity_model = ElectricPermittivityModel::constant;
          electric_permittivity_real =
            prm.get_double("electric permittivity real part");
          electric_permittivity_imag =
            prm.get_double("electric permittivity imag part");
        }
      else if (op == "polynomial")
        {
          electric_permittivity_model = ElectricPermittivityModel::polynomial;

          const std::string electric_permittivity_real_coefficients_string =
            prm.get("electric permittivity real part polynomial coefficients");
          const std::vector<std::string>
            electric_permittivity_real_coefficients_vec =
              Utilities::split_string_list(
                electric_permittivity_real_coefficients_string);

          for (const std::string &coefficient :
               electric_permittivity_real_coefficients_vec)
            {
              electric_permittivity_real_polynomial_coefficients.push_back(
                std::stod(coefficient));
            }

          const std::string electric_permittivity_imag_coefficients_string =
            prm.get("electric permittivity imag part polynomial coefficients");
          const std::vector<std::string>
            electric_permittivity_imag_coefficients_vec =
              Utilities::split_string_list(
                electric_permittivity_imag_coefficients_string);

          for (const std::string &coefficient :
               electric_permittivity_imag_coefficients_vec)
            {
              electric_permittivity_imag_polynomial_coefficients.push_back(
                std::stod(coefficient));
            }
        }
      op = prm.get("magnetic permeability model");
      if (op == "constant")
        {
          magnetic_permeability_model = MagneticPermeabilityModel::constant;
          magnetic_permeability_real =
            prm.get_double("magnetic permeability real part");
          magnetic_permeability_imag =
            prm.get_double("magnetic permeability imag part");
        }
      else if (op == "polynomial")
        {
          magnetic_permeability_model = MagneticPermeabilityModel::polynomial;

          const std::string magnetic_permeability_real_coefficients_string =
            prm.get("magnetic permeability real part polynomial coefficients");
          const std::vector<std::string>
            magnetic_permeability_real_coefficients_vec =
              Utilities::split_string_list(
                magnetic_permeability_real_coefficients_string);
          for (const std::string &coefficient :
               magnetic_permeability_real_coefficients_vec)
            {
              magnetic_permeability_real_polynomial_coefficients.push_back(
                std::stod(coefficient));
            }

          const std::string magnetic_permeability_imag_coefficients_string =
            prm.get("magnetic permeability imag part polynomial coefficients");
          const std::vector<std::string>
            magnetic_permeability_imag_coefficients_vec =
              Utilities::split_string_list(
                magnetic_permeability_imag_coefficients_string);
          for (const std::string &coefficient :
               magnetic_permeability_imag_coefficients_vec)
            {
              magnetic_permeability_imag_polynomial_coefficients.push_back(
                std::stod(coefficient));
            }
        }
    }
    prm.leave_subsection();
  }

  namespace
  {
    std::string
    to_string(const MaterialInteractions::MaterialInteractionsType type)
    {
      switch (type)
        {
          case MaterialInteractions::MaterialInteractionsType::fluid_fluid:
            return "fluid-fluid";
          case MaterialInteractions::MaterialInteractionsType::fluid_solid:
            return "fluid-solid";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const MaterialInteractions::SurfaceTensionModel model)
    {
      switch (model)
        {
          case MaterialInteractions::SurfaceTensionModel::constant:
            return "constant";
          case MaterialInteractions::SurfaceTensionModel::linear:
            return "linear";
          case MaterialInteractions::SurfaceTensionModel::phase_change:
            return "phase change";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const MaterialInteractions::MobilityCahnHilliardModel model)
    {
      switch (model)
        {
          case MaterialInteractions::MobilityCahnHilliardModel::constant:
            return "constant";
          case MaterialInteractions::MobilityCahnHilliardModel::quartic:
            return "quartic";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  MaterialInteractions::declare_parameters(ParameterHandler  &prm,
                                           const unsigned int id) const
  {
    prm.enter_subsection("material interaction " +
                         Utilities::int_to_string(id, 1));
    {
      prm.declare_entry(
        "type",
        to_string(material_interaction_type),
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
          to_string(surface_tension_model),
          Patterns::Selection("constant|linear|phase change"),
          "Model used for the calculation of the surface tension coefficient\n"
          "The choices are <constant|linear|phase change>.");
        surface_tension_parameters.declare_parameters(prm);

        // Cahn-Hilliard mobility
        prm.declare_entry(
          "cahn hilliard mobility model",
          to_string(mobility_cahn_hilliard_model),
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
          to_string(surface_tension_model),
          Patterns::Selection("constant|linear|phase change"),
          "Model used for the calculation of the surface tension coefficient\n"
          "The choices are <constant|linear|phase change>.");
        surface_tension_parameters.declare_parameters(prm);

        // Cahn-Hilliard mobility
        prm.declare_entry(
          "cahn hilliard mobility model",
          to_string(mobility_cahn_hilliard_model),
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
    const FEM defaults;
    prm.enter_subsection("FEM");
    {
      prm.declare_entry("velocity degree",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.velocity_degree),
                        Patterns::Integer(0),
                        "interpolation degree for velocity");
      prm.declare_alias("velocity degree", "velocity order", true);
      prm.declare_entry("pressure degree",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.pressure_degree),
                        Patterns::Integer(0),
                        "interpolation degree for pressure");
      prm.declare_alias("pressure degree", "pressure order", true);
      prm.declare_entry("void fraction degree",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.void_fraction_degree),
                        Patterns::Integer(0),
                        "interpolation degree for void fraction");
      prm.declare_alias("void fraction degree", "void fraction order", true);
      prm.declare_entry("temperature degree",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.temperature_degree),
                        Patterns::Integer(0),
                        "interpolation degree for temperature");
      prm.declare_alias("temperature degree", "temperature order", true);
      prm.declare_entry("tracer degree",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.tracer_degree),
                        Patterns::Integer(0),
                        "interpolation degree for tracer");
      prm.declare_alias("tracer degree", "tracer order", true);
      prm.declare_entry("cls degree",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.CLS_degree),
                        Patterns::Integer(0),
                        "interpolation degree for cls");
      prm.declare_alias("cls degree", "VOF degree");
      prm.declare_alias("cls degree", "cls order", true);
      prm.declare_alias("cls degree", "VOF order", true);
      prm.declare_entry(
        "phase cahn hilliard degree",
        Patterns::Tools::Convert<unsigned int>::to_string(
          defaults.phase_cahn_hilliard_degree),
        Patterns::Integer(),
        "interpolation degree phase parameter for the Cahn-Hilliard equations");
      prm.declare_alias("phase cahn hilliard degree",
                        "phase cahn hilliard order",
                        true);
      prm.declare_entry(
        "potential cahn hilliard degree",
        Patterns::Tools::Convert<unsigned int>::to_string(
          defaults.potential_cahn_hilliard_degree),
        Patterns::Integer(),
        "interpolation degree chemical potential for the Cahn-Hilliard equations");
      prm.declare_alias("potential cahn hilliard degree",
                        "potential cahn hilliard order",
                        true);
      prm.declare_entry(
        "electromagnetics trial degree",
        Patterns::Tools::Convert<unsigned int>::to_string(
          defaults.electromagnetics_trial_degree),
        Patterns::Integer(),
        "interpolation degree for the trial space of the electromagnetics physics (time-harmonic Maxwell equations).");
      prm.declare_alias("electromagnetics trial degree",
                        "electromagnetics trial order",
                        true);
      prm.declare_entry(
        "electromagnetics test degree",
        Patterns::Tools::Convert<unsigned int>::to_string(
          defaults.electromagnetics_test_degree),
        Patterns::Integer(),
        "interpolation degree for the test space of the electromagnetics physics (time-harmonic Maxwell equations).");
      prm.declare_alias("electromagnetics test degree",
                        "electromagnetics test order",
                        true);

      prm.declare_entry(
        "tracer uses dg",
        Patterns::Tools::Convert<bool>::to_string(defaults.tracer_uses_dg),
        Patterns::Bool(),
        "Switch tracer to Discontinuous Galerkin (DG) formulation");

      prm.declare_entry(
        "cls uses dg",
        Patterns::Tools::Convert<bool>::to_string(defaults.CLS_uses_dg),
        Patterns::Bool(),
        "Switch CLS to Discontinuous Galerkin (DG) formulation");
      prm.declare_alias("cls uses dg", "VOF uses dg", true);

      prm.declare_entry(
        "enable bubble function velocity",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.enable_bubble_function_velocity),
        Patterns::Bool(),
        "Enable bubble enrichment function for the velocity field");

      prm.declare_entry(
        "enable bubble function pressure",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.enable_bubble_function_pressure),
        Patterns::Bool(),
        "Enable bubble enrichment function for the pressure field");
    }
    prm.leave_subsection();
  }

  void
  FEM::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("FEM");
    {
      velocity_degree      = prm.get_integer("velocity degree");
      pressure_degree      = prm.get_integer("pressure degree");
      void_fraction_degree = prm.get_integer("void fraction degree");
      temperature_degree   = prm.get_integer("temperature degree");
      tracer_degree        = prm.get_integer("tracer degree");
      tracer_uses_dg       = prm.get_bool("tracer uses dg");
      CLS_degree           = prm.get_integer("cls degree");
      CLS_uses_dg          = prm.get_bool("cls uses dg");
      phase_cahn_hilliard_degree =
        prm.get_integer("phase cahn hilliard degree");
      potential_cahn_hilliard_degree =
        prm.get_integer("potential cahn hilliard degree");
      electromagnetics_trial_degree =
        prm.get_integer("electromagnetics trial degree");
      electromagnetics_test_degree =
        prm.get_integer("electromagnetics test degree");
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
    const Forces defaults;
    prm.enter_subsection("forces");
    {
      prm.declare_entry(
        "verbosity",
        to_string(defaults.verbosity),
        Patterns::Selection("quiet|verbose"),
        "State whether from the non-linear solver should be printed "
        "Choices are <quiet|verbose>.");
      prm.declare_entry("calculate force",
                        Patterns::Tools::Convert<bool>::to_string(
                          defaults.calculate_force),
                        Patterns::Bool(),
                        "Enable calculation of force");
      prm.declare_entry("calculate torque",
                        Patterns::Tools::Convert<bool>::to_string(
                          defaults.calculate_torque),
                        Patterns::Bool(),
                        "Enable calculation of torque");
      prm.declare_entry("force name",
                        defaults.force_output_name,
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("torque name",
                        defaults.torque_output_name,
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("output precision",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.output_precision),
                        Patterns::Integer(),
                        "Precision of the values outputted.");
      prm.declare_entry("calculation frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.calculation_frequency),
                        Patterns::Integer(),
                        "Calculation frequency");
      prm.declare_entry("output frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.output_frequency),
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
    const Laser_FreeSurfaceRadiation defaults;
    prm.enter_subsection("free surface radiation");
    {
      prm.declare_entry(
        "enable",
        Patterns::Tools::Convert<bool>::to_string(defaults.enable_radiation),
        Patterns::Bool(),
        "Enable radiation at the free surface (air/metal interface) <true|false>");
      prm.declare_entry("Stefan-Boltzmann constant",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.Stefan_Boltzmann_constant),
                        Patterns::Double(),
                        "Stefan-Boltzmann constant");
      prm.declare_entry("emissivity",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.emissivity),
                        Patterns::Double(),
                        "Emissivity of the free surface (air/metal interface)");
      prm.declare_entry(
        "Tinf",
        Patterns::Tools::Convert<double>::to_string(defaults.Tinf),
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

  namespace
  {
    template <int dim>
    std::string
    to_string(const typename Laser<dim>::LaserType type)
    {
      switch (type)
        {
          case Laser<dim>::LaserType::exponential_decay:
            return "exponential_decay";
          case Laser<dim>::LaserType::gaussian_heat_flux_cls_interface:
            return "gaussian_heat_flux_cls_interface";
          case Laser<dim>::LaserType::uniform_heat_flux_cls_interface:
            return "uniform_heat_flux_cls_interface";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  template <int dim>
  void
  Laser<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("laser parameters");
    {
      prm.declare_entry("enable",
                        Patterns::Tools::Convert<bool>::to_string(
                          activate_laser),
                        Patterns::Bool(),
                        "Activate laser");
      prm.declare_entry(
        "type",
        to_string<dim>(laser_type),
        Patterns::Selection(
          "exponential_decay|gaussian_heat_flux_cls_interface|uniform_heat_flux_cls_interface"),
        "Type of laser model used."
        "Choices are <exponential_decay|gaussian_heat_flux_cls_interface|uniform_heat_flux_cls_interface>.");
      prm.declare_entry(
        "enable angle of incidence dependence",
        Patterns::Tools::Convert<bool>::to_string(
          enable_angle_of_incidence_dependence),
        Patterns::Bool(),
        "Enable the multiplication of the laser heat flux by the cosine of the angle of incidence of the laser with respect to the surface.");
      prm.declare_entry("concentration factor",
                        Patterns::Tools::Convert<double>::to_string(
                          concentration_factor),
                        Patterns::Double(),
                        "Concentration factor");
      prm.declare_entry("power",
                        Patterns::Tools::Convert<double>::to_string(
                          laser_power),
                        Patterns::Double(),
                        "Laser power");
      prm.declare_entry("absorptivity",
                        Patterns::Tools::Convert<double>::to_string(
                          laser_absorptivity),
                        Patterns::Double(),
                        "Laser absorptivity");
      prm.declare_entry("penetration depth",
                        Patterns::Tools::Convert<double>::to_string(
                          penetration_depth),
                        Patterns::Double(),
                        "Penetration depth");
      prm.declare_entry("beam radius",
                        Patterns::Tools::Convert<double>::to_string(
                          beam_radius),
                        Patterns::Double(),
                        "Laser beam radius");
      radiation.declare_parameters(prm);

      prm.enter_subsection("path");
      laser_scan_path = std::make_shared<Functions::ParsedFunction<dim>>(dim);
      laser_scan_path->declare_parameters(prm, dim);
      prm.leave_subsection();

      prm.declare_entry("start time",
                        Patterns::Tools::Convert<double>::to_string(start_time),
                        Patterns::Double(),
                        "Start time of laser");
      prm.declare_entry("end time",
                        Patterns::Tools::Convert<double>::to_string(end_time),
                        Patterns::Double(),
                        "End time of laser");

      // Not derived from a member default: see the comments on
      // beam_orientation and rotation_axis in the header (dim-dependent).
      prm.declare_entry("beam orientation",
                        "z-",
                        Patterns::Selection("x+|x-|y+|y-|z+|z-"),
                        "Laser beam orientation "
                        "Choices are <x+|x-|y+|y-|z+|z->.");

      prm.declare_entry(
        "beam rotation angle",
        Patterns::Tools::Convert<double>::to_string(rotation_angle),
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
      else if (type_string == "gaussian_heat_flux_cls_interface")
        laser_type = LaserType::gaussian_heat_flux_cls_interface;
      else
        laser_type = LaserType::uniform_heat_flux_cls_interface;
      enable_angle_of_incidence_dependence =
        prm.get_bool("enable angle of incidence dependence");
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

  template <int dim>
  void
  PostProcessing<dim>::IsocontourBoundingBoxes::declare_parameters(
    ParameterHandler &prm)
  {
    prm.enter_subsection("isocontour bounding box");
    {
      prm.declare_entry("number of isocontour bounding boxes",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          number_of_isocontour_bounding_boxes),
                        Patterns::Integer(0),
                        "Number of monitored isocontours");
      prm.declare_entry(
        "variable",
        "temperature",
        Patterns::List(Patterns::Selection("temperature|phase")),
        "Variable(s) of monitored isocontour(s). "
        "Choices are <temperature|phase>. "
        "Each entry corresponds to a different isocontour. "
        "When multiple isocontours are defined, the different variables must "
        "be separated by commas (e.g., ``set variable = phase, temperature, temperature``) "
        "and follow the same order as ``isovalue`` and ``bounding box filename``.");
      prm.declare_entry(
        "isovalue",
        "0.0",
        Patterns::List(Patterns::Double()),
        "Isovalue(s) of monitored isocontour(s). "
        "When multiple isocontours are defined, the different isovalues "
        "must be separated by commas (e.g., ``set isovalue = 0.5, 300, 500``) "
        "and follow the same order as ``variable`` and ``bounding box filename``.");
      prm.declare_entry(
        "bounding box filename",
        "isocontour_bounding_box",
        Patterns::List(Patterns::FileName()),
        "Filename(s) for outputted isocontour(s). "
        "When multiple isocontours are defined, the different filenames must be "
        "separated by commas (e.g., ``set bounding box filename = interface_bounding_box, solidus_bounding_box, liquidus_bounding_box``) "
        "and follow the same order as ``variable`` and ``isovalue``.");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  PostProcessing<dim>::IsocontourBoundingBoxes::parse_parameters(
    ParameterHandler &prm)
  {
    prm.enter_subsection("isocontour bounding box");
    {
      number_of_isocontour_bounding_boxes =
        prm.get_integer("number of isocontour bounding boxes");

      // Get lists and split them into a vector
      const std::string        variables_list = prm.get("variable");
      std::vector<std::string> variables_vec =
        Utilities::split_string_list(variables_list);
      const std::string        isovalue_list = prm.get("isovalue");
      std::vector<std::string> isovalue_vec =
        Utilities::split_string_list(isovalue_list);
      const std::string        filename_list = prm.get("bounding box filename");
      std::vector<std::string> filename_vec =
        Utilities::split_string_list(filename_list);

      // Check that sizes are coherent with each other
      if (number_of_isocontour_bounding_boxes > 0)
        AssertThrow(variables_vec.size() == number_of_isocontour_bounding_boxes,
                    ListSizeIncoherentWithDeclaredNumber(
                      "number of isocontour bounding boxes",
                      number_of_isocontour_bounding_boxes,
                      "variable",
                      variables_vec.size()));
      AssertThrow(variables_vec.size() == isovalue_vec.size(),
                  ListsSizeMismatch("variable",
                                    "isovalue",
                                    variables_vec.size(),
                                    isovalue_vec.size()));
      AssertThrow(variables_vec.size() == filename_vec.size(),
                  ListsSizeMismatch("variable",
                                    "bounding box filename",
                                    variables_vec.size(),
                                    filename_vec.size()));

      if (number_of_isocontour_bounding_boxes > 0)
        {
          // Build map of isocontours
          for (unsigned int i = 0; i < variables_vec.size(); ++i)
            {
              // Initialize Isocontour
              Isocontour isocontour;

              // Get isovalue and output filename
              isocontour.isovalue =
                Utilities::string_to_double(isovalue_vec[i]);
              isocontour.output_name = filename_vec[i];

              // Parse variables
              if (variables_vec[i] == "temperature")
                {
                  ids_and_isocontours_per_variable.insert(
                    {Variable::temperature, std::make_pair(i, isocontour)});
                }
              else if (variables_vec[i] == "phase")
                {
                  ids_and_isocontours_per_variable.insert(
                    {Variable::phase, std::make_pair(i, isocontour)});
                }
              else
                throw std::invalid_argument(
                  "Error, the only valid variables for an 'isocontour bounding box' are: 'temperature' or 'phase'.");
            }
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  PostProcessing<dim>::ProbingPoints::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("probing points");
    {
      prm.declare_entry("number of probing points",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          number_of_probing_points),
                        Patterns::Integer(0),
                        "Number of probing points");

      const std::string default_point_entry_string =
        (dim == 2) ? "0., 0." : "0., 0., 0.";

      for (unsigned int i = 0; i < max_number_of_probing_points; ++i)
        {
          prm.enter_subsection("probe " + Utilities::int_to_string(i));
          {
            prm.declare_entry(
              "location",
              default_point_entry_string,
              Patterns::List(Patterns::Double(), dim, dim),
              "Probe point location in the mesh's reference frame. "
              "The different components of the point must be separated by commas (e.g., ``set location = 0.0, 0.0, 0.0``). ");
            prm.declare_entry(
              "variable",
              "velocity",
              Patterns::List(
                Patterns::Selection("velocity|pressure|phase|temperature")),
              "Variable(s) evaluated at the probing point. "
              "Choices are <velocity|pressure|phase|temperature>. "
              "When multiple variables are defined, the different variables must "
              "be separated by commas (e.g., ``set variable = velocity, pressure, temperature``).");
            prm.declare_entry("probing point filename",
                              "probe_" + Utilities::int_to_string(i, 2),
                              Patterns::FileName(),
                              "Filename for outputted probing point values.");
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  PostProcessing<dim>::ProbingPoints::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("probing points");
    {
      number_of_probing_points = prm.get_integer("number of probing points");
      AssertThrow(number_of_probing_points <= max_number_of_probing_points,
                  ExcMessage(
                    "The current maximum of probing points allowed is " +
                    Utilities::int_to_string(max_number_of_probing_points) +
                    ". Please adjust 'number of probing points'."));

      for (unsigned int id = 0; id < number_of_probing_points; ++id)
        {
          prm.enter_subsection("probe " + Utilities::int_to_string(id));
          {
            const std::string        variables_list = prm.get("variable");
            std::vector<std::string> variables_vec =
              Utilities::split_string_list(variables_list);

            const Point<dim> point(
              value_string_to_tensor<dim>(prm.get("location")));

            for (const std::string &variable : variables_vec)
              {
                if (variable == "velocity")
                  {
                    add_probing_point(Variable::velocity, id, point);
                  }
                else if (variable == "pressure")
                  {
                    add_probing_point(Variable::pressure, id, point);
                  }
                else if (variable == "phase")
                  {
                    add_probing_point(Variable::phase, id, point);
                  }
                else if (variable == "temperature")
                  {
                    add_probing_point(Variable::temperature, id, point);
                  }
                else
                  throw std::invalid_argument(
                    "Error, the only valid variables for a 'probing point' are: "
                    "'velocity', 'pressure','phase' and 'temperature'.");
              }

            const std::string filename = prm.get("probing point filename");
            probing_points_output_names.emplace_back(filename);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  PostProcessing<dim>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("post-processing");
    {
      prm.declare_entry(
        "verbosity",
        to_string(verbosity),
        Patterns::Selection("quiet|verbose"),
        "State whether from the post-processing values should be printed "
        "Choices are <quiet|verbose>.");

      prm.declare_entry(
        "calculate kinetic energy",
        Patterns::Tools::Convert<bool>::to_string(calculate_kinetic_energy),
        Patterns::Bool(),
        "Enable calculation of total kinetic energy. The total kinetic "
        "energy is calculated from the volumetric integral of the kinetic energy over the domain.");

      prm.declare_entry(
        "calculate enstrophy",
        Patterns::Tools::Convert<bool>::to_string(calculate_enstrophy),
        Patterns::Bool(),
        "Enable calculation of total enstrophy. The total enstrophy "
        "is calculated from the volumetric integral of the enstrophy over the domain.");

      prm.declare_entry(
        "calculate pressure power",
        Patterns::Tools::Convert<bool>::to_string(calculate_pressure_power),
        Patterns::Bool(),
        "Enable calculation of the pressure power. The pressure power "
        "is calculated from the volumetric integral of u.grad(p) over the domain.");

      prm.declare_entry(
        "calculate viscous dissipation",
        Patterns::Tools::Convert<bool>::to_string(
          calculate_viscous_dissipation),
        Patterns::Bool(),
        "Enable calculation of the viscous dissipation. The viscous dissipation "
        "is calculated from the volumetric integral of grad(u).tau over the domain.");

      prm.declare_entry("calculate apparent viscosity",
                        Patterns::Tools::Convert<bool>::to_string(
                          calculate_apparent_viscosity),
                        Patterns::Bool(),
                        "Enable calculation of apparent viscosity");

      prm.declare_entry("calculate average velocities",
                        Patterns::Tools::Convert<bool>::to_string(
                          calculate_average_velocities),
                        Patterns::Bool(),
                        "Enable calculation of average velocities.");

      prm.declare_entry(
        "calculate average temperature and heat flux",
        Patterns::Tools::Convert<bool>::to_string(
          calculate_average_temp_and_hf),
        Patterns::Bool(),
        "Enable calculation of time average temperature and time average heat flux");

      prm.declare_entry(
        "calculate pressure drop",
        Patterns::Tools::Convert<bool>::to_string(calculate_pressure_drop),
        Patterns::Bool(),
        "Enable calculation of pressure drop between two boundaries.");

      prm.declare_entry("inlet boundary id",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          inlet_boundary_id),
                        Patterns::Integer(),
                        "Inlet boundary ID for pressure drop calculation");

      prm.declare_entry("outlet boundary id",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          outlet_boundary_id),
                        Patterns::Integer(),
                        "Outlet boundary ID for pressure drop calculation");

      prm.declare_entry("calculate flow rate",
                        Patterns::Tools::Convert<bool>::to_string(
                          calculate_flow_rate),
                        Patterns::Bool(),
                        "Enable calculation of flow rate at boundaries.");

      prm.declare_entry(
        "calculate tracer flow rate",
        Patterns::Tools::Convert<bool>::to_string(calculate_tracer_flow_rate),
        Patterns::Bool(),
        "Enable calculation of tracer flow rate at boundaries.");

      prm.declare_entry(
        "initial time for average velocity",
        Patterns::Tools::Convert<double>::to_string(
          initial_time_for_average_velocities),
        Patterns::Double(),
        "Initial time to start calculations for average velocities");

      prm.declare_entry(
        "initial time for average temperature and heat flux",
        Patterns::Tools::Convert<double>::to_string(
          initial_time_for_average_temp_and_hf),
        Patterns::Double(),
        "Initial time to start calculations for average temperature");

      prm.declare_entry("kinetic energy name",
                        kinetic_energy_output_name,
                        Patterns::FileName(),
                        "File output kinetic energy");

      prm.declare_entry("pressure drop name",
                        pressure_drop_output_name,
                        Patterns::FileName(),
                        "File output pressure drop");

      prm.declare_entry("flow rate name",
                        flow_rate_output_name,
                        Patterns::FileName(),
                        "File output volumetric flux");

      prm.declare_entry("tracer flow rate name",
                        tracer_flow_rate_output_name,
                        Patterns::FileName(),
                        "Output file name for tracer flow rate");

      prm.declare_entry("enstrophy name",
                        enstrophy_output_name,
                        Patterns::FileName(),
                        "File output enstrophy");

      prm.declare_entry("pressure power name",
                        pressure_power_output_name,
                        Patterns::FileName(),
                        "File output pressure power");

      prm.declare_entry("viscous dissipation name",
                        viscous_dissipation_output_name,
                        Patterns::FileName(),
                        "File output viscous dissipation");

      prm.declare_entry("apparent viscosity name",
                        apparent_viscosity_output_name,
                        Patterns::FileName(),
                        "File output apparent viscosity");

      prm.declare_entry("output frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          output_frequency),
                        Patterns::Integer(),
                        "Output frequency");

      prm.declare_entry("calculation frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          calculation_frequency),
                        Patterns::Integer(),
                        "Calculation frequency of post-processed quantities.");

      prm.declare_entry("calculate tracer statistics",
                        Patterns::Tools::Convert<bool>::to_string(
                          calculate_tracer_statistics),
                        Patterns::Bool(),
                        "Enable calculation of tracer statistics.");

      prm.declare_entry("tracer statistics name",
                        tracer_output_name,
                        Patterns::FileName(),
                        "File name output tracer statistics");

      prm.declare_entry(
        "calculate phase statistics",
        Patterns::Tools::Convert<bool>::to_string(calculate_phase_statistics),
        Patterns::Bool(),
        "Enable calculation of phase statistics: maximum, minimum, average and integral over the domain (Cahn-Hilliard).");

      prm.declare_entry("phase statistics name",
                        phase_output_name,
                        Patterns::FileName(),
                        "File name output phase statistics (Cahn-Hilliard)");

      prm.declare_entry("calculate temperature statistics",
                        Patterns::Tools::Convert<bool>::to_string(
                          calculate_temperature_statistics),
                        Patterns::Bool(),
                        "Enable calculation of temperature statistics.");

      prm.declare_entry("temperature statistics name",
                        temperature_output_name,
                        Patterns::FileName(),
                        "File name output temperature statistics");

      prm.declare_entry("monitored fluid with phase change",
                        to_string(monitored_fluid_with_phase_change),
                        Patterns::Selection("fluid 0|fluid 1"),
                        "Fluid with phase change properties <fluid 0|fluid 1>");

      prm.declare_entry(
        "calculate algebraic melt volume",
        Patterns::Tools::Convert<bool>::to_string(
          calculate_algebraic_melt_volume),
        Patterns::Bool(),
        "Enable calculation of the algebraic (phase indicator and liquid fraction weighted) melt volume in the domain. In the case of CLS simulations, the fluid of interest is selected with 'monitored fluid with phase change.'");

      prm.declare_entry("algebraic melt volume name",
                        algebraic_melt_volume_output_name,
                        Patterns::FileName(),
                        "Filename of the algebraic melt volume output file");

      prm.declare_entry(
        "calculate geometric melt volume",
        Patterns::Tools::Convert<bool>::to_string(
          calculate_geometric_melt_volume),
        Patterns::Bool(),
        "Enable calculation of the geometric melt volume. "
        "The melt volume is computed as the volume of fluid over the 'melting temperature'."
        "In the case of CLS simulations, the volume is the geometrical volume within the 'monitored fluid with phase change.'");

      prm.declare_entry("geometric melt volume name",
                        geometric_melt_volume_output_name,
                        Patterns::FileName(),
                        "Filename of the geometric melt volume output file");

      prm.declare_entry(
        "melting temperature",
        Patterns::Tools::Convert<double>::to_string(melting_temperature),
        Patterns::Double(0),
        "Temperature used to define the melting point of the fluid for volume calculation.");

      prm.declare_entry("calculate heat flux",
                        Patterns::Tools::Convert<bool>::to_string(
                          calculate_heat_flux),
                        Patterns::Bool(),
                        "Enable calculation of heat flux.");

      prm.declare_entry("heat flux name",
                        heat_flux_output_name,
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
                        to_string(postprocessed_fluid),
                        Patterns::Selection("fluid 0|fluid 1|both"),
                        "Fluid domain used for thermal postprocesses "
                        "in the heat equation <fluid 0|fluid 1|both>");

      prm.declare_entry(
        "calculate barycenter",
        Patterns::Tools::Convert<bool>::to_string(calculate_barycenter),
        Patterns::Bool(),
        "Enable calculation of the barycenter location and velocity of fluid 1 in CLS and Cahn-Hilliard simulations.");

      prm.declare_entry(
        "barycenter name",
        barycenter_output_name,
        Patterns::FileName(),
        "Name of barycenter information output file in CLS or Cahn-Hilliard simulations");

      prm.declare_entry(
        "calculate mass conservation",
        Patterns::Tools::Convert<bool>::to_string(calculate_mass_conservation),
        Patterns::Bool(),
        "Enable calculation of the mass and momentum of both fluids in CLS simulations.");

      prm.declare_entry(
        "mass conservation name",
        mass_conservation_output_name,
        Patterns::FileName(),
        "Name of mass conservation output file in CLS simulations");

      prm.declare_entry(
        "calculate phase energy",
        Patterns::Tools::Convert<bool>::to_string(calculate_phase_energy),
        Patterns::Bool(),
        "Enable calculation of phase energies, including: total energy, bulk energy, and interface energy");

      prm.declare_entry(
        "phase energy name",
        phase_energy_output_name,
        Patterns::FileName(),
        "Name of energy output file in Cahn-Hilliard simulations. The file is stored in the output folder specified in the simulation control subsection");

      prm.declare_entry(
        "calculate phase volumes",
        Patterns::Tools::Convert<bool>::to_string(calculate_phase_volumes),
        Patterns::Bool(),
        "Enable calculation of total volume each phases in cfd-dem simulation, including: total volume of fluid, and total volume of particles");

      prm.declare_entry(
        "phase volumes name",
        phase_volumes_output_name,
        Patterns::FileName(),
        "Name of phases volume output file in cfd-dem simulations. The file is stored in the output folder specified in the simulation control subsection");

      prm.declare_entry("output qcriterion",
                        Patterns::Tools::Convert<bool>::to_string(
                          output_q_criterion),
                        Patterns::Bool(),
                        "Enable output of Q-criterion field <true|false>");

      prm.declare_entry("output vorticity",
                        Patterns::Tools::Convert<bool>::to_string(
                          output_vorticity),
                        Patterns::Bool(),
                        "Enable output of vorticity field <true|false>");

      prm.declare_entry(
        "output velocity gradient",
        Patterns::Tools::Convert<bool>::to_string(output_velocity_gradient),
        Patterns::Bool(),
        "Enable output of velocity gradient field <true|false>");

      isocontour_bounding_boxes.declare_parameters(prm);
      probing_points.declare_parameters(prm);
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  PostProcessing<dim>::parse_parameters(ParameterHandler &prm)
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
      calculation_frequency = prm.get_integer("calculation frequency");

      AssertThrow(
        output_frequency % calculation_frequency == 0,
        ExcMessage(
          "The post-processing 'output frequency' must be a multiple of the 'calculation frequency'."));

      calculate_tracer_statistics = prm.get_bool("calculate tracer statistics");
      tracer_output_name          = prm.get("tracer statistics name");
      calculate_phase_statistics  = prm.get_bool("calculate phase statistics");
      phase_output_name           = prm.get("phase statistics name");
      calculate_temperature_statistics =
        prm.get_bool("calculate temperature statistics");
      const std::string melt_fluid =
        prm.get("monitored fluid with phase change");
      if (melt_fluid == "fluid 0")
        monitored_fluid_with_phase_change = Parameters::FluidIndicator::fluid0;
      else if (melt_fluid == "fluid 1")
        monitored_fluid_with_phase_change = Parameters::FluidIndicator::fluid1;
      else
        throw(std::invalid_argument("Invalid fluid. "
                                    "Options are 'fluid 0' or 'fluid 1'."));
      calculate_algebraic_melt_volume =
        prm.get_bool("calculate algebraic melt volume");
      algebraic_melt_volume_output_name = prm.get("algebraic melt volume name");
      calculate_geometric_melt_volume =
        prm.get_bool("calculate geometric melt volume");
      geometric_melt_volume_output_name = prm.get("geometric melt volume name");
      melting_temperature               = prm.get_double("melting temperature");
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
      output_q_criterion            = prm.get_bool("output qcriterion");
      output_vorticity              = prm.get_bool("output vorticity");
      output_velocity_gradient      = prm.get_bool("output velocity gradient");

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

      isocontour_bounding_boxes.parse_parameters(prm);
      probing_points.parse_parameters(prm);
    }
    prm.leave_subsection();
  }

  namespace
  {
    std::string
    to_string(const NonLinearSolver::SolverType type)
    {
      switch (type)
        {
          case NonLinearSolver::SolverType::newton:
            return "newton";
          case NonLinearSolver::SolverType::inexact_newton:
            return "inexact_newton";
          case NonLinearSolver::SolverType::kinsol_newton:
            return "kinsol_newton";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const NonLinearSolver::KinsolStrategy strategy)
    {
      switch (strategy)
        {
          case NonLinearSolver::KinsolStrategy::normal_newton:
            return "normal_newton";
          case NonLinearSolver::KinsolStrategy::line_search:
            return "line_search";
          case NonLinearSolver::KinsolStrategy::picard:
            return "picard";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  NonLinearSolver::declare_parameters(ParameterHandler  &prm,
                                      const std::string &physics_name)
  {
    const NonLinearSolver defaults;
    prm.enter_subsection("non-linear solver");
    {
      prm.enter_subsection(physics_name);
      {
        prm.declare_entry(
          "verbosity",
          to_string(defaults.verbosity),
          Patterns::Selection("quiet|verbose"),
          "State whether the outputs from the non-linear solver should be printed. "
          "Choices are <quiet|verbose>.");

        prm.declare_entry(
          "solver",
          to_string(defaults.solver),
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
          to_string(defaults.kinsol_strategy),
          Patterns::Selection("normal_newton|line_search|fixed_point|picard"),
          "Strategy that will be used by the kinsol newton solver");

        prm.declare_entry("tolerance",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.tolerance),
                          Patterns::Double(),
                          "Newton solver tolerance");
        prm.declare_entry("max iterations",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.max_iterations),
                          Patterns::Integer(),
                          "Maximum number of Newton Iterations");

        prm.declare_entry(
          "step tolerance",
          Patterns::Tools::Convert<double>::to_string(defaults.step_tolerance),
          Patterns::Double(),
          "Newton solver relative tolerance between steps."
          " If a newton iteration leads to a residual > step tolerance"
          " * previous residual then the theta relaxation"
          " is applied until this criteria is satisfied");

        prm.declare_entry(
          "matrix tolerance",
          Patterns::Tools::Convert<double>::to_string(
            defaults.matrix_tolerance),
          Patterns::Double(),
          "This parameter controls the frequency at which the matrix is refreshed in the inexact Newton solvers"
          "If the residual after a newton step < previous residual * matrix tolerance, the matrix is not re-assembled");

        prm.declare_entry(
          "force rhs calculation",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.force_rhs_calculation),
          Patterns::Bool(),
          "This is required if there is a fixed point component to the non-linear"
          "solver that is changed at the beginning of every newton iteration."
          "This is notably the case of the sharp edge method."
          "The default value of this parameter is false.");


        prm.declare_entry(
          "reuse matrix",
          Patterns::Tools::Convert<bool>::to_string(defaults.reuse_matrix),
          Patterns::Bool(),
          "Reuse the last jacobian matrix for the next non-linear problem solution");

        prm.declare_entry(
          "reuse preconditioner",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.reuse_preconditioner),
          Patterns::Bool(),
          "Reuse the last preconditioner for the next non-linear problem solution");

        prm.declare_entry(
          "abort at convergence failure",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.abort_at_convergence_failure),
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
  namespace
  {
    std::string
    to_string(const Mesh::Type type)
    {
      switch (type)
        {
          case Mesh::Type::gmsh:
            return "gmsh";
          case Mesh::Type::dealii:
            return "dealii";
          case Mesh::Type::lethe:
            return "lethe";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const Tensor<1, 3> &tensor)
    {
      return Patterns::Tools::Convert<double>::to_string(tensor[0]) + ", " +
             Patterns::Tools::Convert<double>::to_string(tensor[1]) + ", " +
             Patterns::Tools::Convert<double>::to_string(tensor[2]);
    }
  } // namespace

  void
  Mesh::declare_parameters(ParameterHandler &prm)
  {
    const Mesh defaults;
    prm.enter_subsection("mesh");
    {
      prm.declare_entry("type",
                        to_string(defaults.type),
                        Patterns::Selection("gmsh|dealii|lethe"),
                        "Type of mesh "
                        "Choices are <gmsh|dealii|lethe>.");

      prm.declare_entry("file name",
                        defaults.file_name,
                        Patterns::FileName(),
                        "GMSH file name");

      prm.declare_entry("initial refinement",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.initial_refinement),
                        Patterns::Integer(),
                        "Initial refinement of the mesh");

      prm.declare_entry(
        "initial boundary refinement",
        Patterns::Tools::Convert<unsigned int>::to_string(
          defaults.initial_refinement_at_boundaries),
        Patterns::Integer(),
        "Initial refinement of the mesh at the boundaries specified by the user");

      prm.declare_entry(
        "boundaries refined",
        "",
        Patterns::List(Patterns::Integer()),
        "Boundary ids of the boundaries to be initially refined");

      prm.declare_entry(
        "enable target size",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.refine_until_target_size),
        Patterns::Bool(),
        "Enable initial refinement until target size is reached.");

      prm.declare_entry(
        "simplex",
        Patterns::Tools::Convert<bool>::to_string(defaults.simplex),
        Patterns::Bool(),
        "Indicates that the mesh used is a mesh made of only simplex elements.");

      prm.declare_entry(
        "check diamond cells",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.check_for_diamond_cells),
        Patterns::Bool(),
        "Enables checking the input grid for diamond-shaped cells.");

      prm.declare_entry(
        "expand particle-wall contact search",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.expand_particle_wall_contact_search),
        Patterns::Bool(),
        "Enables adding the boundary neighbor cells of boundary cells to the"
        "particle-wall contact search list. This feature should only be "
        "activated in geometries with concave boundaries. (For example, for "
        "particles flow inside a cylinder or sphere). In geometries with "
        "convex boundaries, this feature MUST NOT be activated");

      prm.declare_entry("target size",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.target_size),
                        Patterns::Double(),
                        "Target size of the initial refinement");

      prm.declare_entry("grid type", defaults.grid_type);
      prm.declare_entry("grid arguments", defaults.grid_arguments);

      prm.declare_entry(
        "initial translation",
        to_string(defaults.translation),
        Patterns::List(Patterns::Double()),
        "Component of the desired translation of the mesh at initialization. \n"
        "In 2D, the third value (z-component) is ignored.");

      prm.declare_entry(
        "initial rotation axis",
        to_string(defaults.rotation_axis),
        Patterns::List(Patterns::Double()),
        "Component of the desired rotation of the mesh at initialization.\n"
        "In 2D, this has no effect, only a counter-clockwise rotation around the origin \n "
        "of the coordinate system is applied.");

      prm.declare_entry(
        "initial rotation angle",
        Patterns::Tools::Convert<double>::to_string(defaults.rotation_angle),
        Patterns::Double(),
        "Angle of rotation of the mesh at initialization around the axis in radian");

      prm.declare_entry("scale",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.scale),
                        Patterns::Double(0),
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
        else if (op == "lethe")
          type = Type::lethe;
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
      AssertThrow(
        !(check_for_diamond_cells && expand_particle_wall_contact_search),
        ExcMessage(
          "Lethe DEM does not currently support the simultaneous use of the 'check diamond "
          "cells' and 'expand particle-wall contact search' options. Please set one of "
          "these two options to false."));


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
      prm.declare_entry("number of refinement boxes",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          number_of_refinement_boxes),
                        Patterns::Integer(0),
                        "Number of refinement boxes specified");
      for (unsigned int i_box = 0; i_box < max_number_of_refinement_boxes;
           ++i_box)
        {
          prm.enter_subsection("box " + Utilities::int_to_string(i_box, 1));
          {
            (*refinement_boxes_meshes)[i_box].declare_parameters(prm);
            prm.declare_entry(
              "additional refinement",
              "0",
              Patterns::Integer(0),
              "Additional refinements of the principal mesh within the area delimited by the 'box' mesh.");
            prm.declare_alias("additional refinement",
                              "initial refinement",
                              true);
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  void
  MeshBoxRefinement::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("box refinement");
    {
      number_of_refinement_boxes =
        prm.get_integer("number of refinement boxes");
      AssertThrow(
        number_of_refinement_boxes <= max_number_of_refinement_boxes,
        ExcMessage(
          "The current implementation limits the number of refinement boxes up to " +
          Utilities::int_to_string(max_number_of_refinement_boxes) +
          " refinement boxes.\n You have declared " +
          Utilities::int_to_string(number_of_refinement_boxes) +
          " refinement boxes.\n Please reduce the number of refinement boxes."));

      for (unsigned int i_box = 0; i_box < number_of_refinement_boxes; ++i_box)
        {
          prm.enter_subsection("box " + Utilities::int_to_string(i_box, 1));
          {
            (*refinement_boxes_meshes)[i_box].parse_parameters(prm);
            box_additional_refinements[i_box] =
              prm.get_integer("additional refinement");
          }
          prm.leave_subsection();
        }
    }
    prm.leave_subsection();
  }

  namespace
  {
    std::string
    to_string(const LinearSolver::SolverType type)
    {
      switch (type)
        {
          case LinearSolver::SolverType::gmres:
            return "gmres";
          case LinearSolver::SolverType::bicgstab:
            return "bicgstab";
          case LinearSolver::SolverType::direct:
            return "direct";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const LinearSolver::PreconditionerType type)
    {
      switch (type)
        {
          case LinearSolver::PreconditionerType::ilu:
            return "ilu";
          case LinearSolver::PreconditionerType::amg:
            return "amg";
          case LinearSolver::PreconditionerType::lsmg:
            return "lsmg";
          case LinearSolver::PreconditionerType::gcmg:
            return "gcmg";
          case LinearSolver::PreconditionerType::none:
            return "none";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const LinearSolver::MultigridCoarseningSequenceType type)
    {
      switch (type)
        {
          case LinearSolver::MultigridCoarseningSequenceType::h:
            return "h";
          case LinearSolver::MultigridCoarseningSequenceType::p:
            return "p";
          case LinearSolver::MultigridCoarseningSequenceType::hp:
            return "hp";
          case LinearSolver::MultigridCoarseningSequenceType::ph:
            return "ph";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(
      const MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType
        type)
    {
      switch (type)
        {
          case MGTransferGlobalCoarseningTools::
            PolynomialCoarseningSequenceType::decrease_by_one:
            return "decrease by one";
          case MGTransferGlobalCoarseningTools::
            PolynomialCoarseningSequenceType::bisect:
            return "bisect";
          case MGTransferGlobalCoarseningTools::
            PolynomialCoarseningSequenceType::go_to_one:
            return "go to one";
          default:
            Assert(false, ExcInternalError());
            return "";
        }
    }

    std::string
    to_string(const LinearSolver::MultigridSmootherPreconditionerType type)
    {
      switch (type)
        {
          case LinearSolver::MultigridSmootherPreconditionerType::
            InverseDiagonal:
            return "inverse diagonal";
          case LinearSolver::MultigridSmootherPreconditionerType::
            AdditiveSchwarzMethod:
            return "additive schwarz method";
          case LinearSolver::MultigridSmootherPreconditionerType::Chebyshev:
            return "chebyshev";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const LinearSolver::CoarseGridSolverType type)
    {
      switch (type)
        {
          case LinearSolver::CoarseGridSolverType::gmres:
            return "gmres";
          case LinearSolver::CoarseGridSolverType::amg:
            return "amg";
          case LinearSolver::CoarseGridSolverType::ilu:
            return "ilu";
          case LinearSolver::CoarseGridSolverType::direct:
            return "direct";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  LinearSolver::declare_parameters(ParameterHandler  &prm,
                                   const std::string &physics_name)
  {
    const LinearSolver defaults;
    prm.enter_subsection("linear solver");
    {
      prm.enter_subsection(physics_name);
      {
        prm.declare_entry(
          "verbosity",
          to_string(defaults.verbosity),
          Patterns::Selection("quiet|verbose|extra verbose"),
          "State whether output from solver runs should be printed. "
          "Choices are <quiet|verbose|extra verbose>.");
        prm.declare_entry(
          "method",
          to_string(defaults.solver),
          Patterns::Selection("gmres|bicgstab|direct"),
          "The iterative solver for the linear system of equations. "
          "Choices are <gmres|bicgstab|direct>.");

        prm.declare_entry(
          "rescale residual",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.rescale_residual_by_volume),
          Patterns::Bool(),
          "Rescale the residual by the square root of the volume of the triangulation");
        prm.declare_entry("relative residual",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.relative_residual),
                          Patterns::Double(),
                          "Linear solver residual");
        prm.declare_entry("minimum residual",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.minimum_residual),
                          Patterns::Double(),
                          "Linear solver minimum residual");
        prm.declare_entry("max iters",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.max_iterations),
                          Patterns::Integer(),
                          "Maximum solver iterations");

        prm.declare_entry("max krylov vectors",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.max_krylov_vectors),
                          Patterns::Integer(),
                          "Maximum number of krylov vectors for GMRES");

        prm.declare_entry(
          "enable hessians in jacobian",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.enable_hessians_jacobian),
          Patterns::Bool(),
          "Turns off the terms involving the hessian in the Jacobian");

        prm.declare_entry(
          "enable hessians in residual",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.enable_hessians_residual),
          Patterns::Bool(),
          "Turns off the terms involving the hessian in the rhs");

        prm.declare_entry("preconditioner",
                          to_string(defaults.preconditioner),
                          Patterns::Selection("amg|ilu|lsmg|gcmg|none"),
                          "The preconditioner for the linear solver."
                          "Choices are <amg|ilu|lsmg|gcmg|none>.");


        prm.declare_entry("ilu preconditioner fill",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.ilu_precond_fill),
                          Patterns::Double(),
                          "Ilu preconditioner fill");

        prm.declare_entry("ilu preconditioner absolute tolerance",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.ilu_precond_atol),
                          Patterns::Double(),
                          "Ilu preconditioner tolerance");

        prm.declare_entry("ilu preconditioner relative tolerance",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.ilu_precond_rtol),
                          Patterns::Double(),
                          "Ilu relative tolerance");

        prm.declare_entry("amg preconditioner ilu fill",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.amg_precond_ilu_fill),
                          Patterns::Double(),
                          "amg preconditioner ilu smoother fill");

        prm.declare_entry("amg preconditioner ilu absolute tolerance",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.amg_precond_ilu_atol),
                          Patterns::Double(),
                          "amg preconditioner ilu smoother absolute tolerance");

        prm.declare_entry("amg preconditioner ilu relative tolerance",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.amg_precond_ilu_rtol),
                          Patterns::Double(),
                          "amg preconditioner ilu smoother relative tolerance");

        prm.declare_entry("amg aggregation threshold",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.amg_aggregation_threshold),
                          Patterns::Double(),
                          "amg aggregation threshold");
        prm.declare_entry("amg n cycles",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.amg_n_cycles),
                          Patterns::Integer(),
                          "amg number of cycles");
        prm.declare_entry("amg w cycles",
                          Patterns::Tools::Convert<bool>::to_string(
                            defaults.amg_w_cycles),
                          Patterns::Bool(),
                          "amg w cycling. If this is set to true, W cycling is "
                          "used. Otherwise, V cycling is used.");
        prm.declare_entry("amg smoother sweeps",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.amg_smoother_sweeps),
                          Patterns::Integer(),
                          "amg smoother sweeps");
        prm.declare_entry("amg smoother overlap",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.amg_smoother_overlap),
                          Patterns::Integer(),
                          "amg smoother overlap");
        prm.declare_entry(
          "force linear solver continuation",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.force_linear_solver_continuation),
          Patterns::Bool(),
          "A boolean that will force the linear solver to continue even if it fails");

        prm.declare_entry("mg min level",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.mg_min_level),
                          Patterns::Integer(),
                          "mg min level");

        prm.declare_entry("mg level min cells",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.mg_level_min_cells),
                          Patterns::Integer(),
                          "mg minimum number of cells for coarse level");

        prm.declare_entry("mg int level",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.mg_int_level),
                          Patterns::Integer(),
                          "mg int level");

        prm.declare_entry(
          "mg enable hessians in jacobian",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.mg_enable_hessians_jacobian),
          Patterns::Bool(),
          "Turns off the terms involving the hessian in the Jacobian of mg operators");

        prm.declare_entry("mg smoother iterations",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.mg_smoother_iterations),
                          Patterns::Integer(),
                          "mg smoother iterations for lsmg or gcmg");

        prm.declare_entry("mg smoother relaxation",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.mg_smoother_relaxation),
                          Patterns::Double(),
                          "mg smoother relaxation for lsmg or gcmg");

        prm.declare_entry(
          "mg smoother preconditioner type",
          to_string(defaults.mg_smoother_preconditioner_type),
          Patterns::Selection(
            "inverse diagonal|additive schwarz method|chebyshev"),
          "Preconditioner of smoother");

        prm.declare_entry("mg smoother chebyshev degree",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.mg_smoother_chebyshev_degree),
                          Patterns::Integer(0),
                          "polynomial degree of the Chebyshev smoother");

        prm.declare_entry(
          "mg smoother chebyshev smoothing range",
          Patterns::Tools::Convert<double>::to_string(
            defaults.mg_smoother_chebyshev_smoothing_range),
          Patterns::Double(1.0),
          "smoothing range (lambda_max/lambda_min) of the Chebyshev smoother");

        prm.declare_entry(
          "mg smoother chebyshev eig cg n iterations",
          Patterns::Tools::Convert<unsigned int>::to_string(
            defaults.mg_smoother_chebyshev_eig_cg_n_iterations),
          Patterns::Integer(1),
          "cg/Lanczos iterations to estimate the maximum eigenvalue for the "
          "Chebyshev smoother");

        prm.declare_entry("mg smoother eig estimation",
                          Patterns::Tools::Convert<bool>::to_string(
                            defaults.mg_smoother_eig_estimation),
                          Patterns::Bool(),
                          "estimate eigenvalues for relaxation parameter");

        prm.declare_entry("eig estimation smoothing range",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.eig_estimation_smoothing_range),
                          Patterns::Double(),
                          "sets range between largest and smallest eig");

        prm.declare_entry("eig estimation cg n iterations",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.eig_estimation_cg_n_iterations),
                          Patterns::Integer(),
                          "cg iterations performed to find eigenvalue");

        prm.declare_entry("eig estimation verbosity",
                          to_string(defaults.eig_estimation_verbose),
                          Patterns::Selection("quiet|verbose"),
                          "State whether MG should print max and min eigenvalue"
                          "Choices are <quiet|verbose>.");

        prm.declare_entry("mg coarse grid solver",
                          to_string(defaults.mg_coarse_grid_solver),
                          Patterns::Selection("gmres|amg|ilu|direct"),
                          "The coarse grid solver for lsmg or gcmg"
                          "Choices are <gmres|amg|ilu|direct>.");

        prm.declare_entry(
          "mg coarse grid use fe q iso q1",
          Patterns::Tools::Convert<bool>::to_string(
            defaults.mg_use_fe_q_iso_q1),
          Patterns::Bool(),
          "use elements with linear interpolation for coarse grid");

        prm.declare_entry("mg coarsening type",
                          to_string(defaults.mg_coarsening_type),
                          Patterns::Selection("h|p|hp|ph"),
                          "mg coarsening type for gcmg");

        prm.declare_entry("mg p coarsening type",
                          to_string(defaults.mg_p_coarsening_type),
                          Patterns::Selection(
                            "decrease by one|bisect|go to one"),
                          "mg p coarsening type for gcmg");

        prm.declare_entry("mg p min coarsening degree",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            defaults.mg_p_min_coarsening_degree),
                          Patterns::Integer(),
                          "mg p minimum coarsening degree for gcmg");

        prm.declare_entry("mg gmres max iterations",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.mg_gmres_max_iterations),
                          Patterns::Integer(),
                          "mg gmres iterations for lsmg or gcmg");

        prm.declare_entry("mg gmres tolerance",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.mg_gmres_tolerance),
                          Patterns::Double(),
                          "mg gmres tolerance n for lsmg or gcmg");

        prm.declare_entry("mg gmres reduce",
                          Patterns::Tools::Convert<double>::to_string(
                            defaults.mg_gmres_reduce),
                          Patterns::Double(),
                          "mg gmres reduce for lsmg or gcmg");

        prm.declare_entry("mg gmres max krylov vectors",
                          Patterns::Tools::Convert<int>::to_string(
                            defaults.mg_gmres_max_krylov_vectors),
                          Patterns::Integer(),
                          "mg gmres max krylov vectors for lsmg or gcmg");

        prm.declare_entry("mg gmres preconditioner",
                          to_string(defaults.mg_gmres_preconditioner),
                          Patterns::Selection("amg|ilu"),
                          "The preconditioner for the mg gmres solver"
                          "Choices are <amg|ilu>.");

        prm.declare_entry("mg amg use default parameters",
                          Patterns::Tools::Convert<bool>::to_string(
                            defaults.mg_amg_use_default_parameters),
                          Patterns::Bool(),
                          "Use default parameters for Trilinos AMG");

        prm.declare_entry(
          "mg verbosity",
          to_string(defaults.mg_verbosity),
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
        else if (precond == "none")
          preconditioner = PreconditionerType::none;
        else
          throw std::logic_error(
            "Error, invalid preconditioner type. Choices are amg, ilu, lsmg, gcmg or none.");


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
        else if (mg_smoother_preconditioner_type_str == "chebyshev")
          this->mg_smoother_preconditioner_type =
            MultigridSmootherPreconditionerType::Chebyshev;
        else
          AssertThrow(false, ExcNotImplemented());

        mg_smoother_chebyshev_degree =
          prm.get_integer("mg smoother chebyshev degree");
        mg_smoother_chebyshev_smoothing_range =
          prm.get_double("mg smoother chebyshev smoothing range");
        mg_smoother_chebyshev_eig_cg_n_iterations =
          prm.get_integer("mg smoother chebyshev eig cg n iterations");

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

        mg_p_min_coarsening_degree =
          prm.get_integer("mg p min coarsening degree");

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

  namespace
  {
    std::string
    to_string(const MeshAdaptation::Type type)
    {
      switch (type)
        {
          case MeshAdaptation::Type::none:
            return "none";
          case MeshAdaptation::Type::uniform:
            return "uniform";
          case MeshAdaptation::Type::adaptive:
            return "adaptive";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const MultipleAdaptationParameters::ErrorEstimator estimator)
    {
      switch (estimator)
        {
          case MultipleAdaptationParameters::ErrorEstimator::kelly:
            return "kelly";
          case MultipleAdaptationParameters::ErrorEstimator::dpg:
            return "dpg";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  MeshAdaptation::declare_parameters(ParameterHandler &prm)
  {
    const MeshAdaptation                defaults;
    const MultipleAdaptationParameters &variable_defaults =
      defaults.var_adaptation_param;
    prm.enter_subsection("mesh adaptation");
    {
      prm.declare_entry("initial refinement steps",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.initial_refinement),
                        Patterns::Integer(),
                        "Number of pre-solve adaptive mesh refinement steps");

      prm.declare_entry("type",
                        to_string(defaults.type),
                        Patterns::Selection("none|uniform|adaptive"),
                        "Type of mesh adaptation"
                        "Choices are <none|uniform|adaptive>.");

      prm.declare_entry(
        "error estimator",
        to_string(variable_defaults.error_estimator),
        Patterns::List(Patterns::Selection("kelly|dpg")),
        "Error estimator for adaptive mesh refinement. For multi-variables refinement, separate the different strategies with a comma. They should follow the same order as what is specified in the variable parameter."
        "Choices are <kelly|dpg>.");

      prm.declare_entry(
        "fraction refinement",
        Patterns::Tools::Convert<double>::to_string(
          variable_defaults.refinement_fraction),
        Patterns::List(Patterns::Double()),
        "Fraction of refined elements"
        "For multi-variables refinement, separate the different fractions with a comma "
        "(ex/ 'set fraction refinement = 0.1,0.1')");

      prm.declare_entry(
        "fraction coarsening",
        Patterns::Tools::Convert<double>::to_string(
          variable_defaults.coarsening_fraction),
        Patterns::List(Patterns::Double()),
        "Fraction of coarsened elements"
        "For multi-variables refinement, separate the different fractions with a comma "
        "(ex/ 'set fraction coarsening = 0.05,0.05')");

      prm.declare_entry(
        "variable",
        get_variable_string(defaults.vars),
        Patterns::List(Patterns::Selection(
          "velocity|pressure|phase|temperature|phase_cahn_hilliard|chemical_potential_cahn_hilliard|tracer|electric field|magnetic field|electromagnetic fields")),
        "Variable(s) for error estimation"
        "Choices are <velocity|pressure|phase|temperature|phase_cahn_hilliard|chemical_potential_cahn_hilliard|tracer|electric field|magnetic field|electromagnetic_fields>."
        "For multi-variables refinement, separate the different variables with a comma "
        "(ex/ 'set variable = velocity,temperature')");

      prm.declare_entry(
        "fraction type",
        defaults.fractionType == FractionType::number ? "number" : "fraction",
        Patterns::Selection("number|fraction"),
        "How the fraction of refinement/coarsening are interpreted"
        "Choices are <number|fraction>.");
      prm.declare_entry("max number elements",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.maximum_number_elements),
                        Patterns::Integer(),
                        "Maximum number of elements");
      prm.declare_entry("max refinement level",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.maximum_refinement_level),
                        Patterns::Integer(),
                        "Maximum refinement level");
      prm.declare_entry("min refinement level",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.minimum_refinement_level),
                        Patterns::Integer(),
                        "Minimum refinement level");
      prm.declare_entry("frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.frequency),
                        Patterns::Integer(),
                        "Frequency of the mesh refinement");
      prm.declare_entry(
        "mesh refinement controller",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.mesh_controller_is_enabled),
        Patterns::Bool(),
        "Fraction of refined elements"
        "Enable a controller that will target a specific number of elements in the mesh equal to the maximum number of elements");
      prm.declare_entry("fix boundary refinement",
                        Patterns::Tools::Convert<bool>::to_string(
                          defaults.is_boundary_refinement_fixed),
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
      else if (op == "uniform")
        type = Type::uniform;
      else if (op == "adaptive")
        type = Type::adaptive;
      else
        throw std::logic_error(
          "Error, invalid mesh adaptation type. Choices are <none|uniform|adaptive>.");

      // Getting multivariables refinement parameters
      const std::string        var_op  = prm.get("variable");
      std::vector<std::string> var_vec = Utilities::split_string_list(var_op);
      const std::string        strategy_op = prm.get("error estimator");
      std::vector<std::string> strategy_vec =
        Utilities::split_string_list(strategy_op);
      const std::string        coars_op = prm.get("fraction coarsening");
      std::vector<std::string> coars_vec =
        Utilities::split_string_list(coars_op);
      const std::string        refin_op = prm.get("fraction refinement");
      std::vector<std::string> refin_vec =
        Utilities::split_string_list(refin_op);

      // Checking that the sizes are coherent
      Assert(strategy_vec.size() == var_vec.size(),
             MultipleAdaptationSizeError("error estimator",
                                         strategy_vec.size(),
                                         var_vec.size()));
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
          else if (var_vec[i] == "electric field")
            vars = Variable::electric_field;
          else if (var_vec[i] == "magnetic field")
            vars = Variable::magnetic_field;
          else if (var_vec[i] == "electromagnetic fields")
            vars = Variable::electromagnetic_fields;
          else
            throw std::logic_error(
              "Error, invalid mesh adaptation variable. Choices are velocity, pressure, phase, temperature, phase_cahn_hilliard, chemical_potential_cahn_hilliard, electric field, magnetic field or electromagnetic fields. Note that <electric field> or <magnetic field> and <electromagnetic fields> are mutually exclusive.");

          // Parsing strategy for this variable
          if (strategy_vec[i] == "kelly")
            var_adaptation_param.error_estimator =
              MultipleAdaptationParameters::ErrorEstimator::kelly;
          else if (strategy_vec[i] == "dpg")
            var_adaptation_param.error_estimator =
              MultipleAdaptationParameters::ErrorEstimator::dpg;
          else
            throw std::logic_error(
              "Error, invalid mesh adaptation error estimator. Choices are kelly or dpg");

          var_adaptation_param.coarsening_fraction =
            Utilities::string_to_double(coars_vec[i]);
          var_adaptation_param.refinement_fraction =
            Utilities::string_to_double(refin_vec[i]);

          // defining adaptation map for this variable
          variables[vars] = var_adaptation_param;
        }
      // Verify that the user did not specify both electric_field or
      // magnetic_field with electromagnetic_fields
      const bool has_em =
        variables.find(Variable::electromagnetic_fields) != variables.end();

      const bool has_e =
        variables.find(Variable::electric_field) != variables.end();

      const bool has_h =
        variables.find(Variable::magnetic_field) != variables.end();

      AssertThrow(!(has_em && (has_e || has_h)),
                  ExcMessage(
                    "Invalid mesh adaptation configuration: "
                    "electromagnetic_fields is mutually exclusive with "
                    "electric_field and magnetic_field."));

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

      // Fixing the boundary refinement is a no-op when no boundary id is
      // given, since no cell would ever be matched. Warn the user instead of
      // silently ignoring the request.
      AssertThrow(
        !is_boundary_refinement_fixed || !boundaries_to_fix.empty(),
        ExcMessage(
          "The parameter 'fix boundary refinement' is enabled, but the "
          "parameter 'boundaries fixed' is empty. No boundary would be fixed. "
          "Specify the boundary IDs to fix using 'boundaries fixed'."));
    }
    prm.leave_subsection();
  }

  namespace
  {
    std::string
    to_string(const Testing::TestType type)
    {
      switch (type)
        {
          case Testing::TestType::particles:
            return "particles";
          case Testing::TestType::mobility_status:
            return "mobility_status";
          case Testing::TestType::subdomain:
            return "subdomain";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  Testing::declare_parameters(ParameterHandler &prm)
  {
    const Testing defaults;
    prm.enter_subsection("test");
    {
      prm.declare_entry(
        "enable",
        Patterns::Tools::Convert<bool>::to_string(defaults.enabled),
        Patterns::Bool(),
        "Enable testing mode of a solver. Some solvers have a specific"
        "testing mode which enables the output of debug variables. This"
        "testing mode is generally used only for the automatic testing bench using ctest.");
      prm.declare_entry(
        "type",
        to_string(defaults.test_type),
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
    const Restart defaults;
    prm.enter_subsection("restart");
    {
      prm.declare_entry("filename",
                        defaults.filename,
                        Patterns::FileName(),
                        "Prefix for the filename of checkpoints");
      prm.declare_entry("restart",
                        Patterns::Tools::Convert<bool>::to_string(
                          defaults.restart),
                        Patterns::Bool(),
                        "Frequency for checkpointing");
      prm.declare_entry(
        "checkpoint",
        Patterns::Tools::Convert<bool>::to_string(defaults.checkpoint),
        Patterns::Bool(),
        "Enable checkpointing. Checkpointing creates a restart"
        "point from which the simulation can be restarted from.");

      prm.declare_entry("frequency",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.frequency),
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

  namespace
  {
    std::string
    to_string(const VelocitySource::RotatingFrameType type)
    {
      switch (type)
        {
          case VelocitySource::RotatingFrameType::none:
            return "none";
          case VelocitySource::RotatingFrameType::srf:
            return "srf";
        }
      Assert(false, ExcInternalError());
      return "";
    }

    std::string
    to_string(const VelocitySource::PermeabilityModel model)
    {
      switch (model)
        {
          case VelocitySource::PermeabilityModel::none:
            return "none";
          case VelocitySource::PermeabilityModel::darcy_phase_change:
            return "darcy phase change";
          case VelocitySource::PermeabilityModel::carman_kozeny_phase_change:
            return "carman-kozeny phase change";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  VelocitySource::declare_parameters(ParameterHandler &prm)
  {
    const VelocitySource defaults;
    prm.enter_subsection("velocity source");
    {
      prm.declare_entry(
        "rotating frame type",
        to_string(defaults.rotating_frame_type),
        Patterns::Selection("none|srf"),
        "Rotating frame velocity-dependent source terms"
        "Choices are <none|srf>. The srf stands"
        "for single rotating frame and adds"
        "the coriolis and the centrifugal force to the Navier-Stokes equations");

      prm.declare_entry(
        "omega_x",
        Patterns::Tools::Convert<double>::to_string(defaults.omega_x),
        Patterns::Double(),
        "X component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "omega_y",
        Patterns::Tools::Convert<double>::to_string(defaults.omega_y),
        Patterns::Double(),
        "Y component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "omega_z",
        Patterns::Tools::Convert<double>::to_string(defaults.omega_z),
        Patterns::Double(),
        "Z component of the angular velocity vector of the frame of reference");

      prm.declare_entry(
        "enable Darcy multiply by density",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.enable_darcy_multiply_by_density),
        Patterns::Bool(),
        "Enable the multiplication of the Darcy force term by the density for dimensional consistency when solving the pressure rather than the kinematic pressure in the momentum balance.");

      prm.declare_entry(
        "permeability model",
        to_string(defaults.permeability_model),
        Patterns::Selection(
          "none|darcy phase change|carman-kozeny phase change"),
        "Permeability models for phase change modelling."
        "Choices are <none|darcy phase change|carman-kozeny phase change>.");

      prm.declare_entry("Carman-Kozeny fluid with phase change",
                        to_string(defaults.fluid_with_phase_change),
                        Patterns::Selection("fluid 0|fluid 1|both"),
                        "Select which fluids have phase change"
                        "Choices are <fluid 0|fluid 1|both>.");

      prm.declare_entry(
        "Carman-Kozeny division tolerance",
        "1e-3",
        Patterns::List(Patterns::Double()),
        "This tolerance avoids a division by zero in the Carman-Kozeny source term. For multiple fluids with phase change, separate values with a comma.");

      prm.declare_entry(
        "Carman-Kozeny permeability area",
        "1e-3",
        Patterns::List(Patterns::Double()),
        "This represents the permeability area of the pseudo-porous bed in the Carman-Kozeny source term. For multiple fluids with phase change, separate values with a comma.");
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

      const std::string permeability_model_str = prm.get("permeability model");
      if (permeability_model_str == "none")
        permeability_model = PermeabilityModel::none;
      else if (permeability_model_str == "darcy phase change")
        permeability_model = PermeabilityModel::darcy_phase_change;
      else if (permeability_model_str == "carman-kozeny phase change")
        permeability_model = PermeabilityModel::carman_kozeny_phase_change;
      else
        throw std::logic_error(
          "Error, invalid permeability model. Options are <none|darcy phase change|carman-kozeny phase change>.");

      const std::string fluid_indicator_str =
        prm.get("Carman-Kozeny fluid with phase change");
      if (fluid_indicator_str == "fluid 0")
        fluid_with_phase_change = FluidIndicator::fluid0;
      else if (fluid_indicator_str == "fluid 1")
        fluid_with_phase_change = FluidIndicator::fluid1;
      else if (fluid_indicator_str == "both")
        fluid_with_phase_change = FluidIndicator::both;
      else
        throw std::logic_error(
          "Error, invalid fluid with phase change. Options are <fluid 0|fluid 1|both>.");

      carman_kozeny_permeability_area.resize(2);
      carman_kozeny_tolerance.resize(2);
      const std::string carman_kozeny_permeability_area_list =
        prm.get("Carman-Kozeny permeability area");
      std::vector<std::string> carman_kozeny_permeability_area_vec =
        Utilities::split_string_list(carman_kozeny_permeability_area_list);

      const std::string carman_kozeny_tolerance_list =
        prm.get("Carman-Kozeny division tolerance");
      std::vector<std::string> carman_kozeny_tolerance_vec =
        Utilities::split_string_list(carman_kozeny_tolerance_list);

      // Check that dimensions agree
      AssertThrow(carman_kozeny_permeability_area_vec.size() ==
                    carman_kozeny_tolerance_vec.size(),
                  dealii::ExcDimensionMismatch(
                    (carman_kozeny_permeability_area_vec.size()),
                    (carman_kozeny_tolerance_vec.size())));

      if (fluid_with_phase_change == FluidIndicator::fluid0)
        {
          AssertThrow(
            carman_kozeny_permeability_area_vec.size() == 1,
            ExcMessage(
              "The expected size of the 'Carman-Kozeny permeability area' list of 1 is not met."));
          AssertThrow(
            carman_kozeny_tolerance_vec.size() == 1,
            ExcMessage(
              "The expected size of the 'Carman-Kozeny division tolerance' list of 1 is not met."));

          carman_kozeny_permeability_area[0] =
            Utilities::string_to_double(carman_kozeny_permeability_area_vec[0]);
          carman_kozeny_tolerance[0] =
            Utilities::string_to_double(carman_kozeny_tolerance_vec[0]);

          // Assign default values to fluid 1
          carman_kozeny_permeability_area[1] = 1e-3;
          carman_kozeny_tolerance[1]         = 1e-3;
        }
      else if (fluid_with_phase_change == FluidIndicator::fluid1)
        {
          AssertThrow(
            carman_kozeny_permeability_area_vec.size() == 1,
            ExcMessage(
              "The expected size of the 'Carman-Kozeny permeability area' list of 1 is not met."));
          AssertThrow(
            carman_kozeny_tolerance_vec.size() == 1,
            ExcMessage(
              "The expected size of the 'Carman-Kozeny division tolerance' list of 1 is not met."));

          carman_kozeny_permeability_area[1] =
            Utilities::string_to_double(carman_kozeny_permeability_area_vec[0]);
          carman_kozeny_tolerance[1] =
            Utilities::string_to_double(carman_kozeny_tolerance_vec[0]);

          // Assign default values to fluid 0
          carman_kozeny_permeability_area[0] = 1e-3;
          carman_kozeny_tolerance[0]         = 1e-3;
        }
      else if (fluid_with_phase_change == FluidIndicator::both)
        {
          AssertThrow(
            carman_kozeny_permeability_area_vec.size() ==
              carman_kozeny_permeability_area.size(),
            ExcMessage(
              "The expected size of the 'Carman-Kozeny permeability area' list of 1 is not met."));
          AssertThrow(
            carman_kozeny_tolerance_vec.size() ==
              carman_kozeny_tolerance.size(),
            ExcMessage(
              "The expected size of the 'Carman-Kozeny division tolerance' list of 2 is not met."));

          for (unsigned int i = 0; i < carman_kozeny_permeability_area.size();
               ++i)
            {
              carman_kozeny_permeability_area[i] = Utilities::string_to_double(
                carman_kozeny_permeability_area_vec[i]);
              carman_kozeny_tolerance[i] =
                Utilities::string_to_double(carman_kozeny_tolerance_vec[i]);
            }
        }

      // Check that all values are strictly positive
      for (unsigned int i = 0; i < carman_kozeny_permeability_area.size(); ++i)
        {
          AssertThrow(
            carman_kozeny_permeability_area[i] > 0,
            ExcMessage(
              "The chosen value of 'Carman-Kozeny permeability area' has to be strictly positive."));

          AssertThrow(
            carman_kozeny_tolerance[i] > 0,
            ExcMessage(
              "The 'Carman-Kozeny division tolerance' should be strictly positive."));
        }

      enable_darcy_multiply_by_density =
        prm.get_bool("enable Darcy multiply by density");

      AssertThrow(
        (!enable_darcy_multiply_by_density ||
         permeability_model == PermeabilityModel::darcy_phase_change),
        ExcMessage(
          "Inconsistency in parameters, 'enable Darcy multiply by density' is set to 'true', but 'permeability model' is not set to 'darcy phase change'."));

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
        Patterns::Tools::Convert<unsigned int>::to_string(nb_particles),
        Patterns::Integer(),
        "The number of particles represented by IB. The maximal number of particles is equal to 10 when defined individually. If particles are loaded from a file, this parameter is overridden, and there is no limit to the number of particles.");

      prm.declare_entry(
        "assemble Navier-Stokes inside particles",
        Patterns::Tools::Convert<bool>::to_string(
          assemble_navier_stokes_inside),
        Patterns::Bool(),
        "Bool to define if you assemble the inside of particles with the NS equation.");

      prm.enter_subsection("extrapolation function");
      {
        prm.declare_entry(
          "stencil degree",
          Patterns::Tools::Convert<unsigned int>::to_string(stencil_degree),
          Patterns::Integer(1),
          "The polynomial degree used in the extrapolation function");
        prm.declare_alias("stencil degree", "stencil order", true);
        prm.declare_entry(
          "length ratio",
          Patterns::Tools::Convert<double>::to_string(length_ratio),
          Patterns::Double(),
          "The length ratio used to define the points for the IB stencil. See definition of epsilon_n in the paper on sharp IB.");
        prm.declare_entry(
          "enable extrapolation",
          Patterns::Tools::Convert<bool>::to_string(enable_extrapolation),
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
          Patterns::Tools::Convert<bool>::to_string(load_particles_from_file),
          Patterns::Bool(),
          "Bool to define if particles are loaded from an external file");
        prm.declare_entry("particles file",
                          particles_file,
                          Patterns::FileName(),
                          "The file name from which we load the particles");
        prm.leave_subsection();
      }

      prm.enter_subsection("local mesh refinement");
      {
        prm.declare_entry(
          "initial refinement",
          Patterns::Tools::Convert<unsigned int>::to_string(initial_refinement),
          Patterns::Integer(),
          "Number of refinements around the particles before the start of the simulation ");
        prm.declare_entry(
          "enable distance based coarsening",
          Patterns::Tools::Convert<bool>::to_string(enable_coarsening),
          Patterns::Bool(),
          "Enable distance-based coarsening away from immersed boundary particles");
        prm.declare_entry(
          "refine mesh inside radius factor",
          Patterns::Tools::Convert<double>::to_string(
            refinement_inside_distance_factor),
          Patterns::Double(),
          "The factor that multiplies the radius to define the inside bound for the refinement of the mesh");
        prm.declare_entry(
          "refine mesh outside radius factor",
          Patterns::Tools::Convert<double>::to_string(
            refinement_outside_distance_factor),
          Patterns::Double(),
          "The factor that multiplies the radius to define the outside bound for the refinement of the mesh");
        prm.declare_entry(
          "coarsen mesh outside radius factor",
          Patterns::Tools::Convert<double>::to_string(
            coarsening_distance_factor),
          Patterns::Double(),
          "The factor that multiplies the radius to define the outside bound beyond which distance-based coarsening is allowed");
        prm.declare_entry(
          "refinement zone extrapolation",
          Patterns::Tools::Convert<bool>::to_string(
            time_extrapolation_of_refinement_zone),
          Patterns::Bool(),
          "This parameter enables the extrapolation in time of the refinement zone. This means that it will try to refine where the particle will be at the end of the time step instead of the initial position.");
        prm.leave_subsection();
      }

      prm.enter_subsection("output");
      {
        prm.declare_entry(
          "calculate force",
          Patterns::Tools::Convert<bool>::to_string(calculate_force_ib),
          Patterns::Bool(),
          "Bool to define if the force is evaluated on each particle ");
        prm.declare_entry(
          "ib force output file",
          ib_force_output_file,
          Patterns::FileName(),
          "The name of the file where the data on the force of each particle is stored");
        prm.declare_entry(
          "ib particles pvd file",
          ib_particles_pvd_file,
          Patterns::FileName(),
          "The output files of the pvd data for the ib particles");
        prm.declare_entry(
          "print DEM",
          Patterns::Tools::Convert<bool>::to_string(print_dem),
          Patterns::Bool(),
          "Bool to define if particles' information are printed on the terminal when particles' time step is finished");
        prm.declare_entry(
          "enable extra sharp interface vtu output field",
          Patterns::Tools::Convert<bool>::to_string(
            enable_extra_sharp_interface_vtu_output_field),
          Patterns::Bool(),
          "This parameter enables the output of more information related to the particles in the vtu file.");
        prm.leave_subsection();
      }

      prm.enter_subsection("DEM");
      {
        prm.declare_entry(
          "contact search radius factor",
          Patterns::Tools::Convert<double>::to_string(
            contact_search_radius_factor),
          Patterns::Double(),
          "The factor that multiplies the radius to define the region of contact search around the particle");
        prm.declare_entry(
          "contact search frequency",
          Patterns::Tools::Convert<int>::to_string(contact_search_frequency),
          Patterns::Integer(),
          "The frequency of update in the contact candidates list");
        prm.declare_entry(
          "particle nonlinear tolerance",
          Patterns::Tools::Convert<double>::to_string(
            particle_nonlinear_tolerance),
          Patterns::Double(),
          "The nonlinear tolerance for the coupling of the particle dynamics and the fluid");
        prm.declare_entry("DEM coupling frequency",
                          Patterns::Tools::Convert<unsigned int>::to_string(
                            coupling_frequency),
                          Patterns::Integer(),
                          "The number of DEM time steps per CFD time step");
        prm.declare_entry("alpha",
                          Patterns::Tools::Convert<double>::to_string(alpha),
                          Patterns::Double(),
                          "relaxation parameter");

        prm.declare_entry("enable lubrication force",
                          Patterns::Tools::Convert<bool>::to_string(
                            enable_lubrication_force),
                          Patterns::Bool(),
                          "Bool to enable or disable the lubrication force");
        prm.declare_entry(
          "lubrication range max",
          Patterns::Tools::Convert<double>::to_string(lubrication_range_max),
          Patterns::Double(),
          "Gap require to consider the lubrication force. This value is multiplied the smallest cell size");
        prm.declare_entry(
          "lubrication range min",
          Patterns::Tools::Convert<double>::to_string(lubrication_range_min),
          Patterns::Double(),
          "Smallest gap considered for the lubrification force calculation. This value is multiplied by the smallest cell size");

        prm.declare_entry(
          "explicit contact impulsion",
          Patterns::Tools::Convert<bool>::to_string(
            explicit_contact_impulsion_calculation),
          Patterns::Bool(),
          "Bool to enable or disable the use of explicit contact impulsion evaluation in the resolution of the coupling of the particle. When it is set to true, this parameter results in the code only performing the DEM calculation once per CFD time step and using the resulting contact impulsion to evaluate all the other Newton's iterations. This reduces the number of times the DEM calculation is made.");

        prm.declare_entry(
          "explicit position integration",
          Patterns::Tools::Convert<bool>::to_string(
            explicit_position_integration_calculation),
          Patterns::Bool(),
          "Bool to enable or disable the explicit position integration. This means that the particle position is obtained directly by the integration of the previous velocities only. This avoids multiple cut cell mapping for each newton iteration. Note that this limits the order of convergence in time to one.");

        prm.declare_entry(
          "approximate radius for contact",
          Patterns::Tools::Convert<bool>::to_string(
            approximate_radius_for_contact),
          Patterns::Bool(),
          "Bool to turn on or off using the approximate radius of the particles during contact. If activated, the radius used in the contact calculation is constant and fixed to the effective radius of the shape. If not, the radius of curvature of the shape at the contact point is evaluated. For some shapes, this can be numerically expensive to evaluate.");



        prm.enter_subsection("wall physical properties");
        {
          prm.declare_entry(
            "wall youngs modulus",
            Patterns::Tools::Convert<double>::to_string(wall_youngs_modulus),
            Patterns::Double(),
            "The wall Young's modulus if IB particles are in contact with it");

          prm.declare_entry(
            "wall poisson ratio",
            Patterns::Tools::Convert<double>::to_string(wall_poisson_ratio),
            Patterns::Double(),
            "The wall poisson ratio if IB particles are in contact with it");

          prm.declare_entry(
            "wall rolling friction coefficient",
            Patterns::Tools::Convert<double>::to_string(
              wall_rolling_friction_coefficient),
            Patterns::Double(),
            "The wall rolling friction coefficient if IB particles are in contact with it");

          prm.declare_entry(
            "wall friction coefficient",
            Patterns::Tools::Convert<double>::to_string(
              wall_friction_coefficient),
            Patterns::Double(),
            "The wall friction coefficient if IB particles are in contact with it");

          prm.declare_entry(
            "wall restitution coefficient",
            Patterns::Tools::Convert<double>::to_string(
              wall_restitution_coefficient),
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
        stencil_degree       = prm.get_integer("stencil degree");
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
        enable_coarsening  = prm.get_bool("enable distance based coarsening");
        refinement_inside_distance_factor =
          prm.get_double("refine mesh inside radius factor");
        refinement_outside_distance_factor =
          prm.get_double("refine mesh outside radius factor");
        coarsening_distance_factor =
          prm.get_double("coarsen mesh outside radius factor");
        time_extrapolation_of_refinement_zone =
          prm.get_bool("refinement zone extrapolation");
        if (enable_coarsening)
          {
            AssertThrow(
              coarsening_distance_factor >= refinement_outside_distance_factor,
              ExcMessage(
                "The parameter 'coarsen mesh outside radius factor' must be greater than or equal to 'refine mesh outside radius factor'."));
          }
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
    const DynamicFlowControl defaults;
    prm.enter_subsection("flow control");
    {
      prm.declare_entry("enable",
                        Patterns::Tools::Convert<bool>::to_string(
                          defaults.enable_flow_control),
                        Patterns::Bool(),
                        "Enable flow rate control");
      prm.declare_entry("average velocity",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.average_velocity_0),
                        Patterns::Double(),
                        "The target average velocity");
      prm.declare_entry("inlet boundary id",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.boundary_flow_id),
                        Patterns::Integer(),
                        "Boundary id of the inlet flow");
      prm.declare_entry("flow direction",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          defaults.flow_direction),
                        Patterns::Integer(),
                        "Flow direction at flow inlet");
      prm.declare_entry("initial beta",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.beta_0),
                        Patterns::Double(),
                        "Beta coefficient value for the first step time");
      prm.declare_entry("alpha",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.alpha),
                        Patterns::Double(),
                        "Relaxation coefficient for flow controller");
      prm.declare_entry("enable beta particle",
                        Patterns::Tools::Convert<bool>::to_string(
                          defaults.enable_beta_particle),
                        Patterns::Bool(),
                        "Enable beta force for particles");
      prm.declare_entry("beta threshold",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.beta_threshold),
                        Patterns::Double(),
                        "Enable beta force for particles");
      prm.declare_entry(
        "verbosity",
        to_string(defaults.verbosity),
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

  namespace
  {
    std::string
    to_string(const Evaporation::EvaporativeMassFluxModelType type)
    {
      switch (type)
        {
          case Evaporation::EvaporativeMassFluxModelType::constant:
            return "constant";
          case Evaporation::EvaporativeMassFluxModelType::temperature_dependent:
            return "temperature_dependent";
        }
      Assert(false, ExcInternalError());
      return "";
    }
  } // namespace

  void
  Evaporation::declare_parameters(dealii::ParameterHandler &prm)
  {
    const Evaporation defaults;
    prm.enter_subsection("evaporation");
    {
      prm.declare_entry(
        "evaporation mass flux model",
        to_string(defaults.evaporative_mass_flux_model_type),
        Patterns::Selection("constant|temperature_dependent"),
        "Model used for the calculation of the evaporative mass flux"
        "Choices are <constant|temperature_dependent>.");
      prm.declare_entry(
        "enable evaporative cooling",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.enable_evaporation_cooling),
        Patterns::Bool(),
        "Enable the evaporative cooling at the free surface (gas/liquid interface) in the energy equation <true|false>");
      prm.declare_entry(
        "enable recoil pressure",
        Patterns::Tools::Convert<bool>::to_string(
          defaults.enable_recoil_pressure),
        Patterns::Bool(),
        "Enable the recoil pressure due to evaporation at the free surface (gas/liquid interface) in the momentum equation <true|false>");
      prm.declare_entry(
        "evaporation mass flux",
        Patterns::Tools::Convert<double>::to_string(
          defaults.evaporation_mass_flux),
        Patterns::Double(),
        "Evaporation mass flux used if the constant evaporation model is selected in M*L^-2*T^-1");
      prm.declare_entry(
        "evaporation coefficient",
        Patterns::Tools::Convert<double>::to_string(
          defaults.evaporation_coefficient),
        Patterns::Double(),
        "Evaporation coefficient corresponding to the ratio between the net mass flux (evaporation-condensation) and the mass flux of evaporation");
      prm.declare_entry(
        "recoil pressure coefficient",
        Patterns::Tools::Convert<double>::to_string(
          defaults.recoil_pressure_coefficient),
        Patterns::Double(),
        "Recoil pressure coefficient corresponding to the factor applied to the saturation pressure to compute the recoil pressure in an out of equilibrium evaporation");
      prm.declare_entry("molar mass",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.molar_mass),
                        Patterns::Double(),
                        "Molar mass of the material in M*N^-1");
      prm.declare_entry("boiling temperature",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.boiling_temperature),
                        Patterns::Double(),
                        "Boiling temperature in Theta");
      prm.declare_entry("evaporation latent heat",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.latent_heat_evaporation),
                        Patterns::Double(),
                        "Latent heat of evaporation in L^2*T^-2");
      prm.declare_entry("ambient pressure",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.ambient_pressure),
                        Patterns::Double(),
                        "Pressure of the ambient gas in M*L^-1*T^-2");
      prm.declare_entry("ambient gas density",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.ambient_gas_density),
                        Patterns::Double(),
                        "Ambient gas density in M*L^-3");
      prm.declare_entry("liquid density",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.liquid_density),
                        Patterns::Double(),
                        "Liquid density in M*L^-3");
      prm.declare_entry("universal gas constant",
                        Patterns::Tools::Convert<double>::to_string(
                          defaults.universal_gas_constant),
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
                        Patterns::Tools::Convert<bool>::to_string(enable),
                        Patterns::Bool(),
                        "Enable mortar interface <true|false>");
      prm.declare_entry("interface type",
                        interface_type == InterfaceType::circular ? "circular" :
                                                                    "linear",
                        Patterns::Selection("circular|linear"),
                        "Type of mortar interface"
                        "Choices are <circular|linear>.");
      rotor_mesh = std::make_shared<Mesh>();
      rotor_mesh->declare_parameters(prm);
      prm.declare_entry("rotor boundary id",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          rotor_boundary_id),
                        Patterns::Integer(),
                        "Rotor boundary ID # of the mortar matching interface");
      prm.declare_entry(
        "stator boundary id",
        Patterns::Tools::Convert<unsigned int>::to_string(stator_boundary_id),
        Patterns::Integer(),
        "Stator boundary ID # of the mortar matching interface");
      std::string default_entry_point = (dim == 2) ? "0., 0." : "0., 0., 0.";
      prm.declare_entry("center of rotation",
                        default_entry_point,
                        Patterns::List(Patterns::Double()),
                        "Center of rotation coordinates of rotor domain");

      prm.declare_entry("rotation axis direction",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          rotation_axis_direction),
                        Patterns::Integer(0, 2),
                        "Direction of the rotation axis. Choices are <0|1|2>.");

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
                        Patterns::Tools::Convert<double>::to_string(sip_factor),
                        Patterns::Double(),
                        "Penalty factor for mortar elements");
      prm.declare_entry("oversampling factor",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          oversampling_factor),
                        Patterns::Integer(),
                        "Oversampling factor for quadrature points");
      prm.declare_entry("radius tolerance",
                        Patterns::Tools::Convert<double>::to_string(
                          radius_tolerance),
                        Patterns::Double(),
                        "Tolerance used for the interface radius computation");
      prm.declare_entry("cell weight",
                        Patterns::Tools::Convert<unsigned int>::to_string(
                          cell_weight),
                        Patterns::Integer(),
                        "Cell weight for load balancing of mortar cells");
      prm.declare_entry(
        "verbosity",
        to_string(verbosity),
        Patterns::Selection("quiet|verbose|extra verbose"),
        "State whether from the mortar information should be printed "
        "Choices are <quiet|verbose|extra verbose>.");
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
      {
        const std::string op = prm.get("interface type");
        if (op == "circular")
          interface_type = InterfaceType::circular;
        else if (op == "linear")
          interface_type = InterfaceType::linear;
        else
          AssertThrow(
            false,
            ExcMessage(
              "Error, invalid mortar interface type. Current choices are <circular|linear>."));
      }
      rotor_mesh->parse_parameters(prm);
      rotor_boundary_id  = prm.get_integer("rotor boundary id");
      stator_boundary_id = prm.get_integer("stator boundary id");
      center_of_rotation =
        value_string_to_tensor<dim>(prm.get("center of rotation"));
      rotation_axis_direction = prm.get_integer("rotation axis direction");
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
      radius_tolerance    = prm.get_double("radius tolerance");
      cell_weight         = prm.get_integer("cell weight");

      // Enable printing of mortar information
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;
      if (op == "extra verbose")
        verbosity = Verbosity::extra_verbose;
    }
    prm.leave_subsection();
  }

  // Explicitly instantiate template classes and structs
  template class Laser<2>;
  template class Laser<3>;
  template class PostProcessing<2>;
  template class PostProcessing<3>;
  template class IBParticles<2>;
  template class IBParticles<3>;
  template struct ConstrainSolidDomain<2>;
  template struct ConstrainSolidDomain<3>;
  template struct Mortar<2>;
  template struct Mortar<3>;

} // namespace Parameters
