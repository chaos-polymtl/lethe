#include <core/parameters_multiphysics.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>

DeclException1(
  SharpeningThresholdError,
  double,
  << "Sharpening threshold : " << arg1 << " is smaller than 0 or larger than 1."
  << std::endl
  << "Interface sharpening model requires a sharpening threshold between 0 and 1.");

DeclException1(
  SharpeningThresholdErrorMaxDeviation,
  double,
  << "Sharpening threshold max deviation : " << arg1
  << " is smaller than 0 or larger than 0.5." << std::endl
  << "Adaptative interface sharpening requires a maximum deviation of the"
  << " sharpening threshold between 0.0 and 0.5. See documentation for further details");

DeclException1(
  SharpeningFrequencyError,
  int,
  << "Sharpening frequency : " << arg1 << " is equal or smaller than 0."
  << std::endl
  << "Interface sharpening model requires an integer sharpening frequency larger than 0.");

DeclException1(
  AdaptativeSharpeningError,
  bool,
  << "Sharpening type is set to 'adaptative' but monitoring is : " << arg1
  << std::endl
  << "Adaptative sharpening requires to set 'monitoring = true', and to define"
  << " the 'fluid monitored' and the 'tolerance' to reach. See documentation for further details.");

DeclException1(
  MonitoringConservationError,
  bool,
  << "Sharpening type is set to 'adaptative' and monitoring is set to " << arg1
  << " but " << std::endl
  << "conservative fluid (choices are <fluid 0, fluid 1, both>) does not include the monitored fluid. "
  << std::endl
  << "Conservation must be solved on the monitored fluid for adaptative sharpening to work properly.");

void
Parameters::Multiphysics::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("multiphysics");
  {
    prm.declare_entry("fluid dynamics",
                      "true",
                      Patterns::Bool(),
                      "Fluid flow calculation <true|false>");

    prm.declare_entry("heat transfer",
                      "false",
                      Patterns::Bool(),
                      "Thermic calculation <true|false>");

    prm.declare_entry("tracer",
                      "false",
                      Patterns::Bool(),
                      "Passive tracer calculation <true|false>");

    prm.declare_entry("VOF",
                      "false",
                      Patterns::Bool(),
                      "VOF calculation <true|false>");

    prm.declare_entry("cahn hilliard",
                      "false",
                      Patterns::Bool(),
                      "Cahn-Hilliard calculation <true|false>");

    prm.declare_entry(
      "use time average velocity field",
      "false",
      Patterns::Bool(),
      "Use the average velocity field in subphysics, instead of the present velocity field <true|false>");

    // subparameters for heat_transfer
    prm.declare_entry("viscous dissipation",
                      "false",
                      Patterns::Bool(),
                      "Viscous dissipation in heat equation <true|false>");

    prm.declare_entry("buoyancy force",
                      "false",
                      Patterns::Bool(),
                      "Buoyant force calculation <true|false>");
  }
  prm.leave_subsection();

  vof_parameters.declare_parameters(prm);
  ch_parameters.declare_parameters(prm);
}

void
Parameters::Multiphysics::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("multiphysics");
  {
    fluid_dynamics = prm.get_bool("fluid dynamics");
    heat_transfer  = prm.get_bool("heat transfer");
    tracer         = prm.get_bool("tracer");
    VOF            = prm.get_bool("VOF");
    cahn_hilliard  = prm.get_bool("cahn hilliard");
    use_time_average_velocity_field =
      prm.get_bool("use time average velocity field");

    // subparameter for heat_transfer
    viscous_dissipation = prm.get_bool("viscous dissipation");
    buoyancy_force      = prm.get_bool("buoyancy force");
  }
  prm.leave_subsection();
  vof_parameters.parse_parameters(prm);
  ch_parameters.parse_parameters(prm);
}

void
Parameters::VOF::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("VOF");
  {
    conservation.declare_parameters(prm);
    sharpening.declare_parameters(prm);
    peeling_wetting.declare_parameters(prm);
    surface_tension_force.declare_parameters(prm);
    phase_filter.declare_parameters(prm);

    prm.declare_entry("viscous dissipative fluid",
                      "fluid 1",
                      Patterns::Selection("fluid 0|fluid 1|both"),
                      "Fluid to which the viscous dissipation is applied "
                      "in the heat equation <fluid 0|fluid 1|both>");

    prm.declare_entry(
      "diffusivity",
      "0",
      Patterns::Double(),
      "Diffusivity (diffusion coefficient in L^2/s) in the VOF transport equation. "
      "Default value is 0 to have pure advection.");
  }
  prm.leave_subsection();
}

void
Parameters::VOF::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("VOF");
  {
    conservation.parse_parameters(prm);
    sharpening.parse_parameters(prm);
    peeling_wetting.parse_parameters(prm);
    surface_tension_force.parse_parameters(prm);
    phase_filter.parse_parameters(prm);

    // Viscous dissipative fluid
    const std::string op = prm.get("viscous dissipative fluid");
    if (op == "fluid 1")
      viscous_dissipative_fluid = Parameters::FluidIndicator::fluid1;
    else if (op == "fluid 0")
      viscous_dissipative_fluid = Parameters::FluidIndicator::fluid0;
    else if (op == "both")
      viscous_dissipative_fluid = Parameters::FluidIndicator::both;
    else
      throw(std::runtime_error("Invalid viscous dissipative fluid. "
                               "Options are 'fluid 0', 'fluid 1' or 'both'."));

    diffusivity = prm.get_double("diffusivity");

    // Error definitions
    if (sharpening.type == Parameters::SharpeningType::adaptative)
      Assert(conservation.monitoring == true,
             AdaptativeSharpeningError(conservation.monitoring));

    if (sharpening.type == Parameters::SharpeningType::adaptative)

      Assert(
        (conservation.conservative_fluid == Parameters::FluidIndicator::both) or
          (conservation.conservative_fluid == conservation.monitored_fluid),
        MonitoringConservationError(conservation.monitoring));
  }
  prm.leave_subsection();
}

void
Parameters::VOF_MassConservation::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("mass conservation");
  {
    prm.declare_entry(
      "monitoring",
      "false",
      Patterns::Bool(),
      "Enable conservation monitoring in free surface calculation <true|false>");

    prm.declare_entry(
      "tolerance",
      "1e-6",
      Patterns::Double(),
      "Tolerance on the mass conservation of the monitored fluid, used with adaptative sharpening");

    prm.declare_entry(
      "conservative fluid",
      "both",
      Patterns::Selection("fluid 0|fluid 1|both"),
      "Fluid for which conservation is solved <fluid 0|fluid 1|both>. "
      "Conservation on one fluid instead of both can be used to improve the wetting mechanism. "
      "See documentation for further details.");

    prm.declare_entry(
      "monitored fluid",
      "fluid 1",
      Patterns::Selection("fluid 0|fluid 1"),
      "Fluid for which conservation is monitored <fluid 0|fluid 1>.");

    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose|extra verbose"),
      "States whether the mass conservation data should be printed "
      "Choices are <quiet|verbose|extra verbose>.");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_MassConservation::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("mass conservation");
  {
    monitoring = prm.get_bool("monitoring");
    tolerance  = prm.get_double("tolerance");

    // Conservative fluid
    const std::string op_cf = prm.get("conservative fluid");
    if (op_cf == "fluid 1")
      conservative_fluid = Parameters::FluidIndicator::fluid1;
    else if (op_cf == "fluid 0")
      conservative_fluid = Parameters::FluidIndicator::fluid0;
    else if (op_cf == "both")
      conservative_fluid = Parameters::FluidIndicator::both;
    else
      throw(std::runtime_error("Invalid conservative fluid. "
                               "Options are 'fluid 0', 'fluid 1' or 'both'."));

    // Monitored fluid
    const std::string op_mf = prm.get("monitored fluid");
    if (op_mf == "fluid 1")
      monitored_fluid = Parameters::FluidIndicator::fluid1;
    else if (op_mf == "fluid 0")
      monitored_fluid = Parameters::FluidIndicator::fluid0;
    else if (op_mf == "both")
      monitored_fluid = Parameters::FluidIndicator::both;
    else
      throw(std::runtime_error("Invalid monitored fluid. "
                               "Options are 'fluid 0' or 'fluid 1'."));

    // Verbosity
    const std::string op_v = prm.get("verbosity");
    if (op_v == "verbose")
      verbosity = Parameters::Verbosity::verbose;
    else if (op_v == "quiet")
      verbosity = Parameters::Verbosity::quiet;
    else if (op_v == "extra verbose")
      verbosity = Parameters::Verbosity::extra_verbose;
    else
      throw(std::runtime_error("Invalid verbosity level"));
  }
  prm.leave_subsection();
}

void
Parameters::VOF_InterfaceSharpening::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("interface sharpening");
  {
    prm.declare_entry("enable",
                      "false",
                      Patterns::Bool(),
                      "Enable interface sharpening <true|false>");

    prm.declare_entry(
      "type",
      "constant",
      Patterns::Selection("constant|adaptative"),
      "VOF interface sharpening type, "
      "if constant the sharpening threshold is the same throughout the simulation, "
      "if adaptative the sharpening threshold is determined by binary search, "
      "to ensure mass conservation of the monitored phase");

    // Parameters for constant sharpening
    prm.declare_entry(
      "threshold",
      "0.5",
      Patterns::Double(),
      "Interface sharpening threshold that represents the phase fraction at which "
      "the interphase is considered located");

    // Parameters for adaptative sharpening
    prm.declare_entry(
      "threshold max deviation",
      "0.20",
      Patterns::Double(),
      "Maximum deviation (from the base value of 0.5) considered in the search "
      "algorithm to ensure mass conservation. "
      "A threshold max deviation of 0.20 results in a search interval from 0.30 to 0.70");

    prm.declare_entry(
      "max iterations",
      "20",
      Patterns::Integer(),
      "Maximum number of iteration in the bissection algorithm that ensures mass conservation");

    // This parameter must be larger than 1 for interface sharpening. Choosing
    // values less than 1 leads to interface smoothing instead of sharpening.
    prm.declare_entry(
      "interface sharpness",
      "2",
      Patterns::Double(),
      "Sharpness of the moving interface (parameter alpha in the interface sharpening model)");

    prm.declare_entry("frequency",
                      "10",
                      Patterns::Integer(),
                      "VOF interface sharpening frequency");

    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose|extra verbose"),
      "States whether the interface sharpening calculations should be printed "
      "Choices are <quiet|verbose>.");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_InterfaceSharpening::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("interface sharpening");
  {
    enable              = prm.get_bool("enable");
    interface_sharpness = prm.get_double("interface sharpness");
    frequency           = prm.get_integer("frequency");

    // Sharpening type
    const std::string t = prm.get("type");
    if (t == "constant")
      type = Parameters::SharpeningType::constant;
    if (t == "adaptative")
      type = Parameters::SharpeningType::adaptative;

    // Parameters for constant sharpening
    threshold = prm.get_double("threshold");

    // Parameters for adaptative sharpening
    threshold_max_deviation = prm.get_double("threshold max deviation");
    max_iterations          = prm.get_double("max iterations");

    // Error definitions
    Assert(threshold > 0.0 && threshold < 1.0,
           SharpeningThresholdError(threshold));

    Assert(threshold_max_deviation > 0.0 && threshold_max_deviation < 0.5,
           SharpeningThresholdErrorMaxDeviation(threshold_max_deviation));

    Assert(frequency > 0, SharpeningFrequencyError(frequency));

    // Verbosity
    const std::string op = prm.get("verbosity");
    if (op == "verbose")
      verbosity = Parameters::Verbosity::verbose;
    else if (op == "quiet")
      verbosity = Parameters::Verbosity::quiet;
    else if (op == "extra verbose")
      verbosity = Parameters::Verbosity::extra_verbose;
    else
      throw(std::runtime_error("Invalid verbosity level"));
  }
  prm.leave_subsection();
}

void
Parameters::VOF_PeelingWetting::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("peeling wetting");
  {
    prm.declare_entry(
      "enable peeling",
      "false",
      Patterns::Bool(),
      "Enable peeling mechanism in free surface simulation <true|false>");

    prm.declare_entry(
      "enable wetting",
      "false",
      Patterns::Bool(),
      "Enable wetting mechanism in free surface simulation <true|false>");

    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose"),
      "State whether from the number of wet and peeled cells should be printed. "
      "Choices are <quiet|verbose>.");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_PeelingWetting::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("peeling wetting");
  {
    enable_peeling = prm.get_bool("enable peeling");
    enable_wetting = prm.get_bool("enable wetting");

    const std::string op = prm.get("verbosity");
    if (op == "verbose")
      verbosity = Parameters::Verbosity::verbose;
    else if (op == "quiet")
      verbosity = Parameters::Verbosity::quiet;
    else
      throw(std::runtime_error("Invalid verbosity level"));
  }
  prm.leave_subsection();
}

void
Parameters::VOF_SurfaceTensionForce::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("surface tension force");
  {
    prm.declare_entry("enable",
                      "false",
                      Patterns::Bool(),
                      "Enable surface tension force calculation <true|false>");

    prm.declare_entry("output auxiliary fields",
                      "false",
                      Patterns::Bool(),
                      "Output the phase fraction gradient and curvature");

    prm.declare_entry(
      "phase fraction gradient diffusion factor",
      "4",
      Patterns::Double(),
      "Factor applied to the filter for phase fraction gradient calculations to damp high-frequency errors");

    prm.declare_entry(
      "curvature diffusion factor",
      "1",
      Patterns::Double(),
      "Factor applied to the filter for curvature calculations to damp high-frequency errors");

    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose"),
      "State whether the output from the surface tension force calculations should be printed "
      "Choices are <quiet|verbose>.");

    prm.enter_subsection("marangoni effect");
    {
      prm.declare_entry("enable",
                        "false",
                        Patterns::Bool(),
                        "Enable marangoni effect calculation <true|false>");

      prm.declare_entry("surface tension gradient",
                        "0.0",
                        Patterns::Double(),
                        "Surface tension gradient with respect to temperature");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}

void
Parameters::VOF_SurfaceTensionForce::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("surface tension force");
  {
    enable = prm.get_bool("enable");
    phase_fraction_gradient_diffusion_factor =
      prm.get_double("phase fraction gradient diffusion factor");
    curvature_diffusion_factor = prm.get_double("curvature diffusion factor");

    output_vof_auxiliary_fields = prm.get_bool("output auxiliary fields");

    const std::string op = prm.get("verbosity");
    if (op == "verbose")
      verbosity = Parameters::Verbosity::verbose;
    else if (op == "quiet")
      verbosity = Parameters::Verbosity::quiet;
    else
      throw(std::runtime_error("Invalid verbosity level"));

    prm.enter_subsection("marangoni effect");
    {
      enable_marangoni_effect = prm.get_bool("enable");

      // Surface tension gradient
      surface_tension_gradient = prm.get_double("surface tension gradient");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}

void
Parameters::VOF_PhaseFilter::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("phase filtration");
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection("none|tanh"),
      "VOF phase filtration type, "
      "if <none> is selected, the phase won't be filtered"
      "if <tanh> is selected, the filtered phase will be a result of the "
      "following function: \\alpha_f = 0.5 \\tanh(\\beta(\\alpha-0.5)) + 0.5; "
      "where \\beta is a parameter influencing the interface thickness that "
      "must be defined");
    prm.declare_entry(
      "beta",
      "20",
      Patterns::Double(),
      "This parameter appears in the tanh filter function. It influence "
      "the thickness and the shape of the interface. For higher values of "
      "beta, a thinner and 'sharper/pixelated' interface will be seen.");
    prm.declare_entry("verbosity",
                      "quiet",
                      Patterns::Selection("quiet|verbose|extra verbose"),
                      "States whether the filtered data should be printed "
                      "Choices are <quiet|verbose>.");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_PhaseFilter::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("phase filtration");
  {
    // filter type
    const std::string t = prm.get("type");
    if (t == "none")
      type = Parameters::FilterType::none;
    else if (t == "tanh")
      type = Parameters::FilterType::tanh;
    else
      throw(std::logic_error(
        "Error, invalid filter type. Choices are 'none' or 'tanh'"));

    // beta
    beta = prm.get_double("beta");

    // Verbosity
    const std::string filter_v = prm.get("verbosity");
    if (filter_v == "verbose")
      verbosity = Parameters::Verbosity::verbose;
    else if (filter_v == "quiet")
      verbosity = Parameters::Verbosity::quiet;
    else
      throw(std::logic_error("Invalid verbosity level"));
  }
  prm.leave_subsection();
}

void
Parameters::CahnHilliard::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("cahn hilliard");
  {
    prm.declare_entry("well height",
                      "1",
                      Patterns::Double(),
                      "Potential height well for the Cahn-Hilliard equations.");

    prm.enter_subsection("epsilon");
    {
      prm.declare_entry(
        "method",
        "automatic",
        Patterns::Selection("automatic|manual"),
        "Epsilon is either set to two times the characteristic length (automatic) of the element or user defined on all the domain (manual)");

      prm.declare_entry(
        "value",
        "1.0",
        Patterns::Double(),
        "Parameter linked to the interface thickness. Should always be bigger than the characteristic size of the smallest element");
    }
    prm.leave_subsection();

    prm.enter_subsection("mobility");
    {
      prm.declare_entry("model",
                        "constant",
                        Patterns::Selection("constant|quadratic|quartic"),
                        "Mobility model for the Cahn-Hilliard equations"
                        "Choices are <constant|quadratic|quartic>.");

      prm.declare_entry("mobility constant",
                        "1",
                        Patterns::Double(),
                        "Mobility constant for the Cahn-Hilliard equations");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}

void
Parameters::CahnHilliard::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("cahn hilliard");
  {
    well_height = prm.get_double("well height");

    prm.enter_subsection("epsilon");
    {
      const std::string op_epsilon = prm.get("method");
      if (op_epsilon == "automatic")
        {
          epsilon_set_method = Parameters::EpsilonSetStrategy::automatic;
        }
      else if (op_epsilon == "manual")
        {
          epsilon_set_method = Parameters::EpsilonSetStrategy::manual;
        }
      else
        throw(std::runtime_error("Invalid epsilon setting strategy. "
                                 "Options are 'automatic' or 'manual'."));

      epsilon = prm.get_double("value");
    }
    prm.leave_subsection();

    prm.enter_subsection("mobility");
    {
      const std::string op = prm.get("model");
      if (op == "constant")
        mobility_model = Parameters::MobilityModel::constant;
      else if (op == "quartic")
        mobility_model = Parameters::MobilityModel::quartic;
      else
        throw(std::runtime_error("Invalid mobility model "
                                 "Options are 'constant' or 'quartic'."));

      mobility_constant = prm.get_double("mobility constant");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}
