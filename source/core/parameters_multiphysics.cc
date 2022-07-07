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

    // subparameter for heat_transfer
    viscous_dissipation = prm.get_bool("viscous dissipation");
    buoyancy_force      = prm.get_bool("buoyancy force");
  }
  prm.leave_subsection();
  vof_parameters.parse_parameters(prm);
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

    // Error definitions
    if (sharpening.type == Parameters::SharpeningType::adaptative)
      Assert(conservation.monitoring == true,
             AdaptativeSharpeningError(conservation.monitoring));
  }
  prm.leave_subsection();
}

void
Parameters::VOF_MassConservation::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("mass conservation");
  {
    prm.declare_entry(
      "skip mass conservation in fluid 0",
      "false",
      Patterns::Bool(),
      "Enable skipping mass conservation in fluid 0 <true|false>."
      "Can be used to improve the wetting mechanism, along with a small time step."
      "See documentation for further details.");

    prm.declare_entry(
      "skip mass conservation in fluid 1",
      "false",
      Patterns::Bool(),
      "Enable skipping mass conservation in fluid 1 <true|false>."
      "Can be used to improve the wetting mechanism, along with a small time step."
      "See documentation for further details.");

    prm.declare_entry(
      "monitoring",
      "false",
      Patterns::Bool(),
      "Enable conservation monitoring in free surface calculation <true|false>");

    prm.declare_entry(
      "fluid monitored",
      "1",
      Patterns::Integer(),
      "Index of the fluid which conservation is monitored <0|1>");

    prm.declare_entry(
      "tolerance",
      "1e-2",
      Patterns::Double(),
      "Tolerance on the mass conservation of the monitored fluid, used with adaptative sharpening");

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
    skip_mass_conservation_fluid_0 =
      prm.get_bool("skip mass conservation in fluid 0");
    skip_mass_conservation_fluid_1 =
      prm.get_bool("skip mass conservation in fluid 1");
    monitoring         = prm.get_bool("monitoring");
    id_fluid_monitored = prm.get_integer("fluid monitored");
    tolerance          = prm.get_double("tolerance");

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
      "5",
      Patterns::Integer(),
      "Maximum number of iteration in the binary search algorithm");

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
      "enable",
      "false",
      Patterns::Bool(),
      "Enable peeling/wetting mechanism in free surface simulation <true|false>");

    prm.declare_entry("peeling pressure value",
                      "-0.05",
                      Patterns::Double(),
                      "Value (Double) for pressure value at bc below which "
                      "peeling of the higher density fluid can occur.");

    prm.declare_entry("peeling pressure gradient",
                      "-1e-3",
                      Patterns::Double(),
                      "Value (Double) for pressure gradient at bc below which "
                      "peeling of the higher density fluid can occur.");

    prm.declare_entry("wetting pressure value",
                      "0.05",
                      Patterns::Double(),
                      "Value (Double) for pressure value at bc above which "
                      "wetting of the lower density fluid can occur.");

    prm.declare_entry(
      "wetting phase distance",
      "0",
      Patterns::Double(),
      "Value (Double) for wetting distance at bc, "
      "distance (on the phase value) from the interface above which wetting can occur. "
      "For wetting phase distance>0, the wetting area is larger than "
      "the area occupied by the higher density fluid.");

    prm.declare_entry(
      "diffusivity",
      "0",
      Patterns::Double(),
      "Diffusivity (diffusion coefficient in L^2/s) in the VOF transport equation. "
      "Default value is 0 to have pure advection. Use this parameter, "
      "along with interface sharpening, to improve the wetting mechanism. "
      "See documentation for more details.");

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
    enable                 = prm.get_bool("enable");
    peeling_p_value        = prm.get_double("peeling pressure value");
    peeling_grad_p         = prm.get_double("peeling pressure gradient");
    wetting_p_value        = prm.get_double("wetting pressure value");
    wetting_phase_distance = prm.get_double("wetting phase distance");
    diffusivity            = prm.get_double("diffusivity");

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

    prm.declare_entry("surface tension coefficient",
                      "0.0",
                      Patterns::Double(),
                      "Surface tension coefficient");

    prm.declare_entry("output auxiliary fields",
                      "false",
                      Patterns::Bool(),
                      "Output the phase fraction gradient and curvature");

    prm.declare_entry(
      "phase fraction gradient filter",
      "0.5",
      Patterns::Double(),
      "The filter value for phase fraction gradient calculations to damp high-frequency errors");

    prm.declare_entry(
      "curvature filter",
      "0.5",
      Patterns::Double(),
      "The filter value for curvature calculations to damp high-frequency errors");

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
    // Surface tension coefficient
    surface_tension_coef = prm.get_double("surface tension coefficient");
    phase_fraction_gradient_filter_value =
      prm.get_double("phase fraction gradient filter");
    curvature_filter_value = prm.get_double("curvature filter");

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
