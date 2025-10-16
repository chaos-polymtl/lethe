// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_multiphysics.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>

DeclException1(
  SharpeningThresholdError,
  double,
  << "Sharpening threshold : " << arg1 << " is smaller than 0 or larger than 1."
  << std::endl
  << "Projection-based interface sharpening model requires a sharpening threshold between 0 and 1.");

DeclException1(
  SharpeningThresholdErrorMaxDeviation,
  double,
  << "Sharpening threshold max deviation : " << arg1
  << " is smaller than 0 or larger than 0.5." << std::endl
  << "Adaptive projection-based interface sharpening requires a maximum deviation of the"
  << " sharpening threshold between 0.0 and 0.5. See documentation for further details");

DeclException1(
  RegularizationMethodFrequencyError,
  int,
  << "Regularization method frequency : " << arg1
  << " is equal or smaller than 0." << std::endl
  << "Interface regularization method requires an frequency larger than 0.");

void
Parameters::Multiphysics::declare_parameters(ParameterHandler &prm) const
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
      "electromagnetics",
      "false",
      Patterns::Bool(),
      "Time harmonic electromagnetics calculation <true|false>");

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
  cahn_hilliard_parameters.declare_parameters(prm);
}

void
Parameters::Multiphysics::parse_parameters(ParameterHandler     &prm,
                                           const Dimensionality &dimensions)
{
  prm.enter_subsection("multiphysics");
  {
    fluid_dynamics   = prm.get_bool("fluid dynamics");
    heat_transfer    = prm.get_bool("heat transfer");
    tracer           = prm.get_bool("tracer");
    VOF              = prm.get_bool("VOF");
    cahn_hilliard    = prm.get_bool("cahn hilliard");
    electromagnetics = prm.get_bool("electromagnetics");

    // subparameters for heat_transfer
    viscous_dissipation = prm.get_bool("viscous dissipation");
    buoyancy_force      = prm.get_bool("buoyancy force");
  }
  prm.leave_subsection();
  vof_parameters.parse_parameters(prm);
  cahn_hilliard_parameters.parse_parameters(prm, dimensions);
}

void
Parameters::VOF::declare_parameters(ParameterHandler &prm) const
{
  prm.enter_subsection("VOF");
  {
    regularization_method.declare_parameters(prm);
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

    prm.declare_entry(
      "compressible",
      "false",
      Patterns::Bool(),
      "Enable phase compressibility in the VOF equation. This leads to the inclusion of the phase * div(u) term in the VOF conservation equation. "
      "It should be set to false when the phases are incompressible");
  }
  prm.leave_subsection();
}

void
Parameters::VOF::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("VOF");
  {
    regularization_method.parse_parameters(prm);
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

    compressible = prm.get_bool("compressible");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_RegularizationMethod::declare_parameters(
  ParameterHandler &prm) const
{
  prm.enter_subsection("interface regularization method");
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection(
        "none|projection-based interface sharpening|algebraic interface reinitialization|geometric interface reinitialization"),
      "VOF interface regularization method");

    prm.declare_entry(
      "frequency",
      "10",
      Patterns::Integer(),
      "Reinitialization frequency (number of time-steps) at which the "
      "interface regularization process will be applied to the VOF "
      "phase fraction field.");
    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose|extra verbose"),
      "States whether the output from the interface regularization method "
      "should be printed."
      "Choices are <quiet|verbose|extra verbose>.");

    sharpening.declare_parameters(prm);
    algebraic_interface_reinitialization.declare_parameters(prm);
    geometric_interface_reinitialization.declare_parameters(prm);
  }
  prm.leave_subsection();
}

void
Parameters::VOF_RegularizationMethod::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("interface regularization method");
  {
    const std::string t = prm.get("type");
    if (t == "none")
      this->regularization_method_type =
        Parameters::RegularizationMethodType::none;
    else if (t == "projection-based interface sharpening")
      {
        regularization_method_type =
          Parameters::RegularizationMethodType::sharpening;
        sharpening.enable = true;
      }
    else if (t == "algebraic interface reinitialization")
      {
        this->regularization_method_type =
          Parameters::RegularizationMethodType::algebraic;
        algebraic_interface_reinitialization.enable = true;
      }
    else if (t == "geometric interface reinitialization")
      {
        this->regularization_method_type =
          Parameters::RegularizationMethodType::geometric;
        geometric_interface_reinitialization.enable = true;
      }
    else
      throw(
        std::runtime_error("Invalid interface regularization method type!"));

    this->frequency = prm.get_integer("frequency");
    Assert(this->frequency > 0,
           RegularizationMethodFrequencyError(this->frequency));

    const std::string op2 = prm.get("verbosity");
    if (op2 == "quiet")
      this->verbosity = Parameters::Verbosity::quiet;
    else if (op2 == "verbose")
      this->verbosity = Parameters::Verbosity::verbose;
    else if (op2 == "extra verbose")
      this->verbosity = Parameters::Verbosity::extra_verbose;
    else
      throw(std::invalid_argument("Invalid verbosity level\n "
                                  "Options are: \n"
                                  " <quiet|verbose|extra verbose>"));

    this->sharpening.parse_parameters(prm);
    this->algebraic_interface_reinitialization.parse_parameters(prm);
    this->geometric_interface_reinitialization.parse_parameters(prm);
  }
  prm.leave_subsection();
}

void
Parameters::VOF_InterfaceSharpening::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("projection-based interface sharpening");
  {
    prm.declare_entry(
      "type",
      "constant",
      Patterns::Selection("constant|adaptive"),
      "VOF interface sharpening type, "
      "if constant the sharpening threshold is the same throughout the simulation, "
      "if adaptive the sharpening threshold is determined by binary search, "
      "to ensure mass conservation of the monitored phase");

    // Parameters for constant sharpening
    prm.declare_entry(
      "threshold",
      "0.5",
      Patterns::Double(),
      "Interface sharpening threshold that represents the phase fraction at which "
      "the interphase is considered located");

    // Parameters for adaptive sharpening
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

    prm.declare_entry(
      "monitoring",
      "false",
      Patterns::Bool(),
      "Enable conservation monitoring in multiphase fluid simulations <true|false>");

    prm.declare_entry(
      "tolerance",
      "1e-6",
      Patterns::Double(),
      "Tolerance on the mass conservation of the monitored fluid, used with adaptive sharpening");

    prm.declare_entry(
      "monitored fluid",
      "fluid 1",
      Patterns::Selection("fluid 0|fluid 1"),
      "Fluid for which conservation is monitored <fluid 0|fluid 1>, used with adaptive sharpening.");

    // This parameter must be larger than 1 for interface sharpening. Choosing
    // values less than 1 leads to interface smoothing instead of sharpening.
    prm.declare_entry(
      "interface sharpness",
      "2",
      Patterns::Double(),
      "Sharpness of the moving interface (parameter alpha in the interface sharpening model)");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_InterfaceSharpening::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("projection-based interface sharpening");
  {
    interface_sharpness = prm.get_double("interface sharpness");

    // Sharpening type
    const std::string t = prm.get("type");
    if (t == "constant")
      type = Parameters::SharpeningType::constant;
    if (t == "adaptive")
      type = Parameters::SharpeningType::adaptive;

    // Parameters for constant sharpening
    threshold = prm.get_double("threshold");

    // Parameters for adaptive sharpening
    threshold_max_deviation = prm.get_double("threshold max deviation");
    max_iterations          = prm.get_integer("max iterations");
    monitoring              = prm.get_bool("monitoring");
    tolerance               = prm.get_double("tolerance");

    // Monitored fluid
    const std::string op_mf = prm.get("monitored fluid");
    if (op_mf == "fluid 1")
      monitored_fluid = Parameters::FluidIndicator::fluid1;
    else if (op_mf == "fluid 0")
      monitored_fluid = Parameters::FluidIndicator::fluid0;
    else
      throw(std::runtime_error("Invalid monitored fluid. "
                               "Options are 'fluid 0' or 'fluid 1'."));

    // Error definitions
    Assert(threshold > 0.0 && threshold < 1.0,
           SharpeningThresholdError(threshold));

    Assert(threshold_max_deviation > 0.0 && threshold_max_deviation < 0.5,
           SharpeningThresholdErrorMaxDeviation(threshold_max_deviation));
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


    prm.declare_entry("enable marangoni effect",
                      "false",
                      Patterns::Bool(),
                      "Enable marangoni effect calculation <true|false>");
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

    enable_marangoni_effect = prm.get_bool("enable marangoni effect");
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
Parameters::VOF_AlgebraicInterfaceReinitialization::declare_parameters(
  dealii::ParameterHandler &prm)
{
  prm.enter_subsection("algebraic interface reinitialization");
  {
    prm.declare_entry(
      "output reinitialization steps",
      "false",
      Patterns::Bool(),
      "Enables pvtu format outputs of the algebraic interface reinitialization "
      "steps <true|false>");
    prm.declare_entry(
      "diffusivity multiplier",
      "1.",
      Patterns::Double(),
      "Factor that multiplies the mesh-size in the mesh-dependant diffusion "
      "coefficient of the algebraic interface reinitialization.");
    prm.declare_entry(
      "diffusivity power",
      "1.",
      Patterns::Double(),
      "Power value applied to the mesh-size in the mesh-dependant diffusion "
      "coefficient of the algebraic interface reinitialization.");
    prm.declare_entry("steady-state criterion",
                      "1e-2",
                      Patterns::Double(),
                      "Tolerance for the pseudo-time-stepping scheme.");
    prm.declare_entry("max steps number",
                      "5",
                      Patterns::Integer(),
                      "Maximum number of reinitialization steps.");
    prm.declare_entry(
      "reinitialization CFL",
      "0.5",
      Patterns::Double(),
      "CFL value for pseudo-time-step calculation purposes in the algebraic "
      "interface reinitialization.");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_AlgebraicInterfaceReinitialization::parse_parameters(
  dealii::ParameterHandler &prm)
{
  prm.enter_subsection("algebraic interface reinitialization");
  {
    this->output_reinitialization_steps =
      prm.get_bool("output reinitialization steps");
    this->diffusivity_multiplier = prm.get_double("diffusivity multiplier");
    this->diffusivity_power      = prm.get_double("diffusivity power");
    this->reinitialization_cfl   = prm.get_double("reinitialization CFL");
    this->steady_state_criterion = prm.get_double("steady-state criterion");
    this->max_steps_number       = prm.get_integer("max steps number");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_GeometricInterfaceReinitialization::declare_parameters(
  dealii::ParameterHandler &prm)
{
  prm.enter_subsection("geometric interface reinitialization");
  {
    prm.declare_entry(
      "output signed distance",
      "false",
      Patterns::Bool(),
      "Enables pvtu format outputs of the geometric interface reinitialization "
      "steps <true|false>");
    prm.declare_entry("max reinitialization distance",
                      "1.",
                      Patterns::Double(),
                      "Maximum reinitialization distance value");
    prm.declare_entry(
      "transformation type",
      "tanh",
      Patterns::Selection("tanh|piecewise polynomial"),
      "Transformation function used to get the phase indicator from the signed "
      "distance");
    prm.declare_entry("tanh thickness",
                      "1.",
                      Patterns::Double(),
                      "Interface thickness for the tanh transformation");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_GeometricInterfaceReinitialization::parse_parameters(
  dealii::ParameterHandler &prm)
{
  prm.enter_subsection("geometric interface reinitialization");
  {
    this->output_signed_distance = prm.get_bool("output signed distance");
    this->max_reinitialization_distance =
      prm.get_double("max reinitialization distance");
    const std::string t = prm.get("transformation type");
    if (t == "tanh")
      this->transformation_type =
        Parameters::RedistanciationTransformationType::tanh;
    else if (t == "piecewise polynomial")
      {
        this->transformation_type =
          Parameters::RedistanciationTransformationType::piecewise_polynomial;
      }
    else
      throw(std::runtime_error(
        "Invalid transformation type for the geometric interface reinitialization method!"));
    this->tanh_thickness = prm.get_double("tanh thickness");
  }
  prm.leave_subsection();
}

void
Parameters::CahnHilliard_PhaseFilter::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("phase filtration");
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection("none|clip|tanh"),
      "CahnHilliard phase filtration type, "
      "if <none> is selected, the phase won't be filtered"
      "if <clip> is selected, the phase order values above 1 (respectively below -1) will be brought back to 1 (respectively -1)"
      "if <tanh> is selected, the filtered phase will be a result of the "
      "following function: \\alpha_f = \\tanh(\\beta\\alpha); "
      "where beta is a parameter influencing the interface thickness that "
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
Parameters::CahnHilliard_PhaseFilter::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("phase filtration");
  {
    // filter type
    const std::string t = prm.get("type");
    if (t == "none")
      {
        type = Parameters::FilterType::none;
      }
    else if (t == "clip")
      {
        type = Parameters::FilterType::clip;
      }
    else if (t == "tanh")
      {
        type = Parameters::FilterType::tanh;
      }
    else
      throw(std::logic_error(
        "Error, invalid filter type. Choices are 'none', 'clip' or 'tanh'"));

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
Parameters::CahnHilliard::declare_parameters(ParameterHandler &prm) const
{
  prm.enter_subsection("cahn hilliard");
  {
    cahn_hilliard_phase_filter.declare_parameters(prm);

    prm.declare_entry(
      "potential smoothing coefficient",
      "1",
      Patterns::Double(),
      "Smoothing coefficient for the chemical potential in the Cahn-Hilliard equations.");

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

      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "Display the value of epsilon for each time iteration if set to verbose");
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}

void
Parameters::CahnHilliard::parse_parameters(ParameterHandler     &prm,
                                           const Dimensionality &dimensions)
{
  prm.enter_subsection("cahn hilliard");
  {
    cahn_hilliard_phase_filter.parse_parameters(prm);

    CahnHilliard::potential_smoothing_coefficient =
      prm.get_double("potential smoothing coefficient");

    prm.enter_subsection("epsilon");
    {
      const std::string op_epsilon = prm.get("method");
      if (op_epsilon == "automatic")
        {
          CahnHilliard::epsilon_set_method =
            Parameters::EpsilonSetMethod::automatic;
        }
      else if (op_epsilon == "manual")
        {
          CahnHilliard::epsilon_set_method =
            Parameters::EpsilonSetMethod::manual;
        }
      else
        throw(std::runtime_error("Invalid epsilon setting strategy. "
                                 "Options are 'automatic' or 'manual'."));

      const std::string op_epsilon_verbosity = prm.get("verbosity");
      if (op_epsilon_verbosity == "quiet")
        {
          CahnHilliard::epsilon_verbosity = Parameters::EpsilonVerbosity::quiet;
        }
      else if (op_epsilon_verbosity == "verbose")
        {
          CahnHilliard::epsilon_verbosity =
            Parameters::EpsilonVerbosity::verbose;
        }
      else
        AssertThrow(false,
                    ExcMessage("Invalid epsilon verbosity. "
                               "Options are 'quiet' or 'verbose'."));

      epsilon = prm.get_double("value");
      epsilon *= dimensions.cahn_hilliard_epsilon_scaling;
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();
}
