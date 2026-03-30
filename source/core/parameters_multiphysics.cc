// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
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
  << "Reinitialization method frequency : " << arg1
  << " is equal or smaller than 0." << std::endl
  << "Interface reinitialization method requires an frequency larger than 0.");

template <int dim>
void
Parameters::Multiphysics<dim>::declare_parameters(ParameterHandler &prm) const
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

    prm.declare_entry("cls",
                      "false",
                      Patterns::Bool(),
                      "CLS calculation <true|false>");
    prm.declare_alias("cls", "VOF", true);

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

    prm.declare_entry("thermal buoyancy force",
                      "false",
                      Patterns::Bool(),
                      "Thermal buoyancy force calculation <true|false>");
    prm.declare_alias("thermal buoyancy force",
                      "buoyancy force",
                      true); // temporary alias for backward compatibility
  }
  prm.leave_subsection();

  cls_parameters.declare_parameters(prm);
  cahn_hilliard_parameters.declare_parameters(prm);
  time_harmonic_maxwell_parameters.declare_parameters(prm);
}

template <int dim>
void
Parameters::Multiphysics<dim>::parse_parameters(
  ParameterHandler     &prm,
  const Dimensionality &dimensions)
{
  prm.enter_subsection("multiphysics");
  {
    fluid_dynamics   = prm.get_bool("fluid dynamics");
    heat_transfer    = prm.get_bool("heat transfer");
    tracer           = prm.get_bool("tracer");
    CLS              = prm.get_bool("cls");
    cahn_hilliard    = prm.get_bool("cahn hilliard");
    electromagnetics = prm.get_bool("electromagnetics");

    // subparameters for heat_transfer
    viscous_dissipation    = prm.get_bool("viscous dissipation");
    thermal_buoyancy_force = prm.get_bool("thermal buoyancy force");
  }
  prm.leave_subsection();
  cls_parameters.parse_parameters(prm);
  cahn_hilliard_parameters.parse_parameters(prm, dimensions);
  time_harmonic_maxwell_parameters.parse_parameters(prm, dimensions);
}

void
Parameters::CLS::declare_parameters(ParameterHandler &prm) const
{
  prm.enter_subsection("CLS");
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
      "Diffusivity (diffusion coefficient in L^2/s) in the phase indicator transport equation. "
      "Default value is 0 to have pure advection.");

    prm.declare_entry(
      "compressible",
      "false",
      Patterns::Bool(),
      "Enable phase compressibility in the CLS equation. This leads to the inclusion of the phase * div(u) term in the CLS equation. "
      "It should be set to false when the phases are incompressible");
  }
  prm.leave_subsection();
}

void
Parameters::CLS::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("CLS");
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
Parameters::CLS_RegularizationMethod::declare_parameters(
  ParameterHandler &prm) const
{
  prm.enter_subsection("interface reinitialization method");
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection(
        "none|projection-based interface sharpening|pde-based interface reinitialization|geometric interface reinitialization"),
      "CLS interface reinitialization method");

    prm.declare_entry(
      "frequency",
      "10",
      Patterns::Integer(),
      "Reinitialization frequency (number of time steps) at which the "
      "interface reinitialization process will be applied to the CLS "
      "phase indicator field.");
    prm.declare_entry(
      "verbosity",
      "quiet",
      Patterns::Selection("quiet|verbose|extra verbose"),
      "States whether the output from the interface reinitialization method "
      "should be printed."
      "Choices are <quiet|verbose|extra verbose>.");

    sharpening.declare_parameters(prm);
    algebraic_interface_reinitialization.declare_parameters(prm);
    geometric_interface_reinitialization.declare_parameters(prm);
  }
  prm.leave_subsection();
}

void
Parameters::CLS_RegularizationMethod::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("interface reinitialization method");
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
    else if (t == "pde-based interface reinitialization")
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
        std::runtime_error("Invalid interface reinitialization method type!"));

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
Parameters::CLS_InterfaceSharpening::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("projection-based interface sharpening");
  {
    prm.declare_entry(
      "type",
      "constant",
      Patterns::Selection("constant|adaptive"),
      "CLS interface sharpening type, "
      "if constant the sharpening threshold is the same throughout the simulation, "
      "if adaptive the sharpening threshold is determined by binary search, "
      "to ensure mass conservation of the monitored phase");

    // Parameters for constant sharpening
    prm.declare_entry(
      "threshold",
      "0.5",
      Patterns::Double(),
      "Interface sharpening threshold that represents the phase indicator at which "
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
Parameters::CLS_InterfaceSharpening::parse_parameters(ParameterHandler &prm)
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
Parameters::CLS_SurfaceTensionForce::declare_parameters(ParameterHandler &prm)
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
                      "Output the phase indicator gradient and curvature");

    prm.declare_entry(
      "phase indicator gradient diffusion factor",
      "4",
      Patterns::Double(),
      "Factor applied to the filter for phase indicator gradient calculations to damp high-frequency errors");
    prm.declare_alias("phase indicator gradient diffusion factor",
                      "phase fraction gradient diffusion factor",
                      true);

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
Parameters::CLS_SurfaceTensionForce::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("surface tension force");
  {
    enable = prm.get_bool("enable");
    phase_fraction_gradient_diffusion_factor =
      prm.get_double("phase indicator gradient diffusion factor");
    curvature_diffusion_factor = prm.get_double("curvature diffusion factor");

    output_cls_auxiliary_fields = prm.get_bool("output auxiliary fields");

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
Parameters::CLS_PhaseFilter::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("phase filtration");
  {
    prm.declare_entry(
      "type",
      "none",
      Patterns::Selection("none|tanh"),
      "CLS phase indicator filtration type, "
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
Parameters::CLS_PhaseFilter::parse_parameters(ParameterHandler &prm)
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
Parameters::CLS_AlgebraicInterfaceReinitialization::declare_parameters(
  dealii::ParameterHandler &prm)
{
  prm.enter_subsection("PDE-based interface reinitialization");
  {
    prm.declare_entry(
      "output reinitialization steps",
      "false",
      Patterns::Bool(),
      "Enables pvtu format outputs of the PDE-based interface reinitialization "
      "steps <true|false>");
    prm.declare_entry(
      "diffusivity multiplier",
      "1.",
      Patterns::Double(),
      "Factor that multiplies the mesh-size in the mesh-dependant diffusion "
      "coefficient of the PDE-based interface reinitialization.");
    prm.declare_entry(
      "diffusivity power",
      "1.",
      Patterns::Double(),
      "Power value applied to the mesh-size in the mesh-dependant diffusion "
      "coefficient of the PDE-based interface reinitialization.");
    prm.declare_entry("steady-state criterion",
                      "1e-4",
                      Patterns::Double(),
                      "Tolerance for the artificial time-stepping scheme.");
    prm.declare_entry("max steps number",
                      "10000",
                      Patterns::Integer(),
                      "Maximum number of reinitialization steps.");
    prm.declare_entry(
      "artificial time-step factor",
      "0.25",
      Patterns::Double(),
      "Factor multiplying the artificial time step in the PDE-based "
      "interface reinitialization.");
    prm.declare_alias("artificial time-step factor",
                      "reinitialization CFL",
                      true);
  }
  prm.leave_subsection();
}

void
Parameters::CLS_AlgebraicInterfaceReinitialization::parse_parameters(
  dealii::ParameterHandler &prm)
{
  prm.enter_subsection("PDE-based interface reinitialization");
  {
    this->output_reinitialization_steps =
      prm.get_bool("output reinitialization steps");
    this->diffusivity_multiplier = prm.get_double("diffusivity multiplier");
    this->diffusivity_power      = prm.get_double("diffusivity power");
    this->dtau_factor = prm.get_double("artificial time-step factor");
    this->steady_state_criterion = prm.get_double("steady-state criterion");
    this->max_steps_number       = prm.get_integer("max steps number");
  }
  prm.leave_subsection();
}

void
Parameters::CLS_GeometricInterfaceReinitialization::declare_parameters(
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
Parameters::CLS_GeometricInterfaceReinitialization::parse_parameters(
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

template <int dim>
void
Parameters::TimeHarmonicMaxwell<dim>::declare_parameters(
  ParameterHandler &prm) const
{
  prm.enter_subsection("time harmonic maxwell");
  {
    prm.declare_entry(
      "electromagnetic frequency",
      "1",
      Patterns::Double(0),
      "Frequency of the time harmonic electromagnetic wave excitation (in Hz).");

    prm.declare_entry(
      "electromagnetic scaling type",
      "none",
      Patterns::Selection("none|electric field|magnetic field|power"),
      "The type of electromagnetic scaling to apply to the solution of the time-harmonic Maxwell solver after solving the linear system. This is relevant when the user wants to recover the physical solution in dimensional units instead of the dimensionless solution used for better conditioning of the linear system.");

    prm.declare_entry(
      "electric field amplitude",
      "1",
      Patterns::Double(0),
      "The amplitude of the electric field used for the normalization of the solution in [V/m].");

    prm.declare_entry(
      "magnetic field amplitude",
      "1",
      Patterns::Double(0),
      "The amplitude of the magnetic field used for the normalization of the solution in [A/m].");

    prm.declare_entry("number of waveguide inlets",
                      "0",
                      Patterns::Integer(0),
                      "Number of waveguide inlets in the simulation.");

    // Declare a fixed maximum number of waveguide inlets.
    // Only the ones specified by "number of waveguide inlets" will be parsed.
    // This is necessary because declare_parameters runs before the file is
    // read, so we can't know the actual number of inlets at declaration time.
    constexpr unsigned int max_waveguide_inlets = 10;

    for (unsigned int inlet = 0; inlet < max_waveguide_inlets; ++inlet)
      {
        prm.enter_subsection("waveguide inlet " + std::to_string(inlet));
        {
          prm.declare_entry(
            "port boundary id",
            "0",
            Patterns::Integer(0),
            "The boundary id where the waveguide inlet is applied.");

          prm.declare_entry(
            "waveguide power",
            "1",
            Patterns::Double(0),
            "The power of the waveguide mode excitation in Watts. This is used to compute the amplitude of the electromagnetic wave at the inlet in dimensional units.");

          prm.enter_subsection("waveguide mode");
          {
            prm.declare_entry(
              "mode type",
              "TE",
              Patterns::Selection("TE|TM"),
              "The waveguide mode excitation for a rectangular waveguide can be either Transverse Electric (TE) or Transverse Magnetic (TM).");

            prm.declare_entry(
              "mode order m",
              "1",
              Patterns::Integer(0),
              "The mode order m in the first transverse direction of the rectangular waveguide.");

            prm.declare_entry(
              "mode order n",
              "0",
              Patterns::Integer(0),
              "The mode order n in the second transverse direction of the rectangular waveguide.");
          }
          prm.leave_subsection();

          const unsigned int num_corners = Utilities::fixed_power<2>(dim - 1);
          for (unsigned int corner = 0; corner < num_corners; ++corner)
            {
              if constexpr (dim == 2)
                {
                  prm.declare_entry("corner " + std::to_string(corner),
                                    "0, 0",
                                    Patterns::List(Patterns::Double(), 2, 2),
                                    "Coordinates of corner " +
                                      std::to_string(corner));
                }
              else if constexpr (dim == 3)
                {
                  prm.declare_entry("corner " + std::to_string(corner),
                                    "0, 0, 0",
                                    Patterns::List(Patterns::Double(), 3, 3),
                                    "Coordinates of corner " +
                                      std::to_string(corner));
                }
            }
        }
        prm.leave_subsection();
      }
  }
  prm.leave_subsection();
}

template <int dim>
void
Parameters::TimeHarmonicMaxwell<dim>::parse_parameters(
  ParameterHandler     &prm,
  const Dimensionality &dimensions)
{
  prm.enter_subsection("time harmonic maxwell");
  {
    TimeHarmonicMaxwell::electromagnetic_frequency =
      prm.get_double("electromagnetic frequency") *
      dimensions.electromagnetic_frequency_scaling;

    const std::string op_scaling_type = prm.get("electromagnetic scaling type");
    if (op_scaling_type == "none")
      TimeHarmonicMaxwell::electromagnetic_scaling_type =
        Parameters::ElectromagneticScalingType::none;
    else if (op_scaling_type == "electric field")
      TimeHarmonicMaxwell::electromagnetic_scaling_type =
        Parameters::ElectromagneticScalingType::electric_field;
    else if (op_scaling_type == "magnetic field")
      TimeHarmonicMaxwell::electromagnetic_scaling_type =
        Parameters::ElectromagneticScalingType::magnetic_field;
    else if (op_scaling_type == "power")
      TimeHarmonicMaxwell::electromagnetic_scaling_type =
        Parameters::ElectromagneticScalingType::power;
    else
      throw(std::runtime_error(
        "Invalid electromagnetic scaling type. "
        "Options are <none|electric field|magnetic field|power>."));

    // By default, the electric field dimensionality is in V/m, but if the user
    // changed the dimensionality of the problem, we need to change the
    // dimensionality of the electric field accordingly to ensure that the
    // correct physical solution is obtained in dimensional units. This needs to
    // be apply to the power later on as well if the user chooses power-based
    // scaling which we cannot make dimensionless yet to recover the amplitude
    // scaling.
    TimeHarmonicMaxwell::electric_field_dimensionality =
      dimensions.electric_amplitude_scaling;
    // The same applies for the magnetic field, which is in A/m by default.
    TimeHarmonicMaxwell::magnetic_field_dimensionality =
      dimensions.magnetic_amplitude_scaling;

    // If the choose other scaling, we can already apply the dimensionality to
    // the amplitude provided.
    TimeHarmonicMaxwell::electric_field_amplitude =
      prm.get_double("electric field amplitude") *
      TimeHarmonicMaxwell::electric_field_dimensionality;

    TimeHarmonicMaxwell::magnetic_field_amplitude =
      prm.get_double("magnetic field amplitude") *
      TimeHarmonicMaxwell::magnetic_field_dimensionality;

    TimeHarmonicMaxwell::number_of_waveguide_inlets =
      prm.get_integer("number of waveguide inlets");

    // Ensure that the number of waveguide inlets is smaller than the maximum
    // declared in declare_parameters.
    AssertThrow(
      TimeHarmonicMaxwell::number_of_waveguide_inlets <= 10,
      ExcMessage(
        "The number of waveguide inlets specified exceeds "
        "the maximum allowed. Please increase the value "
        "of max_waveguide_inlets in the code declare_parameters section of the TimeHarmonicMaxwell class."));

    // Resize vectors to hold the correct number of inlets
    TimeHarmonicMaxwell::waveguide_mode.resize(number_of_waveguide_inlets);
    TimeHarmonicMaxwell::mode_order_m.resize(number_of_waveguide_inlets);
    TimeHarmonicMaxwell::mode_order_n.resize(number_of_waveguide_inlets);
    TimeHarmonicMaxwell::waveguide_boundary_ids.resize(
      number_of_waveguide_inlets);
    TimeHarmonicMaxwell::waveguide_power.resize(number_of_waveguide_inlets);

    for (unsigned int inlet = 0; inlet < number_of_waveguide_inlets; ++inlet)
      {
        prm.enter_subsection("waveguide inlet " + std::to_string(inlet));
        {
          TimeHarmonicMaxwell::waveguide_boundary_ids[inlet] =
            prm.get_integer("port boundary id");

          TimeHarmonicMaxwell::waveguide_power[inlet] =
            prm.get_double("waveguide power");

          // Check that the waveguide power is not zero which would imply no
          // excitation at the inlet.
          AssertThrow(
            TimeHarmonicMaxwell::waveguide_power[inlet] > 0,
            ExcMessage(
              "The waveguide power for the inlet " + std::to_string(inlet) +
              " is zero. Please check the waveguide power parameter in the input prm file for this inlet. If you really want to have no excitation but an open boundary condition at this inlet, you should change the boundary condition type to the impedance boundary condition with zero excitation and the desired admittance."));

          prm.enter_subsection("waveguide mode");
          {
            const std::string op_mode_type = prm.get("mode type");
            if (op_mode_type == "TE")
              {
                TimeHarmonicMaxwell::waveguide_mode[inlet] =
                  Parameters::WaveguideMode::TE;
              }
            else if (op_mode_type == "TM")
              {
                TimeHarmonicMaxwell::waveguide_mode[inlet] =
                  Parameters::WaveguideMode::TM;
              }
            else
              throw(std::runtime_error("Invalid waveguide mode type. "
                                       "Options are 'TE' or 'TM'."));

            TimeHarmonicMaxwell::mode_order_m[inlet] =
              prm.get_integer("mode order m");

            TimeHarmonicMaxwell::mode_order_n[inlet] =
              prm.get_integer("mode order n");

            // Check that the choice of mode orders is valid
            if (TimeHarmonicMaxwell::mode_order_m[inlet] == 0 &&
                TimeHarmonicMaxwell::mode_order_n[inlet] == 0)
              {
                throw(std::runtime_error(
                  "Invalid waveguide mode orders m and n. "
                  "At least one of the two mode orders must be non-zero."));
              }
          }
          prm.leave_subsection();

          const unsigned int num_corners = Utilities::fixed_power<2>(dim - 1);
          std::array<Tensor<1, dim>, num_corners> tmp_corners;

          for (unsigned int corner = 0; corner < num_corners; ++corner)
            {
              const std::string corner_name =
                "corner " + std::to_string(corner);
              const Tensor<1, dim> corner_value =
                value_string_to_tensor<dim>(prm.get(corner_name));

              tmp_corners[corner] = corner_value;
            }

          // In 3D, verify that the provided corners define a coplanar
          // quadrilateral. This is not useful in 2D as any 2 points are always
          // collinear.
          if constexpr (dim == 3)
            {
              const Tensor<1, dim> vec1 =
                tmp_corners[1] -
                tmp_corners[0]; // Vector from corner 1 to corner 2
              const Tensor<1, dim> vec2 =
                tmp_corners[2] -
                tmp_corners[0]; // Vector from corner 1 to corner 3
              const Tensor<1, dim> vec3 =
                tmp_corners[3] -
                tmp_corners[0]; // Vector from corner 1 to corner 4

              // The triple scalar product (determinant) scales as the volume of
              // the parallelepiped defined by the three vectors. So we make the
              // it dimensionless by a characteristic volume (product of the
              // norms of the three vectors). This way we can set a tolerance
              // that is independent of the actual size of the waveguide.
              const double determinant =
                scalar_product(vec3, cross_product_3d(vec1, vec2)) /
                (vec1.norm() * vec2.norm() * vec3.norm());

              const double tolerance = 1e-10;

              AssertThrow(std::abs(determinant) < tolerance,
                          ExcMessage(
                            "The provided corners for waveguide inlet " +
                            std::to_string(inlet) +
                            " do not define a coplanar quadrilateral. "
                            "Please check the corner coordinates."));
            }

          TimeHarmonicMaxwell::waveguide_corners.push_back(tmp_corners);
        }
        prm.leave_subsection();
      }

    // Check that there is at least one waveguide inlet if the user use the
    // `power` electromagnetic scaling type.
    AssertThrow(
      !((TimeHarmonicMaxwell::electromagnetic_scaling_type ==
         Parameters::ElectromagneticScalingType::power) &&
        TimeHarmonicMaxwell::number_of_waveguide_inlets == 0),
      ExcMessage(
        "The power-based electromagnetic scaling type requires at least one waveguide inlet to be defined. Please check the number of waveguide inlets specified in the input prm file."));
  }
  prm.leave_subsection();
}

template struct Parameters::TimeHarmonicMaxwell<2>;
template struct Parameters::TimeHarmonicMaxwell<3>;
template struct Parameters::Multiphysics<2>;
template struct Parameters::Multiphysics<3>;
