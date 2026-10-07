// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_multiphysics.h>
#include <core/utilities.h>

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
  ReinitializationMethodFrequencyError,
  int,
  << "Reinitialization method frequency : " << arg1
  << " is equal or smaller than 0." << std::endl
  << "Interface reinitialization method requires an frequency larger than 0.");

namespace Parameters
{
  namespace
  {
    /// Deprecated strings of the "transformation type" parameter
    const DeprecatedEnumNames<Parameters::RedistanciationTransformationType>
      deprecated_transformation_type_names = {
        {"piecewise polynomial",
         Parameters::RedistanciationTransformationType::piecewise_polynomial}};

    /// Deprecated strings of the interface reinitialization method "type"
    /// parameter
    const DeprecatedEnumNames<Parameters::ReinitializationMethodType>
      deprecated_reinitialization_method_type_names = {
        {"projection-based interface sharpening",
         Parameters::ReinitializationMethodType::projection_based_sharpening},
        {"pde-based interface reinitialization",
         Parameters::ReinitializationMethodType::pde_based},
        {"geometric interface reinitialization",
         Parameters::ReinitializationMethodType::geometric}};

    /// Deprecated strings of the "electromagnetic scaling type" parameter
    const DeprecatedEnumNames<Parameters::ElectromagneticScalingType>
      deprecated_electromagnetic_scaling_type_names = {
        {"electric field",
         Parameters::ElectromagneticScalingType::electric_field},
        {"magnetic field",
         Parameters::ElectromagneticScalingType::magnetic_field}};
  } // namespace
} // namespace Parameters

template <int dim>
void
Parameters::Multiphysics<dim>::declare_parameters(ParameterHandler &prm) const
{
  prm.enter_subsection("multiphysics");
  {
    prm.declare_entry("fluid dynamics",
                      Patterns::Tools::Convert<bool>::to_string(fluid_dynamics),
                      Patterns::Bool(),
                      "Fluid flow calculation <true|false>");

    prm.declare_entry("heat transfer",
                      Patterns::Tools::Convert<bool>::to_string(heat_transfer),
                      Patterns::Bool(),
                      "Thermic calculation <true|false>");

    prm.declare_entry("tracer",
                      Patterns::Tools::Convert<bool>::to_string(tracer),
                      Patterns::Bool(),
                      "Passive tracer calculation <true|false>");

    prm.declare_entry("cls",
                      Patterns::Tools::Convert<bool>::to_string(CLS),
                      Patterns::Bool(),
                      "CLS calculation <true|false>");

    prm.declare_entry("cahn hilliard",
                      Patterns::Tools::Convert<bool>::to_string(cahn_hilliard),
                      Patterns::Bool(),
                      "Cahn-Hilliard calculation <true|false>");

    prm.declare_entry(
      "electromagnetics",
      Patterns::Tools::Convert<bool>::to_string(electromagnetics),
      Patterns::Bool(),
      "Time harmonic electromagnetics calculation <true|false>");

    // subparameters for heat_transfer
    prm.declare_entry("viscous dissipation",
                      Patterns::Tools::Convert<bool>::to_string(
                        viscous_dissipation),
                      Patterns::Bool(),
                      "Viscous dissipation in heat equation <true|false>");

    prm.declare_entry("thermal buoyancy force",
                      Patterns::Tools::Convert<bool>::to_string(
                        thermal_buoyancy_force),
                      Patterns::Bool(),
                      "Thermal buoyancy force calculation <true|false>");
    prm.declare_entry("microwave heating",
                      Patterns::Tools::Convert<bool>::to_string(
                        microwave_heating),
                      Patterns::Bool(),
                      "Microwave heating calculation <true|false>");
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
    microwave_heating      = prm.get_bool("microwave heating");
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
    reinitialization_method.declare_parameters(prm);
    surface_tension_force.declare_parameters(prm);
    phase_filter.declare_parameters(prm);

    prm.declare_entry("viscous dissipative fluid",
                      enum_to_string(viscous_dissipative_fluid),
                      Patterns::Selection(enum_to_selection<FluidIndicator>(
                        deprecated_fluid_indicator_names())),
                      "Fluid to which the viscous dissipation is applied "
                      "in the heat equation <fluid0|fluid1|both>");

    prm.declare_entry(
      "diffusivity",
      Patterns::Tools::Convert<double>::to_string(diffusivity),
      Patterns::Double(),
      "Diffusivity (diffusion coefficient in L^2/s) in the phase indicator transport equation. "
      "Default value is 0 to have pure advection.");

    prm.declare_entry(
      "compressible",
      Patterns::Tools::Convert<bool>::to_string(compressible),
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
    reinitialization_method.parse_parameters(prm);
    surface_tension_force.parse_parameters(prm);
    phase_filter.parse_parameters(prm);

    // Viscous dissipative fluid
    viscous_dissipative_fluid =
      string_to_enum<FluidIndicator>(prm.get("viscous dissipative fluid"),
                                     deprecated_fluid_indicator_names(),
                                     "viscous dissipative fluid");

    diffusivity = prm.get_double("diffusivity");

    compressible = prm.get_bool("compressible");
  }
  prm.leave_subsection();
}

void
Parameters::CLS_ReinitializationMethod::declare_parameters(
  ParameterHandler &prm) const
{
  prm.enter_subsection("interface reinitialization method");
  {
    prm.declare_entry(
      "type",
      enum_to_string(reinitialization_method_type),
      Patterns::Selection(
        enum_to_selection<Parameters::ReinitializationMethodType>(
          deprecated_reinitialization_method_type_names)),
      "CLS interface reinitialization method. "
      "Choices are <none|projection_based_sharpening|pde_based|geometric>.");

    prm.declare_entry(
      "frequency",
      Patterns::Tools::Convert<int>::to_string(frequency),
      Patterns::Integer(0),
      "Reinitialization frequency (number of time steps) at which the "
      "interface reinitialization process will be applied to the CLS "
      "phase indicator field.");
    prm.declare_entry(
      "verbosity",
      enum_to_string(verbosity),
      Patterns::Selection(
        enum_to_selection<Verbosity>(deprecated_verbosity_names())),
      "States whether the output from the interface reinitialization method "
      "should be printed. "
      "Choices are <quiet|verbose|extra_verbose>.");

    sharpening.declare_parameters(prm);
    pde_based_interface_reinitialization.declare_parameters(prm);
    geometric_interface_reinitialization.declare_parameters(prm);
  }
  prm.leave_subsection();
}

void
Parameters::CLS_ReinitializationMethod::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("interface reinitialization method");
  {
    this->reinitialization_method_type =
      string_to_enum<Parameters::ReinitializationMethodType>(
        prm.get("type"), deprecated_reinitialization_method_type_names, "type");
    if (this->reinitialization_method_type ==
        Parameters::ReinitializationMethodType::projection_based_sharpening)
      sharpening.enable = true;
    else if (this->reinitialization_method_type ==
             Parameters::ReinitializationMethodType::pde_based)
      pde_based_interface_reinitialization.enable = true;
    else if (this->reinitialization_method_type ==
             Parameters::ReinitializationMethodType::geometric)
      geometric_interface_reinitialization.enable = true;

    this->frequency = prm.get_integer("frequency");
    Assert(this->frequency > 0,
           ReinitializationMethodFrequencyError(this->frequency));

    this->verbosity = string_to_enum<Verbosity>(prm.get("verbosity"),
                                                deprecated_verbosity_names(),
                                                "verbosity");

    this->sharpening.parse_parameters(prm);
    this->pde_based_interface_reinitialization.parse_parameters(prm);
    this->geometric_interface_reinitialization.parse_parameters(prm);
  }
  prm.leave_subsection();
}

void
Parameters::CLS_InterfaceSharpening::declare_parameters(ParameterHandler &prm)
{
  const CLS_InterfaceSharpening defaults;
  prm.enter_subsection("projection-based interface sharpening");
  {
    prm.declare_entry(
      "type",
      enum_to_string(defaults.type),
      Patterns::Selection(enum_to_selection<Parameters::SharpeningType>()),
      "CLS interface sharpening type, "
      "if constant the sharpening threshold is the same throughout the simulation, "
      "if adaptive the sharpening threshold is determined by binary search, "
      "to ensure mass conservation of the monitored phase");

    // Parameters for constant sharpening
    prm.declare_entry(
      "threshold",
      Patterns::Tools::Convert<double>::to_string(defaults.threshold),
      Patterns::Double(),
      "Interface sharpening threshold that represents the phase indicator at which "
      "the interphase is considered located");

    // Parameters for adaptive sharpening
    prm.declare_entry(
      "threshold max deviation",
      Patterns::Tools::Convert<double>::to_string(
        defaults.threshold_max_deviation),
      Patterns::Double(),
      "Maximum deviation (from the base value of 0.5) considered in the search "
      "algorithm to ensure mass conservation. "
      "A threshold max deviation of 0.20 results in a search interval from 0.30 to 0.70");

    prm.declare_entry(
      "max iterations",
      Patterns::Tools::Convert<int>::to_string(defaults.max_iterations),
      Patterns::Integer(),
      "Maximum number of iteration in the bissection algorithm that ensures mass conservation");

    prm.declare_entry(
      "monitoring",
      Patterns::Tools::Convert<bool>::to_string(defaults.monitoring),
      Patterns::Bool(),
      "Enable conservation monitoring in multiphase fluid simulations <true|false>");

    prm.declare_entry(
      "tolerance",
      Patterns::Tools::Convert<double>::to_string(defaults.tolerance),
      Patterns::Double(),
      "Tolerance on the mass conservation of the monitored fluid, used with adaptive sharpening");

    prm.declare_entry(
      "monitored fluid",
      enum_to_string(defaults.monitored_fluid),
      Patterns::Selection(enum_to_selection<FluidIndicator>(
        deprecated_fluid_indicator_names(), single_fluid_indicators())),
      "Fluid for which conservation is monitored <fluid0|fluid1>, used with adaptive sharpening.");

    // This parameter must be larger than 1 for interface sharpening. Choosing
    // values less than 1 leads to interface smoothing instead of sharpening.
    prm.declare_entry(
      "interface sharpness",
      Patterns::Tools::Convert<double>::to_string(defaults.interface_sharpness),
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
    type = string_to_enum<Parameters::SharpeningType>(prm.get("type"));

    // Parameters for constant sharpening
    threshold = prm.get_double("threshold");

    // Parameters for adaptive sharpening
    threshold_max_deviation = prm.get_double("threshold max deviation");
    max_iterations          = prm.get_integer("max iterations");
    monitoring              = prm.get_bool("monitoring");
    tolerance               = prm.get_double("tolerance");

    // Monitored fluid
    monitored_fluid =
      string_to_enum<FluidIndicator>(prm.get("monitored fluid"),
                                     deprecated_fluid_indicator_names(),
                                     "monitored fluid");

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
  const CLS_SurfaceTensionForce defaults;
  prm.enter_subsection("surface tension force");
  {
    prm.declare_entry("enable",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.enable),
                      Patterns::Bool(),
                      "Enable surface tension force calculation <true|false>");

    prm.declare_entry("output auxiliary fields",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.output_cls_auxiliary_fields),
                      Patterns::Bool(),
                      "Output the phase indicator gradient and curvature");

    prm.declare_entry(
      "phase indicator gradient diffusion factor",
      Patterns::Tools::Convert<double>::to_string(
        defaults.phase_indicator_gradient_diffusion_factor),
      Patterns::Double(),
      "Factor applied to the filter for phase indicator gradient calculations to damp high-frequency errors");

    prm.declare_entry(
      "curvature diffusion factor",
      Patterns::Tools::Convert<double>::to_string(
        defaults.curvature_diffusion_factor),
      Patterns::Double(),
      "Factor applied to the filter for curvature calculations to damp high-frequency errors");

    prm.declare_entry(
      "verbosity",
      enum_to_string(defaults.verbosity),
      Patterns::Selection(enum_to_selection<Verbosity>(
        deprecated_verbosity_names(), quiet_or_verbose())),
      "State whether the output from the surface tension force calculations should be printed "
      "Choices are <quiet|verbose>.");


    prm.declare_entry("enable marangoni effect",
                      Patterns::Tools::Convert<bool>::to_string(
                        defaults.enable_marangoni_effect),
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
    phase_indicator_gradient_diffusion_factor =
      prm.get_double("phase indicator gradient diffusion factor");
    curvature_diffusion_factor = prm.get_double("curvature diffusion factor");

    output_cls_auxiliary_fields = prm.get_bool("output auxiliary fields");

    verbosity = string_to_enum<Verbosity>(prm.get("verbosity"),
                                          deprecated_verbosity_names(),
                                          "verbosity");

    enable_marangoni_effect = prm.get_bool("enable marangoni effect");
  }
  prm.leave_subsection();
}

void
Parameters::CLS_PhaseFilter::declare_parameters(ParameterHandler &prm)
{
  const CLS_PhaseFilter defaults;
  prm.enter_subsection("phase filtration");
  {
    prm.declare_entry(
      "type",
      enum_to_string(defaults.type),
      Patterns::Selection(enum_to_selection<Parameters::FilterType>(
        {}, {Parameters::FilterType::none, Parameters::FilterType::tanh})),
      "CLS phase indicator filtration type, "
      "if <none> is selected, the phase won't be filtered; "
      "if <tanh> is selected, the filtered phase will be a result of the "
      "following function: \\alpha_f = 0.5 \\tanh(\\beta(\\alpha-0.5)) + 0.5; "
      "where \\beta is a parameter influencing the interface thickness that "
      "must be defined");
    prm.declare_entry(
      "beta",
      Patterns::Tools::Convert<double>::to_string(defaults.beta),
      Patterns::Double(),
      "This parameter appears in the tanh filter function. It influence "
      "the thickness and the shape of the interface. For higher values of "
      "beta, a thinner and 'sharper/pixelated' interface will be seen.");
    prm.declare_entry("verbosity",
                      enum_to_string(defaults.verbosity),
                      Patterns::Selection(enum_to_selection<Verbosity>(
                        deprecated_verbosity_names(), quiet_or_verbose())),
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
    type = string_to_enum<Parameters::FilterType>(prm.get("type"));

    // beta
    beta = prm.get_double("beta");

    // Verbosity
    verbosity = string_to_enum<Verbosity>(prm.get("verbosity"),
                                          deprecated_verbosity_names(),
                                          "verbosity");
  }
  prm.leave_subsection();
}

void
Parameters::CLS_PDEBasedInterfaceReinitialization::declare_parameters(
  dealii::ParameterHandler &prm)
{
  const CLS_PDEBasedInterfaceReinitialization defaults;
  prm.enter_subsection("PDE-based interface reinitialization");
  {
    prm.declare_entry(
      "output reinitialization steps",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.output_reinitialization_steps),
      Patterns::Bool(),
      "Enables pvtu format outputs of the PDE-based interface reinitialization "
      "steps <true|false>");
    prm.declare_entry(
      "diffusivity multiplier",
      Patterns::Tools::Convert<double>::to_string(
        defaults.diffusivity_multiplier),
      Patterns::Double(),
      "Factor that multiplies the mesh-size in the mesh-dependant diffusion "
      "coefficient of the PDE-based interface reinitialization.");
    prm.declare_entry(
      "diffusivity power",
      Patterns::Tools::Convert<double>::to_string(defaults.diffusivity_power),
      Patterns::Double(),
      "Power value applied to the mesh-size in the mesh-dependant diffusion "
      "coefficient of the PDE-based interface reinitialization.");
    prm.declare_entry("steady-state criterion",
                      Patterns::Tools::Convert<double>::to_string(
                        defaults.steady_state_criterion),
                      Patterns::Double(),
                      "Tolerance for the artificial time-stepping scheme.");
    prm.declare_entry("max steps number",
                      Patterns::Tools::Convert<double>::to_string(
                        defaults.max_steps_number),
                      Patterns::Integer(),
                      "Maximum number of reinitialization steps.");
    prm.declare_entry(
      "artificial time-step factor",
      Patterns::Tools::Convert<double>::to_string(defaults.dtau_factor),
      Patterns::Double(),
      "Factor multiplying the artificial time step in the PDE-based "
      "interface reinitialization.");
  }
  prm.leave_subsection();
}

void
Parameters::CLS_PDEBasedInterfaceReinitialization::parse_parameters(
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
  const CLS_GeometricInterfaceReinitialization defaults;
  prm.enter_subsection("geometric interface reinitialization");
  {
    prm.declare_entry(
      "output signed distance",
      Patterns::Tools::Convert<bool>::to_string(
        defaults.output_signed_distance),
      Patterns::Bool(),
      "Enables pvtu format outputs of the geometric interface reinitialization "
      "steps <true|false>");
    prm.declare_entry("max reinitialization distance",
                      Patterns::Tools::Convert<double>::to_string(
                        defaults.max_reinitialization_distance),
                      Patterns::Double(),
                      "Maximum reinitialization distance value");
    prm.declare_entry(
      "transformation type",
      enum_to_string(defaults.transformation_type),
      Patterns::Selection(
        enum_to_selection<Parameters::RedistanciationTransformationType>(
          deprecated_transformation_type_names)),
      "Transformation function used to get the phase indicator from the signed "
      "distance");
    prm.declare_entry("tanh thickness",
                      Patterns::Tools::Convert<double>::to_string(
                        defaults.tanh_thickness),
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
    this->transformation_type =
      string_to_enum<Parameters::RedistanciationTransformationType>(
        prm.get("transformation type"),
        deprecated_transformation_type_names,
        "transformation type");
    this->tanh_thickness = prm.get_double("tanh thickness");
  }
  prm.leave_subsection();
}

void
Parameters::CahnHilliard_PhaseFilter::declare_parameters(ParameterHandler &prm)
{
  const CahnHilliard_PhaseFilter defaults;
  prm.enter_subsection("phase filtration");
  {
    prm.declare_entry(
      "type",
      enum_to_string(defaults.type),
      Patterns::Selection(enum_to_selection<Parameters::FilterType>()),
      "CahnHilliard phase filtration type, "
      "if <none> is selected, the phase won't be filtered; "
      "if <clip> is selected, the phase order values above 1 (respectively below -1) will be brought back to 1 (respectively -1); "
      "if <tanh> is selected, the filtered phase will be a result of the "
      "following function: \\alpha_f = \\tanh(\\beta\\alpha); "
      "where beta is a parameter influencing the interface thickness that "
      "must be defined");
    prm.declare_entry(
      "beta",
      Patterns::Tools::Convert<double>::to_string(defaults.beta),
      Patterns::Double(),
      "This parameter appears in the tanh filter function. It influence "
      "the thickness and the shape of the interface. For higher values of "
      "beta, a thinner and 'sharper/pixelated' interface will be seen.");
    prm.declare_entry("verbosity",
                      enum_to_string(defaults.verbosity),
                      Patterns::Selection(enum_to_selection<Verbosity>(
                        deprecated_verbosity_names(), quiet_or_verbose())),
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
    type = string_to_enum<Parameters::FilterType>(prm.get("type"));

    // beta
    beta = prm.get_double("beta");

    // Verbosity
    verbosity = string_to_enum<Verbosity>(prm.get("verbosity"),
                                          deprecated_verbosity_names(),
                                          "verbosity");
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
      Patterns::Tools::Convert<double>::to_string(
        potential_smoothing_coefficient),
      Patterns::Double(),
      "Smoothing coefficient for the chemical potential in the Cahn-Hilliard equations.");

    prm.enter_subsection("epsilon");
    {
      prm.declare_entry(
        "method",
        enum_to_string(epsilon_set_method),
        Patterns::Selection(enum_to_selection<Parameters::EpsilonSetMethod>()),
        "Epsilon is either set to two times the characteristic length (automatic) of the element or user defined on all the domain (manual)");

      prm.declare_entry(
        "value",
        Patterns::Tools::Convert<double>::to_string(epsilon),
        Patterns::Double(),
        "Parameter linked to the interface thickness. Should always be bigger than the characteristic size of the smallest element");

      prm.declare_entry(
        "verbosity",
        enum_to_string(epsilon_verbosity),
        Patterns::Selection(enum_to_selection<Parameters::EpsilonVerbosity>()),
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
      CahnHilliard::epsilon_set_method =
        string_to_enum<Parameters::EpsilonSetMethod>(prm.get("method"));

      CahnHilliard::epsilon_verbosity =
        string_to_enum<Parameters::EpsilonVerbosity>(prm.get("verbosity"));

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
    prm.enter_subsection("time coupling strategy");
    {
      prm.declare_entry(
        "type",
        enum_to_string(time_coupling_strategy),
        Patterns::Selection(
          enum_to_selection<Parameters::TimeHarmonicMaxwellCouplingStrategy>()),
        "The type of time coupling strategy to use.");

      prm.declare_entry(
        "coupling iteration",
        Patterns::Tools::Convert<unsigned int>::to_string(coupling_iteration),
        Patterns::Integer(1),
        "Coupling parameter for the time coupling strategy based on "
        "iteration, this parameter represents the number of time iterations "
        "between two consecutive resolutions of the electromagnetic fields.");

      prm.declare_entry(
        "coupling time",
        Patterns::Tools::Convert<double>::to_string(coupling_time),
        Patterns::Double(DBL_MIN),
        "Coupling parameter for the time coupling strategy based on time. "
        "This parameter represents the real time interval between two "
        "consecutive resolutions of the electromagnetic fields.");

      prm.declare_entry(
        "coupling threshold",
        Patterns::Tools::Convert<double>::to_string(coupling_threshold),
        Patterns::Double(0),
        "Coupling parameter for the time coupling strategy based on a "
        "threshold. This parameter represents the change in the "
        "electromagnetic properties of the medium (i.e., permittivity, "
        "permeability or conductivity) that triggers the recomputation of "
        "the electromagnetic fields.)");
    }
    prm.leave_subsection();


    prm.declare_entry(
      "electromagnetic frequency",
      Patterns::Tools::Convert<double>::to_string(electromagnetic_frequency),
      Patterns::Double(0),
      "Frequency of the time harmonic electromagnetic wave excitation (in Hz).");

    prm.declare_entry(
      "electromagnetic scaling type",
      enum_to_string(electromagnetic_scaling_type),
      Patterns::Selection(
        enum_to_selection<Parameters::ElectromagneticScalingType>(
          deprecated_electromagnetic_scaling_type_names)),
      "The type of electromagnetic scaling to apply to the solution of the time-harmonic Maxwell solver after solving the linear system. This is relevant when the user wants to recover the physical solution in dimensional units instead of the dimensionless solution used for better conditioning of the linear system.");

    prm.declare_entry(
      "electric field amplitude",
      Patterns::Tools::Convert<double>::to_string(electric_field_amplitude),
      Patterns::Double(0),
      "The amplitude of the electric field used for the normalization of the solution in [V/m].");

    prm.declare_entry(
      "magnetic field amplitude",
      Patterns::Tools::Convert<double>::to_string(magnetic_field_amplitude),
      Patterns::Double(0),
      "The amplitude of the magnetic field used for the normalization of the solution in [A/m].");

    prm.declare_entry("number of waveguide inlets",
                      Patterns::Tools::Convert<unsigned int>::to_string(
                        number_of_waveguide_inlets),
                      Patterns::Integer(0),
                      "Number of waveguide inlets in the simulation.");

    // Declare a fixed maximum number of waveguide inlets.
    // Only the ones specified by "number of waveguide inlets" will be
    // parsed. This is necessary because declare_parameters runs before the
    // file is read, so we can't know the actual number of inlets at
    // declaration time.
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
              enum_to_string(Parameters::WaveguideMode::TE),
              Patterns::Selection(
                enum_to_selection<Parameters::WaveguideMode>()),
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
    prm.enter_subsection("time coupling strategy");
    {
      time_coupling_strategy =
        string_to_enum<Parameters::TimeHarmonicMaxwellCouplingStrategy>(
          prm.get("type"));

      TimeHarmonicMaxwell::coupling_iteration =
        prm.get_integer("coupling iteration");
      TimeHarmonicMaxwell::coupling_time = prm.get_double("coupling time");
      TimeHarmonicMaxwell::coupling_threshold =
        prm.get_double("coupling threshold");
    }
    prm.leave_subsection();

    TimeHarmonicMaxwell::electromagnetic_frequency =
      prm.get_double("electromagnetic frequency") *
      dimensions.electromagnetic_frequency_scaling;

    TimeHarmonicMaxwell::electromagnetic_scaling_type =
      string_to_enum<Parameters::ElectromagneticScalingType>(
        prm.get("electromagnetic scaling type"),
        deprecated_electromagnetic_scaling_type_names,
        "electromagnetic scaling type");

    // The user always provides the electric and magnetic field amplitudes
    // in dimensional units (V/m and A/m respectively). The scaling factors
    // will be applied in the time_harmonic_maxwell class.
    TimeHarmonicMaxwell::electric_field_amplitude =
      prm.get_double("electric field amplitude");

    TimeHarmonicMaxwell::magnetic_field_amplitude =
      prm.get_double("magnetic field amplitude");

    TimeHarmonicMaxwell::number_of_waveguide_inlets =
      prm.get_integer("number of waveguide inlets");

    // Ensure that the number of waveguide inlets is smaller than the
    // maximum declared in declare_parameters.
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
            TimeHarmonicMaxwell::waveguide_mode[inlet] =
              string_to_enum<Parameters::WaveguideMode>(prm.get("mode type"));

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
          // quadrilateral. This is not useful in 2D as any 2 points are
          // always collinear.
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

              // The triple scalar product (determinant) scales as the
              // volume of the parallelepiped defined by the three vectors.
              // So we make the it dimensionless by a characteristic volume
              // (product of the norms of the three vectors). This way we
              // can set a tolerance that is independent of the actual size
              // of the waveguide.
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
