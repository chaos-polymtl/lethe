
#include <core/parameters_multiphysics.h>

#include <deal.II/base/parameter_handler.h>



void
Parameters::VOF_Monitoring::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("monitoring");
  {
    prm.declare_entry(
      "conservation monitoring",
      "false",
      Patterns::Bool(),
      "Conservation monitoring in free surface calculation <true|false>");

    prm.declare_entry(
      "fluid monitored",
      "1",
      Patterns::Integer(),
      "Index of the fluid which conservation is monitored <0|1>");
  }
  prm.leave_subsection();
}

void
Parameters::VOF_Monitoring::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("monitoring");
  {
    conservation_monitoring = prm.get_bool("conservation monitoring");
    id_fluid_monitored      = prm.get_integer("fluid monitored");
  }
  prm.leave_subsection();
}

void
Parameters::VOF::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("VOF");
  {
    prm.declare_entry("interface sharpening",
                      "false",
                      Patterns::Bool(),
                      "Interface sharpening <true|false>");

    prm.declare_entry("continuum surface force",
                      "false",
                      Patterns::Bool(),
                      "Continuum surface force calculation <true|false>");

    prm.declare_entry(
      "peeling wetting",
      "false",
      Patterns::Bool(),
      "Enable peeling/wetting in free surface calculation <true|false>");

    prm.declare_entry(
      "skip mass conservation in fluid 0",
      "false",
      Patterns::Bool(),
      "Enable skipping mass conservation in fluid 0 <true|false>");

    prm.declare_entry(
      "skip mass conservation in fluid 1",
      "false",
      Patterns::Bool(),
      "Enable skipping mass conservation in fluid 1 <true|false>");

    monitoring.declare_parameters(prm);
  }
  prm.leave_subsection();
}

void
Parameters::VOF::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("VOF");
  {
    interface_sharpening    = prm.get_bool("interface sharpening");
    continuum_surface_force = prm.get_bool("continuum surface force");
    peeling_wetting         = prm.get_bool("peeling wetting");
    skip_mass_conservation_fluid_0 =
      prm.get_bool("skip mass conservation in fluid 0");
    skip_mass_conservation_fluid_1 =
      prm.get_bool("skip mass conservation in fluid 1");

    monitoring.parse_parameters(prm);
  }
  prm.leave_subsection();
}

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
