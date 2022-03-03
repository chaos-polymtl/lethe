
#include <core/parameters_multiphysics.h>

#include <deal.II/base/parameter_handler.h>



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

    // subparameters for VOF
    prm.declare_entry("interface sharpening",
                      "false",
                      Patterns::Bool(),
                      "Interface sharpening <true|false>");

    prm.declare_entry("surface tension force",
                      "false",
                      Patterns::Bool(),
                      "Surface tension force calculation <true|false>");

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

    prm.declare_entry(
      "peeling wetting",
      "false",
      Patterns::Bool(),
      "Enable peeling/wetting in free surface calculation <true|false>");
  }
  prm.leave_subsection();
}

void
Parameters::Multiphysics::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("multiphysics");
  {
    fluid_dynamics        = prm.get_bool("fluid dynamics");
    heat_transfer         = prm.get_bool("heat transfer");
    tracer                = prm.get_bool("tracer");
    VOF                   = prm.get_bool("VOF");
    interface_sharpening  = prm.get_bool("interface sharpening");
    buoyancy_force        = prm.get_bool("buoyancy force");
    surface_tension_force = prm.get_bool("surface tension force");

    // subparameter for heat_transfer
    viscous_dissipation = prm.get_bool("viscous dissipation");
    buoyancy_force      = prm.get_bool("buoyancy force");

    // subparameters for VOF
    interface_sharpening    = prm.get_bool("interface sharpening");
    surface_tension_force = prm.get_bool("surface tension force");
    conservation_monitoring = prm.get_bool("conservation monitoring");
    id_fluid_monitored      = prm.get_integer("fluid monitored");
    peeling_wetting         = prm.get_bool("peeling wetting");
  }
  prm.leave_subsection();
}
