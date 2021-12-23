
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

    prm.declare_entry("interface sharpening",
                      "false",
                      Patterns::Bool(),
                      "Interface sharpening <true|false>");

    prm.declare_entry("buoyancy force",
                      "false",
                      Patterns::Bool(),
                      "Buoyant force calculation <true|false>");

    // subparameter for heat_transfer
    prm.declare_entry("viscous dissipation",
                      "false",
                      Patterns::Bool(),
                      "Viscous dissipation in heat equation <true|false>");

    // subparameter for free_surface
    prm.declare_entry(
      "conservation monitoring",
      "false",
      Patterns::Bool(),
      "Conservation monitoring in free surface calculation <true|false>");

    prm.declare_entry(
      "fluid index",
      "0",
      Patterns::Integer(),
      "Index of the fluid which conservation is monitored <0|1>");
  }
  prm.leave_subsection();
}

void
Parameters::Multiphysics::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("multiphysics");
  {
    fluid_dynamics       = prm.get_bool("fluid dynamics");
    heat_transfer        = prm.get_bool("heat transfer");
    tracer               = prm.get_bool("tracer");
    VOF                  = prm.get_bool("VOF");
    interface_sharpening = prm.get_bool("interface sharpening");
    buoyancy_force       = prm.get_bool("buoyancy force");

    // subparameter for heat_transfer
    viscous_dissipation = prm.get_bool("viscous dissipation");

    // subparameter for free_surface
    conservation_monitoring = prm.get_bool("conservation monitoring");
    fluid_index             = prm.get_integer("fluid index");
  }
  prm.leave_subsection();
}
