#include <core/dimensionality.h>


using namespace dealii;
namespace Parameters
{
  void
  Dimensionality::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("dimensionality");
    {
      prm.declare_entry(
        "length",
        "1",
        Patterns::Double(),
        "Length scale. The default value assumed is meters. If the simulation is carried out in centimeters,"
        " a length of 0.01 should be specified. If the simulation is carried out in kilometers, a length of 1000 should be specified.");

      prm.declare_entry(
        "time",
        "1",
        Patterns::Double(),
        "Time scale. The default value assumed is seconds. If the simulation is carried out in milliseconds,"
        " a time of 0.001 should be specified.");

      prm.declare_entry(
        "mass",
        "1",
        Patterns::Double(),
        "Mass scale. The default value assumed is kilogram. If the simulation is carried out in grams,"
        " a mass of 0.001 should be specified.");

      prm.declare_entry(
        "temperature",
        "1",
        Patterns::Double(),
        "Temperature scale. The default value assumed is Kelvin. If the simulation is carried out in kiloKelvin,"
        " a temperature of 1000 should be specified.");
    }
    prm.leave_subsection();
  }


  void
  Dimensionality::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("dimensionality");
    {
      length      = prm.get_double("length");
      time        = prm.get_double("time");
      mass        = prm.get_double("mass");
      temperature = prm.get_double("temperature");
    }
    prm.leave_subsection();

    // Define all scales using these dimensions
    define_all_scales();
  }
} // namespace Parameters
