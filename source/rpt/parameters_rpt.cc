#include <rpt/parameters_rpt.h>
#include <time.h>

void
Parameters::RPTParameters::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("rpt parameters");
  {
    prm.declare_entry("particle positions file",
                      "none",
                      Patterns::FileName(),
                      "Particle positions file name");

    prm.declare_entry(
      "export counts",
      "false",
      Patterns::Bool(),
      "Enable to export counts result in a .csv file <true|false>");

    prm.declare_entry("monte carlo iteration",
                      "1",
                      Patterns::Integer(),
                      "Number of Monte Carlo iteration");

    prm.declare_entry("random number seed",
                      "auto",
                      Patterns::Anything(),
                      "Seed of the random number generator <auto|number>");

    prm.declare_entry("reactor radius",
                      "1",
                      Patterns::Double(),
                      "Radius of the reactor or tank");

    prm.declare_entry("peak-to-total ratio",
                      "1",
                      Patterns::Double(),
                      "Peak-to-total ratio");

    prm.declare_entry("sampling time",
                      "1",
                      Patterns::Double(),
                      "Sampling time of the counting");
  }
  prm.leave_subsection();
}

void
Parameters::RPTParameters::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("rpt parameters");
  {
    particle_positions_file = prm.get("particle positions file");
    export_counts           = prm.get_bool("export counts");
    n_monte_carlo_iteration = prm.get_integer("monte carlo iteration");
    reactor_radius          = prm.get_double("reactor radius");
    peak_to_total_ratio     = prm.get_double("peak-to-total ratio");
    sampling_time           = prm.get_double("sampling time");

    seed = (prm.get("random number seed") == "auto") ?
             time(NULL) :
             prm.get_integer("random number seed");
  }
  prm.leave_subsection();
}


void
Parameters::InitialRPTParameters::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("parameter tuning");
  {
    prm.declare_entry("tuning",
                      "false",
                      Patterns::Bool(),
                      "Enable parameter tuning <true|false>");

    prm.declare_entry("experimental data file",
                      "none",
                      Patterns::FileName(),
                      "Experimental counts data file name");

    prm.declare_entry("dead time",
                      "1",
                      Patterns::Double(),
                      "Dead time of the detector per accepted pulse");

    prm.declare_entry("activity",
                      "1",
                      Patterns::Double(),
                      "Activity of the tracer");

    prm.declare_entry("gamma-rays emitted",
                      "1",
                      Patterns::Double(),
                      "Number of gamma-rays emitted by each disintegration");

    prm.declare_entry("attenuation coefficient reactor",
                      "1",
                      Patterns::Double(),
                      "Total linear attenuation coefficient of the medium");

    prm.declare_entry("attenuation coefficient detector",
                      "1",
                      Patterns::Double(),
                      "Total linear attenuation coefficient of the detector");
  }
  prm.leave_subsection();
}

void
Parameters::InitialRPTParameters::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("parameter tuning");
  {
    tuning             = prm.get_bool("tuning");
    experimental_file  = prm.get("experimental data file");
    dead_time          = prm.get_double("dead time");
    activity           = prm.get_double("activity");
    gamma_rays_emitted = prm.get_double("gamma-rays emitted");
    attenuation_coefficient_reactor =
      prm.get_double("attenuation coefficient reactor");
    attenuation_coefficient_detector =
      prm.get_double("attenuation coefficient detector");
  }
  prm.leave_subsection();
}

void
Parameters::DetectorParameters::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("detector parameters");
  {
    prm.declare_entry("number", "1", Patterns::Integer(), "Number of detector");

    prm.declare_entry("radius",
                      "1",
                      Patterns::Double(),
                      "Radius of all detectors");

    prm.declare_entry("length",
                      "1",
                      Patterns::Double(),
                      "Length of all detectors");

    prm.declare_entry("detector positions file",
                      "none",
                      Patterns::FileName(),
                      "Detector positions file name");
  }
  prm.leave_subsection();
}

void
Parameters::DetectorParameters::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("detector parameters");
  {
    radius                  = prm.get_double("radius");
    length                  = prm.get_double("length");
    detector_positions_file = prm.get("detector positions file");
  }
  prm.leave_subsection();
}
