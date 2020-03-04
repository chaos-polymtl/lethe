#include "core/parameters_lagrangian.h"

namespace Parameters
{
  namespace Lagrangian
  {
    void
    SimulationControl::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("simulation control");
      {
        prm.declare_entry("time step",
                          "1.",
                          Patterns::Double(),
                          "Time step value");
        prm.declare_entry("step end", "1", Patterns::Integer(), "End step");
        prm.declare_entry("n total",
                          "1",
                          Patterns::Integer(),
                          "Total number of particles");
        prm.declare_entry("write frequency",
                          "1",
                          Patterns::Integer(),
                          "Write frequency");
      }
      prm.leave_subsection();
    }

    void
    SimulationControl::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("simulation control");
      {
        dt                    = prm.get_double("time step");
        final_time_step       = prm.get_integer("step end");
        total_particle_number = prm.get_integer("n total");
        write_frequency       = prm.get_integer("write frequency");
      }
      prm.leave_subsection();
    }

    void
    PhysicalProperties::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("physical properties");
      {
        prm.declare_entry("gx",
                          "1.",
                          Patterns::Double(),
                          "Gravitational acceleration in x direction");
        prm.declare_entry("gy",
                          "1.",
                          Patterns::Double(),
                          "Gravitational acceleration in y direction");
        prm.declare_entry("gz",
                          "1.",
                          Patterns::Double(),
                          "Gravitational acceleration in z direction");
        prm.declare_entry("diameter",
                          "1.",
                          Patterns::Double(),
                          "Particle diameter");
        prm.declare_entry("density",
                          "1.",
                          Patterns::Double(),
                          "Particle density");
        prm.declare_entry("Youngs_modulus_particle",
                          "1.",
                          Patterns::Double(),
                          "Young's modulus of particle");
        prm.declare_entry("Youngs_modulus_wall",
                          "1.",
                          Patterns::Double(),
                          "Young's modulus of wall");
        prm.declare_entry("Poisson_ratio_particle",
                          "1.",
                          Patterns::Double(),
                          "Poisson's ratio of particle");
        prm.declare_entry("Poisson_ratio_wall",
                          "1.",
                          Patterns::Double(),
                          "Poisson's ratio of wall");
        prm.declare_entry("restitution_coefficient_particle",
                          "1.",
                          Patterns::Double(),
                          "Coefficient of restitution of particle");
        prm.declare_entry("restitution_coefficient_wall",
                          "1.",
                          Patterns::Double(),
                          "Coefficient of restitution of wall");
        prm.declare_entry("friction_coefficient_particle",
                          "1.",
                          Patterns::Double(),
                          "Friction coefficient of particle");
        prm.declare_entry("friction_coefficient_wall",
                          "1.",
                          Patterns::Double(),
                          "Friction coefficient of wall");
        prm.declare_entry("rolling_friction_particle",
                          "1.",
                          Patterns::Double(),
                          "Rolling friction coefficient of particle");
        prm.declare_entry("rolling_friction_wall",
                          "1.",
                          Patterns::Double(),
                          "Rolling friction coefficient of wall");
      }
      prm.leave_subsection();
    }

    void
    PhysicalProperties::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("physical properties");
      {
        gx                      = prm.get_double("gx");
        gy                      = prm.get_double("gy");
        gz                      = prm.get_double("gz");
        diameter                = prm.get_double("diameter");
        density                 = prm.get_double("density");
        Youngs_modulus_particle = prm.get_integer("Youngs_modulus_particle");
        Youngs_modulus_wall     = prm.get_integer("Youngs_modulus_wall");
        Poisson_ratio_particle  = prm.get_double("Poisson_ratio_particle");
        Poisson_ratio_wall      = prm.get_double("Poisson_ratio_wall");
        restitution_coefficient_particle =
          prm.get_double("restitution_coefficient_particle");
        restitution_coefficient_wall =
          prm.get_double("restitution_coefficient_wall");
        friction_coefficient_particle =
          prm.get_double("friction_coefficient_particle");
        friction_coefficient_wall = prm.get_double("friction_coefficient_wall");
        rolling_friction_particle = prm.get_double("rolling_friction_particle");
        rolling_friction_wall     = prm.get_double("rolling_friction_wall");
      }
      prm.leave_subsection();
    }

    void
    InsertionInfo::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        prm.declare_entry("Insertion time step",
                          "1",
                          Patterns::Integer(),
                          "Insertion time steps");
        prm.declare_entry("Inserted number of particles at each time step",
                          "1",
                          Patterns::Integer(),
                          "Inserted number of particles at each time step");
        prm.declare_entry("Insertion frequency",
                          "1",
                          Patterns::Integer(),
                          "Insertion frequncy");
        prm.declare_entry("Insertion box minimum x",
                          "1",
                          Patterns::Double(),
                          "Insertion x min");
        prm.declare_entry("Insertion box minimum y",
                          "1",
                          Patterns::Double(),
                          "Insertion y min");
        prm.declare_entry("Insertion box minimum z",
                          "1",
                          Patterns::Double(),
                          "Insertion z min");
        prm.declare_entry("Insertion box maximum x",
                          "1",
                          Patterns::Double(),
                          "Insertion x max");
        prm.declare_entry("Insertion box maximum y",
                          "1",
                          Patterns::Double(),
                          "Insertion y max");
        prm.declare_entry("Insertion box maximum z",
                          "1",
                          Patterns::Double(),
                          "Insertion z max");
        prm.declare_entry("Insertion distance threshold",
                          "1",
                          Patterns::Double(),
                          "Distance threshold");
        prm.declare_entry("Insertion random number range",
                          "1",
                          Patterns::Double(),
                          "Random number bin");
      }
      prm.leave_subsection();
    }

    void
    InsertionInfo::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        insertion_steps_number = prm.get_integer("Insertion time step");
        inserted_this_step =
          prm.get_integer("Inserted number of particles at each time step");
        insertion_frequency = prm.get_integer("Insertion frequency");
        x_min               = prm.get_double("Insertion box minimum x");
        y_min               = prm.get_double("Insertion box minimum y");
        z_min               = prm.get_double("Insertion box minimum z");
        x_max               = prm.get_double("Insertion box maximum x");
        y_max               = prm.get_double("Insertion box maximum y");
        z_max               = prm.get_double("Insertion box maximum z");
        distance_threshold  = prm.get_double("Insertion distance threshold");
        random_number_bin   = prm.get_double("Insertion random number range");
      }
      prm.leave_subsection();
    }

    void
    OutputProperties::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("output properties");
      {
        prm.declare_entry("Number of properties",
                          "1",
                          Patterns::Integer(),
                          "Number of properties for visualization");

        prm.declare_entry("Output directory",
                          "1",
                          Patterns::FileName(),
                          "Visulization output folder");

        prm.declare_entry("General information file prefix",
                          "1",
                          Patterns::FileName(),
                          "General information file (.pvtu) prefix");

        prm.declare_entry("Result name prefix",
                          "1",
                          Patterns::FileName(),
                          "Result (.vtu) name prefix");
      }
      prm.leave_subsection();
    }

    void
    OutputProperties::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("output properties");
      {
        properties_number   = prm.get_integer("Number of properties");
        output_folder       = prm.get("Output directory");
        general_file_prefix = prm.get("General information file prefix");
        result_prefix       = prm.get("Result name prefix");
      }

      prm.leave_subsection();
    }

    void
    ModelParameters::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        prm.declare_entry("pp_broad_search_frequency",
                          "1",
                          Patterns::Integer(),
                          "Particle-particle broad search frequency");

        prm.declare_entry("pw_broad_search_frequency",
                          "1",
                          Patterns::Integer(),
                          "Particle-wall broad search frequency");

        prm.declare_entry("pp_contact_force_method",
                          "pp_nonlinear",
                          Patterns::Selection("pp_linear|pp_nonlinear"),
                          "Choosing particle-particle contact force model. "
                          "Choices are <pp_linear|pp_nonlinear>.");

        prm.declare_entry("pw_contact_force_method",
                          "pw_nonlinear",
                          Patterns::Selection("pw_linear|pw_nonlinear"),
                          "Choosing particle-wall contact force model. "
                          "Choices are <pw_linear|pw_nonlinear>.");
      }
      prm.leave_subsection();
    }

    void
    ModelParameters::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        pp_broad_search_frequency =
          prm.get_integer("pp_broad_search_frequency");
        pw_broad_search_frequency =
          prm.get_integer("pw_broad_search_frequency");

        const std::string ppcf = prm.get("pp_contact_force_method");
        if (ppcf == "pp_linear")
          pp_contact_force_method = PPContactForceModel::pp_linear;
        else if (ppcf == "pp_nonlinear")
          pp_contact_force_method = PPContactForceModel::pp_nonlinear;

        else
          {
            std::runtime_error("Invalid particle-particle contact force model");
          }

        const std::string pwcf = prm.get("pw_contact_force_method");
        if (pwcf == "pw_linear")
          pw_contact_force_method = PWContactForceModel::pw_linear;
        else if (pwcf == "pw_nonlinear")
          pw_contact_force_method = PWContactForceModel::pw_nonlinear;

        else
          {
            std::runtime_error("Invalid particle-wall contact force model");
          }
      }
      prm.leave_subsection();
    }

  } // namespace Lagrangian

} // namespace Parameters
