#include "core/parameters_lagrangian.h"

namespace Parameters
{
  namespace Lagrangian
  {
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
        prm.declare_entry("insertion_method",
                          "non_uniform",
                          Patterns::Selection("uniform|non_uniform"),
                          "Choosing insertion method. "
                          "Choices are <uniform|non_uniform>.");
        prm.declare_entry("n total",
                          "1",
                          Patterns::Integer(),
                          "Total number of particles");
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
                          "Random number range");
        prm.declare_entry("Insertion random number seed",
                          "1",
                          Patterns::Integer(),
                          "Random number seed");
      }
      prm.leave_subsection();
    }

    void
    InsertionInfo::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        const std::string insertion = prm.get("insertion_method");
        if (insertion == "uniform")
          insertion_method = InsertionMethod::uniform;
        else if (insertion == "non_uniform")
          insertion_method = InsertionMethod::non_uniform;
        else
          {
            std::runtime_error("Invalid insertion method");
          }
        total_particle_number = prm.get_integer("n total");
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
        random_number_range = prm.get_double("Insertion random number range");
        random_number_seed  = prm.get_double("Insertion random number seed");
      }
      prm.leave_subsection();
    }

    void
    ModelParameters::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        prm.declare_entry("contact_detection_frequency",
                          "1",
                          Patterns::Integer(),
                          "Particle-particle broad search frequency");

        prm.declare_entry(
          "neighborhood_threshold",
          "1",
          Patterns::Double(),
          "Contact search zone diameter to particle diameter ratio");

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

        prm.declare_entry("integration_method",
                          "velocity_verlet",
                          Patterns::Selection("velocity_verlet|explicit_euler"),
                          "Choosing integration method. "
                          "Choices are <velocity_verlet|explicit_euler>.");
      }
      prm.leave_subsection();
    }

    void
    ModelParameters::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        contact_detection_frequency =
          prm.get_integer("contact_detection_frequency");
        neighborhood_threshold = prm.get_double("neighborhood_threshold");

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

        const std::string integration = prm.get("integration_method");
        if (integration == "velocity_verlet")
          integration_method = IntegrationMethod::velocity_verlet;
        else if (integration == "explicit_euler")
          integration_method = IntegrationMethod::explicit_euler;
        else
          {
            std::runtime_error("Invalid integration method");
          }
      }
      prm.leave_subsection();
    }

  } // namespace Lagrangian

} // namespace Parameters
