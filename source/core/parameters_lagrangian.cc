#include "core/parameters_lagrangian.h"

#include <deal.II/grid/grid_in.h>

namespace Parameters
{
  namespace Lagrangian
  {
    void
    LagrangianPhysicalProperties::declareDefaultEntry(ParameterHandler &prm)
    {
      prm.declare_entry("size distribution type",
                        "uniform",
                        Patterns::Selection("uniform|normal"),
                        "Particle size distribution"
                        "Choices are <uniform|normall>.");
      prm.declare_entry("diameter",
                        "0.001",
                        Patterns::Double(),
                        "Particle diameter");
      prm.declare_entry("average diameter",
                        "0.001",
                        Patterns::Double(),
                        "Average particle diameter");
      prm.declare_entry("standard deviation",
                        "0",
                        Patterns::Double(),
                        "Particle size standard deviation");
      prm.declare_entry("number",
                        "0",
                        Patterns::Integer(),
                        "Number of particles of this type");
      prm.declare_entry("density particles",
                        "1000",
                        Patterns::Double(),
                        "Particle density");
      prm.declare_entry("young modulus particles",
                        "1000000",
                        Patterns::Double(),
                        "Particle Young's modulus");
      prm.declare_entry("poisson ratio particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle Poisson ratio");
      prm.declare_entry("restitution coefficient particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle restitution coefficient");
      prm.declare_entry("friction coefficient particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle friction coefficient");
      prm.declare_entry("rolling friction particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle rolling friction");
    }

    void
    LagrangianPhysicalProperties::parse_particle_properties(
      const unsigned int &particle_type,
      ParameterHandler &  prm)
    {
      const std::string size_distribution_type =
        prm.get("size distribution type");
      if (size_distribution_type == "uniform")
        {
          particle_average_diameter.at(particle_type) =
            prm.get_double("diameter");
        }
      else if (size_distribution_type == "normal")
        {
          particle_average_diameter.at(particle_type) =
            prm.get_double("average diameter");
          particle_size_std.at(particle_type) =
            prm.get_double("standard deviation");
        }
      else
        {
          throw(std::runtime_error("Invalid size distribution type "));
        }
      number.at(particle_type)           = prm.get_integer("number");
      density_particle.at(particle_type) = prm.get_double("density particles");
      youngs_modulus_particle.at(particle_type) =
        prm.get_double("young modulus particles");
      poisson_ratio_particle.at(particle_type) =
        prm.get_double("poisson ratio particles");
      restitution_coefficient_particle.at(particle_type) =
        prm.get_double("restitution coefficient particles");
      friction_coefficient_particle.at(particle_type) =
        prm.get_double("friction coefficient particles");
      rolling_friction_coefficient_particle.at(particle_type) =
        prm.get_double("rolling friction particles");
    }

    void
    LagrangianPhysicalProperties::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("lagrangian physical properties");
      {
        prm.declare_entry("gx",
                          "0.",
                          Patterns::Double(),
                          "Gravitational acceleration in x direction");
        prm.declare_entry("gy",
                          "0.",
                          Patterns::Double(),
                          "Gravitational acceleration in y direction");
        prm.declare_entry("gz",
                          "0.",
                          Patterns::Double(),
                          "Gravitational acceleration in z direction");


        prm.declare_entry("number of particle types",
                          "1",
                          Patterns::Integer(),
                          "Number of particle types");

        for (unsigned int counter = 0; counter < particle_type_maximum_number;
             ++counter)
          {
            prm.enter_subsection("particle type " +
                                 Utilities::int_to_string(counter, 1));
            {
              declareDefaultEntry(prm);
            }
            prm.leave_subsection();
          }

        prm.declare_entry("young modulus wall",
                          "1000000.",
                          Patterns::Double(),
                          "Young's modulus of wall");
        prm.declare_entry("poisson ratio wall",
                          "1000000.",
                          Patterns::Double(),
                          "Poisson's ratio of wall");
        prm.declare_entry("restitution coefficient wall",
                          "0.1",
                          Patterns::Double(),
                          "Coefficient of restitution of wall");
        prm.declare_entry("friction coefficient wall",
                          "0.1",
                          Patterns::Double(),
                          "Friction coefficient of wall");
        prm.declare_entry("rolling friction wall",
                          "0.1",
                          Patterns::Double(),
                          "Rolling friction coefficient of wall");
      }
      prm.leave_subsection();
    }

    void
    LagrangianPhysicalProperties::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("lagrangian physical properties");
      initialize_containers(particle_average_diameter,
                            particle_size_std,
                            number,
                            density_particle,
                            youngs_modulus_particle,
                            poisson_ratio_particle,
                            restitution_coefficient_particle,
                            friction_coefficient_particle,
                            rolling_friction_coefficient_particle);
      {
        g[0] = prm.get_double("gx");
        g[1] = prm.get_double("gy");
        g[2] = prm.get_double("gz");

        particle_type_number = prm.get_integer("number of particle types");

        if (particle_type_number >= 1)
          {
            prm.enter_subsection("particle type 0");
            {
              parse_particle_properties(0, prm);
            }
            prm.leave_subsection();
          }
        if (particle_type_number >= 2)
          {
            prm.enter_subsection("particle type 1");
            {
              parse_particle_properties(1, prm);
            }
            prm.leave_subsection();
          }
        if (particle_type_number >= 3)
          {
            prm.enter_subsection("particle type 2");
            {
              parse_particle_properties(2, prm);
            }
            prm.leave_subsection();
          }
        if (particle_type_number >= 4)
          {
            prm.enter_subsection("particle type 3");
            {
              parse_particle_properties(3, prm);
            }
            prm.leave_subsection();
          }
        if (particle_type_number >= 5)
          {
            prm.enter_subsection("particle type 4");
            {
              parse_particle_properties(4, prm);
            }
            prm.leave_subsection();
          }


        youngs_modulus_wall = prm.get_double("young modulus wall");
        poisson_ratio_wall  = prm.get_double("poisson ratio wall");
        restitution_coefficient_wall =
          prm.get_double("restitution coefficient wall");
        friction_coefficient_wall = prm.get_double("friction coefficient wall");
        rolling_friction_wall     = prm.get_double("rolling friction wall");
      }
      prm.leave_subsection();
    }

    void
    LagrangianPhysicalProperties::initialize_containers(
      std::unordered_map<unsigned int, double> &particle_average_diameter,
      std::unordered_map<unsigned int, double> &particle_size_std,
      std::unordered_map<unsigned int, int> &   number,
      std::unordered_map<unsigned int, double> &density_particle,
      std::unordered_map<unsigned int, double> &youngs_modulus_particle,
      std::unordered_map<unsigned int, double> &poisson_ratio_particle,
      std::unordered_map<unsigned int, double>
        &restitution_coefficient_particle,
      std::unordered_map<unsigned int, double> &friction_coefficient_particle,
      std::unordered_map<unsigned int, double>
        &rolling_friction_coefficient_particle)
    {
      for (unsigned int counter = 0; counter < particle_type_maximum_number;
           ++counter)
        {
          particle_average_diameter.insert({counter, 0.});
          particle_size_std.insert({counter, 0.});
          number.insert({counter, 0.});
          density_particle.insert({counter, 0.});
          youngs_modulus_particle.insert({counter, 0.});
          poisson_ratio_particle.insert({counter, 0.});
          restitution_coefficient_particle.insert({counter, 0.});
          friction_coefficient_particle.insert({counter, 0.});
          rolling_friction_coefficient_particle.insert({counter, 0.});
        }
    }

    void
    InsertionInfo::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        prm.declare_entry("insertion method",
                          "non_uniform",
                          Patterns::Selection("uniform|non_uniform|list"),
                          "Choosing insertion method. "
                          "Choices are <uniform|non_uniform|list>.");
        prm.declare_entry("inserted number of particles at each time step",
                          "1",
                          Patterns::Integer(),
                          "Inserted number of particles at each time step");
        prm.declare_entry("insertion frequency",
                          "1",
                          Patterns::Integer(),
                          "Insertion frequency");

        prm.declare_entry("insertion first direction",
                          "0",
                          Patterns::Integer(),
                          "First direction of particle insertion");
        prm.declare_entry("insertion second direction",
                          "1",
                          Patterns::Integer(),
                          "Second direction of particle insertion");
        prm.declare_entry("insertion third direction",
                          "2",
                          Patterns::Integer(),
                          "Third direction of particle insertion");

        prm.declare_entry("insertion box minimum x",
                          "1",
                          Patterns::Double(),
                          "Insertion x min");
        prm.declare_entry("insertion box minimum y",
                          "1",
                          Patterns::Double(),
                          "Insertion y min");
        prm.declare_entry("insertion box minimum z",
                          "1",
                          Patterns::Double(),
                          "Insertion z min");
        prm.declare_entry("insertion box maximum x",
                          "1",
                          Patterns::Double(),
                          "Insertion x max");
        prm.declare_entry("insertion box maximum y",
                          "1",
                          Patterns::Double(),
                          "Insertion y max");
        prm.declare_entry("insertion box maximum z",
                          "1",
                          Patterns::Double(),
                          "Insertion z max");
        prm.declare_entry("insertion distance threshold",
                          "1",
                          Patterns::Double(),
                          "Distance threshold");
        prm.declare_entry("insertion random number range",
                          "1",
                          Patterns::Double(),
                          "Random number range");
        prm.declare_entry("insertion random number seed",
                          "1",
                          Patterns::Integer(),
                          "Random number seed");

        prm.declare_entry("list x",
                          "0",
                          Patterns::List(Patterns::Double()),
                          "List of particles x positions");
        prm.declare_entry("list y",
                          "0",
                          Patterns::List(Patterns::Double()),
                          "List of particles y positions");
        prm.declare_entry("list z",
                          "0",
                          Patterns::List(Patterns::Double()),
                          "List of particles z positions");

        prm.declare_entry("velocity x",
                          "0.0",
                          Patterns::Double(),
                          "Initial velocity x");
        prm.declare_entry("velocity y",
                          "0.0",
                          Patterns::Double(),
                          "Initial velocity y");
        prm.declare_entry("velocity z",
                          "0.0",
                          Patterns::Double(),
                          "Initial velocity z");
        prm.declare_entry("omega x",
                          "0.0",
                          Patterns::Double(),
                          "Initial omega x");
        prm.declare_entry("omega y",
                          "0.0",
                          Patterns::Double(),
                          "Initial omega y");
        prm.declare_entry("omega z",
                          "0.0",
                          Patterns::Double(),
                          "Initial omega z");
      }
      prm.leave_subsection();
    }

    void
    InsertionInfo::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        const std::string insertion = prm.get("insertion method");
        if (insertion == "uniform")
          insertion_method = InsertionMethod::uniform;
        else if (insertion == "non_uniform")
          insertion_method = InsertionMethod::non_uniform;
        else if (insertion == "list")
          insertion_method = InsertionMethod::list;
        else
          {
            throw(std::runtime_error("Invalid insertion method "));
          }
        inserted_this_step =
          prm.get_integer("inserted number of particles at each time step");
        insertion_frequency = prm.get_integer("insertion frequency");
        axis_0              = prm.get_integer("insertion first direction");
        axis_1              = prm.get_integer("insertion second direction");
        axis_2              = prm.get_integer("insertion third direction");
        x_min               = prm.get_double("insertion box minimum x");
        y_min               = prm.get_double("insertion box minimum y");
        z_min               = prm.get_double("insertion box minimum z");
        x_max               = prm.get_double("insertion box maximum x");
        y_max               = prm.get_double("insertion box maximum y");
        z_max               = prm.get_double("insertion box maximum z");
        distance_threshold  = prm.get_double("insertion distance threshold");
        random_number_range = prm.get_double("insertion random number range");
        random_number_seed  = prm.get_double("insertion random number seed");

        vel_x   = prm.get_double("velocity x");
        vel_y   = prm.get_double("velocity y");
        vel_z   = prm.get_double("velocity z");
        omega_x = prm.get_double("omega x");
        omega_y = prm.get_double("omega y");
        omega_z = prm.get_double("omega z");

        // Read x, y and z list as a single string
        std::string x_str = prm.get("list x");
        std::string y_str = prm.get("list y");
        std::string z_str = prm.get("list z");

        // Convert x,y and z string to vector of strings
        std::vector<std::string> x_str_list(
          Utilities::split_string_list(x_str));
        std::vector<std::string> y_str_list(
          Utilities::split_string_list(y_str));
        std::vector<std::string> z_str_list(
          Utilities::split_string_list(z_str));

        // Convert x,y and z string vector to double vector
        list_x = Utilities::string_to_double(x_str_list);
        list_y = Utilities::string_to_double(y_str_list);
        list_z = Utilities::string_to_double(z_str_list);
      }
      prm.leave_subsection();
    }

    void
    ModelParameters::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        prm.declare_entry(
          "load balance method",
          "none",
          Patterns::Selection(
            "none|once|frequent|dynamic|dynamic_with_disabling_contacts"),
          "Choosing load-balance method"
          "Choices are <none|once|frequent|dynamic|dynamic_with_disabling_contacts>.");

        prm.declare_entry(
          "load balance step",
          "1000000000",
          Patterns::Integer(),
          "Step at which the triangulation is repartitioned "
          "and load is balanced for single-step load-balancing");

        prm.declare_entry(
          "load balance frequency",
          "1000000000",
          Patterns::Integer(),
          "Frequency at which the triangulation is repartitioned "
          "and load is balanced for frequent load-balancing");

        prm.declare_entry("load balance threshold",
                          "1",
                          Patterns::Double(),
                          "Threshold for dynamic load-balancing");

        prm.declare_entry("dynamic load balance check frequency",
                          "100000",
                          Patterns::Integer(),
                          "Checking frequency for dynamic load-balancing");

        prm.declare_entry(
          "load balance particle weight",
          "10000",
          Patterns::Integer(),
          "The particle weight based on a default cell weight of 1000");

        prm.declare_entry(
          "load balance active weight factor",
          "1.0",
          Patterns::Double(),
          "Factor applied on the particle weight in load balancing if the cell is active");

        prm.declare_entry(
          "load balance inactive weight factor",
          "1.0",
          Patterns::Double(),
          "Factor applied on the particle weight in load balancing if the cell is inactive");

        prm.declare_entry("contact detection method",
                          "dynamic",
                          Patterns::Selection("constant|dynamic"),
                          "Choosing contact detection method"
                          "Choices are <constant|dynamic>.");

        prm.declare_entry("contact detection frequency",
                          "1",
                          Patterns::Integer(),
                          "Particle-particle contact list");

        prm.declare_entry("dynamic contact search size coefficient",
                          "0.8",
                          Patterns::Double(),
                          "Security coefficient for dynamic contact detection");

        prm.declare_entry(
          "neighborhood threshold",
          "1",
          Patterns::Double(),
          "Contact search zone diameter to particle diameter ratio");

        prm.declare_entry(
          "particle particle contact force method",
          "hertz_mindlin_limit_overlap",
          Patterns::Selection(
            "linear|hertz_mindlin_limit_force|hertz_mindlin_limit_overlap|hertz"),
          "Choosing particle-particle contact force model"
          "Choices are <linear|hertz_mindlin_limit_force|hertz_mindlin_limit_overlap|hertz>.");

        prm.declare_entry("particle wall contact force method",
                          "nonlinear",
                          Patterns::Selection("linear|nonlinear"),
                          "Choosing particle-wall contact force model"
                          "Choices are <linear|nonlinear>.");

        prm.declare_entry(
          "rolling resistance torque method",
          "constant_resistance",
          Patterns::Selection(
            "no_resistance|constant_resistance|viscous_resistance"),
          "Choosing rolling resistance torque model"
          "Choices are <no_resistance|constant_resistance|viscous_resistance>.");

        prm.declare_entry(
          "integration method",
          "velocity_verlet",
          Patterns::Selection("velocity_verlet|explicit_euler|gear3"),
          "Choosing integration method"
          "Choices are <velocity_verlet|explicit_euler|gear3>.");

        prm.enter_subsection("dynamic disabling contacts");
        {
          prm.declare_entry(
            "enable dynamic disabling contacts",
            "false",
            Patterns::Selection("true|false"),
            "Enable the dynamic search for disabling particle contacts"
            "Choices are <true|false>.");

          prm.declare_entry(
            "enable particle advection",
            "false",
            Patterns::Selection("true|false"),
            "Enable the advection of particles with hydrodynamic forces"
            "Choices are <true|false>.");

          prm.declare_entry(
            "granular temperature threshold",
            "1e-4",
            Patterns::Double(),
            "Minimum granular temperature where particle contacts are considered");

          prm.declare_entry(
            "solid fraction threshold",
            "0.4",
            Patterns::Double(),
            "Maximum solid fraction where particle contacts are considered "
            "no matter the granular temperature");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    void
    ModelParameters::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        const std::string load_balance = prm.get("load balance method");

        prm.enter_subsection("dynamic disabling contacts");
        {
          disable_particle_contacts =
            prm.get_bool("enable dynamic disabling contacts");
          advect_particles = prm.get_bool("enable particle advection");

          // Thresholds for disabling contacts
          granular_temperature_threshold =
            prm.get_double("granular temperature threshold");
          solid_fraction_threshold = prm.get_double("solid fraction threshold");
        }
        prm.leave_subsection();

        load_balance_particle_weight =
          prm.get_integer("load balance particle weight");

        if (load_balance == "once")
          {
            load_balance_method = LoadBalanceMethod::once;

            load_balance_step = prm.get_integer("load balance step");
          }
        else if (load_balance == "frequent")
          {
            load_balance_method    = LoadBalanceMethod::frequent;
            load_balance_frequency = prm.get_integer("load balance frequency");
          }
        else if (load_balance == "dynamic")
          {
            load_balance_method    = LoadBalanceMethod::dynamic;
            load_balance_threshold = prm.get_double("load balance threshold");
            dynamic_load_balance_check_frequency =
              prm.get_integer("dynamic load balance check frequency");
          }
        else if (load_balance == "dynamic_with_disabling_contacts")
          {
            // Check if dynamic disabling contacts is enabled, otherwise throw
            // an error message indicating that the user should use dynamic
            // load balancing instead or enable dynamic disabling contacts
            if (disable_particle_contacts)
              {
                load_balance_method =
                  LoadBalanceMethod::dynamic_with_disabling_contacts;
                load_balance_threshold =
                  prm.get_double("load balance threshold");
                dynamic_load_balance_check_frequency =
                  prm.get_integer("dynamic load balance check frequency");

                // Weights for load balancing of active and inactive cells
                active_load_balancing_factor =
                  prm.get_double("load balance active weight factor");
                inactive_load_balancing_factor =
                  prm.get_double("load balance inactive weight factor");
              }
            else
              {
                throw(std::runtime_error(
                  "Invalid contact detection method: dynamic disabling contacts is not enabled "
                  "while dynamic_with_disabling_contacts is selected, use dynamic instead"));
              }
          }
        else if (load_balance == "none")
          {
            load_balance_method = LoadBalanceMethod::none;
          }
        else
          {
            throw(std::runtime_error("Invalid load-balance method "));
          }

        const std::string contact_search = prm.get("contact detection method");
        if (contact_search == "constant")
          {
            contact_detection_method = ContactDetectionMethod::constant;
            contact_detection_frequency =
              prm.get_integer("contact detection frequency");
          }
        else if (contact_search == "dynamic")
          {
            contact_detection_method = ContactDetectionMethod::dynamic;
            dynamic_contact_search_factor =
              prm.get_double("dynamic contact search size coefficient");
          }
        else
          {
            throw(std::runtime_error("Invalid contact detection method "));
          }

        neighborhood_threshold = prm.get_double("neighborhood threshold");

        const std::string ppcf =
          prm.get("particle particle contact force method");
        if (ppcf == "linear")
          particle_particle_contact_force_model =
            ParticleParticleContactForceModel::linear;
        else if (ppcf == "hertz_mindlin_limit_force")
          particle_particle_contact_force_model =
            ParticleParticleContactForceModel::hertz_mindlin_limit_force;
        else if (ppcf == "hertz_mindlin_limit_overlap")
          particle_particle_contact_force_model =
            ParticleParticleContactForceModel::hertz_mindlin_limit_overlap;
        else if (ppcf == "hertz")
          particle_particle_contact_force_model =
            ParticleParticleContactForceModel::hertz;
        else
          {
            throw(std::runtime_error(
              "Invalid particle-particle contact force model "));
          }

        const std::string pwcf = prm.get("particle wall contact force method");
        if (pwcf == "linear")
          particle_wall_contact_force_method =
            ParticleWallContactForceModel::linear;
        else if (pwcf == "nonlinear")
          particle_wall_contact_force_method =
            ParticleWallContactForceModel::nonlinear;
        else
          {
            throw(
              std::runtime_error("Invalid particle-wall contact force model "));
          }

        const std::string rolling_resistance_torque =
          prm.get("rolling resistance torque method");
        if (rolling_resistance_torque == "no_resistance")
          {
            rolling_resistance_method = RollingResistanceMethod::no_resistance;
          }
        else if (rolling_resistance_torque == "constant_resistance")
          {
            rolling_resistance_method =
              RollingResistanceMethod::constant_resistance;
          }
        else if (rolling_resistance_torque == "viscous_resistance")
          {
            rolling_resistance_method =
              RollingResistanceMethod::viscous_resistance;
          }
        else
          {
            throw(
              std::runtime_error("Invalid rolling resistance torque method "));
          }

        const std::string integration = prm.get("integration method");
        if (integration == "velocity_verlet")
          integration_method = IntegrationMethod::velocity_verlet;
        else if (integration == "explicit_euler")
          integration_method = IntegrationMethod::explicit_euler;
        else if (integration == "gear3")
          integration_method = IntegrationMethod::gear3;
        else
          {
            throw(std::runtime_error("Invalid integration method "));
          }
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ForceTorqueOnWall<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("boundary forces");
      prm.declare_entry("calculation",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of forces");
      prm.declare_entry(
        "verbosity",
        "verbose",
        Patterns::Selection("quiet|verbose"),
        "State whether output from solver runs should be printed. "
        "Choices are <quiet|verbose>.");
      prm.declare_entry("filename",
                        "force",
                        Patterns::FileName(),
                        "File output force prefix");
      prm.declare_entry("output frequency",
                        "1",
                        Patterns::Integer(),
                        "Output frequency");
      prm.enter_subsection("center of mass coordinate");
      prm.declare_entry("x",
                        "0",
                        Patterns::Double(),
                        "X coordinate of center of mass");
      prm.declare_entry("y",
                        "0",
                        Patterns::Double(),
                        "Y coordinate of center of mass");
      prm.declare_entry("z",
                        "0",
                        Patterns::Double(),
                        "Z coordinate of center of mass");
      prm.leave_subsection();

      prm.leave_subsection();
    }

    template <int dim>
    void
    ForceTorqueOnWall<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("boundary forces");
      calculate_force_torque    = prm.get_bool("calculation");
      const std::string verbose = prm.get("verbosity");
      if (verbose == "quiet")
        force_torque_verbosity = Parameters::Verbosity::quiet;
      else if (verbose == "verbose")
        force_torque_verbosity = Parameters::Verbosity::verbose;
      else
        {
          throw(std::runtime_error("Invalid verbosity choice "));
        }
      force_torque_output_name = prm.get("filename");
      output_frequency         = prm.get_integer("output frequency");
      prm.enter_subsection("center of mass coordinate");
      point_center_mass[0] = prm.get_double("x");
      point_center_mass[1] = prm.get_double("y");
      point_center_mass[2] = prm.get_double("z");
      prm.leave_subsection();
      prm.leave_subsection();
    }

    template <int dim>
    void
    FloatingWalls<dim>::declareDefaultEntry(ParameterHandler &prm)
    {
      prm.enter_subsection("point on wall");
      prm.declare_entry("x", "0.", Patterns::Double(), "X Point on wall");
      prm.declare_entry("y", "0.", Patterns::Double(), "Y Point on wall");
      prm.declare_entry("z", "0.", Patterns::Double(), "Z Point on wall");
      prm.leave_subsection();

      prm.enter_subsection("normal vector");
      prm.declare_entry("nx", "0.", Patterns::Double(), "X Normal vector wall");
      prm.declare_entry("ny", "0.", Patterns::Double(), "Y Normal vector wall");
      prm.declare_entry("nz", "0.", Patterns::Double(), "Z Normal vector wall");
      prm.leave_subsection();

      prm.declare_entry("start time", "0.", Patterns::Double(), "Start time");

      prm.declare_entry("end time", "0.", Patterns::Double(), "End time");
    }

    template <int dim>
    void
    FloatingWalls<dim>::parse_floating_wall(ParameterHandler &prm)
    {
      prm.enter_subsection("point on wall");
      Point<dim> wall_point;
      wall_point[0] = prm.get_double("x");
      wall_point[1] = prm.get_double("y");
      if (dim == 3)
        wall_point[2] = prm.get_double("z");
      this->points_on_walls.push_back(wall_point);
      prm.leave_subsection();

      prm.enter_subsection("normal vector");
      Tensor<1, dim> wall_normal;
      wall_normal[0] = prm.get_double("nx");
      wall_normal[1] = prm.get_double("ny");
      if (dim == 3)
        wall_normal[2] = prm.get_double("nz");
      this->floating_walls_normal_vectors.push_back(wall_normal);
      prm.leave_subsection();

      time_start.push_back(prm.get_double("start time"));
      time_end.push_back(prm.get_double("end time"));
    }

    template <int dim>
    void
    FloatingWalls<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("floating walls");
      {
        prm.declare_entry("number of floating walls",
                          "0",
                          Patterns::Integer(),
                          "Number of floating walls");

        for (unsigned int counter = 0; counter < max_number_floating_walls;
             ++counter)
          {
            prm.enter_subsection("wall " +
                                 Utilities::int_to_string(counter, 1));
            {
              declareDefaultEntry(prm);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    FloatingWalls<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("floating walls");
      {
        floating_walls_number = prm.get_integer("number of floating walls");

        if (floating_walls_number >= 1)
          {
            prm.enter_subsection("wall 0");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 2)
          {
            prm.enter_subsection("wall 1");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 3)
          {
            prm.enter_subsection("wall 2");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 4)
          {
            prm.enter_subsection("wall 3");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 5)
          {
            prm.enter_subsection("wall 4");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 6)
          {
            prm.enter_subsection("wall 5");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 7)
          {
            prm.enter_subsection("wall 6");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 8)
          {
            prm.enter_subsection("wall 7");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 9)
          {
            prm.enter_subsection("wall 8");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
        if (floating_walls_number >= 10)
          {
            prm.enter_subsection("wall 9");
            {
              parse_floating_wall(prm);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();
    }

    void
    BCDEM::declareDefaultEntry(ParameterHandler &prm)
    {
      prm.declare_entry("boundary id", "0", Patterns::Integer(), "Boundary ID");
      prm.declare_entry(
        "type",
        "fixed_wall",
        Patterns::Selection(
          "fixed_wall|outlet|translational|rotational|periodic"),
        "Type of boundary condition"
        "Choices are <fixed_wall|outlet|translational|rotational|periodic>.");

      prm.declare_entry("speed x",
                        "0.",
                        Patterns::Double(),
                        "Translational boundary speed in x direction");
      prm.declare_entry("speed y",
                        "0.",
                        Patterns::Double(),
                        "Translational boundary speed in y direction");
      prm.declare_entry("speed z",
                        "0.",
                        Patterns::Double(),
                        "Translational boundary speed in z direction");

      prm.declare_entry("rotational speed",
                        "0.",
                        Patterns::Double(),
                        "Rotational boundary speed");
      prm.declare_entry("rotational vector x",
                        "0.",
                        Patterns::Double(),
                        "Rotational vector element in x direction");
      prm.declare_entry("rotational vector y",
                        "0.",
                        Patterns::Double(),
                        "Rotational vector element in y direction");
      prm.declare_entry("rotational vector z",
                        "0.",
                        Patterns::Double(),
                        "Rotational vector element in z direction");
      prm.declare_entry("periodic id 0",
                        "0",
                        Patterns::Integer(),
                        "Periodic boundary ID 0");
      prm.declare_entry("periodic id 1",
                        "0",
                        Patterns::Integer(),
                        "Periodic boundary ID 1");
      prm.declare_entry(
        "periodic direction",
        "0",
        Patterns::Integer(),
        "Periodic direction or normal direction of periodic boundary");
    }

    void
    BCDEM::parse_boundary_conditions(ParameterHandler &prm)
    {
      const unsigned int boundary_id   = prm.get_integer("boundary id");
      const std::string  boundary_type = prm.get("type");

      if (boundary_type == "outlet")
        {
          BC_type = BoundaryType::outlet;
          this->outlet_boundaries.push_back(boundary_id);
        }
      else if (boundary_type == "translational")
        {
          BC_type = BoundaryType::translational;
          Tensor<1, 3> translational_velocity;
          translational_velocity[0] = prm.get_double("speed x");
          translational_velocity[1] = prm.get_double("speed y");
          translational_velocity[2] = prm.get_double("speed z");

          this->boundary_translational_velocity.at(boundary_id) =
            translational_velocity;
        }
      else if (boundary_type == "rotational")
        {
          BC_type                       = BoundaryType::rotational;
          double       rotational_speed = prm.get_double("rotational speed");
          Tensor<1, 3> rotational_vector;

          rotational_vector[0] = prm.get_double("rotational vector x");
          rotational_vector[1] = prm.get_double("rotational vector y");
          rotational_vector[2] = prm.get_double("rotational vector z");

          this->boundary_rotational_speed.at(boundary_id)  = rotational_speed;
          this->boundary_rotational_vector.at(boundary_id) = rotational_vector;
        }
      else if (boundary_type == "fixed_wall")
        {
          BC_type = BoundaryType::fixed_wall;
        }
      else if (boundary_type == "periodic")
        {
          BC_type                   = BoundaryType::periodic;
          this->periodic_boundary_0 = prm.get_integer("periodic id 0");
          this->periodic_boundary_1 = prm.get_integer("periodic id 1");
          this->periodic_direction  = prm.get_integer("periodic direction");
        }
      else
        {
          throw(std::runtime_error("Invalid boundary condition type "));
        }
    }

    void
    BCDEM::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("DEM boundary conditions");
      {
        prm.declare_entry("number of boundary conditions",
                          "0",
                          Patterns::Integer(),
                          "Number of boundary conditions");

        for (unsigned int counter = 0; counter < DEM_BC_number_max; ++counter)
          { // Example: "boundary condition 0"
            prm.enter_subsection("boundary condition " +
                                 Utilities::int_to_string(counter, 1));
            {
              declareDefaultEntry(prm);
            }
            prm.leave_subsection();
          }
      }
      prm.leave_subsection();
    }

    void
    BCDEM::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("DEM boundary conditions");

      DEM_BC_number = prm.get_integer("number of boundary conditions");

      initialize_containers(boundary_translational_velocity,
                            boundary_rotational_speed,
                            boundary_rotational_vector,
                            outlet_boundaries);

      for (unsigned int counter = 0; counter < DEM_BC_number; ++counter)
        {
          prm.enter_subsection("boundary condition " +
                               Utilities::int_to_string(counter, 1));
          {
            parse_boundary_conditions(prm);
          }
          prm.leave_subsection();
        }
      prm.leave_subsection();
    }

    void
    BCDEM::initialize_containers(
      std::unordered_map<unsigned int, Tensor<1, 3>>
        &                                       boundary_translational_velocity,
      std::unordered_map<unsigned int, double> &boundary_rotational_speed,
      std::unordered_map<unsigned int, Tensor<1, 3>>
        &                        boundary_rotational_vector,
      std::vector<unsigned int> &outlet_boundaries)
    {
      Tensor<1, 3> zero_tensor({0.0, 0.0, 0.0});

      for (unsigned int counter = 0; counter < DEM_BC_number_max; ++counter)
        {
          boundary_translational_velocity.insert({counter, zero_tensor});
          boundary_rotational_speed.insert({counter, 0});
          boundary_rotational_vector.insert({counter, zero_tensor});
        }

      outlet_boundaries.reserve(DEM_BC_number);
    }

    template <int dim>
    void
    GridMotion<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("grid motion");
      {
        prm.declare_entry(
          "motion type",
          "none",
          Patterns::Selection(
            "none|translational|rotational|translational_rotational"),
          "Choosing grid motion type. "
          "Choices are <none|translational|rotational|translational_rotational>.");

        prm.declare_entry("grid translational velocity x",
                          "0",
                          Patterns::Double(),
                          "grid translational velocity x");
        prm.declare_entry("grid translational velocity y",
                          "0",
                          Patterns::Double(),
                          "grid translational velocity y");
        prm.declare_entry("grid translational velocity z",
                          "0",
                          Patterns::Double(),
                          "grid translational velocity z");

        prm.declare_entry("grid rotational speed",
                          "0",
                          Patterns::Double(),
                          "grid rotational speed");

        prm.declare_entry("grid rotational axis",
                          "0",
                          Patterns::Integer(),
                          "grid rotational axis");
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    GridMotion<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("grid motion");
      {
        const std::string motion = prm.get("motion type");
        if (motion == "rotational")
          {
            motion_type           = MotionType::rotational;
            grid_rotational_speed = prm.get_double("grid rotational speed");
            grid_rotational_axis  = prm.get_integer("grid rotational axis");
          }
        else if (motion == "translational")
          {
            motion_type = MotionType::translational;
            grid_translational_velocity[0] =
              prm.get_double("grid translational velocity x");
            grid_translational_velocity[1] =
              prm.get_double("grid translational velocity y");
            if (dim == 3)
              grid_translational_velocity[2] =
                prm.get_double("grid translational velocity z");
          }
        else if (motion == "none")
          {
            motion_type = MotionType::none;
          }
        else
          {
            throw(std::runtime_error("Invalid grid motion "));
          }
      }
      prm.leave_subsection();
    }

    void
    LagrangianPostProcessing::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("post-processing");
      {
        prm.declare_entry(
          "Lagrangian post-processing",
          "false",
          Patterns::Bool(),
          "State whether Lagrangian post-processing should be performed.");
      }
      prm.leave_subsection();
    }

    void
    LagrangianPostProcessing::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("post-processing");
      {
        Lagrangian_post_processing = prm.get_bool("Lagrangian post-processing");
      }
      prm.leave_subsection();
    }
    template <int dim>
    void
    FloatingGrid<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("floating grid");
      {
        mesh.declare_parameters(prm);
        motion.declare_parameters(prm);
        prm.declare_entry("start time", "0.", Patterns::Double(), "Start time");
        prm.declare_entry("end time", "0.", Patterns::Double(), "End time");
      }

      prm.leave_subsection();
    }

    template <int dim>
    void
    FloatingGrid<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("floating grid");
      {
        mesh.parse_parameters(prm);
        motion.parse_parameters(prm);
        time_start = prm.get_double("start time");
        time_end   = prm.get_double("end time");
      }

      prm.leave_subsection();
    }

    template class ForceTorqueOnWall<2>;
    template class ForceTorqueOnWall<3>;
    template class FloatingWalls<2>;
    template class FloatingWalls<3>;
    template class FloatingGrid<2>;
    template class FloatingGrid<3>;
    template class GridMotion<2>;
    template class GridMotion<3>;

  } // namespace Lagrangian

} // namespace Parameters
