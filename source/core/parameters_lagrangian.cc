// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_lagrangian.h>

#include <deal.II/grid/grid_in.h>

namespace Parameters
{
  namespace Lagrangian
  {
    void
    LagrangianPhysicalProperties::declare_parameters(
      ParameterHandler &prm) const
    {
      prm.enter_subsection("lagrangian physical properties");
      {
        // Parameter <g> is a list of values, its deprecated version are
        // individual parameters <gx>, <gy> and <gz>
        prm.declare_entry("g",
                          "0., 0., 0.",
                          Patterns::List(Patterns::Double()),
                          "Gravitational acceleration vector");
        prm.declare_alias("g", "gx", true);
        prm.declare_entry(
          "gy",
          "0.",
          Patterns::Double(),
          "Gravitational acceleration in y direction (deprecated, use <g> as vector)");
        prm.declare_entry(
          "gz",
          "0.",
          Patterns::Double(),
          "Gravitational acceleration in z direction (deprecated, use <g> as vector)");

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
                          "0.3",
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
        prm.declare_entry("rolling viscous damping wall",
                          "0.1",
                          Patterns::Double(),
                          "Rolling viscous damping wall");
        prm.declare_entry("surface energy wall",
                          "0.0",
                          Patterns::Double(),
                          "Surface energy of wall");
        prm.declare_entry("hamaker constant wall",
                          "4.e-19",
                          Patterns::Double(),
                          "Hamaker constant of wall");
        prm.declare_entry("thermal conductivity wall",
                          "100",
                          Patterns::Double(),
                          "Thermal conductivity of wall");
        prm.declare_entry("microhardness wall",
                          "1.e9",
                          Patterns::Double(),
                          "Microhardness of wall");
        prm.declare_entry("surface slope wall",
                          "0.1",
                          Patterns::Double(),
                          "Surface slope of wall");
        prm.declare_entry("surface roughness wall",
                          "1.e-10",
                          Patterns::Double(),
                          "Surface roughness of wall");
        prm.declare_entry("thermal accommodation wall",
                          "0.7",
                          Patterns::Double(),
                          "Thermal accommodation of wall");
        prm.declare_entry("real young modulus wall",
                          "0.",
                          Patterns::Double(),
                          "Real Young's modulus of wall");

        prm.declare_entry("thermal conductivity gas",
                          "0.01",
                          Patterns::Double(),
                          "Thermal conductivity of interstitial gas");
        prm.declare_entry("specific heat gas",
                          "1000",
                          Patterns::Double(),
                          "Specific heat of interstitial gas");
        prm.declare_entry("dynamic viscosity gas",
                          "1.e-5",
                          Patterns::Double(),
                          "Dynamic viscosity of interstitial gas");
        prm.declare_entry("specific heats ratio gas",
                          "1",
                          Patterns::Double(),
                          "Specific heats ratio of interstitial gas");
        prm.declare_entry("molecular mean free path gas",
                          "68.e-9",
                          Patterns::Double(),
                          "Molecular mean free path of interstitial gas");
      }
      prm.leave_subsection();
    }

    void
    LagrangianPhysicalProperties::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("lagrangian physical properties");
      initialize_containers(particle_average_diameter,
                            particle_size_std,
                            distribution_type,
                            particle_custom_diameter,
                            particle_custom_probability,
                            seed_for_distributions,
                            diameter_min_cutoff,
                            diameter_max_cutoff,
                            number,
                            density_particle,
                            youngs_modulus_particle,
                            poisson_ratio_particle,
                            restitution_coefficient_particle,
                            friction_coefficient_particle,
                            rolling_viscous_damping_coefficient_particle,
                            rolling_friction_coefficient_particle,
                            surface_energy_particle,
                            hamaker_constant_particle,
                            thermal_conductivity_particle,
                            specific_heat_particle,
                            microhardness_particle,
                            surface_slope_particle,
                            surface_roughness_particle,
                            thermal_accommodation_particle,
                            real_youngs_modulus_particle);

      // Deprecated parameter handling
      // <g> used to be 3 parameters: <gx>, <gy> and <gz>
      // If <gx> is in the input file, it will be used as the value for <g>
      // as an alias. This way, the parameter <g> is not a tensor and allows the
      // parsing of deprecated parameters.
      g = value_string_to_tensor<3>(prm.get("g"),
                                    prm.get_double("gy"),
                                    prm.get_double("gz"));

      particle_type_number = prm.get_integer("number of particle types");

      for (unsigned int id = 0; id < particle_type_number; ++id)
        {
          prm.enter_subsection("particle type " +
                               Utilities::int_to_string(id, 1));
          parse_particle_properties(id, prm);
          prm.leave_subsection();
        }

      youngs_modulus_wall = prm.get_double("young modulus wall");
      poisson_ratio_wall  = prm.get_double("poisson ratio wall");
      restitution_coefficient_wall =
        prm.get_double("restitution coefficient wall");
      friction_coefficient_wall = prm.get_double("friction coefficient wall");
      rolling_friction_wall     = prm.get_double("rolling friction wall");
      rolling_viscous_damping_wall =
        prm.get_double("rolling viscous damping wall");
      surface_energy_wall        = prm.get_double("surface energy wall");
      hamaker_constant_wall      = prm.get_double("hamaker constant wall");
      thermal_conductivity_wall  = prm.get_double("thermal conductivity wall");
      microhardness_wall         = prm.get_double("microhardness wall");
      surface_slope_wall         = prm.get_double("surface slope wall");
      surface_roughness_wall     = prm.get_double("surface roughness wall");
      thermal_accommodation_wall = prm.get_double("thermal accommodation wall");
      real_youngs_modulus_wall   = prm.get_double("real young modulus wall");
      // Only use the real Young's modulus if it is higher than the Young's
      // modulus
      if (real_youngs_modulus_wall < youngs_modulus_wall)
        real_youngs_modulus_wall = youngs_modulus_wall;

      thermal_conductivity_gas = prm.get_double("thermal conductivity gas");
      specific_heat_gas        = prm.get_double("specific heat gas");
      dynamic_viscosity_gas    = prm.get_double("dynamic viscosity gas");
      specific_heats_ratio_gas = prm.get_double("specific heats ratio gas");
      molecular_mean_free_path_gas =
        prm.get_double("molecular mean free path gas");

      prm.leave_subsection();
    }

    void
    LagrangianPhysicalProperties::declareDefaultEntry(ParameterHandler &prm)
    {
      prm.declare_entry("size distribution type",
                        "uniform",
                        Patterns::Selection("uniform|normal|lognormal|custom"),
                        "Particle size distribution"
                        "Choices are <uniform|normal|lognormal|custom>.");
      prm.declare_entry("average diameter",
                        "0.001",
                        Patterns::Double(),
                        "Particle diameter");
      prm.declare_alias("average diameter", "diameter", false);
      prm.declare_entry("standard deviation",
                        "0",
                        Patterns::Double(),
                        "Particle size standard deviation");
      prm.declare_entry("custom diameters",
                        "0.001 , 0.0005",
                        Patterns::List(Patterns::Double()),
                        "Diameter values for a custom distribution");
      prm.declare_entry("custom volume fractions",
                        "0.6 , 0.4",
                        Patterns::List(Patterns::Double()),
                        "Probabilities of each diameter of the custom"
                        " distribution based on the volume fraction");
      prm.declare_entry("random seed distribution",
                        "1",
                        Patterns::Integer(),
                        "Seed for generation of random numbers"
                        " for the size distribution");
      prm.declare_entry("minimum diameter cutoff",
                        "-1.",
                        Patterns::Double(),
                        "Cutoff values used when the log-normal distribution "
                        "is used. If equal to -1., the cut of will be fixed at "
                        "0.1% of the cumulative density function of the "
                        "log-normal distribution");
      prm.declare_entry("maximum diameter cutoff",
                        "-1.",
                        Patterns::Double(),
                        "Cutoff values used when the log-normal distribution "
                        "is used. If equal to -1., the cut of will be fixed at "
                        "99.9% of the cumulative density function of the "
                        "log-normal distribution");
      prm.declare_entry("number of particles",
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
                        "0.3",
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
      prm.declare_entry("rolling viscous damping particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle rolling viscous damping");
      prm.declare_entry("rolling friction particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle rolling friction");
      prm.declare_entry("surface energy particles",
                        "0.0",
                        Patterns::Double(),
                        "Particle surface energy");
      prm.declare_entry("hamaker constant particles",
                        "4.e-19",
                        Patterns::Double(),
                        "Material Hamaker constant");
      prm.declare_entry("thermal conductivity particles",
                        "1",
                        Patterns::Double(),
                        "Particle thermal conductivity");
      prm.declare_entry("specific heat particles",
                        "1000",
                        Patterns::Double(),
                        "Particle specific heat");
      prm.declare_entry("microhardness particles",
                        "1.e9",
                        Patterns::Double(),
                        "Particle microhardness");
      prm.declare_entry("surface slope particles",
                        "0.1",
                        Patterns::Double(),
                        "Particle surface slope");
      prm.declare_entry("surface roughness particles",
                        "1.e-9",
                        Patterns::Double(),
                        "Particle surface roughness");
      prm.declare_entry("thermal accommodation particles",
                        "0.7",
                        Patterns::Double(),
                        "Particle thermal accommodation");
      prm.declare_entry("real young modulus particles",
                        "0.",
                        Patterns::Double(),
                        "Particle real Young's modulus");
    }

    void
    LagrangianPhysicalProperties::parse_particle_properties(
      const unsigned int     &particle_type,
      const ParameterHandler &prm)
    {
      // unordered maps
      particle_average_diameter.at(particle_type) =
        prm.get_double("average diameter");
      particle_size_std.at(particle_type) =
        prm.get_double("standard deviation");
      particle_custom_diameter.at(particle_type) =
        convert_string_to_vector<double>(prm, "custom diameters");
      particle_custom_probability.at(particle_type) =
        convert_string_to_vector<double>(prm, "custom volume fractions");

      // vectors
      seed_for_distributions.push_back(
        prm.get_integer("random seed distribution"));
      diameter_min_cutoff.push_back(prm.get_double("minimum diameter cutoff"));
      diameter_max_cutoff.push_back(prm.get_double("maximum diameter cutoff"));

      double probability_sum =
        std::reduce(particle_custom_probability.at(particle_type).begin(),
                    particle_custom_probability.at(particle_type).end());

      // We make sure that the cumulative probability is equal to 1.
      if (std::abs(probability_sum - 1.0) > 1.e-5)
        {
          throw(std::runtime_error(
            "Invalid custom volume fraction. The sum of volume fractions should be equal to 1.0 "));
        }
      const std::string size_distribution_type_str =
        prm.get("size distribution type");
      if (size_distribution_type_str == "uniform")
        {
          distribution_type.at(particle_type) = SizeDistributionType::uniform;
        }
      else if (size_distribution_type_str == "normal")
        {
          distribution_type.at(particle_type) = SizeDistributionType::normal;
        }
      else if (size_distribution_type_str == "lognormal")
        {
          distribution_type.at(particle_type) = SizeDistributionType::lognormal;
        }
      else if (size_distribution_type_str == "custom")
        {
          distribution_type.at(particle_type) = SizeDistributionType::custom;
        }
      else
        {
          throw(std::runtime_error(
            "Invalid size distribution type. Choices are <uniform|normal|custom>."));
        }
      number.at(particle_type) = prm.get_integer("number of particles");
      density_particle.at(particle_type) = prm.get_double("density particles");
      youngs_modulus_particle.at(particle_type) =
        prm.get_double("young modulus particles");
      poisson_ratio_particle.at(particle_type) =
        prm.get_double("poisson ratio particles");
      restitution_coefficient_particle.at(particle_type) =
        prm.get_double("restitution coefficient particles");
      friction_coefficient_particle.at(particle_type) =
        prm.get_double("friction coefficient particles");
      rolling_viscous_damping_coefficient_particle.at(particle_type) =
        prm.get_double("rolling viscous damping particles");
      rolling_friction_coefficient_particle.at(particle_type) =
        prm.get_double("rolling friction particles");
      surface_energy_particle.at(particle_type) =
        prm.get_double("surface energy particles");
      hamaker_constant_particle.at(particle_type) =
        prm.get_double("hamaker constant particles");
      thermal_conductivity_particle.at(particle_type) =
        prm.get_double("thermal conductivity particles");
      specific_heat_particle.at(particle_type) =
        prm.get_double("specific heat particles");
      microhardness_particle.at(particle_type) =
        prm.get_double("microhardness particles");
      surface_slope_particle.at(particle_type) =
        prm.get_double("surface slope particles");
      surface_roughness_particle.at(particle_type) =
        prm.get_double("surface roughness particles");
      thermal_accommodation_particle.at(particle_type) =
        prm.get_double("thermal accommodation particles");
      real_youngs_modulus_particle.at(particle_type) =
        prm.get_double("real young modulus particles");
      // Only use the real Young's modulus if it is higher than the Young's
      // modulus
      if (real_youngs_modulus_particle.at(particle_type) <
          youngs_modulus_particle.at(particle_type))
        {
          real_youngs_modulus_particle.at(particle_type) =
            youngs_modulus_particle.at(particle_type);
        }
    }

    void
    LagrangianPhysicalProperties::initialize_containers(
      std::unordered_map<unsigned int, double>              &p_average_diameter,
      std::unordered_map<unsigned int, double>              &p_size_std,
      std::vector<SizeDistributionType>                     &dist_type,
      std::unordered_map<unsigned int, std::vector<double>> &p_custom_diameter,
      std::unordered_map<unsigned int, std::vector<double>>
                                               &p_custom_probability,
      std::vector<unsigned int>                &seed_for_dist,
      std::vector<double>                      &d_min_cutoff,
      std::vector<double>                      &d_max_cutoff,
      std::unordered_map<unsigned int, int>    &p_number,
      std::unordered_map<unsigned int, double> &p_density,
      std::unordered_map<unsigned int, double> &p_youngs_modulus,
      std::unordered_map<unsigned int, double> &p_poisson_ratio,
      std::unordered_map<unsigned int, double> &p_restitution_coefficient,
      std::unordered_map<unsigned int, double> &p_friction_coefficient,
      std::unordered_map<unsigned int, double>
        &p_rolling_viscous_damping_coefficient,
      std::unordered_map<unsigned int, double> &p_rolling_friction_coefficient,
      std::unordered_map<unsigned int, double> &p_surface_energy,
      std::unordered_map<unsigned int, double> &hamaker_constant_p,
      std::unordered_map<unsigned int, double> &thermal_conductivity_p,
      std::unordered_map<unsigned int, double> &specific_heat_p,
      std::unordered_map<unsigned int, double> &microhardness_p,
      std::unordered_map<unsigned int, double> &surface_slope_p,
      std::unordered_map<unsigned int, double> &surface_roughness_p,
      std::unordered_map<unsigned int, double> &thermal_accommodation_p,
      std::unordered_map<unsigned int, double> &real_youngs_modulus_p) const
    {
      for (unsigned int counter = 0; counter < particle_type_maximum_number;
           ++counter)
        {
          p_average_diameter.insert({counter, 0.});
          p_size_std.insert({counter, 0.});
          dist_type.push_back(SizeDistributionType::uniform);
          p_custom_diameter.insert({counter, {0.}});
          p_custom_probability.insert({counter, {1.}});
          p_number.insert({counter, 0});
          p_density.insert({counter, 0.});
          p_youngs_modulus.insert({counter, 0.});
          p_poisson_ratio.insert({counter, 0.});
          p_restitution_coefficient.insert({counter, 0.});
          p_friction_coefficient.insert({counter, 0.});
          p_rolling_viscous_damping_coefficient.insert({counter, 0.});
          p_rolling_friction_coefficient.insert({counter, 0.});
          p_surface_energy.insert({counter, 0.});
          hamaker_constant_p.insert({counter, 0.});
          thermal_conductivity_p.insert({counter, 0.});
          specific_heat_p.insert({counter, 0.});
          microhardness_p.insert({counter, 0.});
          surface_slope_p.insert({counter, 0.});
          surface_roughness_p.insert({counter, 0.});
          thermal_accommodation_p.insert({counter, 0.});
          real_youngs_modulus_p.insert({counter, 0.});
        }
      seed_for_dist.reserve(particle_type_maximum_number);
      d_min_cutoff.reserve(particle_type_maximum_number);
      d_max_cutoff.reserve(particle_type_maximum_number);
    }

    template <int dim>
    void
    InsertionInfo<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        prm.declare_entry("insertion method",
                          "volume",
                          Patterns::Selection("file|list|plane|volume"),
                          "Choosing insertion method. "
                          "Choices are <file|plane|list|volume>.");
        prm.declare_entry("inserted number of particles at each time step",
                          "0",
                          Patterns::Integer(),
                          "Inserted number of particles at each time step");
        prm.declare_entry("insertion frequency",
                          "0",
                          Patterns::Integer(),
                          "Insertion frequency");

        // Removal box:
        prm.declare_entry(
          "remove particles",
          "false",
          Patterns::Bool(),
          "State whether particles should be cleared on insertion.");

        prm.declare_entry(
          "removal box points coordinates",
          "0. , 0. , 0. : 1. , 1. , 1.",
          Patterns::List(
            Patterns::List(Patterns::Double(), 2, 3, ","), 2, 2, ":"),
          "Coordinates of two points for the removal box (x1, y1, z1 : x2, y2, z2)");

        // File:
        prm.declare_entry("list of input files",
                          "particles.input",
                          Patterns::List(Patterns::FileName()),
                          "The file name from which we load the particles");

        // Plane:
        prm.declare_entry("insertion plane point",
                          "0., 0., 0.",
                          Patterns::List(Patterns::Double()),
                          "Insertion plane point location");
        prm.declare_entry("insertion plane normal vector",
                          "1., 0., 0.",
                          Patterns::List(Patterns::Double()),
                          "Insertion plane normal vector");
        prm.declare_entry(
          "insertion plane threshold distance",
          "0.",
          Patterns::Double(),
          "If all the vertices of a cell are closer or equal to this value, than this cell is in the plane");

        // List:
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
        prm.declare_entry("list velocity x",
                          "0",
                          Patterns::List(Patterns::Double()),
                          "List of initial velocities x");
        prm.declare_entry("list velocity y",
                          "0",
                          Patterns::List(Patterns::Double()),
                          "List of initial velocities y");
        prm.declare_entry("list velocity z",
                          "0",
                          Patterns::List(Patterns::Double()),
                          "List of initial velocities z");
        prm.declare_entry("list omega x",
                          "0.",
                          Patterns::List(Patterns::Double()),
                          "List of initial omega x");
        prm.declare_entry("list omega y",
                          "0.",
                          Patterns::List(Patterns::Double()),
                          "List of initial omega y");
        prm.declare_entry("list omega z",
                          "0.",
                          Patterns::List(Patterns::Double()),
                          "List of initial omega z");
        prm.declare_entry("list diameters",
                          "-1.0",
                          Patterns::List(Patterns::Double()),
                          "List of diameters");
        prm.declare_entry("list temperatures",
                          "0.",
                          Patterns::List(Patterns::Double()),
                          "List of initial temperatures");
        // Volume:
        prm.declare_entry(
          "insertion direction sequence",
          "0,1,2",
          Patterns::List(Patterns::Integer(), 2, 3),
          "Direction of particle insertion for the volume insertion method");
        prm.declare_entry(
          "insertion box points coordinates",
          "0. , 0. , 0. : 1. , 1. , 1.",
          Patterns::List(
            Patterns::List(Patterns::Double(), 2, 3, ","), 2, 2, ":"),
          "Coordinates of two points for the insertion box (x1, y1, z1 : x2, y2, z2)");
        prm.declare_entry("insertion distance threshold",
                          "1.",
                          Patterns::Double(),
                          "Distance threshold");

        // Volume or plane:
        prm.declare_entry(
          "insertion maximum offset",
          "1.",
          Patterns::Double(),
          "Maximum position offset when insertion of particles");
        prm.declare_entry(
          "insertion prn seed",
          "1",
          Patterns::Integer(),
          "Pseudo-random number seed used to generate the position offsets");
        prm.declare_entry("initial velocity",
                          "0.0, 0.0, 0.0",
                          Patterns::List(Patterns::Double()),
                          "Initial velocity (x, y, z)");
        prm.declare_entry("initial angular velocity",
                          "0.0, 0.0, 0.0",
                          Patterns::List(Patterns::Double()),
                          "Initial angular velocity (x, y, z)");
        auto initial_temperature_function_parsed =
          std::make_shared<Functions::ParsedFunction<dim>>(1);
        prm.enter_subsection("initial temperature function");
        {
          initial_temperature_function_parsed->declare_parameters(prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    InsertionInfo<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        const std::string insertion = prm.get("insertion method");
        if (insertion == "file")
          insertion_method = InsertionMethod::file;
        else if (insertion == "plane")
          insertion_method = InsertionMethod::plane;
        else if (insertion == "list")
          insertion_method = InsertionMethod::list;
        else if (insertion == "volume")
          insertion_method = InsertionMethod::volume;
        else
          {
            throw(std::runtime_error("Invalid insertion method "));
          }
        inserted_this_step =
          prm.get_integer("inserted number of particles at each time step");
        insertion_frequency = prm.get_integer("insertion frequency");

        // Clear:
        removing_particles_in_region = prm.get_bool("remove particles");

        const std::vector<std::string> removal_box_point_coordinates_list(
          Utilities::split_string_list(
            prm.get("removal box points coordinates"), ":"));

        std::vector<double> removal_point_coord_temp_1 =
          Utilities::string_to_double(Utilities::split_string_list(
            removal_box_point_coordinates_list.at(0)));
        std::vector<double> removal_point_coord_temp_2 =
          Utilities::string_to_double(Utilities::split_string_list(
            removal_box_point_coordinates_list.at(1)));

        if (removal_point_coord_temp_1.size() == 2 &&
            removal_point_coord_temp_2.size() == 2)
          {
            removal_point_coord_temp_1.push_back(0.);
            removal_point_coord_temp_2.push_back(0.);
          }

        for (unsigned int i = 0; i < 3; ++i)
          {
            clear_box_point_1[i] = removal_point_coord_temp_1.at(i);
            clear_box_point_2[i] = removal_point_coord_temp_2.at(i);
          }

        // File:
        // File for the insertion
        list_of_input_files =
          convert_string_to_vector<std::string>(prm, "list of input files");

        // Plane:
        // Insertion plane normal vector
        insertion_plane_normal_vector =
          value_string_to_tensor<3>(prm.get("insertion plane normal vector"));

        // Insertion plane point
        insertion_plane_point =
          value_string_to_tensor<3>(prm.get("insertion plane point"));

        // List:
        // Read x, y and z lists
        list_x = convert_string_to_vector<double>(prm, "list x");
        list_y = convert_string_to_vector<double>(prm, "list y");
        list_z = convert_string_to_vector<double>(prm, "list z");

        // Find which vector is the longest
        int max_size = std::max({list_x.size(), list_y.size(), list_z.size()});

        // Resize the vectors in the case that one was longer
        list_x.resize(max_size);
        list_y.resize(max_size);
        list_z.resize(max_size);

        // Read vx, vv and vz
        list_vx = convert_string_to_vector<double>(prm, "list velocity x");
        list_vy = convert_string_to_vector<double>(prm, "list velocity y");
        list_vz = convert_string_to_vector<double>(prm, "list velocity z");

        // Fill the velocity vectors with zeros to match the size of list_x
        if (list_vx != list_x)
          list_vx.resize(max_size);
        if (list_vy != list_x)
          list_vy.resize(max_size);
        if (list_vz != list_x)
          list_vz.resize(max_size);

        // Read wx, wy and wz
        list_wx = convert_string_to_vector<double>(prm, "list omega x");
        list_wy = convert_string_to_vector<double>(prm, "list omega y");
        list_wz = convert_string_to_vector<double>(prm, "list omega z");

        // Fill the angular velocity vectors with zeros to match the size of
        // list_x
        if (list_wx != list_x)
          list_wx.resize(max_size);
        if (list_wy != list_x)
          list_wy.resize(max_size);
        if (list_wz != list_x)
          list_wz.resize(max_size);

        // Read the diameters list
        list_d = convert_string_to_vector<double>(prm, "list diameters");

        // Read the temperatures list
        list_T = convert_string_to_vector<double>(prm, "list temperatures");

        // Volume:
        std::vector<int> axis_order =
          convert_string_to_vector<int>(prm, "insertion direction sequence");
        if (axis_order.size() == 2)
          axis_order.push_back(0.);

        direction_sequence.reserve(3);
        direction_sequence.push_back(axis_order[0]);
        direction_sequence.push_back(axis_order[1]);
        direction_sequence.push_back(axis_order[2]);

        const std::vector<std::string> point_coordinates_list(
          Utilities::split_string_list(
            prm.get("insertion box points coordinates"), ":"));

        std::vector<double> point_coord_temp_1 = Utilities::string_to_double(
          Utilities::split_string_list(point_coordinates_list.at(0)));
        std::vector<double> point_coord_temp_2 = Utilities::string_to_double(
          Utilities::split_string_list(point_coordinates_list.at(1)));

        if (point_coord_temp_1.size() == 2 && point_coord_temp_2.size() == 2)
          {
            point_coord_temp_1.push_back(0.);
            point_coord_temp_2.push_back(0.);
          }
        for (int i = 0; i < 3; ++i)
          {
            insertion_box_point_1[i] = point_coord_temp_1.at(i);
            insertion_box_point_2[i] = point_coord_temp_2.at(i);
          }

        distance_threshold = prm.get_double("insertion distance threshold");
        insertion_maximum_offset = prm.get_double("insertion maximum offset");
        seed_for_insertion       = prm.get_integer("insertion prn seed");

        initial_vel = value_string_to_tensor<3>(prm.get("initial velocity"));
        initial_omega =
          value_string_to_tensor<3>(prm.get("initial angular velocity"));

        auto initial_temperature_function_parsed =
          std::make_shared<Functions::ParsedFunction<dim>>(1);
        prm.enter_subsection("initial temperature function");
        {
          initial_temperature_function_parsed->parse_parameters(prm);
          initial_temperature_function = initial_temperature_function_parsed;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ModelParameters<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        prm.enter_subsection("load balancing");
        {
          prm.declare_entry(
            "load balance method",
            "none",
            Patterns::Selection(
              "none|once|frequent|dynamic|dynamic_with_sparse_contacts"),
            "Choosing load-balance method"
            "Choices are <none|once|frequent|dynamic|dynamic_with_sparse_contacts>.");

          prm.declare_entry(
            "step",
            "100000",
            Patterns::Integer(),
            "Step at which the triangulation is repartitioned "
            "and load is balanced for single-step load-balancing");

          prm.declare_entry(
            "frequency",
            "100000",
            Patterns::Integer(),
            "Frequency at which the triangulation is repartitioned "
            "and load is balanced for frequent load-balancing");

          prm.declare_entry("threshold",
                            "0.5",
                            Patterns::Double(),
                            "Threshold for dynamic load-balancing");

          prm.declare_entry("dynamic check frequency",
                            "10000",
                            Patterns::Integer(),
                            "Checking frequency for dynamic load-balancing");


          auto cell_weight_function_parsed =
            std::make_shared<Functions::ParsedFunction<dim>>(1);
          prm.enter_subsection("cell weight function");
          {
#if DEAL_II_VERSION_GTE(9, 7, 0)
            cell_weight_function_parsed->declare_parameters(prm, 1, "1000");

#else
            cell_weight_function_parsed->declare_parameters(prm, 1);

#endif
          }
          prm.leave_subsection();

          prm.declare_entry(
            "particle weight",
            "2000",
            Patterns::Integer(),
            "The particle weight based on a default cell weight of 1000");

          prm.declare_entry(
            "active weight factor",
            "1.0",
            Patterns::Double(),
            "Factor applied on the particle weight in load balancing if the cell is active");

          prm.declare_entry(
            "inactive weight factor",
            "1.0",
            Patterns::Double(),
            "Factor applied on the particle weight in load balancing if the cell is inactive");
        }
        prm.leave_subsection();

        prm.enter_subsection("contact detection");
        {
          prm.declare_entry("contact detection method",
                            "dynamic",
                            Patterns::Selection("constant|dynamic"),
                            "Choosing contact detection method"
                            "Choices are <constant|dynamic>.");

          prm.declare_entry("frequency",
                            "1",
                            Patterns::Integer(),
                            "Particle-particle contact list");

          prm.declare_entry(
            "dynamic contact search size coefficient",
            "0.8",
            Patterns::Double(),
            "Security coefficient for dynamic contact detection");

          prm.declare_entry(
            "neighborhood threshold",
            "1.3",
            Patterns::Double(),
            "Contact search zone diameter to particle diameter ratio");
        }
        prm.leave_subsection();

        prm.declare_entry(
          "particle particle contact force method",
          "hertz_mindlin_limit_overlap",
          Patterns::Selection(
            "linear|hertz_mindlin_limit_force|hertz_mindlin_limit_overlap|hertz|hertz_JKR|DMT"),
          "Choosing particle-particle contact force model"
          "Choices are <linear|hertz_mindlin_limit_force|hertz_mindlin_limit_overlap|hertz|hertz_JKR|DMT>.");

        prm.declare_entry("particle wall contact force method",
                          "nonlinear",
                          Patterns::Selection("linear|nonlinear|JKR|DMT"),
                          "Choosing particle-wall contact force model"
                          "Choices are <linear|nonlinear|JKR|DMT>.");

        prm.declare_entry(
          "dmt cut-off threshold",
          "0.1",
          Patterns::Double(),
          "Cut-off threshold above which the Van der Waal forces are "
          "ignored for the DMT model relative to the pull-off force");

        prm.declare_entry(
          "rolling resistance torque method",
          "constant",
          Patterns::Selection(
            "none|no_resistance|constant|constant_resistance|viscous|viscous_resistance|epsd|epsd_resistance"),
          "Choosing rolling resistance torque model"
          "Choices are <no_resistance|constant_resistance|viscous_resistance|epsd_resistance>.");

        prm.declare_entry(
          "f coefficient",
          "0.",
          Patterns::Double(),
          "Model parameter for the EPSD rolling resistance model.");

        prm.declare_entry("integration method",
                          "velocity_verlet",
                          Patterns::Selection("velocity_verlet|explicit_euler"),
                          "Choosing integration method"
                          "Choices are <velocity_verlet|explicit_euler>.");

        prm.declare_entry("solver type",
                          "dem",
                          Patterns::Selection("dem|cfd_dem|dem_mp"),
                          "Choosing solver type"
                          "Choices are <dem|cfd_dem|dem_mp>.");

        prm.enter_subsection("adaptive sparse contacts");
        {
          prm.declare_entry(
            "enable adaptive sparse contacts",
            "false",
            Patterns::Selection("true|false"),
            "Enable the dynamic search for sparse particle contacts"
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

        prm.declare_entry("disable position integration",
                          "false",
                          Patterns::Selection("true|false"),
                          "Disable the integration of position and velocity"
                          "Choices are <true|false>.");
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ModelParameters<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("model parameters");
      {
        prm.enter_subsection("adaptive sparse contacts");
        {
          sparse_particle_contacts =
            prm.get_bool("enable adaptive sparse contacts");
          advect_particles = prm.get_bool("enable particle advection");

          // Thresholds for disabling contacts
          granular_temperature_threshold =
            prm.get_double("granular temperature threshold");
          solid_fraction_threshold = prm.get_double("solid fraction threshold");
        }
        prm.leave_subsection();

        prm.enter_subsection("load balancing");
        {
          const std::string load_balance = prm.get("load balance method");

          if (load_balance == "once")
            {
              load_balance_method = LoadBalanceMethod::once;
              load_balance_step   = prm.get_integer("step");
            }
          else if (load_balance == "frequent")
            {
              load_balance_method    = LoadBalanceMethod::frequent;
              load_balance_frequency = prm.get_integer("frequency");
            }
          else if (load_balance == "dynamic")
            {
              load_balance_method    = LoadBalanceMethod::dynamic;
              load_balance_threshold = prm.get_double("threshold");
              dynamic_load_balance_check_frequency =
                prm.get_integer("dynamic check frequency");
            }
          else if (load_balance == "dynamic_with_sparse_contacts")
            {
              // Check if adaptive sparse contacts is enabled, otherwise
              // throw an error message indicating that the user should use
              // dynamic load balancing instead or enable adaptive sparse
              // contacts
              if (sparse_particle_contacts)
                {
                  load_balance_method =
                    LoadBalanceMethod::dynamic_with_sparse_contacts;
                  load_balance_threshold = prm.get_double("threshold");
                  dynamic_load_balance_check_frequency =
                    prm.get_integer("dynamic check frequency");

                  // Weights for load balancing of active and inactive cells
                  active_load_balancing_factor =
                    prm.get_double("active weight factor");
                  inactive_load_balancing_factor =
                    prm.get_double("inactive weight factor");
                }
              else
                {
                  throw(std::runtime_error(
                    "Invalid contact detection method: adaptive sparse contacts is not enabled "
                    "while dynamic_with_sparse_contacts is selected, use dynamic instead"));
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
          auto cell_weight_function_parsed =
            std::make_shared<Functions::ParsedFunction<dim>>(1);

          prm.enter_subsection("cell weight function");
          {
            cell_weight_function_parsed->parse_parameters(prm);

            cell_weight_function = cell_weight_function_parsed;
          }
          prm.leave_subsection();
          load_balance_particle_weight = prm.get_integer("particle weight");
        }
        prm.leave_subsection();

        prm.enter_subsection("contact detection");
        {
          contact_detection_frequency = prm.get_integer("frequency");
          dynamic_contact_search_factor =
            prm.get_double("dynamic contact search size coefficient");
          neighborhood_threshold = prm.get_double("neighborhood threshold");

          const std::string contact_search =
            prm.get("contact detection method");

          if (contact_search == "constant")
            contact_detection_method = ContactDetectionMethod::constant;
          else if (contact_search == "dynamic")
            contact_detection_method = ContactDetectionMethod::dynamic;
          else
            throw(std::runtime_error("Invalid contact detection method "));
        }
        prm.leave_subsection();

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
        else if (ppcf == "hertz_JKR")
          particle_particle_contact_force_model =
            ParticleParticleContactForceModel::hertz_JKR;
        else if (ppcf == "DMT")
          particle_particle_contact_force_model =
            ParticleParticleContactForceModel::DMT;
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
        else if (pwcf == "JKR")
          particle_wall_contact_force_method =
            ParticleWallContactForceModel::JKR;
        else if (pwcf == "DMT")
          particle_wall_contact_force_method =
            ParticleWallContactForceModel::DMT;
        else
          {
            throw(
              std::runtime_error("Invalid particle-wall contact force model "));
          }

        dmt_cut_off_threshold = prm.get_double("dmt cut-off threshold");

        const std::string rolling_resistance_torque =
          prm.get("rolling resistance torque method");

        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          {
            if (rolling_resistance_torque == "no_resistance")
              std::cout
                << "Warning, the \"no_resistance\" entry to the \"rolling "
                   "resistance torque method\" parameter will be deprecated. "
                   "Please use \"none\" instead."
                << std::endl;

            if (rolling_resistance_torque == "constant_resistance")
              std::cout
                << "Warning, the \"constant_resistance\" entry to the \"rolling"
                   " resistance torque method\" parameter will be deprecated."
                   " Please use \"constant\" instead."
                << std::endl;

            if (rolling_resistance_torque == "viscous_resistance")
              std::cout
                << "Warning, the \"viscous_resistance\" entry to the \"rolling"
                   " resistance torque method\" parameter will be deprecated. "
                   "Please use \"viscous\" instead."
                << std::endl;
            if (rolling_resistance_torque == "epsd_resistance")
              std::cout
                << "Warning, the \"epsd_resistance\" entry to the \"rolling "
                   "resistance torque method\" parameter will be deprecated. "
                   "Please use \"epsd\" instead."
                << std::endl;
          }
        if (rolling_resistance_torque == "no_resistance" ||
            rolling_resistance_torque == "none")
          {
            rolling_resistance_method = RollingResistanceMethod::none;
          }
        else if (rolling_resistance_torque == "constant_resistance" ||
                 rolling_resistance_torque == "constant")
          {
            rolling_resistance_method = RollingResistanceMethod::constant;
          }
        else if (rolling_resistance_torque == "viscous_resistance" ||
                 rolling_resistance_torque == "viscous")
          {
            rolling_resistance_method = RollingResistanceMethod::viscous;
          }
        else if (rolling_resistance_torque == "epsd_resistance" ||
                 rolling_resistance_torque == "epsd")
          {
            rolling_resistance_method = RollingResistanceMethod::epsd;
          }


        // Model parameter for the EPSD rolling resistance model
        f_coefficient_epsd = prm.get_double("f coefficient");

        const std::string integration = prm.get("integration method");
        if (integration == "velocity_verlet")
          integration_method = IntegrationMethod::velocity_verlet;
        else if (integration == "explicit_euler")
          integration_method = IntegrationMethod::explicit_euler;
        else
          {
            throw(std::runtime_error("Invalid integration method "));
          }

        const std::string solver_type_str = prm.get("solver type");
        if (solver_type_str == "dem")
          solver_type = DEM::SolverType::dem;
        else if (solver_type_str == "cfd_dem")
          solver_type = DEM::SolverType::cfd_dem;
        else if (solver_type_str == "dem_mp")
          solver_type = DEM::SolverType::dem_mp;
        else
          {
            throw(std::runtime_error("Invalid solver type"));
          }

        disable_position_integration =
          prm.get_bool("disable position integration");
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
    FloatingWalls<dim>::declare_parameters(ParameterHandler &prm) const
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

    template <int dim>
    void
    FloatingWalls<dim>::declareDefaultEntry(ParameterHandler &prm)
    {
      prm.declare_entry("point on wall",
                        "0., 0., 0.",
                        Patterns::List(Patterns::Double(), 2, 3),
                        "Point on wall");
      prm.declare_entry("normal vector",
                        "1., 0., 0.",
                        Patterns::List(Patterns::Double(), 2, 3),
                        "Point on wall");

      prm.enter_subsection("point on wall");
      prm.declare_entry("x", "0.", Patterns::Double(), "X Point on wall");
      prm.declare_entry("y", "0.", Patterns::Double(), "Y Point on wall");
      prm.declare_entry("z", "0.", Patterns::Double(), "Z Point on wall");
      prm.mark_as_deprecated("x");
      prm.mark_as_deprecated("y");
      prm.mark_as_deprecated("z");
      prm.leave_subsection();

      prm.enter_subsection("normal vector");
      prm.declare_entry("nx", "0.", Patterns::Double(), "X Normal vector wall");
      prm.declare_entry("ny", "0.", Patterns::Double(), "Y Normal vector wall");
      prm.declare_entry("nz", "0.", Patterns::Double(), "Z Normal vector wall");
      prm.mark_as_deprecated("nx");
      prm.mark_as_deprecated("ny");
      prm.mark_as_deprecated("nz");
      prm.leave_subsection();

      prm.declare_entry("start time", "0.", Patterns::Double(), "Start time");

      prm.declare_entry("end time", "0.", Patterns::Double(), "End time");
    }

    template <int dim>
    void
    FloatingWalls<dim>::parse_floating_wall(ParameterHandler &prm)
    {
      // Position
      Point<dim> wall_point(
        value_string_to_tensor<dim>(prm.get("point on wall")));
      this->points_on_walls.push_back(wall_point);

      // Normal
      Tensor<1, dim> wall_normal(
        value_string_to_tensor<dim>(prm.get("normal vector")));
      wall_normal = wall_normal / wall_normal.norm();
      this->floating_walls_normal_vectors.push_back(wall_normal);

      // Time
      time_start.push_back(prm.get_double("start time"));
      time_end.push_back(prm.get_double("end time"));
    }

    void
    BCDEM::declare_parameters(ParameterHandler &prm) const
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
                            point_on_rotation_axis,
                            outlet_boundaries,
                            bc_types);

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

      prm.declare_entry("rotational vector",
                        "1.,0.,0.",
                        Patterns::List(Patterns::Double()),
                        "Rotational vector elements");

      prm.declare_entry("periodic id 0",
                        "0",
                        Patterns::Integer(),
                        "Periodic boundary ID 0");
      prm.declare_entry("periodic id 1",
                        "0",
                        Patterns::Integer(),
                        "Periodic boundary ID 1");

      prm.declare_entry("point on rotational vector",
                        "0, 0, 0",
                        Patterns::List(Patterns::Double()),
                        "Point on the rotational vector");
      prm.declare_entry(
        "periodic direction",
        "0",
        Patterns::Integer(),
        "Periodic direction or normal direction of periodic boundary");
    }

    void
    BCDEM::parse_boundary_conditions(const ParameterHandler &prm)
    {
      const unsigned int boundary_id   = prm.get_integer("boundary id");
      const std::string  boundary_type = prm.get("type");

      if (boundary_type == "outlet")
        {
          bc_types.push_back(BoundaryType::outlet);
          this->outlet_boundaries.push_back(boundary_id);
        }
      else if (boundary_type == "translational")
        {
          bc_types.push_back(BoundaryType::translational);
          Tensor<1, 3> translational_velocity;
          translational_velocity[0] = prm.get_double("speed x");
          translational_velocity[1] = prm.get_double("speed y");
          translational_velocity[2] = prm.get_double("speed z");

          this->boundary_translational_velocity.at(boundary_id) =
            translational_velocity;
        }
      else if (boundary_type == "rotational")
        {
          bc_types.push_back(BoundaryType::rotational);
          double rotational_speed = prm.get_double("rotational speed");

          // Read the rotational vector from a list of doubles
          Tensor<1, 3> rotational_vector =
            value_string_to_tensor<3>(prm.get("rotational vector"));
          if (rotational_vector.norm() == 0.)
            {
              throw(std::runtime_error(
                "Invalid rotational vector. Its norm cannot be equal to zero."));
            }
          // Read the point from a list of doubles
          Tensor<1, 3> point_on_rotation_axis_tensor =
            value_string_to_tensor<3>(prm.get("point on rotational vector"));

          this->boundary_rotational_speed.at(boundary_id) = rotational_speed;
          this->boundary_rotational_vector.at(boundary_id) =
            rotational_vector / rotational_vector.norm();
          this->point_on_rotation_axis.at(boundary_id) =
            point_on_rotation_axis_tensor;
        }
      else if (boundary_type == "fixed_wall")
        {
          bc_types.push_back(BoundaryType::fixed_wall);
        }
      else if (boundary_type == "periodic")
        {
          bc_types.push_back(BoundaryType::periodic);
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
    BCDEM::initialize_containers(
      std::unordered_map<unsigned int, Tensor<1, 3>> &boundary_trans_velocity,
      std::unordered_map<unsigned int, double>       &boundary_rot_speed,
      std::unordered_map<unsigned int, Tensor<1, 3>> &boundary_rot_vector,
      std::unordered_map<unsigned int, Point<3>>     &point_on_rot_axis,
      std::vector<unsigned int>                      &outlet_boundaries_id,
      std::vector<BoundaryType>                      &boundaries_types) const
    {
      Tensor<1, 3> zero_tensor({0.0, 0.0, 0.0});

      for (unsigned int counter = 0; counter < DEM_BC_number_max; ++counter)
        {
          boundary_trans_velocity.insert({counter, zero_tensor});
          boundary_rot_speed.insert({counter, 0});
          boundary_rot_vector.insert({counter, zero_tensor});
          point_on_rot_axis.insert({counter, Point<3>(zero_tensor)});
        }

      // NOTE This first vector should not be initialized this big.
      outlet_boundaries_id.reserve(DEM_BC_number);
      boundaries_types.reserve(DEM_BC_number);
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
          "lagrangian post-processing",
          "false",
          Patterns::Bool(),
          "State whether lagrangian post-processing should be performed.");
        prm.declare_entry(
          "force chains",
          "false",
          Patterns::Bool(),
          "State whether force chains visualization should be performed.");
        prm.enter_subsection("particle wall collision statistics");
        {
          prm.declare_entry(
            "enable particle wall collision statistics",
            "false",
            Patterns::Selection("true|false"),
            "Enable the logging of particle-wall collision statistics"
            "Choices are <true|false>.");
          prm.declare_entry(
            "log collisions with all walls",
            "true",
            Patterns::Selection("true|false"),
            "State whether collisions with all walls should be logged"
            "Choices are <true|false>.");
          prm.declare_entry("wall boundary ids",
                            "0",
                            Patterns::List(Patterns::Integer()),
                            "Boundary ids of the walls to log collisions with");
          prm.declare_entry(
            "collision statistics file",
            "collision_statistics.csv",
            Patterns::FileName(),
            "Exported particle-wall collision results filename");

          prm.declare_entry(
            "verbosity",
            "quiet",
            Patterns::Selection("quiet|verbose"),
            "State whether collision starts and ends should be printed. "
            "Choices are <quiet|verbose>.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    void
    LagrangianPostProcessing::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("post-processing");
      {
        lagrangian_post_processing_enabled =
          prm.get_bool("lagrangian post-processing");
        force_chains = prm.get_bool("force chains");
        prm.enter_subsection("particle wall collision statistics");
        {
          particle_wall_collision_statistics =
            prm.get_bool("enable particle wall collision statistics");
          log_collisions_with_all_walls =
            prm.get_bool("log collisions with all walls");
          particle_wall_collision_boundary_ids =
            convert_string_to_vector<int>(prm, "wall boundary ids");
          collision_stats_file_name = prm.get("collision statistics file");
          const std::string verbose = prm.get("verbosity");
          if (verbose == "quiet")
            collision_verbosity = Parameters::Verbosity::quiet;
          else if (verbose == "verbose")
            collision_verbosity = Parameters::Verbosity::verbose;
          else
            {
              throw(std::runtime_error("Invalid verbosity choice "));
            }
        }
        prm.leave_subsection();
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


    template <int dim>
    void
    ParticleRayTracing<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("particle ray tracing");
      {
        // Location of the first photon to be inserted
        prm.declare_entry("starting photon insertion position",
                          "0.,0.,0.",
                          Patterns::List(Patterns::Double()),
                          "Location of the first photon being inserted.");

        // In which direction will the photons be inserted  relative to the
        // first photon.
        prm.declare_entry(
          "insertion unit tensors",
          "1.,0.,0. : 0., 1., 0. : 0., 0., 1.",
          Patterns::List(
            Patterns::List(Patterns::Double(), 3, 3, ","), 3, 3, ":"),
          "Directions used to insert photons.");

        // How many photon will be inserted in each of those directions.
        prm.declare_entry("number of inserted photons per direction",
                          "1 : 1 : 1",
                          Patterns::List(Patterns::Integer(), 3, 3, ":"),
                          "Number of inserted photon in each direction.");

        // What is the distance between each photon in each of those directions
        // considering an offset equal to 0.
        prm.declare_entry("distance between photons on insertion per direction",
                          "1. :  1. : 1.",
                          Patterns::List(Patterns::Double(), 3, 3, ":"),
                          "Number of inserted photon in each direction.");

        // In which direction photons will move considering a photon maximum
        // angle offset equal to 0.
        prm.declare_entry("reference displacement vector",
                          "0.,0.,1.",
                          Patterns::List(Patterns::Double()),
                          "Reference displacement vector of each photons.");

        // Insertion location
        prm.declare_entry(
          "photon insertion maximum offset",
          "0.",
          Patterns::Double(),
          "Set the maximum offset applied on each photon position during their "
          "insertion. If set to 0., photons will be perfectly aligned."
          "respectively to the insertion unit tensors ");

        prm.declare_entry(
          "photon insertion prn seed",
          "0",
          Patterns::Integer(),
          "Pseudo random seed used to generate the offset for each photon "
          "insertion location.");

        // Displacement unit vector
        prm.declare_entry(
          "photon maximum angular offset",
          "0.",
          Patterns::Double(),
          "Used to introduce randomness in the displacement direction of each "
          "photon. This parameter defines the maximum angle between a given "
          "photon displacement vector and the prescribed displacement vector "
          "parameter. If set to zero, every photon will move in the same "
          "direction defined by the prescribed displacement vector parameter."
          "Otherwise, the offset is applied in a random orientation relative to "
          "the reference displacement vector parameter.");

        prm.declare_entry("photon angular offset prn seed",
                          "1",
                          Patterns::Integer(),
                          "Pseudo random seed used to generated the angle "
                          "offset and the random orientation.");
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ParticleRayTracing<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("particle ray tracing");
      {
        // Location of the first photon to be inserted
        starting_point = point_nd_to_3d(Point<dim>(value_string_to_tensor<dim>(
          prm.get("starting photon insertion position"))));

        // Extract the strings related to the three std::vector
        const std::vector<std::string> insertion_unit_tensors_string(
          Utilities::split_string_list(prm.get("insertion unit tensors"), ":"));
        const std::vector<std::string> n_photons_per_directions_strings(
          Utilities::split_string_list(
            prm.get("number of inserted photons per direction"), ":"));
        const std::vector<std::string>
          step_between_photon_per_direction_strings(
            Utilities::split_string_list(
              prm.get("distance between photons on insertion per direction"),
              ":"));

        // We always use 3d tensor even in a dim=2 simulation. Still, we want to
        // make sure the user write 2d tensor in the prm when the simulation is
        // in 2d. Those tensors will be but in 3d afterward.
        AssertThrow(
          insertion_unit_tensors_string.size() == dim &&
            step_between_photon_per_direction_strings.size() == dim &&
            step_between_photon_per_direction_strings.size() == dim,
          dealii::ExcMessage(
            "The \"insertion unit tensors\", \"distance between photons on"
            " insertion per direction\" and \"number of inserted photons per "
            "direction\" all need to have a number of dimension equal to the "
            "\"dimension\" parameter."));

        // Always of size 3.
        insertion_directions_units_vector.reserve(3);
        n_photons_each_directions.reserve(3);
        step_between_photons_each_directions.reserve(3);
        for (int i = 0; i < dim; i++)
          {
            Tensor<1, dim> temp_direction_tensor =
              value_string_to_tensor<dim>(insertion_unit_tensors_string.at(i));
            temp_direction_tensor =
              temp_direction_tensor / temp_direction_tensor.norm();
            insertion_directions_units_vector.emplace_back(
              tensor_nd_to_3d(temp_direction_tensor));

            n_photons_each_directions.emplace_back(Utilities::string_to_double(
              n_photons_per_directions_strings.at(i)));

            step_between_photons_each_directions.emplace_back(
              Utilities::string_to_double(
                step_between_photon_per_direction_strings.at(i)));
          }
        // Reference displacement unit tensor
        ref_displacement_tensor_unit =
          tensor_nd_to_3d(value_string_to_tensor<dim>(
            prm.get("reference displacement vector")));
        ref_displacement_tensor_unit =
          ref_displacement_tensor_unit / ref_displacement_tensor_unit.norm();

        // Insertion offset
        max_insertion_offset =
          prm.get_double("photon insertion maximum offset");
        prn_seed_photon_insertion =
          prm.get_integer("photon insertion prn seed");

        // Displacement direction offset
        max_angular_offset = prm.get_double("photon maximum angular offset");
        prn_seed_photon_displacement =
          prm.get_integer("photon insertion prn seed");
      }
      prm.leave_subsection();
    }

    template class ForceTorqueOnWall<2>;
    template class ForceTorqueOnWall<3>;
    template class FloatingWalls<2>;
    template class FloatingWalls<3>;
    template class ModelParameters<2>;
    template class ModelParameters<3>;
    template class FloatingGrid<2>;
    template class FloatingGrid<3>;
    template class GridMotion<2>;
    template class GridMotion<3>;
    template class InsertionInfo<2>;
    template class InsertionInfo<3>;
    template struct ParticleRayTracing<2>;
    template struct ParticleRayTracing<3>;

  } // namespace Lagrangian
} // namespace Parameters
