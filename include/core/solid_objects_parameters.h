// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_solid_objects_parameters_h
#define lethe_solid_objects_parameters_h

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  enum ThermalBoundaryType
  {
    adiabatic,
    isothermal
  };

  template <int dim>
  class NitscheObject
  {
  public:
    NitscheObject()
      : solid_velocity(dim)
    {}

    void
    declare_parameters(ParameterHandler &prm, unsigned int id);
    void
    parse_parameters(ParameterHandler &prm, unsigned int id);

    // Solid mesh
    Parameters::Mesh solid_mesh;
    unsigned int     number_quadrature_points;

    // Solid velocity
    Functions::ParsedFunction<dim> solid_velocity;
    Functions::ParsedFunction<dim> solid_temperature;
    bool                           enable_particles_motion;
    bool                           enable_heat_bc;

    // Penalization term
    double beta;
    double beta_heat;

    // Particle motion integration parameters
    unsigned int particles_sub_iterations;
    bool         stop_particles_lost;

    // Information for force calculation
    Point<dim>
      center_of_rotation; // Center of rotation used for torque calculation

    bool        calculate_force_on_solid;
    bool        calculate_torque_on_solid;
    Point<dim>  cor; // Center of rotation used for torque calculation
    std::string force_output_name;
    std::string torque_output_name;
  };

  template <int dim>
  void
  NitscheObject<dim>::declare_parameters(ParameterHandler &prm, unsigned int id)
  {
    prm.enter_subsection("nitsche solid " + Utilities::int_to_string(id, 1));
    {
      solid_mesh.declare_parameters(prm);

      prm.enter_subsection("solid velocity");
      solid_velocity.declare_parameters(prm, dim);
      prm.leave_subsection();



      prm.enter_subsection("solid temperature");
      solid_temperature.declare_parameters(prm);
      prm.leave_subsection();

      prm.declare_entry("enable particles motion",
                        "false",
                        Patterns::Bool(),
                        "Condition on the motion of particles");
      prm.declare_entry(
        "enable heat boundary condition",
        "false",
        Patterns::Bool(),
        "controls if a heat boundary condition is imposed on the Nitsche immersed boundary");

      prm.declare_entry("beta",
                        "10",
                        Patterns::Double(),
                        "Penalization term for Nitsche method");
      prm.declare_entry(
        "beta heat",
        "10",
        Patterns::Double(),
        "Penalization term for Nitsche method applied to the heat equation");

      prm.declare_entry(
        "particles sub iterations",
        "1",
        Patterns::Integer(),
        "Number of sub iterations for the motion of the particles. This parameter"
        "enables the uses of a higher CFL condition for the Nitsche solver while preventing the loss of particles");

      prm.declare_entry(
        "stop if particles lost",
        "true",
        Patterns::Bool(),
        "Enable stopping the simulation if Nitsche particles have been lost");

      prm.declare_entry(
        "number quadrature points",
        "2",
        Patterns::Integer(),
        "Number of Nitsche (quadrature) points to insert in a 1D cell. The number of"
        "inserted points will be higher for higher dimensions. Increasing this"
        "number will lead to a higher points density inside the solid.");

      if constexpr (dim == 2)
        {
          prm.declare_entry("center of rotation",
                            "0., 0.",
                            Patterns::List(Patterns::Double()),
                            "Solid object center of rotation");
        }
      if constexpr (dim == 3)
        {
          prm.declare_entry("center of rotation",
                            "0., 0., 0.",
                            Patterns::List(Patterns::Double()),
                            "Solid object center of rotation");
        }

      prm.declare_entry("calculate force on solid",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of forces on solid");
      prm.declare_entry("calculate torque on solid",
                        "false",
                        Patterns::Bool(),
                        "Enable calculation of torques on solid");
      prm.declare_entry("solid force name",
                        "force_solid",
                        Patterns::FileName(),
                        "File output solid force prefix");
      prm.declare_entry("solid torque name",
                        "torque_solid",
                        Patterns::FileName(),
                        "File output solid torque prefix");
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  NitscheObject<dim>::parse_parameters(ParameterHandler &prm, unsigned int id)
  {
    prm.enter_subsection("nitsche solid " + Utilities::int_to_string(id, 1));
    {
      solid_mesh.parse_parameters(prm);
      prm.enter_subsection("solid velocity");
      solid_velocity.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("solid temperature");
      solid_temperature.parse_parameters(prm);
      prm.leave_subsection();

      enable_particles_motion  = prm.get_bool("enable particles motion");
      enable_heat_bc           = prm.get_bool("enable heat boundary condition");
      beta                     = prm.get_double("beta");
      beta_heat                = prm.get_double("beta heat");
      particles_sub_iterations = prm.get_integer("particles sub iterations");
      stop_particles_lost      = prm.get_bool("stop if particles lost");
      number_quadrature_points = prm.get_integer("number quadrature points");

      const std::vector<double> temp =
        convert_string_to_vector<double>(prm, "center of rotation");

      AssertThrow(temp.size() == dim,
                  ExcMessage("Invalid center of rotation. This should be a " +
                             Utilities::int_to_string(dim) +
                             " dimensional point."));

      for (unsigned int i = 0; i < dim; ++i)
        {
          center_of_rotation[i] = temp.at(i);
        }

      calculate_force_on_solid  = prm.get_bool("calculate force on solid");
      calculate_torque_on_solid = prm.get_bool("calculate torque on solid");
      force_output_name         = prm.get("solid force name");
      torque_output_name        = prm.get("solid torque name");
    }
    prm.leave_subsection();
  }

  template <int dim>
  class Nitsche
  {
  public:
    Nitsche()
    {}

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    // Calculate forces
    Verbosity verbosity;

    // Nitsche solid objects
    std::vector<std::shared_ptr<NitscheObject<dim>>> nitsche_solids;
    unsigned int                                     number_solids;
    static const unsigned int                        max_nitsche_solids = 10;
  };

  template <int dim>
  void
  Nitsche<dim>::declare_parameters(ParameterHandler &prm)
  {
    nitsche_solids.resize(max_nitsche_solids);

    prm.enter_subsection("nitsche");
    {
      prm.declare_entry(
        "verbosity",
        "quiet",
        Patterns::Selection("quiet|verbose"),
        "State whether the force on the solid should be printed "
        "Choices are <quiet|verbose>.");

      prm.declare_entry("number of solids",
                        "0",
                        Patterns::Integer(),
                        "Number of solid object");

      for (unsigned int i_solid = 0; i_solid < max_nitsche_solids; ++i_solid)
        {
          nitsche_solids[i_solid] = std::make_shared<NitscheObject<dim>>();
          nitsche_solids[i_solid]->declare_parameters(prm, i_solid);
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Nitsche<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("nitsche");
    {
      const std::string op = prm.get("verbosity");
      if (op == "verbose")
        verbosity = Verbosity::verbose;
      if (op == "quiet")
        verbosity = Verbosity::quiet;

      number_solids = prm.get_integer("number of solids");
      nitsche_solids.resize(number_solids);

      for (unsigned int i_solid = 0; i_solid < number_solids; ++i_solid)
        nitsche_solids[i_solid]->parse_parameters(prm, i_solid);
    }
    prm.leave_subsection();
  }


  template <int dim>
  class RigidSolidObject
  {
  public:
    RigidSolidObject()
    {}

    void
    declare_parameters(ParameterHandler &prm, unsigned int id);
    void
    parse_parameters(ParameterHandler &prm, unsigned int id);

    // Solid mesh
    Parameters::Mesh solid_mesh;

    // Output management
    bool output_bool;

    // Solid object velocity
    // Velocities are std::shared_ptr<Function<dim>> type, this way it is
    // possible to define a velocity using any type of function, like
    // ParsedFunction (which is used in the declared and parsing functions)
    // or ConstantFunctions (which is useful for unit tests).
    std::shared_ptr<Function<dim>> translational_velocity;
    std::shared_ptr<Function<dim>> angular_velocity;
    Point<dim>
      center_of_rotation; // Center of rotation used to locate the center of the
                          // object and also used to rotate the object
    // Temperature function for solid object
    std::shared_ptr<Function<dim>> solid_temperature;
    // Thermal boundary type for solid object
    ThermalBoundaryType thermal_boundary_type;
  };


  template <int dim>
  void
  RigidSolidObject<dim>::declare_parameters(ParameterHandler &prm,
                                            unsigned int      id)
  {
    // Use ParsedFunction<dim> during parameter declaration. We need to do this
    // to use the declare_parameters function, which is not possible with
    // std::make_shared<Function<dim>> type. Please refer to the comment before
    // the translational_velocity and angular_velocity attributes declaration.
    auto translational_velocity_parsed =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    auto angular_velocity_parsed =
      std::make_shared<Functions::ParsedFunction<dim>>(3);
    auto solid_temperature_parsed =
      std::make_shared<Functions::ParsedFunction<dim>>(1);

    prm.enter_subsection("solid object " + Utilities::int_to_string(id, 1));
    {
      solid_mesh.declare_parameters(prm);

      prm.enter_subsection("translational velocity");
      translational_velocity_parsed->declare_parameters(prm, dim);
      prm.leave_subsection();

      prm.enter_subsection("angular velocity");
      angular_velocity_parsed->declare_parameters(prm, 3);
      prm.leave_subsection();

      if constexpr (dim == 2)
        {
          prm.declare_entry("center of rotation",
                            "0., 0.",
                            Patterns::List(Patterns::Double()),
                            "Solid object center of rotation");
        }
      if constexpr (dim == 3)
        {
          prm.declare_entry("center of rotation",
                            "0., 0., 0.",
                            Patterns::List(Patterns::Double()),
                            "Solid object center of rotation");
        }

      prm.declare_entry("thermal boundary type",
                        "adiabatic",
                        Patterns::Selection("adiabatic|isothermal"),
                        "Choosing thermal boundary type"
                        "Choices are <adiabatic|isothermal>.");

      // Isothermal boundary
      prm.enter_subsection("temperature");
      solid_temperature_parsed->declare_parameters(prm, 1);
      prm.leave_subsection();

      prm.declare_entry("output solid object",
                        "true",
                        Patterns::Bool(),
                        "Controls the generation of output files");
    }
    prm.leave_subsection();

    // Cast to std::shared_ptr<Function<dim>> after parameter declaration.
    translational_velocity = translational_velocity_parsed;
    angular_velocity       = angular_velocity_parsed;
    solid_temperature      = solid_temperature_parsed;
  }

  template <int dim>
  void
  RigidSolidObject<dim>::parse_parameters(ParameterHandler &prm,
                                          unsigned int      id)
  {
    // Use ParsedFunction<dim> during parameter declaration. We need to do this
    // to use the parse_parameters function, which is not possible with
    // std::make_shared<Function<dim>> type. Please refer to the comment before
    // the translational_velocity and angular_velocity attributes declaration.
    auto translational_velocity_parsed =
      std::make_shared<Functions::ParsedFunction<dim>>(dim);
    auto angular_velocity_parsed =
      std::make_shared<Functions::ParsedFunction<dim>>(3);
    auto solid_temperature_parsed =
      std::make_shared<Functions::ParsedFunction<dim>>(1);

    prm.enter_subsection("solid object " + Utilities::int_to_string(id, 1));
    {
      solid_mesh.parse_parameters(prm);
      prm.enter_subsection("translational velocity");
      translational_velocity_parsed->parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("angular velocity");
      angular_velocity_parsed->parse_parameters(prm);
      prm.leave_subsection();

      const std::vector<double> temp =
        convert_string_to_vector<double>(prm, "center of rotation");

      AssertThrow(temp.size() == dim,
                  ExcMessage("Invalid center of rotation. This should be a " +
                             Utilities::int_to_string(dim) +
                             " dimensional point."));

      for (unsigned int i = 0; i < dim; ++i)
        {
          center_of_rotation[i] = temp.at(i);
        }

      output_bool = prm.get_bool("output solid object");

      const std::string thermal_type = prm.get("thermal boundary type");
      if (thermal_type == "adiabatic")
        thermal_boundary_type = ThermalBoundaryType::adiabatic;
      else if (thermal_type == "isothermal")
        thermal_boundary_type = ThermalBoundaryType::isothermal;
      else
        {
          throw(std::runtime_error("Invalid thermal boundary type"));
        }

      // Isothermal boundary
      prm.enter_subsection("temperature");
      solid_temperature_parsed->parse_parameters(prm);
      prm.leave_subsection();
    }
    prm.leave_subsection();
    translational_velocity = translational_velocity_parsed;
    angular_velocity       = angular_velocity_parsed;
    solid_temperature      = solid_temperature_parsed;
  }


  template <int dim>
  class DEMSolidObjects
  {
  public:
    DEMSolidObjects()
    {}

    void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);

    // Calculate forces
    Verbosity verbosity;

    // DEM solid objects
    std::vector<std::shared_ptr<RigidSolidObject<dim>>> solid_surfaces;
    std::vector<std::shared_ptr<RigidSolidObject<dim>>> solid_volumes;
    unsigned int                                        number_solid_surfaces;
    unsigned int                                        number_solid_volumes;

    static const unsigned int max_number_of_solids = 50;
  };

  template <int dim>
  void
  DEMSolidObjects<dim>::declare_parameters(ParameterHandler &prm)
  {
    solid_surfaces.resize(max_number_of_solids);
    solid_volumes.resize(max_number_of_solids);
    number_solid_surfaces = 0;
    number_solid_volumes  = 0;

    prm.enter_subsection("solid objects");
    {
      prm.enter_subsection("solid surfaces");
      {
        prm.declare_entry("number of solids",
                          "0",
                          Patterns::Integer(),
                          "Number of solid surfaces");

        for (unsigned int i_solid = 0; i_solid < max_number_of_solids;
             ++i_solid)
          {
            solid_surfaces[i_solid] = std::make_shared<RigidSolidObject<dim>>();
            solid_surfaces[i_solid]->declare_parameters(prm, i_solid);
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("solid volumes");
      {
        prm.declare_entry("number of solids",
                          "0",
                          Patterns::Integer(),
                          "Number of solid volumes");

        for (unsigned int i_solid = 0; i_solid < max_number_of_solids;
             ++i_solid)
          {
            solid_volumes[i_solid] = std::make_shared<RigidSolidObject<dim>>();
            solid_volumes[i_solid]->declare_parameters(prm, i_solid);
          }
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  DEMSolidObjects<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("solid objects");
    {
      prm.enter_subsection("solid surfaces");
      {
        number_solid_surfaces = prm.get_integer("number of solids");
        for (unsigned int i_solid = 0; i_solid < number_solid_surfaces;
             ++i_solid)
          {
            solid_surfaces[i_solid]->parse_parameters(prm, i_solid);
          }
      }
      prm.leave_subsection();

      prm.enter_subsection("solid volumes");
      {
        number_solid_volumes = prm.get_integer("number of solids");
        for (unsigned int i_solid = 0; i_solid < number_solid_volumes;
             ++i_solid)
          {
            solid_volumes[i_solid]->parse_parameters(prm, i_solid);
          }
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
} // namespace Parameters

#endif
