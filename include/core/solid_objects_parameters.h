/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Carole-Anne Daunais, Polytechnique Montreal, 2020 -
 */


#ifndef nitsche_h
#define nitsche_h

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
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
      solid_temperature.declare_parameters(prm, 1);
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

      prm.enter_subsection("center of rotation");
      prm.declare_entry("x", "0", Patterns::Double(), "X COR");
      prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
      prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
      prm.leave_subsection();

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

      prm.enter_subsection("center of rotation");
      center_of_rotation[0] = prm.get_double("x");
      center_of_rotation[1] = prm.get_double("y");
      if (dim == 3)
        center_of_rotation[2] = prm.get_double("z");
      prm.leave_subsection();

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
      : translational_velocity(dim)
      , angular_velocity(3)
    {}

    void
    declare_parameters(ParameterHandler &prm, unsigned int id);
    void
    parse_parameters(ParameterHandler &prm, unsigned int id);

    // Solid mesh
    Parameters::Mesh solid_mesh;

    // Solid velocity
    Functions::ParsedFunction<dim> translational_velocity;
    Functions::ParsedFunction<dim> angular_velocity;
    Point<dim>
      center_of_rotation; // Center of rotation used to locate the center of the
                          // object and also used to rotate the object
  };


  template <int dim>
  void
  RigidSolidObject<dim>::declare_parameters(ParameterHandler &prm,
                                            unsigned int      id)
  {
    prm.enter_subsection("solid object " + Utilities::int_to_string(id, 1));
    {
      solid_mesh.declare_parameters(prm);

      prm.enter_subsection("translational velocity");
      translational_velocity.declare_parameters(prm, dim);
      prm.leave_subsection();

      prm.enter_subsection("angular velocity");
      angular_velocity.declare_parameters(prm, 3);
      prm.leave_subsection();


      prm.enter_subsection("center of rotation");
      prm.declare_entry("x", "0", Patterns::Double(), "X COR");
      prm.declare_entry("y", "0", Patterns::Double(), "Y COR");
      prm.declare_entry("z", "0", Patterns::Double(), "Z COR");
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  RigidSolidObject<dim>::parse_parameters(ParameterHandler &prm,
                                          unsigned int      id)
  {
    prm.enter_subsection("solid object " + Utilities::int_to_string(id, 1));
    {
      solid_mesh.parse_parameters(prm);
      prm.enter_subsection("translational velocity");
      translational_velocity.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("angular velocity");
      angular_velocity.parse_parameters(prm);
      prm.leave_subsection();

      prm.enter_subsection("center of rotation");
      center_of_rotation[0] = prm.get_double("x");
      center_of_rotation[1] = prm.get_double("y");
      if (dim == 3)
        center_of_rotation[2] = prm.get_double("z");
      prm.leave_subsection();
    }
    prm.leave_subsection();
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
    std::vector<std::shared_ptr<RigidSolidObject<dim>>> solids;
    unsigned int                                        number_solids;
    static const unsigned int max_number_of_solids = 10;
  };

  template <int dim>
  void
  DEMSolidObjects<dim>::declare_parameters(ParameterHandler &prm)
  {
    solids.resize(max_number_of_solids);
    number_solids = 0;

    prm.enter_subsection("solid objects");
    {
      prm.declare_entry("number of solids",
                        "0",
                        Patterns::Integer(),
                        "Number of solid object");

      for (unsigned int i_solid = 0; i_solid < max_number_of_solids; ++i_solid)
        {
          solids[i_solid] = std::make_shared<RigidSolidObject<dim>>();
          solids[i_solid]->declare_parameters(prm, i_solid);
        }
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  DEMSolidObjects<dim>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("solid objects");
    {
      number_solids = prm.get_integer("number of solids");
      for (unsigned int i_solid = 0; i_solid < number_solids; ++i_solid)
        {
          solids[i_solid]->parse_parameters(prm, i_solid);
        }
    }
    prm.leave_subsection();
  }

} // namespace Parameters

#endif
