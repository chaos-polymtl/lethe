// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file solid_objects_parameters.h
 * @brief Parameter classes for solid objects used in Nitsche immersed boundary
 * and DEM simulations.
 *
 * This file defines the parameter structures for configuring solid objects
 * within Lethe simulations. It provides three levels of solid object
 * management:
 * - NitscheObject: individual solid for the Nitsche immersed boundary method.
 * - Nitsche: container managing multiple NitscheObject instances.
 * - RigidSolidObject: individual rigid solid for DEM simulations.
 * - DEMSolidObjects: container managing multiple RigidSolidObject instances
 *   as surfaces and volumes.
 */

#ifndef lethe_solid_objects_parameters_h
#define lethe_solid_objects_parameters_h

#include <core/parameters.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

using namespace dealii;

namespace Parameters
{
  /**
   * @brief Thermal boundary condition type applied on a solid object surface.
   */
  enum ThermalBoundaryType
  {
    /// No heat flux through the solid boundary.
    adiabatic,
    /// Fixed temperature imposed on the solid boundary.
    isothermal
  };

  /**
   * @brief Parameters for a single solid object used with the Nitsche immersed
   * boundary method.
   *
   * This class stores all the parameters needed to define and control an
   * individual Nitsche solid, including its mesh, velocity, temperature,
   * penalization coefficients, particle motion settings, and force/torque
   * calculation options.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class NitscheObject
  {
  public:
    NitscheObject()
      : solid_velocity(dim)
    {}

    /**
     * @brief Declare the parameters in the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     * @param[in] id Index of the solid object.
     */
    void
    declare_parameters(ParameterHandler &prm, unsigned int id);

    /**
     * @brief Parse the parameters from the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     * @param[in] id Index of the solid object.
     */
    void
    parse_parameters(ParameterHandler &prm, unsigned int id);

    /// Mesh parameters for the solid object.
    Parameters::Mesh solid_mesh;

    /// Number of quadrature points per 1D cell used to represent the solid.
    unsigned int number_quadrature_points;

    /// Velocity function imposed on the solid.
    Functions::ParsedFunction<dim> solid_velocity;

    /// Temperature function imposed on the solid boundary.
    Functions::ParsedFunction<dim> solid_temperature;

    /// Enable the motion of Nitsche particles.
    bool enable_particles_motion;

    /// Enable heat boundary condition on the immersed boundary.
    bool enable_heat_bc;

    /// Penalization parameter for the Nitsche method (momentum equation).
    double beta;

    /// Penalization parameter for the Nitsche method (heat equation).
    double beta_heat;

    /// Number of sub-iterations for particle motion integration.
    unsigned int particles_sub_iterations;

    /// Stop the simulation if Nitsche particles are lost.
    bool stop_particles_lost;

    /// Center of rotation used for torque calculation.
    Point<dim> center_of_rotation;

    /// Enable calculation of forces on the solid.
    bool calculate_force_on_solid;

    /// Enable calculation of torques on the solid.
    bool calculate_torque_on_solid;

    /// Center of rotation used for torque calculation.
    Point<dim> cor;

    /// File name prefix for the force output.
    std::string force_output_name;

    /// File name prefix for the torque output.
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

  /**
   * @brief Container class managing multiple Nitsche solid objects.
   *
   * This class holds a collection of NitscheObject instances and provides
   * parameter declaration and parsing for the Nitsche immersed boundary method
   * section of the parameter file.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class Nitsche
  {
  public:
    Nitsche()
    {}

    /**
     * @brief Declare the parameters in the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters from the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);

    /// Verbosity level for force/torque output.
    Verbosity verbosity;

    /// Collection of Nitsche solid objects.
    std::vector<std::shared_ptr<NitscheObject<dim>>> nitsche_solids;

    /// Number of active Nitsche solid objects.
    unsigned int number_solids;

    /// Maximum number of Nitsche solid objects allowed.
    static const unsigned int max_nitsche_solids = 10;
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


  /**
   * @brief Parameters for a single rigid solid object used in DEM simulations.
   *
   * This class stores the mesh, motion, temperature, and output parameters for
   * a rigid solid. The velocity fields are stored as
   * std::shared_ptr<Function<dim>> so that any function type can be used
   * (e.g., ParsedFunction for parameter-file input or ConstantFunction for
   * unit tests).
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class RigidSolidObject
  {
  public:
    RigidSolidObject()
    {}

    /**
     * @brief Declare the parameters in the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     * @param[in] id Index of the solid object.
     */
    void
    declare_parameters(ParameterHandler &prm, unsigned int id);

    /**
     * @brief Parse the parameters from the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     * @param[in] id Index of the solid object.
     */
    void
    parse_parameters(ParameterHandler &prm, unsigned int id);

    /// Mesh parameters for the solid object.
    Parameters::Mesh solid_mesh;

    /// Controls the generation of output files for this solid.
    bool output_bool;

    /// Translational velocity function of the solid object.
    std::shared_ptr<Function<dim>> translational_velocity;

    /// Angular velocity function of the solid object.
    std::shared_ptr<Function<dim>> angular_velocity;

    /// Center of rotation used to locate and rotate the solid object.
    Point<dim> center_of_rotation;

    /// Temperature function for the solid object boundary.
    std::shared_ptr<Function<dim>> solid_temperature;

    /// Type of thermal boundary condition applied on the solid surface.
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


  /**
   * @brief Container class managing rigid solid objects for DEM simulations.
   *
   * This class holds two collections of RigidSolidObject instances:
   * solid surfaces (2D boundaries) and solid volumes (3D bodies). It provides
   * parameter declaration and parsing for the "solid objects" section of the
   * parameter file.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class DEMSolidObjects
  {
  public:
    DEMSolidObjects()
    {}

    /**
     * @brief Declare the parameters in the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * @brief Parse the parameters from the parameter handler.
     *
     * @param[in,out] prm The parameter handler.
     */
    void
    parse_parameters(ParameterHandler &prm);

    /// Verbosity level for force/torque output.
    Verbosity verbosity;

    /// Collection of rigid solid surface objects.
    std::vector<std::shared_ptr<RigidSolidObject<dim>>> solid_surfaces;

    /// Collection of rigid solid volume objects.
    std::vector<std::shared_ptr<RigidSolidObject<dim>>> solid_volumes;

    /// Number of active solid surface objects.
    unsigned int number_solid_surfaces;

    /// Number of active solid volume objects.
    unsigned int number_solid_volumes;

    /// Maximum number of solid objects allowed.
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
