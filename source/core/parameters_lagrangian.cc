/*
 * parameters.cpp
 *
 *  Created on: Dec 16, 2019
 *      Author: shahab
 */

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
        dt             = prm.get_double("time step");
        tFinal         = prm.get_integer("step end");
        nTotal         = prm.get_integer("n total");
        writeFrequency = prm.get_integer("write frequency");
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
        prm.declare_entry("Yp",
                          "1.",
                          Patterns::Double(),
                          "Young's modulus of particle");
        prm.declare_entry("Yw",
                          "1.",
                          Patterns::Double(),
                          "Young's modulus of wall");
        prm.declare_entry("vp",
                          "1.",
                          Patterns::Double(),
                          "Poisson's ratio of particle");
        prm.declare_entry("vw",
                          "1.",
                          Patterns::Double(),
                          "Poisson's ratio of wall");
        prm.declare_entry("ep",
                          "1.",
                          Patterns::Double(),
                          "Coefficient of restitution of particle");
        prm.declare_entry("ew",
                          "1.",
                          Patterns::Double(),
                          "Coefficient of restitution of wall");
        prm.declare_entry("mup",
                          "1.",
                          Patterns::Double(),
                          "Friction coefficient of particle");
        prm.declare_entry("muw",
                          "1.",
                          Patterns::Double(),
                          "Friction coefficient of wall");
        prm.declare_entry("murp",
                          "1.",
                          Patterns::Double(),
                          "Rolling friction coefficient of particle");
        prm.declare_entry("murw",
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
        gx       = prm.get_double("gx");
        gy       = prm.get_double("gy");
        gz       = prm.get_double("gz");
        diameter = prm.get_double("diameter");
        density  = prm.get_double("density");
        Yp       = prm.get_double("Yp");
        Yw       = prm.get_double("Yw");
        vp       = prm.get_double("vp");
        vw       = prm.get_double("vw");
        ep       = prm.get_double("ep");
        ew       = prm.get_double("ew");
        mup      = prm.get_double("mup");
        muw      = prm.get_double("muw");
        murp     = prm.get_double("murp");
        murw     = prm.get_double("murw");
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
      }
      prm.leave_subsection();
    }

    void
    InsertionInfo::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("insertion info");
      {
        tInsertion = prm.get_integer("Insertion time step");
        nInsert =
          prm.get_integer("Inserted number of particles at each time step");
        insertFrequency    = prm.get_integer("Insertion frequency");
        x_min              = prm.get_double("Insertion box minimum x");
        y_min              = prm.get_double("Insertion box minimum y");
        z_min              = prm.get_double("Insertion box minimum z");
        x_max              = prm.get_double("Insertion box maximum x");
        y_max              = prm.get_double("Insertion box maximum y");
        z_max              = prm.get_double("Insertion box maximum z");
        distance_threshold = prm.get_double("Insertion distance threshold");
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
        prm.declare_entry("Number of fields",
                          "1",
                          Patterns::Integer(),
                          "Number of fields of properties");
      }
      prm.leave_subsection();
    }

    void
    OutputProperties::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("output properties");
      {
        numProperties = prm.get_integer("Number of properties");
        numFields     = prm.get_integer("Number of fields");
      }
      prm.leave_subsection();
    }

    void
    SimulationModel::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("simulation model");
      {}
      prm.leave_subsection();
    }

    void
    SimulationModel::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("simulation model");
      {}
      prm.leave_subsection();
    }

  } // namespace Lagrangian

} // namespace Parameters
