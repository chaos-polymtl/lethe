/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */
#include <core/parameters.h>
#include <core/parameters_lagrangian.h>
#include <core/simulation_control.h>

#ifndef parameters_DEM_h
#  define parameters_DEM_h

/**
 * Handles all the parameters declared in the parameter handler file
 */

template <int dim>
class DEMSolverParameters
{
public:
  Parameters::Mesh                                mesh;
  Parameters::Testing                             test;
  Parameters::Restart                             restart;
  Parameters::Timer                               timer;
  Parameters::SimulationControl                   simulation_control;
  Parameters::Lagrangian::PhysicalProperties<dim> physical_properties;
  Parameters::Lagrangian::InsertionInfo           insertion_info;
  Parameters::Lagrangian::ModelParameters         model_parameters;
  Parameters::Lagrangian::FloatingWalls<dim>      floating_walls;
  Parameters::Lagrangian::BCDEM<dim>              boundary_conditions;
  Parameters::Lagrangian::ForceTorqueOnWall<dim>  forces_torques;
  Parameters::Lagrangian::GridMotion<dim>         grid_motion;

  void
  declare(ParameterHandler &prm)
  {
    Parameters::SimulationControl::declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    Parameters::Restart::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);
    physical_properties.declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo::declare_parameters(prm);
    Parameters::Lagrangian::ModelParameters::declare_parameters(prm);
    floating_walls.declare_parameters(prm);
    boundary_conditions.declare_parameters(prm);
    forces_torques.declare_parameters(prm);
    grid_motion.declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    mesh.parse_parameters(prm);
    test.parse_parameters(prm);
    restart.parse_parameters(prm);
    timer.parse_parameters(prm);
    physical_properties.parse_parameters(prm);
    insertion_info.parse_parameters(prm);
    model_parameters.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    floating_walls.parse_parameters(prm);
    boundary_conditions.parse_parameters(prm);
    forces_torques.parse_parameters(prm);
    grid_motion.parse_parameters(prm);
  }
};

#endif /* parameters_DEM_h */
