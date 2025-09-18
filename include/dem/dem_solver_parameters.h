// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_solver_parameters_h
#define lethe_dem_solver_parameters_h

#include <core/parameters.h>
#include <core/parameters_lagrangian.h>
#include <core/solid_objects_parameters.h>

/**
 * @brief Handles all the parameters declared in the parameter handler file.
 */
template <int dim>
class DEMSolverParameters
{
public:
  Parameters::Mesh              mesh;
  Parameters::Testing           test;
  Parameters::Restart           restart;
  Parameters::Timer             timer;
  Parameters::SimulationControl simulation_control;
  Parameters::Lagrangian::LagrangianPhysicalProperties
                                                 lagrangian_physical_properties;
  Parameters::Lagrangian::InsertionInfo<dim>     insertion_info;
  Parameters::Lagrangian::ModelParameters<dim>   model_parameters;
  Parameters::Lagrangian::FloatingWalls<dim>     floating_walls;
  Parameters::Lagrangian::BCDEM                  boundary_conditions;
  Parameters::Lagrangian::FloatingGrid<dim>      floating_grid;
  Parameters::Lagrangian::ForceTorqueOnWall<dim> forces_torques;
  Parameters::Lagrangian::GridMotion<dim>        grid_motion;
  Parameters::Lagrangian::LagrangianPostProcessing  post_processing;
  std::shared_ptr<Parameters::DEMSolidObjects<dim>> solid_objects;

  void
  declare(ParameterHandler &prm)
  {
    prm.declare_entry("dimension",
                      "0",
                      Patterns::Integer(),
                      "Dimension of the problem");

    prm.declare_entry("print parameters",
                      "none",
                      Patterns::Selection("none|only changed|all"),
                      "Print all the parameters, or only"
                      "the changed parameters or none");

    Parameters::SimulationControl::declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    Parameters::Restart::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);
    lagrangian_physical_properties.declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo<dim>::declare_parameters(prm);
    Parameters::Lagrangian::ModelParameters<dim>::declare_parameters(prm);
    floating_walls.declare_parameters(prm);
    floating_grid.declare_parameters(prm);
    boundary_conditions.declare_parameters(prm);
    forces_torques.declare_parameters(prm);
    grid_motion.declare_parameters(prm);
    post_processing.declare_parameters(prm);
    solid_objects = std::make_shared<Parameters::DEMSolidObjects<dim>>();
    solid_objects->declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    mesh.parse_parameters(prm);
    test.parse_parameters(prm);
    restart.parse_parameters(prm);
    timer.parse_parameters(prm);
    lagrangian_physical_properties.parse_parameters(prm);
    insertion_info.parse_parameters(prm);
    model_parameters.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    floating_walls.parse_parameters(prm);
    floating_grid.parse_parameters(prm);
    boundary_conditions.parse_parameters(prm);
    forces_torques.parse_parameters(prm);
    grid_motion.parse_parameters(prm);
    post_processing.parse_parameters(prm);
    solid_objects->parse_parameters(prm);
  }
};

#endif
