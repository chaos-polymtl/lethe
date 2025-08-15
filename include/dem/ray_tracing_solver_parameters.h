// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ray_tracing_solver_parameters_h
#define lethe_ray_tracing_solver_parameters_h

#include <core/parameters.h>
#include <core/parameters_lagrangian.h>
#include <core/simulation_control.h>
// #include <core/solid_objects_parameters.h>

/**
 * @brief Handles all the parameters declared in the parameter handler file.
 */
template <int dim>
class RayTracingSolverParameters
{
public:
  Parameters::Mesh              mesh;
  Parameters::Timer             timer;
  Parameters::SimulationControl simulation_control;

  Parameters::Lagrangian::InsertionInfo<dim> particle_insertion_info;

  // Maybe usefull later
  // Parameters::Lagrangian::FloatingWalls<dim>     floating_walls;
  // Parameters::Lagrangian::FloatingGrid<dim>      floating_grid;
  // std::shared_ptr<Parameters::DEMSolidObjects<dim>> solid_objects;

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
    Parameters::Lagrangian::InsertionInfo<dim>::declare_parameters(prm);
    Parameters::Lagrangian::ModelParameters<dim>::declare_parameters(prm);
    // floating_walls.declare_parameters(prm);
    // floating_grid.declare_parameters(prm);
    // solid_objects = std::make_shared<Parameters::DEMSolidObjects<dim>>();
    // solid_objects->declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    mesh.parse_parameters(prm);
    timer.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    // floating_walls.parse_parameters(prm);
    // floating_grid.parse_parameters(prm);
    // solid_objects->parse_parameters(prm);
  }
};

#endif
