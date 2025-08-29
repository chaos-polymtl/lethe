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
  Parameters::Mesh                                mesh;
  Parameters::Testing                             test;
  Parameters::Timer                               timer;
  Parameters::SimulationControl                   simulation_control;
  Parameters::Lagrangian::InsertionInfo<dim>      particle_insertion_info;
  Parameters::Lagrangian::ParticleRayTracing<dim> ray_tracing_info;
  Parameters::Lagrangian::ModelParameters<dim>    model_parameters;

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

    Parameters::Mesh::declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::SimulationControl::declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo<dim>::declare_parameters(prm);
    Parameters::Lagrangian::ParticleRayTracing<dim>::declare_parameters(prm);
    Parameters::Lagrangian::ModelParameters<dim>::declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    mesh.parse_parameters(prm);
    test.parse_parameters(prm);
    timer.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    particle_insertion_info.parse_parameters(prm);
    ray_tracing_info.parse_parameters(prm);
    model_parameters.parse_parameters(prm);
  }
};

#endif
