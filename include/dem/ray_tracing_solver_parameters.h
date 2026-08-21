// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ray_tracing_solver_parameters_h
#define lethe_ray_tracing_solver_parameters_h

#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/parameters_lagrangian.h>
#include <core/simulation_control.h>

/**
 * @brief Handles all the parameters declared in the parameter handler file.
 */
template <int dim>
class RayTracingSolverParameters
{
public:
  Parameters::Mesh                                mesh;
  Parameters::Manifolds                           manifolds_parameters;
  Parameters::Testing                             test;
  Parameters::Timer                               timer;
  Parameters::SimulationControl                   simulation_control;
  Parameters::Lagrangian::InsertionInfo<dim>      particle_insertion_info;
  Parameters::Lagrangian::ParticleRayTracing<dim> ray_tracing_info;
  Parameters::Lagrangian::ModelParameters<dim>    model_parameters;

  /**
   * @brief Declare all the parameters of a ray tracing simulation.
   *
   * @param[in,out] prm The parameter handler.
   *
   * @param[in] size_of_subsections The maximum size of the variable-size
   * subsections of the parameter file, as returned by
   * Parameters::get_size_of_subsections.
   */
  void
  declare(ParameterHandler                   &prm,
          const Parameters::SizeOfSubsections size_of_subsections)
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

    prm.declare_entry(
      "comment message",
      "",
      Patterns::Anything(),
      "Print a comment at the beginning of the console output.");

    Parameters::Mesh::declare_parameters(prm);
    manifolds_parameters.declare_parameters(prm, size_of_subsections.manifolds);
    Parameters::Testing::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::SimulationControl::declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo<dim>::declare_parameters(prm);
    Parameters::Lagrangian::ParticleRayTracing<dim>::declare_parameters(prm);
    Parameters::Lagrangian::ModelParameters<dim>::declare_parameters(prm);
  }

  /**
   * @brief Handles the parsing of all the parameters in the parameter handler file.
   *
   * @param[in,out] prm The parameter handler.
   */
  void
  parse(ParameterHandler &prm)
  {
    mesh.parse_parameters(prm);
    manifolds_parameters.parse_parameters(prm);
    test.parse_parameters(prm);
    timer.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    particle_insertion_info.parse_parameters(prm);
    ray_tracing_info.parse_parameters(prm);
    model_parameters.parse_parameters(prm);
  }
};

#endif
