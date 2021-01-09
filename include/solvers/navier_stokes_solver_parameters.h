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
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef lethe_navier_stokes_solver_parameters_h
#define lethe_navier_stokes_solver_parameters_h

#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/parameters_multiphysics.h>

#include "analytical_solutions.h"
#include "initial_conditions.h"
#include "nitsche.h"
#include "source_terms.h"

template <int dim>
class NavierStokesSolverParameters
{
public:
  Parameters::Testing                            test;
  Parameters::LinearSolver                       linear_solver;
  Parameters::NonLinearSolver                    non_linear_solver;
  Parameters::MeshAdaptation                     mesh_adaptation;
  Parameters::Mesh                               mesh;
  std::shared_ptr<Parameters::Nitsche<dim>>      nitsche;
  Parameters::PhysicalProperties                 physical_properties;
  Parameters::SimulationControl                  simulation_control;
  Parameters::Timer                              timer;
  Parameters::FEM                                fem_parameters;
  Parameters::Forces                             forces_parameters;
  Parameters::PostProcessing                     post_processing;
  Parameters::Restart                            restart_parameters;
  Parameters::Manifolds                          manifolds_parameters;
  BoundaryConditions::NSBoundaryConditions<dim>  boundary_conditions;
  BoundaryConditions::HTBoundaryConditions<dim>  boundary_conditions_ht;
  Parameters::InitialConditions<dim> *           initial_condition;
  AnalyticalSolutions::AnalyticalSolution<dim> * analytical_solution;
  SourceTerms::SourceTerm<dim> *                 sourceTerm;
  Parameters::VelocitySource                     velocitySource;
  Parameters::IBParticles<dim>                   particlesParameters;
  std::shared_ptr<Parameters::VoidFraction<dim>> void_fraction;
  Parameters::DynamicFlowControl                 flow_control;
  Parameters::Multiphysics                       multiphysics;



  void
  declare(ParameterHandler &prm)
  {
    Parameters::SimulationControl::declare_parameters(prm);
    Parameters::PhysicalProperties::declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    nitsche = std::make_shared<Parameters::Nitsche<dim>>();
    nitsche->declare_parameters(prm);
    Parameters::Restart::declare_parameters(prm);
    boundary_conditions.declare_parameters(prm);
    boundary_conditions_ht.declare_parameters(prm);

    initial_condition = new Parameters::InitialConditions<dim>;
    initial_condition->declare_parameters(prm);

    Parameters::FEM::declare_parameters(prm);
    Parameters::Multiphysics::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::Forces::declare_parameters(prm);
    Parameters::MeshAdaptation::declare_parameters(prm);
    Parameters::NonLinearSolver::declare_parameters(prm);
    Parameters::LinearSolver::declare_parameters(prm);
    Parameters::PostProcessing::declare_parameters(prm);
    Parameters::DynamicFlowControl ::declare_parameters(prm);
    particlesParameters.declare_parameters(prm);
    manifolds_parameters.declare_parameters(prm);

    analytical_solution = new AnalyticalSolutions::AnalyticalSolution<dim>;
    analytical_solution->declare_parameters(prm);
    sourceTerm = new SourceTerms::SourceTerm<dim>;
    sourceTerm->declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);

    Parameters::VelocitySource::declare_parameters(prm);

    void_fraction = std::make_shared<Parameters::VoidFraction<dim>>();
    void_fraction->declare_parameters(prm);

    multiphysics.declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    test.parse_parameters(prm);
    linear_solver.parse_parameters(prm);
    non_linear_solver.parse_parameters(prm);
    mesh_adaptation.parse_parameters(prm);
    mesh.parse_parameters(prm);
    nitsche->parse_parameters(prm);
    physical_properties.parse_parameters(prm);
    multiphysics.parse_parameters(prm);
    timer.parse_parameters(prm);
    fem_parameters.parse_parameters(prm);
    forces_parameters.parse_parameters(prm);
    post_processing.parse_parameters(prm);
    flow_control.parse_parameters(prm);
    restart_parameters.parse_parameters(prm);
    boundary_conditions.parse_parameters(prm);
    boundary_conditions_ht.parse_parameters(prm);
    manifolds_parameters.parse_parameters(prm);
    initial_condition->parse_parameters(prm);
    analytical_solution->parse_parameters(prm);
    sourceTerm->parse_parameters(prm);
    simulation_control.parse_parameters(prm);
    velocitySource.parse_parameters(prm);
    particlesParameters.parse_parameters(prm);
    void_fraction->parse_parameters(prm);
    multiphysics.parse_parameters(prm);
  }
};

#endif
