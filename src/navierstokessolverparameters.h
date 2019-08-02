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


#ifndef LETHE_NAVIERSTOKESSOLVERPARAMETERS_H
#define LETHE_NAVIERSTOKESSOLVERPARAMETERS_H

#include "exactsolutions.h"
#include "forcingfunctions.h"
#include "boundaryconditions.h"
#include "initialconditions.h"
#include "parameters.h"
#include "simulationcontrol.h"



template <int dim>
class NavierStokesSolverParameters
{
public:

  Parameters::Testing                   test;
  Parameters::LinearSolver              linearSolver;
  Parameters::NonLinearSolver           nonLinearSolver;
  Parameters::MeshAdaptation            meshAdaptation;
  Parameters::Mesh                      mesh;
  Parameters::PhysicalProperties        physicalProperties;
  Parameters::Timer                     timer;
  Parameters::FEM                       femParameters;
  Parameters::Forces                    forcesParameters;
  Parameters::AnalyticalSolution        analyticalSolution;
  Parameters::Restart                   restartParameters;
  BoundaryConditions::NSBoundaryConditions<dim>   boundaryConditions;
  Parameters::InitialConditions<dim>    *initialCondition;
  SimulationControl                     simulationControl;

  void declare(ParameterHandler &prm)
  {
    initialCondition = new Parameters::InitialConditions<dim>;
    Parameters::SimulationControl::declare_parameters (prm);
    Parameters::PhysicalProperties::declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    Parameters::Restart::declare_parameters(prm);
    boundaryConditions.declare_parameters(prm);
    initialCondition->declare_parameters(prm);
    Parameters::FEM::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::Forces::declare_parameters(prm);
    Parameters::MeshAdaptation::declare_parameters (prm);
    Parameters::NonLinearSolver::declare_parameters (prm);
    Parameters::LinearSolver::declare_parameters (prm);
    Parameters::AnalyticalSolution::declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);
  }

  void parse(ParameterHandler &prm)
  {
    test.parse_parameters(prm);
    linearSolver.parse_parameters (prm);
    nonLinearSolver.parse_parameters (prm);
    meshAdaptation.parse_parameters(prm);
    mesh.parse_parameters(prm);
    physicalProperties.parse_parameters(prm);
    timer.parse_parameters(prm);
    femParameters.parse_parameters(prm);
    forcesParameters.parse_parameters(prm);
    analyticalSolution.parse_parameters(prm);
    restartParameters.parse_parameters(prm);
    boundaryConditions.parse_parameters(prm);
    initialCondition->parse_parameters(prm);
    simulationControl.initialize(prm);
  }
};

#endif
