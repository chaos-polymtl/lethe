#include "exactsolutions.h"
#include "forcingfunctions.h"
#include "boundaryconditions.h"
#include "initialconditions.h"
#include "parameters.h"
#include "simulationcontrol.h"

#ifndef LETHE_NAVIERSTOKESSOLVERPARAMETERS_H
#define LETHE_NAVIERSTOKESSOLVERPARAMETERS_H

template <int dim>
class NavierStokesSolverParameters
{
public:

  Parameters::Testing                   test;
  Parameters::LinearSolver              linearSolver;
  Parameters::NonLinearSolver           nonLinearSolver;
  Parameters::MeshAdaptation            meshAdaptation;
  Parameters::Mesh                      meshParameters;
  Parameters::PhysicalProperties        physicalProperties;
  Parameters::Timer                     clock;
  Parameters::FEM                       femParameters;
  Parameters::Forces                    forcesParameters;
  Parameters::AnalyticalSolution        analyticalSolution;
  BoundaryConditions::NSBoundaryConditions<dim>   boundaryConditions;
  Parameters::InitialConditions<dim>    *initialCondition;
  SimulationControl                     simulationControl;

  void declare(ParameterHandler &prm)
  {
    initialCondition = new Parameters::InitialConditions<dim>;
    Parameters::Testing::declare_parameters(prm);
    Parameters::NonLinearSolver::declare_parameters (prm);
    Parameters::LinearSolver::declare_parameters (prm);
    Parameters::SimulationControl::declare_parameters (prm);
    Parameters::MeshAdaptation::declare_parameters (prm);
    Parameters::Mesh::declare_parameters(prm);
    Parameters::PhysicalProperties::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::FEM::declare_parameters(prm);
    Parameters::Forces::declare_parameters(prm);
    Parameters::AnalyticalSolution::declare_parameters(prm);
    boundaryConditions.declare_parameters(prm);
    initialCondition->declare_parameters(prm);
  }

  void parse(ParameterHandler &prm)
  {
    test.parse_parameters(prm);
    linearSolver.parse_parameters (prm);
    nonLinearSolver.parse_parameters (prm);
    meshAdaptation.parse_parameters(prm);
    meshParameters.parse_parameters(prm);
    physicalProperties.parse_parameters(prm);
    clock.parse_parameters(prm);
    femParameters.parse_parameters(prm);
    analyticalSolution.parse_parameters(prm);
    forcesParameters.parse_parameters(prm);
    initialCondition->parse_parameters(prm);
    boundaryConditions.parse_parameters(prm);
    simulationControl.initialize(prm);
  }
};

#endif
