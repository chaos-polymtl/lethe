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

#ifndef PARAMETERS_DEM_H_
#  define PARAMETERS_DEM_H_

/**
 * Handles all the parameters declared in the parameter handler file
 */

template <int dim>
class DEMSolverParameters
{
public:
  Parameters::Mesh                           mesh;
  Parameters::Testing                        test;
  Parameters::Timer                          timer;
  Parameters::SimulationControl              simulation_control;
  Parameters::Lagrangian::PhysicalProperties physicalProperties;
  Parameters::Lagrangian::InsertionInfo      insertionInfo;
  Parameters::Lagrangian::ModelParameters    model_parmeters;


  void
  declare(ParameterHandler &prm)
  {
    Parameters::SimulationControl::declare_parameters(prm);
    Parameters::Mesh::declare_parameters(prm);
    Parameters::Timer::declare_parameters(prm);
    Parameters::Testing::declare_parameters(prm);
    Parameters::Lagrangian::PhysicalProperties::declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo::declare_parameters(prm);
    Parameters::Lagrangian::ModelParameters::declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    mesh.parse_parameters(prm);
    test.parse_parameters(prm);
    timer.parse_parameters(prm);
    physicalProperties.parse_parameters(prm);
    insertionInfo.parse_parameters(prm);
    model_parmeters.parse_parameters(prm);
    simulation_control.parse_parameters(prm);
  }
};

#endif /* DEMPARAMTERS_H_ */
