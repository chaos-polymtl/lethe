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

#ifndef PARAMETERS_DEM_H_
#define PARAMETERS_DEM_H_

#include <core/parameters_lagrangian.h>

template <int dim>
class ParametersDEM
{
public:
  Parameters::Lagrangian::SimulationControl  simulationControl;
  Parameters::Lagrangian::PhysicalProperties physicalProperties;
  Parameters::Lagrangian::InsertionInfo      insertionInfo;
  Parameters::Lagrangian::OutputProperties   outputProperties;
  Parameters::Lagrangian::SimulationModel    simulationModel;

  void
  declare(ParameterHandler &prm)
  {
    Parameters::Lagrangian::SimulationControl::declare_parameters(prm);
    Parameters::Lagrangian::PhysicalProperties::declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo::declare_parameters(prm);
    Parameters::Lagrangian::OutputProperties::declare_parameters(prm);
    Parameters::Lagrangian::SimulationModel::declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    physicalProperties.parse_parameters(prm);
    insertionInfo.parse_parameters(prm);
    simulationControl.parse_parameters(prm);
    outputProperties.parse_parameters(prm);
    simulationModel.parse_parameters(prm);
  }
};

#endif /* DEMPARAMTERS_H_ */
