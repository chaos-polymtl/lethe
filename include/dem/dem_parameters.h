/*
 * DEMParamters.h
 *
 *  Created on: Dec 16, 2019
 *      Author: shahab
 */

#ifndef DEMPARAMTERS_H_
#define DEMPARAMTERS_H_

//#include "parameters.h"
#include <core/parameters_lagrangian.h>

template <int dim>
class demparameters
{
public:
  Parameters::Lagrangian::SimulationControl  simulationControl;
  Parameters::Lagrangian::PhysicalProperties physicalProperties;
  Parameters::Lagrangian::InsertionInfo      insertionInfo;
  Parameters::Lagrangian::OutputProperties   outputProperties;

  void
  declare(ParameterHandler &prm)
  {
    Parameters::Lagrangian::SimulationControl::declare_parameters(prm);
    Parameters::Lagrangian::PhysicalProperties::declare_parameters(prm);
    Parameters::Lagrangian::InsertionInfo::declare_parameters(prm);
    Parameters::Lagrangian::OutputProperties::declare_parameters(prm);
  }

  void
  parse(ParameterHandler &prm)
  {
    physicalProperties.parse_parameters(prm);
    insertionInfo.parse_parameters(prm);
    simulationControl.parse_parameters(prm);
    outputProperties.parse_parameters(prm);
  }
};

#endif /* DEMPARAMTERS_H_ */
