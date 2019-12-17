/*
 * DEMParamters.h
 *
 *  Created on: Dec 16, 2019
 *      Author: shahab
 */

#ifndef DEMPARAMTERS_H_
#define DEMPARAMTERS_H_

//#include "parameters.h"
#include <core/parameters_dem.h>

template <int dim>
class demparamters {
public:
	Parameters::SimulationControl			simulationControl;
	Parameters::PhysicalProperties			physicalProperties;
	Parameters::InsertionInfo				insertionInfo;
	Parameters::OutputProperties			outputProperties;

	void declare(ParameterHandler &prm)
	  {
		Parameters::SimulationControl::declare_parameters(prm);
		Parameters::PhysicalProperties::declare_parameters(prm);
		Parameters::InsertionInfo::declare_parameters(prm);
		Parameters::OutputProperties::declare_parameters(prm);
	  }

	  void parse(ParameterHandler &prm)
	  {
		  physicalProperties.parse_parameters(prm);
		  insertionInfo.parse_parameters(prm);
		  simulationControl.parse_parameters(prm);
		  outputProperties.parse_parameters(prm);

	  }
};

#endif /* DEMPARAMTERS_H_ */
