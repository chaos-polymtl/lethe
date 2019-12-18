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
class demparameters {
public:
	Parameters::DEM::SimulationControl			simulationControl;
	Parameters::DEM::PhysicalProperties			physicalProperties;
	Parameters::DEM::InsertionInfo				insertionInfo;
	Parameters::DEM::OutputProperties			outputProperties;

	void declare(ParameterHandler &prm)
	  {
		Parameters::DEM::SimulationControl::declare_parameters(prm);
		Parameters::DEM::PhysicalProperties::declare_parameters(prm);
		Parameters::DEM::InsertionInfo::declare_parameters(prm);
		Parameters::DEM::OutputProperties::declare_parameters(prm);
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
