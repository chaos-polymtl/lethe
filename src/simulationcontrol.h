#include "parameters.h"

#ifndef LETHE_SIMULATIONCONTROL_H
#define LETHE_SIMULATIONCONTROL_H

class SimulationControl
{
  // Time step
  double dt;
  // CFL
  double CFL;
  // Maximal CFL condition
  double maxCFL;
  // Time
  double time;
  // Simulation end time
  double endTime;
  // Iteration number
  unsigned int iter;
  // Number of mesh adaptation iteration
  unsigned int nbMeshAdapt;

  // Time stepping method
  Parameters::SimulationControl::TimeSteppingMethod method;

  // Parameters from the parser that do not change during the simulation (names, etc.)
  Parameters::SimulationControl parameterControl;
public:
  void initialize(ParameterHandler &prm);
  void initialize(Parameters::SimulationControl param);

  Parameters::SimulationControl::TimeSteppingMethod getMethod(){return method;}
  void setMethod(Parameters::SimulationControl::TimeSteppingMethod  p_method){ method=p_method;}

  std::string getOuputName(){return parameterControl.output_name;}
  std::string getOutputFolder(){return parameterControl.output_folder;}

  double getTimeStep(){return dt;}
  double getTime(){return time;}
  double getEndTime(){return endTime;}
  unsigned int getIter(){return iter;}
  bool         firstIter(){return iter==1;}
  double       getCFL(){return CFL;}
  void         setCFL(double p_CFL) {CFL=p_CFL;}
  double       getMaxCFL(){return maxCFL;}
  unsigned int getNbMeshAdapt(){return nbMeshAdapt;}
  unsigned int getSubdivision(){return parameterControl.subdivision;}


  bool isOutputIteration() {return (iter%parameterControl.outputFrequency==0);}

  bool integrate();

  void save(std::string filename);
  void read(std::string filename);

};

void printTime(ConditionalOStream pcout, SimulationControl control);

#endif

