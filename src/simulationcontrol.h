
#ifndef LETHE_SIMULATIONCONTROL_H
#define LETHE_SIMULATIONCONTROL_H

#include "parameters.h"


class SimulationControl
{
  // Time step
  std::vector<double> dt;
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

  // number of time steps stored
  static const unsigned int numberTimeStepStored=4;

  // Calculate time step based on either CFL or fixed;
  double calculateTimeStep();

  // Add a time step and stores the previous one in a list
  void  addTimeStep(double p_timestep);


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


  void   setTimeStep(double p_timestep) {addTimeStep(p_timestep);}
  double getCurrentTimeStep(){return dt[0];}
  std::vector<double> getTimeSteps(){return dt;}
  double getTime(){return time;}
  double getEndTime(){return endTime;}


  unsigned int getIter(){return iter;}
  bool         firstIter(){return iter==1;}
  double       getCFL(){return CFL;}
  void         setCFL(double p_CFL) {CFL=p_CFL;}
  double       getMaxCFL(){return maxCFL;}
  Parameters::SimulationControl getParameters() {return parameterControl;}

  unsigned int getNbMeshAdapt(){return nbMeshAdapt;}
  unsigned int getSubdivision(){return parameterControl.subdivision;}


  bool isOutputIteration() {return (iter%parameterControl.outputFrequency==0);}

  bool integrate();

  void save(std::string filename);
  void read(std::string filename);

};

void printTime(ConditionalOStream pcout, SimulationControl control);

#endif

