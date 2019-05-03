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
  void initialize(ParameterHandler &prm)
  {
    parameterControl.parse_parameters (prm);
    method     = parameterControl.method;
    dt         = parameterControl.dt;
    endTime    = parameterControl.timeEnd;
    maxCFL     = parameterControl.maxCFL;
    nbMeshAdapt= parameterControl.nbMeshAdapt;
    time=0;
    iter=0;
    CFL=0;
  }
  Parameters::SimulationControl::TimeSteppingMethod getMethod(){return method;}
  void setMethod(Parameters::SimulationControl::TimeSteppingMethod  p_method){ method=p_method;}

  std::string getOuputName(){return parameterControl.output_name;}
  std::string getOutputFolder(){return parameterControl.output_folder;}

  double getTimeStep(){return dt;}
  double getTime(){return time;}
  unsigned int getIter(){return iter;}
  bool         firstIter(){return iter==1;}
  unsigned int getCFL(){return CFL;}
  unsigned int getNbMeshAdapt(){return nbMeshAdapt;}
  unsigned int getSubdivision(){return parameterControl.subdivision;}


  bool isOutputIteration() {return (iter%parameterControl.outputFrequency==0);}

  bool integrate()
  {
    if ( (parameterControl.method==parameterControl.steady && iter>=(nbMeshAdapt+1)) || (parameterControl.method==parameterControl.backward && time >=(endTime-1e-6*dt))) return false;
    else
    {
      iter++;
      time += dt;
      return true;
    }
  }
};

void printTime(ConditionalOStream pcout, SimulationControl control);

#endif

