#include "simulationcontrol.h"
#include <fstream>


void printTime(ConditionalOStream pcout, SimulationControl control)
{
  if (control.getMethod()==Parameters::SimulationControl::steady)
  {
    pcout<< std::endl;
    pcout<< "***********************************************************************************" << std::endl;
    pcout<< "Steady iteration : "   << std::setw(8) << std::right << control.getIter() << "/" << control.getNbMeshAdapt()+1 << std::endl;
    pcout<< "***********************************************************************************" << std::endl;
  }
  else
  {
    pcout<< std::endl;
    pcout<< "***********************************************************************************" << std::endl;
    pcout<< "Transient iteration : "   << std::setw(8) << std::left << control.getIter()
         << " Time : "       << std::setw(8) << std::left << control.getTime()
         << " Time step : " << std::setw(8) << std::left <<  control.getTimeStep()
         << " CFL : "       << std::setw(8) << std::left <<  control.getCFL() << std::endl;
    pcout<< "***********************************************************************************" << std::endl;
  }
}

bool SimulationControl::integrate()
{
  if ( (parameterControl.method==parameterControl.steady && iter>=(nbMeshAdapt+1))
       || (parameterControl.method==parameterControl.backward && time >=(endTime-1e-6*dt)))
    return false;
  else
  {
    iter++;
    time += dt;
    return true;
  }
}

void SimulationControl::initialize(ParameterHandler &prm)
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
void SimulationControl::initialize(Parameters::SimulationControl param)
{
  parameterControl=param;
  method     = parameterControl.method;
  dt         = parameterControl.dt;
  endTime    = parameterControl.timeEnd;
  maxCFL     = parameterControl.maxCFL;
  nbMeshAdapt= parameterControl.nbMeshAdapt;
  time=0;
  iter=0;
  CFL=0;
}

void SimulationControl::save(std::string prefix)
{
  std::string filename = prefix + ".simulationcontrol";
  std::ofstream output (filename.c_str());
  output << "Simulation control" << std::endl;
  output << "dt   " << dt << std::endl;
  output << "CFL  " << CFL << std::endl;
  output << "Time " << time << std::endl;
  output << "Iter " << iter << std::endl;
}

void SimulationControl::read(std::string prefix)
{
  std::string filename = prefix + ".simulationcontrol";
  std::ifstream input (filename.c_str());
  if (!input) {
       throw("Unable to open file");
   }
  std::string buffer;
  std::getline(input,buffer);
  input >> buffer >> dt          ;
  input >> buffer >> CFL         ;
  input >> buffer >> time        ;
  input >> buffer >> iter        ;
}
