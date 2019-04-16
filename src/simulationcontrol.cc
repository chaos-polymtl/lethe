#include "simulationcontrol.h"

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
