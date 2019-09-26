#include "core/parameters.h"
#include "core/simulationcontrol.h"

#include <fstream>

void
printTime(ConditionalOStream pcout, SimulationControl control)
{
  if (control.getMethod() == Parameters::SimulationControl::steady)
    {
      pcout << std::endl;
      pcout
        << "*****************************************************************"
           "******************"
        << std::endl;
      pcout << "Steady iteration : " << std::setw(8) << std::right
            << control.getIter() << "/" << control.getNbMeshAdapt() + 1
            << std::endl;
      pcout
        << "*****************************************************************"
           "******************"
        << std::endl;
    }
  else
    {
      pcout << std::endl;
      pcout
        << "*****************************************************************"
           "******************"
        << std::endl;
      pcout << "Transient iteration : " << std::setw(8) << std::left
            << control.getIter() << " Time : " << std::setw(8) << std::left
            << control.getTime() << " Time step : " << std::setw(8) << std::left
            << control.getCurrentTimeStep() << " CFL : " << std::setw(8)
            << std::left << control.getCFL() << std::endl;
      pcout
        << "*****************************************************************"
           "******************"
        << std::endl;
    }
}

void
SimulationControl::addTimeStep(double p_timestep)
{
  // Store previous time step in table
  for (unsigned int i_time = dt.size() - 1; i_time > 0; --i_time)
    dt[i_time] = dt[i_time - 1];

  // Calculate time step, right now this is a dummy function
  dt[0] = p_timestep;
}

bool
SimulationControl::integrate()
{
  if ((parameterControl.method == parameterControl.steady &&
       iter >= (nbMeshAdapt + 1)) ||
      (parameterControl.method != parameterControl.steady &&
       time >= (endTime - 1e-6 * dt[0])))
    return false;
  else
    {
      iter++;
      addTimeStep(calculateTimeStep());
      // Increment time
      time += dt[0];
      return true;
    }
}

// Calculate the time step depending on the time stepping control parameters
double
SimulationControl::calculateTimeStep()
{
  return parameterControl.dt;
}

void
SimulationControl::initialize(ParameterHandler &prm)
{
  parameterControl.parse_parameters(prm);
  method = parameterControl.method;
  // Even if high order time stepping schemes are not used
  // dt always contains the necessary information to restart
  // at the highest order method available
  dt.resize(numberTimeStepStored);
  dt[0]       = parameterControl.dt;
  endTime     = parameterControl.timeEnd;
  maxCFL      = parameterControl.maxCFL;
  nbMeshAdapt = parameterControl.nbMeshAdapt;
  time        = 0;
  iter        = 0;
  CFL         = 0;
}

void
SimulationControl::initialize(Parameters::SimulationControl param)
{
  parameterControl = param;
  method           = parameterControl.method;
  dt.resize(numberTimeStepStored);
  dt[0]       = parameterControl.dt;
  endTime     = parameterControl.timeEnd;
  maxCFL      = parameterControl.maxCFL;
  nbMeshAdapt = parameterControl.nbMeshAdapt;
  time        = 0;
  iter        = 0;
  CFL         = 0;
}

void
SimulationControl::save(std::string prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ofstream output(filename.c_str());
  output << "Simulation control" << std::endl;
  for (unsigned int i = 0; i < dt.size(); ++i)
    output << "dt_" << i << " " << dt[i] << std::endl;
  output << "CFL  " << CFL << std::endl;
  output << "Time " << time << std::endl;
  output << "Iter " << iter << std::endl;
}

void
SimulationControl::read(std::string prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ifstream input(filename.c_str());
  if (!input)
    {
      throw("Unable to open file");
    }
  std::string buffer;
  std::getline(input, buffer);
  for (unsigned int i = 0; i < dt.size(); ++i)
    input >> buffer >> dt[i];
  input >> buffer >> CFL;
  input >> buffer >> time;
  input >> buffer >> iter;
}
