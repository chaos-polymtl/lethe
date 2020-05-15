#include "core/simulation_flow_control.h"

#include <fstream>

#include "core/parameters.h"


SimulationFlowControl::SimulationFlowControl(
  Parameters::SimulationControl param,
  double                        p_start_time,
  double                        p_end_time,
  double                        p_step)
  : DiscreteTime(p_start_time, p_end_time, p_step)
{
  time_step_vector.resize(numberTimeStepStored);
  time_step_vector[0] = param.dt;
  end_time            = param.timeEnd;
  max_CFL             = param.maxCFL;
  number_mesh_adapt   = param.nbMeshAdapt;
  time                = 0;
  iter                = 0;
  CFL                 = 0;
}

// void
// printTime(ConditionalOStream pcout, SimulationFlowControl control)
//{
//  if (control.getMethod() ==
//      Parameters::SimulationControl::TimeSteppingMethod::steady)
//    {
//      pcout << std::endl;
//      pcout
//        << "*****************************************************************"
//           "******************"
//        << std::endl;
//      pcout << "Steady iteration : " << std::setw(8) << std::right
//            << control.getIter() << "/" << control.getNbMeshAdapt() + 1
//            << std::endl;
//      pcout
//        << "*****************************************************************"
//           "******************"
//        << std::endl;
//    }
//  else
//    {
//      pcout << std::endl;
//      pcout
//        << "*****************************************************************"
//           "******************"
//        << std::endl;
//      pcout << "Transient iteration : " << std::setw(8) << std::left
//            << control.getIter() << " Time : " << std::setw(8) << std::left
//            << control.getTime() << " Time step : " << std::setw(8) <<
//            std::left
//            << control.getCurrentTimeStep() << " CFL : " << std::setw(8)
//            << std::left << control.getCFL() << std::endl;
//      pcout
//        << "*****************************************************************"
//           "******************"
//        << std::endl;
//    }
//}

void
SimulationFlowControl::addTimeStep(double p_timestep)
{
  // Store previous time step in table
  for (unsigned int i_time = time_step_vector.size() - 1; i_time > 0; --i_time)
    time_step_vector[i_time] = time_step_vector[i_time - 1];

  // Calculate time step, right now this is a dummy function
  time_step_vector[0] = p_timestep;
}

// bool
// SimulationFlowControl::integrate()
//{
//  if ((parameterControl.method ==
//         Parameters::SimulationControl::TimeSteppingMethod::steady &&
//       iter >= (number_mesh_adapt + 1)) ||
//      (parameterControl.method !=
//         Parameters::SimulationControl::TimeSteppingMethod::steady &&
//       time >= (end_time - 1e-6 * dt[0])))
//    return false;
//  else
//    {
//      iter++;
//      addTimeStep(calculateTimeStep());
//      // Increment time
//      time += dt[0];
//      return true;
//    }
//}

void
SimulationFlowControl::save(std::string prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ofstream output(filename.c_str());
  output << "Simulation control" << std::endl;
  for (unsigned int i = 0; i < time_step_vector.size(); ++i)
    output << "dt_" << i << " " << time_step_vector[i] << std::endl;
  output << "CFL  " << CFL << std::endl;
  output << "Time " << time << std::endl;
  output << "Iter " << iter << std::endl;
}

void
SimulationFlowControl::read(std::string prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ifstream input(filename.c_str());
  if (!input)
    {
      throw("Unable to open file");
    }
  std::string buffer;
  std::getline(input, buffer);
  for (unsigned int i = 0; i < time_step_vector.size(); ++i)
    input >> buffer >> time_step_vector[i];
  input >> buffer >> CFL;
  input >> buffer >> time;
  input >> buffer >> iter;
}


SimulationControlTransient::SimulationControlTransient(
  Parameters::SimulationControl param)
  : SimulationFlowControl(param, 0., param.timeEnd, param.dt)
{}

void
SimulationControlTransient::print_progression(ConditionalOStream &pcout)
{
  pcout << std::endl;
  pcout << "*****************************************************************"
           "******************"
        << std::endl;
  pcout << "Transient iteration : " << std::setw(8) << std::left
        << DiscreteTime::get_step_number() << " Time : " << std::setw(8)
        << std::left << DiscreteTime::get_next_time()
        << " Time step : " << std::setw(8) << std::left
        << DiscreteTime::get_next_step_size() << " CFL : " << std::setw(8)
        << std::left << SimulationFlowControl::get_CFL() << std::endl;
  pcout << "*****************************************************************"
           "******************"
        << std::endl;
}

void
SimulationControlSteady::print_progression(ConditionalOStream &pcout)
{
  pcout << std::endl;
  pcout << "*****************************************************************"
           "******************"
        << std::endl;
  pcout << "Steady iteration : " << std::setw(8) << std::right
        << DiscreteTime::get_step_number() << "/" << number_mesh_adapt + 1
        << std::endl;
  pcout << "*****************************************************************"
           "******************"
        << std::endl;
}
