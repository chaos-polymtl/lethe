#include "core/simulation_control.h"

#include <cfloat>
#include <fstream>
#include <iomanip>
#include <sstream>


SimulationControl::SimulationControl(const Parameters::SimulationControl param)
  : method(param.method)
  , assembly_method(param.method)
  , current_time(0)
  , time_step(param.dt)
  , initial_time_step(param.dt)
  , end_time(param.timeEnd)
  , iteration_number(0)
  , number_mesh_adapt(param.number_mesh_adaptation)
  , CFL(0)
  , max_CFL(param.maxCFL)
  , residual(DBL_MAX)
  , stop_tolerance(param.stop_tolerance)
  , output_frequency(param.output_frequency)
  , output_time_frequency(param.output_time)
  , log_frequency(param.log_frequency)
  , log_precision(param.log_precision)
  , subdivision(param.subdivision)
  , group_files(param.group_files)
  , output_name(param.output_name)
  , output_path(param.output_folder)
  , output_boundaries(param.output_boundaries)
  , first_assembly(true)
  , bdf_start_method(param.bdf_startup_method)
  , startup_timestep_scaling(param.startup_timestep_scaling)
{
  time_step_vector.resize(numberTimeStepStored);
  time_step_vector[0] = param.dt;
  time_step_vector[1] = param.dt;
  time_step_vector[2] = param.dt;
  time_step_vector[3] = param.dt;
}

void
SimulationControl::add_time_step(double p_timestep)
{
  time_step = p_timestep;
  // Store previous time step in table
  for (unsigned int i_time = time_step_vector.size() - 1; i_time > 0; --i_time)
    time_step_vector[i_time] = time_step_vector[i_time - 1];

  // Calculate time step, right now this is a dummy function
  time_step_vector[0] = p_timestep;
}

bool
SimulationControl::is_output_iteration()
{
  if (output_frequency == 0)
    return false;
  else
    {
      return (get_step_number() % output_frequency == 0);
    }
}

void
SimulationControl::update_assembly_method()
{
  if (iteration_number <= 1 &&
      method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 &&
      bdf_start_method ==
        Parameters::SimulationControl::BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf1;
      set_current_time_step(initial_time_step * startup_timestep_scaling);
    }
  else if (iteration_number == 2 &&
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 &&
           bdf_start_method == Parameters::SimulationControl::
                                 BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf2;
      set_suggested_time_step(initial_time_step *
                              (1 - startup_timestep_scaling));
    }
  else if (iteration_number == 3 &&
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf2 &&
           bdf_start_method == Parameters::SimulationControl::
                                 BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf2;
      set_suggested_time_step(initial_time_step);
    }
  else if (iteration_number <= 1 &&
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf3 &&
           bdf_start_method == Parameters::SimulationControl::
                                 BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf1;
      set_current_time_step(initial_time_step * startup_timestep_scaling);
    }
  else if (iteration_number == 2 &&
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf3 &&
           bdf_start_method == Parameters::SimulationControl::
                                 BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf2;
      set_suggested_time_step(initial_time_step * startup_timestep_scaling);
    }
  else if (iteration_number == 3 &&
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf3 &&
           bdf_start_method == Parameters::SimulationControl::
                                 BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf3;
      set_suggested_time_step(initial_time_step *
                              (1 - 2 * startup_timestep_scaling));
    }
  else if (iteration_number == 4 &&
           method == Parameters::SimulationControl::TimeSteppingMethod::bdf3 &&
           bdf_start_method == Parameters::SimulationControl::
                                 BDFStartupMethods::multiple_step_bdf)
    {
      assembly_method = Parameters::SimulationControl::TimeSteppingMethod::bdf3;
      set_suggested_time_step(initial_time_step);
    }
  else if (iteration_number <= 1 &&
           (method ==
            Parameters::SimulationControl::TimeSteppingMethod::bdf2) &&
           bdf_start_method ==
             Parameters::SimulationControl::BDFStartupMethods::sdirk_step)
    {
      assembly_method =
        Parameters::SimulationControl::TimeSteppingMethod::sdirk22;
      set_suggested_time_step(initial_time_step);
    }
  else if (iteration_number <= 2 &&
           (method ==
            Parameters::SimulationControl::TimeSteppingMethod::bdf3) &&
           bdf_start_method ==
             Parameters::SimulationControl::BDFStartupMethods::sdirk_step)
    {
      assembly_method =
        Parameters::SimulationControl::TimeSteppingMethod::sdirk33;
      set_suggested_time_step(initial_time_step);
    }
  else
    {
      assembly_method = method;
    }
}



bool
SimulationControl::is_verbose_iteration()
{
  return (get_step_number() % log_frequency == 0);
}

void
SimulationControl::save(std::string prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ofstream output(filename.c_str());
  output << "Simulation control" << std::endl;
  for (unsigned int i = 0; i < time_step_vector.size(); ++i)
    output << "dt_" << i << " " << time_step_vector[i] << std::endl;
  output << "CFL  " << CFL << std::endl;
  output << "Time " << current_time << std::endl;
  output << "Iter " << iteration_number << std::endl;
}

void
SimulationControl::read(std::string prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ifstream input(filename.c_str());
  AssertThrow(input, ExcFileNotOpen(filename));

  std::string buffer;
  std::getline(input, buffer);
  for (unsigned int i = 0; i < time_step_vector.size(); ++i)
    input >> buffer >> time_step_vector[i];
  input >> buffer >> CFL;
  input >> buffer >> current_time;
  input >> buffer >> iteration_number;
}


SimulationControlTransient::SimulationControlTransient(
  Parameters::SimulationControl param)
  : SimulationControl(param)
  , adapt(param.adapt)
  , adaptative_time_step_scaling(param.adaptative_time_step_scaling)
{}

void
SimulationControlTransient::print_progression(const ConditionalOStream &pcout)
{
  if (!is_verbose_iteration())
    return;

  pcout << std::endl;

  std::stringstream ss;

  // Copy information into a string stream
  ss << "Transient iteration: " << std::setw(8) << std::left << iteration_number
     << " Time: " << std::setw(8) << std::left << current_time
     << " Time step: " << std::setw(8) << std::left << time_step
     << " CFL: " << std::setw(8) << std::left << SimulationControl::get_CFL();

  // Announce string
  announce_string(pcout, ss.str(), '*');
}

bool
SimulationControlTransient::integrate()
{
  if (!is_at_end())
    {
      previous_time  = current_time;
      first_assembly = true;
      iteration_number++;
      update_assembly_method();
      add_time_step(calculate_time_step());
      current_time += time_step;

      return true;
    }

  else
    return false;
}



bool
SimulationControlTransient::is_at_end()
{
  double floating_point_margin = std::max(1e-6 * time_step, 1e-12 * end_time);
  return current_time >= (end_time - floating_point_margin);
}

double
SimulationControlTransient::calculate_time_step()
{
  double new_time_step = time_step;

  if (adapt && iteration_number > 1)
    {
      new_time_step = time_step * adaptative_time_step_scaling;
      if (CFL > 0 && max_CFL / CFL < adaptative_time_step_scaling)
        new_time_step = time_step * max_CFL / CFL;
    }
  if (current_time + new_time_step > end_time)
    new_time_step = end_time - current_time;

  return new_time_step;
}

SimulationControlTransientDEM::SimulationControlTransientDEM(
  Parameters::SimulationControl param)
  : SimulationControlTransient(param)
{}

void
SimulationControlTransientDEM::print_progression(
  const ConditionalOStream &pcout)
{
  if (!is_verbose_iteration())
    return;

  pcout << std::endl;
  std::stringstream ss;

  // Copy information into a string stream
  ss << "Transient iteration: " << std::setw(8) << std::left << iteration_number
     << " Time: " << std::setw(8) << std::left << current_time
     << " Time step: " << std::setw(8) << std::left << time_step;

  // Announce string
  announce_string(pcout, ss.str(), '*');
}


SimulationControlTransientDynamicOutput::
  SimulationControlTransientDynamicOutput(Parameters::SimulationControl param)
  : SimulationControlTransient(param)
  , time_step_forced_output(false)
  // To be fixed for restarts
  , last_output_time(0.)
{}

double
SimulationControlTransientDynamicOutput::calculate_time_step()
{
  double new_time_step = time_step;
  if (time_step_forced_output)
    {
      new_time_step           = time_step_vector[1];
      time_step_forced_output = false;
    }
  else if (iteration_number > 1)
    {
      new_time_step = time_step * adaptative_time_step_scaling;
      if (CFL > 0 && max_CFL / CFL < adaptative_time_step_scaling)
        new_time_step = time_step * max_CFL / CFL;
    }

  if (current_time + new_time_step > end_time)
    new_time_step = end_time - current_time;

  if (current_time + new_time_step > last_output_time + output_time_frequency)
    {
      new_time_step = last_output_time + output_time_frequency - current_time;
      time_step_forced_output = true;
    }

  return new_time_step;
}

bool
SimulationControlTransientDynamicOutput::is_output_iteration()
{
  bool is_output_time =
    (current_time - last_output_time) - output_time_frequency >
    -1e-12 * output_time_frequency;
  if (is_output_time)
    last_output_time = current_time;

  return is_output_time;
}


SimulationControlSteady::SimulationControlSteady(
  Parameters::SimulationControl param)
  : SimulationControl(param)
{}

bool
SimulationControlSteady::integrate()
{
  if (!is_at_end())
    {
      previous_time = current_time;
      iteration_number++;
      // Fix the time to the iteration number so that pvd outputs
      // do not lead to confusing results
      current_time = iteration_number;
      return true;
    }

  else
    return false;
}

void
SimulationControlSteady::print_progression(const ConditionalOStream &pcout)
{
  if (!is_verbose_iteration())
    return;

  pcout << std::endl;
  std::stringstream ss;

  // Copy information into a string stream
  ss << "Steady iteration: " << std::setw(8) << std::right << iteration_number
     << "/" << number_mesh_adapt + 1;

  // Announce string
  announce_string(pcout, ss.str(), '*');
}

bool
SimulationControlSteady::is_at_end()
{
  return iteration_number >= (number_mesh_adapt + 1);
}

void
SimulationControlAdjointSteady::print_progression(
  const ConditionalOStream &pcout)
{
  if (!is_verbose_iteration())
    return;

  pcout << std::endl;
  std::stringstream ss;

  // Copy information into a string stream
  ss << "Pseudo steady-state iteration: " << std::setw(8) << std::left
     << iteration_number << " Time: " << std::setw(8) << std::left
     << current_time << " Time step: " << std::setw(8) << std::left << time_step
     << " CFL: " << std::setw(8) << std::left << SimulationControl::get_CFL();

  // Announce string
  announce_string(pcout, ss.str(), '*');
}

bool
SimulationControlAdjointSteady::is_at_end()
{
  return residual <= stop_tolerance;
}

SimulationControlAdjointSteady::SimulationControlAdjointSteady(
  Parameters::SimulationControl param)
  : SimulationControlTransient(param)
{}

double
SimulationControlAdjointSteady::calculate_time_step()
{
  double new_time_step = time_step;

  if (adapt && iteration_number > 1)
    {
      new_time_step = time_step * adaptative_time_step_scaling;
      if (CFL > 0 && max_CFL / CFL < adaptative_time_step_scaling)
        new_time_step = time_step * max_CFL / CFL;
    }

  return new_time_step;
}
