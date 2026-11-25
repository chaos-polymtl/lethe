// SPDX-FileCopyrightText: Copyright (c) 2019-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/simulation_control.h>
#include <core/utilities.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>


SimulationControl::SimulationControl(const Parameters::SimulationControl &param)
  : method(param.method)
  , assembly_method(param.method)
  , current_time(0)
  , previous_time(0)
  , time_step(param.dt)
  , initial_time_step(param.dt)
  , end_time(param.time_end)
  , time_step_independent_of_end_time(param.time_step_independent_of_end_time)
  , iteration_number(0)
  , end_iteration(param.iteration_end)
  , number_mesh_adapt(param.number_mesh_adaptation)
  , CFL(0)
  , max_CFL(param.maxCFL)
  , capillary_time_step_constraint(std::numeric_limits<double>::max())
  , max_capillary_time_step_ratio(param.max_capillary_time_step_ratio)
  , current_capillary_time_step_ratio(0)
  , adapt_with_capillary_time_step_ratio(
      param.adapt_with_capillary_time_step_ratio)
  , residual(DBL_MAX)
  , stop_tolerance(param.stop_tolerance)
  , output_iteration_frequency(param.output_iteration_frequency)
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
  time_step_vector.resize(n_previous_time_steps);
  time_step_vector[0] = param.dt;
  time_step_vector[1] = param.dt;
  time_step_vector[2] = param.dt;
  time_step_vector[3] = param.dt;

  // Resize the bdf_coefficients to ensure they have a default size;
  bdf_coefs.reinit(n_previous_time_steps + 1);
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
  if (output_iteration_frequency == 0)
    return false;
  else
    {
      // Check if the current step number matches the output frequency
      return (get_iteration_number() % output_iteration_frequency == 0);
    }
}

bool
SimulationControl::is_refinement_step(
  const Parameters::MeshAdaptation &mesh_adaptation)
{
  bool is_refinement_step;

  if (mesh_adaptation.refinement_at_frequency)
    is_refinement_step =
      get_iteration_number() % mesh_adaptation.frequency == 0;
  else
    is_refinement_step = get_iteration_number() == 0;

  return is_refinement_step;
}

void
SimulationControl::update_assembly_method()
{
  if (bdf_start_method ==
      Parameters::SimulationControl::BDFStartupMethods::multiple_step_bdf)
    {
      if (iteration_number <= 1 &&
          method == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf1;
          set_current_time_step(initial_time_step * startup_timestep_scaling);
        }
      else if (iteration_number == 2 &&
               method ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf2;
          set_current_time_step(initial_time_step *
                                (1 - startup_timestep_scaling));
        }
      else if (iteration_number == 3 &&
               method ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf2)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf2;
          set_current_time_step(initial_time_step);
        }
      else if (iteration_number <= 1 &&
               method ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf1;
          set_current_time_step(initial_time_step * startup_timestep_scaling);
        }
      else if (iteration_number == 2 &&
               method ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf2;
          set_current_time_step(initial_time_step * startup_timestep_scaling);
        }
      else if (iteration_number == 3 &&
               method ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf3;
          set_current_time_step(initial_time_step *
                                (1 - 2 * startup_timestep_scaling));
        }
      else if (iteration_number == 4 &&
               method ==
                 Parameters::SimulationControl::TimeSteppingMethod::bdf3)
        {
          assembly_method =
            Parameters::SimulationControl::TimeSteppingMethod::bdf3;
          set_current_time_step(initial_time_step);
        }
    }
  else
    {
      assembly_method = method;
    }
}

unsigned int
SimulationControl::get_number_of_stages(
  const Parameters::SimulationControl::TimeSteppingMethod &method)
{
  using Method = Parameters::SimulationControl::TimeSteppingMethod;

  switch (method)
    {
      case Method::steady:
      case Method::steady_bdf:
      case Method::bdf1:
      case Method::bdf2:
      case Method::bdf3:
        return 1;

      // When referencing to SDIRK methods the pattern used is
      // sdirkOrderStage. For instance, sdirk22 indicates SDIRK method with
      // order 2 and 2 stages
      case Method::sdirk22:
        return 2;

      case Method::sdirk33:
        return 3;

      case Method::sdirk45:
        return 5;

      default:
        AssertThrow(
          false,
          ExcMessage(
            "Unknown TimeSteppingMethod passed to "
            "get_number_of_stages(). Please add the number of stages of the new method in core/simulation_control.cc."));
        return 1;
    }
}



bool
SimulationControl::is_verbose_iteration()
{
  return (get_iteration_number() % log_frequency == 0);
}

void
SimulationControl::save(const std::string &prefix)
{
  std::string   filename = prefix + ".simulationcontrol";
  std::ofstream output(filename.c_str());
  output << "Simulation control" << std::endl;

  // We set the precision to 17. If the time-step is significantly smaller then
  // the end time of the simulation, round off errors can become problematic.
  output << std::fixed << std::setprecision(17);

  for (unsigned int i = 0; i < time_step_vector.size(); ++i)
    output << "dt_" << i << " " << time_step_vector[i] << std::endl;
  output << "CFL  " << CFL << std::endl;
  output << "Time " << current_time << std::endl;
  output << "Iter " << iteration_number << std::endl;
}

void
SimulationControl::read(const std::string &prefix)
{
  std::string filename = prefix + ".simulationcontrol";
  check_file_exists(filename,
                    "checkpoint file of the simulation control, given by "
                    "'subsection restart' - 'set filename'");
  std::ifstream input(filename.c_str());

  std::string buffer;
  std::getline(input, buffer);
  for (auto &time_step_value : time_step_vector)
    input >> buffer >> time_step_value;
  input >> buffer >> CFL;
  input >> buffer >> current_time;
  input >> buffer >> iteration_number;

  // The previous time is only known once the next iteration is integrated.
  // Until then, set it to the current time so that queries made before the
  // first iteration after the restart (e.g., is_output_iteration) are defined.
  previous_time = current_time;

  // Fix time step to be the last time_step that was used
  time_step = time_step_vector[0];
}

std::vector<double>
SimulationControl::get_checkpointed_simulation_control_info(
  const std::string &prefix)
{
  std::string filename = prefix + ".simulationcontrol";
  check_file_exists(filename,
                    "checkpoint file of the simulation control, given by "
                    "'subsection restart' - 'set filename'");
  std::ifstream input(filename.c_str());

  // Store the time steps and last checkpointed time without modifying the
  // simulation control information
  std::vector<double> time_steps(time_step_vector.size());
  double              last_checkpointed_time;

  // Returned vector with the last checkpointed time and time step
  std::vector<double> simulation_control_info(2);

  std::string buffer;
  std::getline(input, buffer);
  for (unsigned int i = 0; i < time_step_vector.size(); ++i)
    input >> buffer >> time_steps[i];
  input >> buffer >> buffer;
  input >> buffer >> last_checkpointed_time;
  input >> buffer >> buffer;

  // Fix the current time as final checkpointed time
  simulation_control_info[0] = last_checkpointed_time;

  // Fix time step to be the last time_step that was used
  simulation_control_info[1] = time_steps[0];

  return simulation_control_info;
}

SimulationControlTransient::SimulationControlTransient(
  const Parameters::SimulationControl &param)
  : SimulationControl(param)
  , time_step_adaptation_required(param.time_step_adaptation_required)
  , adapt_with_cfl(param.adapt_with_cfl)
  , adaptative_time_step_scaling(param.adaptative_time_step_scaling)
  , max_dt(param.max_dt)
  , output_time_frequency(param.output_time_frequency)
  , output_times_vector(param.output_times_vector)
  , override_time_step_on_restart(param.override_time_step_on_restart)
  , output_time_interval(param.output_time_interval)
  , end_control(param.end_control)
  , output_control(param.output_control)
{
  AssertThrow(output_control !=
                  Parameters::SimulationControl::OutputControl::time ||
                output_time_frequency == -1 || output_time_frequency > 0,
              ExcMessage("The output time frequency must be strictly positive. "
                         "Set it to -1 to disable it."));
}

void
SimulationControlTransient::print_progression(const ConditionalOStream &pcout)
{
  if (!is_verbose_iteration())
    return;

  pcout << std::endl;

  std::stringstream ss;
  ss << std::setprecision(this->log_precision);

  // Copy information into a string stream
  ss << "Transient iteration: " << std::setw(8) << std::left << iteration_number
     << " Time: " << std::setw(8) << std::left << current_time
     << " Time step: " << std::setw(8) << std::left << time_step
     << " CFL: " << std::setw(8) << std::left << this->get_CFL();

  unsigned int first_line_size = ss.str().size();

  if (adapt_with_capillary_time_step_ratio)
    {
      // Get the length of blank spaces to add
      std::stringstream cfl_ss;
      cfl_ss << std::setprecision(this->log_precision) << std::setw(8)
             << std::left << this->get_CFL();
      const unsigned int cfl_value_length = cfl_ss.str().size();
      const unsigned int blank_space_length =
        first_line_size - std::strlen(" CFL: ") - cfl_value_length;

      // Print CTR value on a new line, under CFL
      ss << '\n'
         << std::setw(blank_space_length) << " "
         << " CTR: " << std::setw(8) << std::left
         << this->get_current_capillary_time_step_ratio();
    }

  // Announce string
  announce_string(pcout, ss.str(), first_line_size, '*');
}

bool
SimulationControlTransient::integrate()
{
  if (!is_at_end())
    {
      // We store the value of the time at the previous iteration
      // before we change the current time.
      previous_time = current_time;

      // Reset the first assembly to true to indicate that we are starting
      // the time step. This variable is used to monitor the initial residual
      // of the set of non-linear equation which can be used for example in
      // steady-bdf methods
      first_assembly = true;

      // We increment the iteration number since its initial value is zero.
      iteration_number++;

      // Update assembly method is used to change the time marching method
      // in the case where the methods are not self-starting (all BDF of orders
      // 2 and above)
      update_assembly_method();
      add_time_step(calculate_time_step());
      current_time += time_step;

      // Calculate CTR to print on the console if capillary time-step constraint
      // is enabled
      if (adapt_with_capillary_time_step_ratio)
        {
          set_current_capillary_time_step_ratio();
        }

      if (is_bdf())
        update_bdf_coefficients();
      return true;
    }

  else
    return false;
}

bool
SimulationControlTransient::is_at_end()
{
  if (end_control == Parameters::SimulationControl::EndControl::iteration)
    return iteration_number >= end_iteration;

  double floating_point_margin = std::max(1e-6 * time_step, 1e-12 * end_time);
  return current_time >= (end_time - floating_point_margin);
}

double
SimulationControlTransient::calculate_time_step()
{
  double new_time_step = time_step;

  if (time_step_adaptation_required && iteration_number > 1)
    {
      new_time_step = time_step * adaptative_time_step_scaling;
      if (adapt_with_cfl && CFL > 0 &&
          max_CFL / CFL < adaptative_time_step_scaling)
        new_time_step = time_step * max_CFL / CFL;

      new_time_step = std::min(new_time_step, max_dt);

      if (adapt_with_capillary_time_step_ratio)
        {
          double capillary_time_step =
            capillary_time_step_constraint * max_capillary_time_step_ratio;
          new_time_step = std::min(new_time_step, capillary_time_step);
        }
    }

  // Ensure that the time step for the last iteration is kept regardless of the
  // end time set. The end time is also irrelevant when the simulation ends at a
  // given iteration, in which case the time step is never clipped.
  if (time_step_independent_of_end_time ||
      end_control == Parameters::SimulationControl::EndControl::iteration)
    return new_time_step;

  // Modify last time step to ensure that the last iteration is exactly the end
  // time specified in the parameter file
  if (current_time + new_time_step > end_time)
    new_time_step = end_time - current_time;

  return new_time_step;
}

bool
SimulationControlTransient::is_output_iteration()
{
  // If iteration control the only options are to not print or at certain
  // frequency of iterations for the time interval given. In this mode, we also
  // always output the last iteration in the simulation.
  if (output_control == Parameters::SimulationControl::OutputControl::iteration)
    {
      if (output_iteration_frequency == 0)
        return false;
      else
        {
          // Check if the current step number matches the following condition:
          // (The step number matches the output frequency OR is the last
          // time step) AND the current time is within the output time interval
          return ((get_iteration_number() % output_iteration_frequency == 0 ||
                   is_at_end()) &&
                  get_current_time() >= output_time_interval[0] &&
                  get_current_time() <= output_time_interval[1]);
        }
    }

  // If time control there are several options:

  // Case 1. A specific output time frequency is given (with or without a
  // specific time interval). The output times are
  // t_k = t_start + k * output_time_frequency, where t_start is the lower bound
  // of the output time interval. This criterion only depends on the current
  // and previous times, so it returns the same result when it is called
  // several times within the same iteration.
  if (output_time_frequency != -1)
    {
      const double start_time = output_time_interval[0];

      // The tolerance prevents round-off errors accumulated in the current time
      // from delaying or skipping an output. It is proportional to the output
      // time frequency, and not to the time step, so that it is the same at
      // every iteration. Consequently, every output time is attributed to
      // exactly one iteration, even when the time step changes.
      const double margin = output_time_tolerance * output_time_frequency;

      // We never output outside of the output time interval. The same
      // tolerance is used so that the bounds of the interval are considered as
      // reached when they are reached up to round-off errors.
      if (current_time < start_time - margin ||
          current_time > output_time_interval[1] + margin)
        return false;

      // As for iteration control, we always output the initial iteration and
      // the last iteration of the simulation, as long as they lie within the
      // output time interval (checked above)
      if (iteration_number == 0 || is_at_end())
        return true;

      // We output if the time step that led to the current time reached or
      // crossed an output time
      const auto output_time_index = [&](const double time) {
        return std::floor((time - start_time + margin) / output_time_frequency);
      };

      return output_time_index(current_time) > output_time_index(previous_time);
    }

  // Time reached at the end of the next time step. It is used by the two
  // remaining cases: since the specific output times and the bounds of the
  // output time interval do not necessarily correspond to the time of an
  // iteration, we also output the iteration that precedes them.
  const double next_time = current_time + calculate_time_step();

  // Case 2. One or several specific output times are given, in any order. We
  // output the iteration whose time step reaches or crosses one of these times,
  // and the iteration that precedes it, i.e., the iteration whose next time
  // step reaches or crosses it. Both conditions amount to an output time lying
  // within (previous_time, next_time]. At the initial iteration, there is no
  // previous time step, so only the output times equal to the current time are
  // considered as reached. The tolerance prevents round-off errors accumulated
  // in the current time from delaying an output by one iteration. This
  // criterion only depends on the current, previous and next times, so it
  // returns the same result when it is called several times within the same
  // iteration.
  if (!output_times_vector.empty() && output_times_vector[0] != -1)
    {
      const double margin      = output_time_tolerance * time_step;
      const double lower_bound = (iteration_number == 0) ?
                                   current_time - margin :
                                   previous_time + margin;
      const double upper_bound = next_time + margin;

      return std::ranges::any_of(output_times_vector,
                                 [&](const double output_time) {
                                   return output_time > lower_bound &&
                                          output_time <= upper_bound;
                                 });
    }

  // Case 3. Only a specific time interval is specified
  const double lower_bound = output_time_interval[0];
  const double upper_bound = output_time_interval[1];

  // We output if the current time is within the interval
  bool is_output_time =
    (current_time >= lower_bound && current_time <= upper_bound);

  // We always write one step before, in case the interval lower bound does not
  // correspond exactly to a time iteration performed (due to time step)
  if (current_time < lower_bound && next_time >= lower_bound)
    is_output_time = true;

  // We always write one step after, in case the interval upper bound does not
  // correspond exactly to a time iteration performed (due to time step)
  if (current_time >= upper_bound && previous_time < upper_bound)
    is_output_time = true;

  return is_output_time;
}

void
SimulationControlTransient::read(const std::string &prefix)
{
  // Checkpoints written by previous versions of Lethe with time output control
  // contain an additional line with the index of the next specific output time.
  // This line follows the information read here and is ignored, since the
  // output criterion no longer depends on it.
  SimulationControl::read(prefix);

  if (override_time_step_on_restart)
    {
      // Fix the time step to the new provided value.
      // We understand that users may wish to override the checkpointed
      // time-step value with another one.
      const double old_CFL       = CFL;
      const double old_time_step = time_step;
      set_current_time_step(initial_time_step);
      set_CFL(old_CFL * initial_time_step / old_time_step);
    }
}

SimulationControlTransientDEM::SimulationControlTransientDEM(
  const Parameters::SimulationControl &param)
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


SimulationControlSteady::SimulationControlSteady(
  const Parameters::SimulationControl &param)
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

SimulationControlAdjointSteady::SimulationControlAdjointSteady(
  const Parameters::SimulationControl &param)
  : SimulationControlTransient(param)
{}

void
SimulationControlAdjointSteady::print_progression(
  const ConditionalOStream &pcout)
{
  if (!is_verbose_iteration())
    return;

  pcout << std::endl;

  std::stringstream ss;
  ss << std::setprecision(this->log_precision);

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

SimulationControlRayTracing::SimulationControlRayTracing(
  const Parameters::SimulationControl &param,
  Particles::ParticleHandler<3>       &input_photon_handler)
  : SimulationControl(param)
  , photon_handler(input_photon_handler)
{}

bool
SimulationControlRayTracing::integrate()
{
  if (!is_at_end())
    {
      iteration_number++;
      return true;
    }

  else
    return false;
}
bool
SimulationControlRayTracing::is_at_end()
{
  return photon_handler.n_global_particles() == 0;
}

void
SimulationControlRayTracing::print_progression(const ConditionalOStream &pcout)
{
  pcout << "Iteration : " << iteration_number << std::endl
        << "Remaining photon : " << photon_handler.n_global_particles()
        << std::endl
        << std::endl;
}
