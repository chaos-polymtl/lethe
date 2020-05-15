/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_simulation_flow_control_h
#define lethe_simulation_flow_control_h

#include <deal.II/base/discrete_time.h>

#include <core/parameters.h>

class SimulationFlowControl : public DiscreteTime
{
protected:
  // Time step
  std::vector<double> time_step_vector;
  // CFL
  double CFL;

  // Maximal CFL condition
  double max_CFL;

  // Time
  double time;

  // Simulation end time
  double end_time;

  // Iteration number
  unsigned int iter;

  // Number of mesh adaptation iteration
  unsigned int number_mesh_adapt;

  // number of time steps stored
  static const unsigned int numberTimeStepStored = 4;

  // Number of parallel file to generate
  unsigned int group_files;

  // Output iteration frequency
  unsigned int output_frequency;

  // Subdivision
  unsigned int subdivision;

  // Output name
  std::string output_name;

  // Output path
  std::string output_path;

  // Add a time step and stores the previous one in a list
  void
  addTimeStep(double p_timestep);

public:
  SimulationFlowControl(Parameters::SimulationControl param,
                        double                        p_start_time,
                        double                        p_end_time,
                        double                        p_step);


  /**
   * @brief print_progress Function that prints the current progress status of the simulation
   * @param pcout the ConditionalOSStream that is use to write
   */
  virtual void
  print_progression(ConditionalOStream &pcout) = 0;


  std::string
  get_output_name()
  {
    return output_name;
  }
  std::string
  get_output_path()
  {
    return output_path;
  }

  unsigned int
  get_group_files()
  {
    return group_files;
  }

  void
  set_time_step(double p_timestep)
  {
    addTimeStep(p_timestep);
  }
  double
  getCurrentTimeStep()
  {
    return time_step_vector[0];
  }
  std::vector<double>
  get_time_steps_vector()
  {
    return time_step_vector;
  }

  double
  get_CFL()
  {
    return CFL;
  }
  void
  set_CFL(double p_CFL)
  {
    CFL = p_CFL;
  }

  double
  get_max_CFL()
  {
    return max_CFL;
  }

  unsigned int
  get_number_subdivision()
  {
    return subdivision;
  }

  bool
  is_output_iteration()
  {
    return (iter % output_frequency == 0);
  }

  void
  save(std::string filename);
  void
  read(std::string filename);
};


class SimulationControlTransient : public SimulationFlowControl
{
public:
  SimulationControlTransient(Parameters::SimulationControl param);

  virtual void
  print_progression(ConditionalOStream &pcout) override;
};

class SimulationControlSteady : public SimulationFlowControl
{
public:
  SimulationControlSteady(Parameters::SimulationControl param);

  virtual void
  print_progression(ConditionalOStream &pcout) override;
};

#endif
