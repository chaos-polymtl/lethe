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

class SimulationFlowControl
{
protected:
  // Time
  double time;

  // Time
  double time_step;

  // Simulation end time
  double end_time;

  // Time step
  std::vector<double> time_step_vector;

  // Iteration
  unsigned int iteration;

  // Number of mesh adaptation iteration
  unsigned int number_mesh_adapt;

  // CFL
  double CFL;

  // Maximal CFL condition
  double max_CFL;

  // number of time steps stored
  static const unsigned int numberTimeStepStored = 4;



  // Output iteration frequency
  unsigned int output_frequency;

  // Subdivision
  unsigned int subdivision;

  // Number of parallel file to generate
  unsigned int group_files;

  // Output name
  std::string output_name;

  // Output path
  std::string output_path;



public:
  SimulationFlowControl(Parameters::SimulationControl param);

  virtual bool
  integrate() = 0;

  virtual bool
  is_at_end() = 0;

  // Add a time step and stores the previous one in a list
  void
  add_time_step(double p_timestep);

  bool
  is_at_start()
  {
    return iteration == 1;
  }

  virtual double
  calculate_time_step()
  {
    return time_step;
  };


  /**
   * @brief print_progress Function that prints the current progress status of the simulation
   * @param pcout the ConditionalOSStream that is use to write
   */
  virtual void
  print_progression(ConditionalOStream &pcout) = 0;

  void
  set_CFL(const double p_CFL)
  {
    CFL = p_CFL;
  }


  void
  set_desired_time_step(const double new_time_step)
  {
    time_step = new_time_step;
  }

  double
  get_time_step() const
  {
    return time_step;
  }

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

  double
  get_current_time()
  {
    return time;
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

  unsigned int
  get_step_number()
  {
    return iteration;
  }

  unsigned int
  get_number_subdivision()
  {
    return subdivision;
  }

  bool
  is_output_iteration()
  {
    return (get_step_number() % output_frequency == 0);
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

  virtual bool
  integrate() override;

  virtual bool
  is_at_end() override;
};

class SimulationControlSteady : public SimulationFlowControl
{
public:
  SimulationControlSteady(Parameters::SimulationControl param);

  virtual void
  print_progression(ConditionalOStream &pcout) override;

  virtual bool
  integrate() override;

  virtual bool
  is_at_end() override;
};

#endif
