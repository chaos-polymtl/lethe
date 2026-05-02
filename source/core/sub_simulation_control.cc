// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/sub_simulation_control.h>


SubSimulationControlDEM::SubSimulationControlDEM(
  const DEMSubIterationLogic iteration_logic,
  const double               time_interval,
  const unsigned int         coupling_frequency,
  const double               rayleigh_characteristic_time,
  const double               fraction_of_rayleigh_characteristic_time)
  : iteration_logic(iteration_logic)
  , time_interval(time_interval)
  , coupling_frequency(coupling_frequency)
  , rayleigh_characteristic_time(rayleigh_characteristic_time)
  , fraction_of_rayleigh_characteristic_time(
      fraction_of_rayleigh_characteristic_time)
  , iteration_number(0)
{
  Assert(this->time_interval > 0,
         ExcMessage("The time interval must be strictly positive"));


  if (this->iteration_logic == DEMSubIterationLogic::fixed_number_of_iterations)
    {
      Assert(
        this->coupling_frequency > 0,
        ExcMessage(
          "The coupling frequency must be positive. A value of zero or a negative value was encountered. The simulation will stop."));

      // A fixed number of iterations is requested, the time step is fixed
      // using the coupling frequency and the time interval.
      time_step = this->time_interval / this->coupling_frequency;
      total_number_of_iterations = this->coupling_frequency;
    }
  else //(iteration_logic ==
       // DEMSubIterationLogic::fixed_fraction_of_rayleigh_time_step)
    {
      Assert(this->rayleigh_characteristic_time > 0,
             ExcMessage(
               "The Rayleigh characteristic time must be strictly positive"));
      Assert(
        this->fraction_of_rayleigh_characteristic_time > 0,
        ExcMessage(
          "The fraction of the Rayleigh characteristic time must be strictly positive"));

      // A fixed fraction of the Rayleigh time step is used. The number of
      // iterations will be manually calculated
      time_step = this->fraction_of_rayleigh_characteristic_time *
                  this->rayleigh_characteristic_time;

      // We calculate the total number of iterations with the proposed time step
      // We then take the ceiling of it to have a integer number of time steps
      total_number_of_iterations =
        static_cast<unsigned int>(std::ceil(this->time_interval / time_step));

      // Since the total number of iterations needs to be an integer, the actual
      // value of the time step may have changed as a result of the previous
      // calculations. We then recalculate the time step.
      time_step = this->time_interval / double(total_number_of_iterations);

      Assert(
        total_number_of_iterations > 0,
        ExcMessage(
          "The total number of iterations for the DEM sub simulation control is below 0."
          " Something went wrong with either the fraction of the Rayleigh characteristic time or the Rayleigh characteristic time provided."));
    }
}

bool
SubSimulationControlDEM::iterate()
{
  if (iteration_number == total_number_of_iterations)
    return false;

  iteration_number++;
  return true;
}
