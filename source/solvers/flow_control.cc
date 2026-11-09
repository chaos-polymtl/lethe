// SPDX-FileCopyrightText: Copyright (c) 2020, 2023-2024, 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <solvers/flow_control.h>

#include <deal.II/base/exceptions.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <utility>

template <int dim>
FlowControl<dim>::FlowControl(
  const Parameters::DynamicFlowControl &flow_control)
  : target_average_velocity(flow_control.average_velocity_0)
  , initial_beta(flow_control.beta_0)
  , alpha(flow_control.alpha)
  , beta_threshold(flow_control.beta_threshold)
  , flow_direction(flow_control.flow_direction)
{
  // The beta force is applied to the component flow_direction of a tensor of
  // dimension dim. Tensor::operator[] only checks its index in Debug mode, so
  // an invalid flow direction would write out of bounds in Release mode.
  AssertThrow(
    !flow_control.enable_flow_control || std::cmp_less(flow_direction, dim),
    ExcMessage("The 'flow direction' of the 'flow control' subsection is " +
               std::to_string(flow_direction) +
               ", but it must be smaller than the dimension of the problem (" +
               std::to_string(dim) + "). Set it to 0 for x, 1 for y" +
               (dim == 3 ? " or 2 for z." : ".")));
}

template <int dim>
double
FlowControl<dim>::proportional_flow_controller(const double &average_velocity_n,
                                               const double &dt) const
{
  return 0.5 * alpha * (target_average_velocity - average_velocity_n) / dt;
}

template <int dim>
double
FlowControl<dim>::main_flow_controller(const double &average_velocity_n,
                                       const double &dt) const
{
  const double new_beta =
    previous_beta + alpha *
                      (target_average_velocity - 2 * average_velocity_n +
                       previous_average_velocity) /
                      dt;

  // If desired average velocity is reached, new beta only maintains the force
  // to keep the flow at the desired value. If calculated beta is
  // negative it is set to 0 to avoid +/- force.
  if (target_average_velocity * new_beta < 0 &&
      !force_free_deceleration_is_allowed)
    return 0.0;

  return new_beta;
}

template <int dim>
void
FlowControl<dim>::calculate_beta(const double       &average_velocity_n,
                                 const double       &dt,
                                 const unsigned int &step_number)
{
  // Check if disabling the force-free deceleration and modify the force-free
  // velocity band (> step time 1):
  // If average velocity target is reached without any beta force applied,
  // the force-free deceleration is disabled. It means that the average
  // velocity has now reached the target value without force OR that slowing
  // down the flow with the only pressure drop (no force applied) is
  // ineffective (too big pressure drop). In case the average velocity is in
  // the initial band but the pressure drop is too small to slow down the
  // flow, the band factor is halved towards 1. Otherwise, the average velocity
  // may take a lot of time to reach the target value only with pressure drop.
  if (std::abs(previous_beta) < zero_beta_tolerance && step_number > 1)
    {
      if (std::abs(average_velocity_n) < std::abs(target_average_velocity))
        force_free_deceleration_is_allowed = false;
      else if (std::abs(previous_average_velocity - average_velocity_n) <
                 std::abs(average_velocity_n - target_average_velocity) &&
               force_free_velocity_band_factor >
                 minimum_force_free_velocity_band_factor)
        force_free_velocity_band_factor =
          0.5 * (1 + force_free_velocity_band_factor);
    }

  double new_beta = 0.0;
  if (step_number == 0)
    {
      // No force applied at first time step.
      // At this moment, calculate_beta is not called at this time step but
      // the condition prevents a beta calculation if the way it's called
      // changes.
      new_beta = 0.0;
    }
  else if (step_number == 1)
    {
      // Apply initial beta force if user defined (> than 0)
      // otherwise, it is calculated by the proportional controller
      new_beta = (std::abs(initial_beta) > zero_beta_tolerance ?
                    initial_beta :
                    proportional_flow_controller(average_velocity_n, dt));
    }
  else if (std::abs(average_velocity_n) > std::abs(target_average_velocity) &&
           std::abs(average_velocity_n) < force_free_velocity_band_factor *
                                            std::abs(target_average_velocity) &&
           force_free_deceleration_is_allowed)
    {
      // If the average velocity is between targeted value and the upper limit
      // of the band (meaning it is slightly over the value) and it has not
      // reached once the value during the simulation (force-free deceleration
      // is allowed), it decreases by itself (pressure drop).
      new_beta = 0.0;
    }
  else if (step_number == 2)
    {
      // Since the main controller takes the previous beta value, it is better
      // to use the proportional controller for this second time step since
      // it won't be affected by the first beta value that may cause overshoot
      // or undershoot.
      new_beta = proportional_flow_controller(average_velocity_n, dt);
    }
  else
    {
      // Main controller
      new_beta = main_flow_controller(average_velocity_n, dt);
    }

  // If the new beta is in the user defined threshold of the old beta, the new
  // beta is kept as the previous beta. This could allow avoiding unnecessary
  // reassembly of the matrix because of the force term when reuse matrix option
  // is enabled for the non-linear solver
  if (std::abs(new_beta - previous_beta) <
      std::abs(beta_threshold * previous_beta))
    new_beta = previous_beta;

  // Setting beta coefficient to the tensor according to the flow direction.
  beta[flow_direction] = new_beta;

  // Assigning values of this time step as previous values prior next
  // calculation.
  previous_beta             = new_beta;
  previous_average_velocity = average_velocity_n;
}

template <int dim>
void
FlowControl<dim>::save(const std::string &prefix) const
{
  const std::string filename = prefix + ".flowcontrol";
  std::ofstream     output(filename);

  // The state is written with enough digits to be read back exactly, so that a
  // restarted simulation continues as if it had not been interrupted.
  output << std::setprecision(std::numeric_limits<double>::max_digits10);
  output << "Flow control" << std::endl;
  output << "Previous_beta " << previous_beta << std::endl;
  output << "Previous_average_velocity " << previous_average_velocity
         << std::endl;
  output << "No_force " << force_free_deceleration_is_allowed << std::endl;
  output << "Threshold_factor " << force_free_velocity_band_factor << std::endl;
}

template <int dim>
void
FlowControl<dim>::read(const std::string &prefix)
{
  const std::string filename = prefix + ".flowcontrol";
  check_file_exists(filename,
                    "checkpoint file of the dynamic flow control, given by "
                    "'subsection restart' - 'set filename'");
  std::ifstream input(filename);

  std::string buffer;
  std::getline(input, buffer);
  input >> buffer >> previous_beta;
  input >> buffer >> previous_average_velocity;
  input >> buffer >> force_free_deceleration_is_allowed;
  input >> buffer >> force_free_velocity_band_factor;

  AssertThrow(!input.fail(),
              ExcMessage(
                "The checkpoint file of the dynamic flow control <" + filename +
                "> could not be read: it is truncated or does not have the "
                "expected format. Restart from a complete checkpoint written "
                "by a simulation with the flow control enabled."));
}

template class FlowControl<2>;
template class FlowControl<3>;
