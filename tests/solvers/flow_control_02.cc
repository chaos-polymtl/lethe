// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the restart of the flow control: a flow controller
 * read from the checkpoint file of another flow controller must calculate the
 * same beta force at every following time step as the uninterrupted one.
 */

// Lethe
#include <core/parameters.h>

#include <solvers/flow_control.h>

// Tests
#include <../tests/tests.h>

// Standard library
#include <cmath>
#include <string>

/**
 * @brief Calculate a fake new average velocity with 1% of pressure drop and an
 * arbitrary velocity increase of 25% of the beta force, as in flow_control_01.
 *
 * @param[in] average_velocity Average velocity of the last time step.
 *
 * @param[in] beta Beta force applied during the time step.
 *
 * @return Average velocity at the end of the time step.
 */
double
advance_fake_flow(const double average_velocity, const double beta)
{
  return 0.99 * average_velocity + 0.25 * beta;
}

void
test()
{
  Parameters::DynamicFlowControl flow_control_parameters;
  flow_control_parameters.enable_flow_control = true;
  flow_control_parameters.flow_direction      = 0;
  flow_control_parameters.average_velocity_0  = 10;
  flow_control_parameters.beta_0              = 20;
  flow_control_parameters.alpha               = 0.5;
  flow_control_parameters.beta_threshold      = 0;

  const double       time_step              = 0.1;
  const unsigned int checkpoint_step_number = 12;
  const unsigned int final_step_number      = 20;
  const unsigned int flow_direction    = flow_control_parameters.flow_direction;
  const std::string  checkpoint_prefix = "flow_control_02";

  // A checkpoint written with 6 significant digits stores the average velocity
  // 10.010968... as 10.011. The main controller multiplies this error by
  // alpha / time_step, which shifts the beta force by about 4e-4 relative to
  // its value, far above this tolerance.
  constexpr double restart_relative_tolerance = 1e-12;

  // Run the flow controller up to the checkpoint. As in the solvers, the first
  // time step is step 1.
  FlowControl<3> flow_control(flow_control_parameters);
  double         average_velocity = 0.0;
  for (unsigned int step_number = 1; step_number <= checkpoint_step_number;
       ++step_number)
    {
      flow_control.calculate_beta(average_velocity, time_step, step_number);
      average_velocity =
        advance_fake_flow(average_velocity,
                          flow_control.get_beta()[flow_direction]);
    }
  flow_control.save(checkpoint_prefix);

  // Restart a second flow controller from the checkpoint
  FlowControl<3> restarted_flow_control(flow_control_parameters);
  restarted_flow_control.read(checkpoint_prefix);

  deallog << "Checkpoint written at step " << checkpoint_step_number
          << std::endl;

  // Both flow controllers see the same average velocity, as a restarted
  // simulation reads the exact solution of its checkpoint.
  bool restart_matches_uninterrupted_run = true;
  for (unsigned int step_number = checkpoint_step_number + 1;
       step_number <= final_step_number;
       ++step_number)
    {
      flow_control.calculate_beta(average_velocity, time_step, step_number);
      restarted_flow_control.calculate_beta(average_velocity,
                                            time_step,
                                            step_number);

      const double beta = flow_control.get_beta()[flow_direction];
      const double restarted_beta =
        restarted_flow_control.get_beta()[flow_direction];

      if (std::abs(restarted_beta - beta) >
          restart_relative_tolerance * std::abs(beta))
        restart_matches_uninterrupted_run = false;

      deallog << "Step number: " << step_number << std::endl;
      deallog << " Beta uninterrupted: " << beta << std::endl;
      deallog << " Beta restarted:     " << restarted_beta << std::endl;

      average_velocity = advance_fake_flow(average_velocity, beta);
    }

  deallog << "Restart matches the uninterrupted run: "
          << (restart_matches_uninterrupted_run ? "true" : "false")
          << std::endl;
}

int
main()
{
  try
    {
      initlog();
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
