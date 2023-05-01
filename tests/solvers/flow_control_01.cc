/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
*/

/**
 * @brief This code tests the flow control algorithm which calculate a beta
 * force coefficient at each time step.
 */

// Lethe
#include <core/parameters.h>

#include <solvers/flow_control.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  Parameters::DynamicFlowControl flow_control_parameters;
  flow_control_parameters.flow_direction     = 0;
  flow_control_parameters.average_velocity_0 = 10;
  flow_control_parameters.beta_0             = 20;
  flow_control_parameters.alpha              = 0.5;

  double average_velocity = 0.0;
  double dt               = 0.1;

  FlowControl<3> flow_control(flow_control_parameters);

  Tensor<1, 3> beta;

  deallog << "*************************************************" << std::endl;
  deallog << " Time step:         " << dt << std::endl;
  deallog << " Initial beta:      " << flow_control_parameters.beta_0
          << std::endl;
  deallog << " Desired average velocity: "
          << flow_control_parameters.average_velocity_0 << std::endl;
  deallog << "*************************************************" << std::endl;

  for (int unsigned step_number = 0; step_number < 25; ++step_number)
    {
      // Calculating beta at each time step
      flow_control.calculate_beta(average_velocity, dt, step_number);
      beta = flow_control.get_beta();

      // Calculating a fake new average velocity with 1% of pressure drop and an
      // arbitrary velocity value of 25% of beta
      average_velocity *= 0.99;
      average_velocity += 0.25 * beta[flow_control_parameters.flow_direction];
      deallog << "" << std::endl;
      deallog << " Step number:  " << step_number << std::endl;
      deallog << " Average velocity:    " << average_velocity << std::endl;
      deallog << " Beta applied: "
              << beta[flow_control_parameters.flow_direction] << std::endl;
    }
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
