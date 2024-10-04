// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the coefficients required for BDF integration
 * of order n.
 */

// Lethe
#include <core/bdf.h>
#include <core/parameters.h>

// Tests (with common definitions)
#include <../tests/tests.h>


void
test()
{
  std::vector<double> dt(5, 0.1);
  dt[1] = 0.2;
  dt[2] = 0.3;
  dt[3] = 0.4;
  dt[4] = 0.5;
  deallog << "Time steps ";
  for (const auto &dt_val : dt)
    {
      deallog << dt_val << " ";
    }

  deallog << std::endl;

  Vector<double> order1_coefficients = calculate_bdf_coefficients(
    Parameters::SimulationControl::TimeSteppingMethod::bdf1, dt);
  deallog << "Order 1 : " << order1_coefficients[0] << " "
          << order1_coefficients[1] << std::endl;

  Vector<double> order2_coefficients = calculate_bdf_coefficients(
    Parameters::SimulationControl::TimeSteppingMethod::bdf2, dt);
  deallog << "Order 2 : " << order2_coefficients[0] << " "
          << order2_coefficients[1] << " " << order2_coefficients[2]
          << std::endl;

  Vector<double> order3_coefficients = calculate_bdf_coefficients(
    Parameters::SimulationControl::TimeSteppingMethod::bdf3, dt);
  deallog << "Order 3 : ";
  for (const auto &coefficient : order3_coefficients)
    {
      deallog << coefficient << " ";
    }

  deallog << std::endl;
}

int
main(int argc, char *argv[])
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
