// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test to check the number of stages returned by the method
 * SimulationControl::get_number_of_stages()
 */

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  using Method = Parameters::SimulationControl::TimeSteppingMethod;

  std::vector<std::pair<Method, unsigned int>> test_cases = {
    {Method::steady, 1},
    {Method::steady_bdf, 1},
    {Method::bdf1, 1},
    {Method::bdf2, 1},
    {Method::bdf3, 1},
    {Method::sdirk22, 2},
    {Method::sdirk33, 3},
    {Method::sdirk43, 3}};

  for (const auto &[method, expected] : test_cases)
    {
      unsigned int result = SimulationControl::get_number_of_stages(method);
      deallog << "Method: " << static_cast<int>(method)
              << " -> Stages: " << result << std::endl;
      AssertThrow(result == expected,
                  ExcMessage(
                    "get_number_of_stages() returned incorrect value"));
    }
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
                << "Aborting!" << std::endl;
      return 1;
    }

  return 0;
}
