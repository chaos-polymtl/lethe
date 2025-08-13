// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/sdirk_stage_data.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  // 12 digits of precision for the output
  deallog << std::setprecision(12) << std::scientific;

  deallog
    << "This test returns the values of the Butcher's table associated with each SDIRK method"
    << std::endl;

  deallog << "Testing SDIRK22 coefficients" << std::endl;

  // Important note : the nomenclature used for the name of the SDIRK methods
  // are sdirkOrderStage sdirk22 means SDIRK with order 2 and 2 stages, sdirk33
  // means SDIRK with order 3 and 3 stages.
  SDIRKTable table22 =
    sdirk_table(Parameters::SimulationControl::TimeSteppingMethod::sdirk22);
  const unsigned int n_stages22 = table22.A.m();

  // Data printed at each stage
  for (unsigned int stage_i = 1; stage_i <= n_stages22; ++stage_i)
    {
      SDIRKStageData data(table22, stage_i);

      deallog << "\nStage " << stage_i << ":" << std::endl;
      deallog << "  a_ij: ";
      for (const auto &a : data.a_ij)
        deallog << std::setw(12) << a << " ";
      deallog << "\n  b_i : " << data.b_i << std::endl;
      deallog << "\n  c_i : " << data.c_i << std::endl;
    }

  deallog << "<----------------------------->\n" << std::endl;

  deallog << "Testing SDIRK33 coefficients" << std::endl;

  SDIRKTable table33 =
    sdirk_table(Parameters::SimulationControl::TimeSteppingMethod::sdirk33);
  const unsigned int n_stages33 = table33.A.m();

  // Data printed at each stage
  for (unsigned int stage_i = 1; stage_i <= n_stages33; ++stage_i)
    {
      SDIRKStageData data(table33, stage_i);

      deallog << "\nStage " << stage_i << ":" << std::endl;
      deallog << "  a_ij: ";
      for (const auto &a : data.a_ij)
        deallog << std::setw(12) << a << " ";
      deallog << "\n  b_i : " << data.b_i << std::endl;
      deallog << "\n  c_i : " << data.c_i << std::endl;
    }

  deallog << "<----------------------------->\n" << std::endl;

  deallog << "Testing SDIRK43 coefficients" << std::endl;

  SDIRKTable table43 =
    sdirk_table(Parameters::SimulationControl::TimeSteppingMethod::sdirk43);
  const unsigned int n_stages43 = table43.A.m();

  // Data printed at each stage
  for (unsigned int stage_i = 1; stage_i <= n_stages43; ++stage_i)
    {
      SDIRKStageData data(table43, stage_i);

      deallog << "\nStage " << stage_i << ":" << std::endl;
      deallog << "  a_ij: ";
      for (const auto &a : data.a_ij)
        deallog << std::setw(12) << a << " ";
      deallog << "\n  b_i : " << data.b_i << std::endl;
      deallog << "\n  c_i : " << data.c_i << std::endl;
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
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
