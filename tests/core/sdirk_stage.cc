// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/sdirk_stage_data.h>
#include <core/sdirk_table.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Testing SDIRK2 coefficients" << std::endl;

  SDIRKTable         table    = sdirk_table("SDIRK3");
  const unsigned int n_stages = table.A.m();

  // Data printed at each stage
  for (unsigned int stage_i = 0; stage_i < n_stages; ++stage_i)
    {
      SDIRKStageData data =
        sdirk_stage_data(table.A, table.c, table.b, stage_i);

      deallog << "\nStage " << stage_i << ":" << std::endl;
      deallog << "  a_ij: ";
      for (const auto &a : data.a_ij)
        deallog << std::setw(12) << a << " ";
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
