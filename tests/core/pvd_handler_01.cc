// SPDX-FileCopyrightText: Copyright (c) 2019-2020 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Check the read and write of simulationcontrol
 */

// Lethe
#include <core/pvd_handler.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  deallog << "Beggining" << std::endl;

  PVDHandler pvdhandlerMaster;

  pvdhandlerMaster.append(0.03, "time_at_step_1");
  pvdhandlerMaster.append(2.00, "time_at_step_2");
  pvdhandlerMaster.append(7.00, "time_at_step_3");
  pvdhandlerMaster.save("restart");

  PVDHandler pvdhandlerWorker;
  pvdhandlerWorker.read("restart");

  if (pvdhandlerWorker.times_and_names.size() !=
      pvdhandlerMaster.times_and_names.size())
    throw std::runtime_error("Size are not equal");
  unsigned int size = pvdhandlerWorker.times_and_names.size();
  for (unsigned int i = 0; i < size; ++i)
    {
      deallog << pvdhandlerWorker.times_and_names[i].first << " "
              << pvdhandlerWorker.times_and_names[i].second << std::endl;
      if (!approximatelyEqual(pvdhandlerMaster.times_and_names[i].first,
                              pvdhandlerWorker.times_and_names[i].first,
                              1e-8))
        throw std::runtime_error("Time not equal");
      if (pvdhandlerMaster.times_and_names[i].second !=
          pvdhandlerWorker.times_and_names[i].second)
        throw std::runtime_error("File not equal");
    }
  deallog << "OK" << std::endl;
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
}
