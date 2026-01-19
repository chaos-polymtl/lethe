// SPDX-FileCopyrightText: Copyright (c) 2019, 2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// check the read and write of simulationcontrol

#include "solvers/simulation_parameters.h"

#include <fstream>

int
main()
{
  try
    {
      // Declare dummy size of subsections to declare a single of each entity
      Parameters::SizeOfSubsections size_of_subsections;
      size_of_subsections.boundary_conditions = 1;

      {
        ParameterHandler        prm;
        SimulationParameters<2> nsparam;

        nsparam.declare(prm, size_of_subsections);
        std::ofstream output_prm("template-2d.prm");
        prm.print_parameters(output_prm, prm.PRM);
      }
      {
        ParameterHandler        prm;
        SimulationParameters<3> nsparam;

        nsparam.declare(prm, size_of_subsections);
        std::ofstream output_prm("template-3d.prm");
        prm.print_parameters(output_prm, prm.PRM);
      }
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
