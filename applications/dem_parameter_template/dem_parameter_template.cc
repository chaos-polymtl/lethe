// SPDX-FileCopyrightText: Copyright (c) 2020-2021, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// check the read and write of simulationcontrol

#include "dem/dem_solver_parameters.h"

#include <fstream>

int
main()
{
  try
    {
      {
        ParameterHandler       prm;
        DEMSolverParameters<2> dem_parameters;
        dem_parameters.declare(prm);
        std::ofstream output_prm("dem-2d.prm");
        prm.print_parameters(output_prm, prm.PRM);
      }
      {
        ParameterHandler       prm;
        DEMSolverParameters<3> dem_parameters;
        dem_parameters.declare(prm);
        std::ofstream output_prm("dem-3d.prm");
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
