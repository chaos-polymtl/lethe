// check the read and write of simulationcontrol

#include <core/parameters.h>

#include "../tests.h"

#include "core/simulation_control.h"
#include "solvers/navier_stokes_solver_parameters.h"

int
main()
{
  try
    {
      initlog();

      Parameters::SimulationControl simulationControlParameters;
      simulationControlParameters.dt     = 0.01;
      simulationControlParameters.adapt  = false;
      simulationControlParameters.maxCFL = 99;
      simulationControlParameters.method =
        Parameters::SimulationControl::TimeSteppingMethod::bdf1;
      simulationControlParameters.timeEnd         = 999;
      simulationControlParameters.nbMeshAdapt     = 9;
      simulationControlParameters.output_name     = "test";
      simulationControlParameters.subdivision     = 7;
      simulationControlParameters.output_folder   = "canard";
      simulationControlParameters.outputFrequency = 8;

      SimulationControl simulationControl(simulationControlParameters);

      for (int i = 0; i < 10; ++i)
        simulationControl.integrate();

      simulationControl.save("testFile");
      simulationControl.read("testFile");

      deallog << "dt                  : "
              << simulationControl.getCurrentTimeStep() << std::endl;
      deallog << "CFL                 : " << simulationControl.getCFL()
              << std::endl;
      deallog << "time                : " << simulationControl.getTime()
              << std::endl;
      deallog << "iter                : " << simulationControl.getIter()
              << std::endl;
      deallog << "OK" << std::endl;
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
