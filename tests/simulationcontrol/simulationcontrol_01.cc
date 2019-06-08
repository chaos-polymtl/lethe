// check the read and write of simulationcontrol

#include "simulationcontrol.h"
#include "parameters.h"
#include "navierstokessolverparameters.h"
#include "../tests.h"

int
main()
{
  try
  {

    initlog();

    Parameters::SimulationControl simulationControlParameters;
    simulationControlParameters.dt=0.01;
    simulationControlParameters.adapt=false;
    simulationControlParameters.maxCFL=99;
    simulationControlParameters.method=simulationControlParameters.backward;
    simulationControlParameters.timeEnd=999;
    simulationControlParameters.nbMeshAdapt=9;
    simulationControlParameters.output_name="test";
    simulationControlParameters.subdivision=7;
    simulationControlParameters.output_folder="canard";
    simulationControlParameters.outputFrequency=8;

    SimulationControl simulationControl;
    simulationControl.initialize(simulationControlParameters);

    for (int i=0 ; i <10 ; ++i)
      simulationControl.integrate();

    if(simulationControl.getIter()!=10)                                 std::runtime_error("Iteration number is wrong - before");
    if(!approximatelyEqual(simulationControl.getTime(),0.1,1e-10))      std::runtime_error("Run time is wrong - before");
    if(!approximatelyEqual(simulationControl.getCurrentTimeStep(),0.01,1e-10)) std::runtime_error("Time step is wrong - before");

    simulationControl.save("testFile");
    simulationControl.read("testFile");

    if(simulationControl.getIter()!=10)                                 std::runtime_error("Iteration number is wrong - after");
    if(!approximatelyEqual(simulationControl.getTime(),0.1,1e-10))      std::runtime_error("Run time is wrong - after");
    if(!approximatelyEqual(simulationControl.getCurrentTimeStep(),0.01,1e-10)) std::runtime_error("Time step is wrong - after");

    deallog << "dt                  : " << simulationControl.getCurrentTimeStep()    << std::endl;
    deallog << "CFL                 : " << simulationControl.getCFL()         << std::endl;
    deallog << "time                : " << simulationControl.getTime()        << std::endl;
    deallog << "iter                : " << simulationControl.getIter()        << std::endl;
    deallog << "OK" << std::endl;
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what()  << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

}
