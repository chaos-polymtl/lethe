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
#if DEAL_II_VERSION_GTE(9, 7, 0)
        prm.print_parameters(output_prm, prm.DefaultStyle);
#else
        prm.print_parameters(output_prm, prm.Text);
#endif
      }
      {
        ParameterHandler        prm;
        SimulationParameters<3> nsparam;

        nsparam.declare(prm, size_of_subsections);
        std::ofstream output_prm("template-3d.prm");
#if DEAL_II_VERSION_GTE(9, 7, 0)
        prm.print_parameters(output_prm, prm.DefaultStyle);
#else
        prm.print_parameters(output_prm, prm.Text);
#endif
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
