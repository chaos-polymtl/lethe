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
        prm.print_parameters(output_prm, prm.Text);

        std::ofstream output_xml("templa"
                                 "te-2d.xml");
        prm.print_parameters(output_xml, prm.XML);
      }
      {
        ParameterHandler        prm;
        SimulationParameters<3> nsparam;

        nsparam.declare(prm, size_of_subsections);
        std::ofstream output_prm("template-3d.prm");
        prm.print_parameters(output_prm, prm.Text);

        std::ofstream output_xml("template-3d.xml");
        prm.print_parameters(output_xml, prm.XML);
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
