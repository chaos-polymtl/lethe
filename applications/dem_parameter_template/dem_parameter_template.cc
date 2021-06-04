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
        prm.print_parameters(output_prm, prm.Text);

        std::ofstream output_xml("dem-2d.xml");
        prm.print_parameters(output_xml, prm.XML);
      }
      {
        ParameterHandler       prm;
        DEMSolverParameters<3> dem_parameters;
        dem_parameters.declare(prm);
        std::ofstream output_prm("dem-3d.prm");
        prm.print_parameters(output_prm, prm.Text);

        std::ofstream output_xml("dem-3d.xml");
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
