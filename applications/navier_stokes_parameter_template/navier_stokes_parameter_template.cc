// check the read and write of simulationcontrol

#include "solvers/navier_stokes_solver_parameters.h"

#include <fstream>

int
main()
{
  try
    {
      {
        ParameterHandler                prm;
        NavierStokesSolverParameters<2> nsparam;
        nsparam.declare(prm);
        std::ofstream output_prm("template-2d.prm");
        prm.print_parameters(output_prm, prm.Text);

        std::ofstream output_xml("templa"
                                 "te-2d.xml");
        prm.print_parameters(output_xml, prm.XML);
      }
      {
        ParameterHandler                prm;
        NavierStokesSolverParameters<3> nsparam;
        nsparam.declare(prm);
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
