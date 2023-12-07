// Executable used to test the validity of parameter files


#include "solvers/simulation_parameters.h"

#include "dem/dem_solver_parameters.h"
#include "fem-dem/cfd_dem_simulation_parameters.h"

#include <fstream>

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

      if (argc != 3)
        {
          std::cout << "Usage:" << argv[0] << " input_file"
                    << "solver_family" << std::endl;
          std::cout
            << "Solver family include: lethe-fluid, lethe-particles, lethe-fluid-particles"
            << std::endl;
          std::exit(1);
        }

      const unsigned int dim = get_dimension(argv[1]);
      std::string        solver_family(argv[2]);



      if (solver_family == "lethe-fluid")
        {
          const Parameters::SizeOfSubsections size_of_subsections =
            Parameters::get_size_of_subsections(argv[1]);
          if (dim == 2)
            {
              ParameterHandler        prm;
              SimulationParameters<2> NSparam;
              NSparam.declare(prm, size_of_subsections);
              prm.parse_input(argv[1]);
            }
          if (dim == 3)
            {
              ParameterHandler        prm;
              SimulationParameters<3> NSparam;
              NSparam.declare(prm, size_of_subsections);
              prm.parse_input(argv[1]);
            }
        }
      else if (solver_family == "lethe-particles")
        {
          if (dim == 2)
            {
              ParameterHandler       prm;
              DEMSolverParameters<2> parameters;
              parameters.declare(prm);
              prm.parse_input(argv[1]);
            }
          if (dim == 3)
            {
              ParameterHandler       prm;
              DEMSolverParameters<3> parameters;
              parameters.declare(prm);
              prm.parse_input(argv[1]);
            }
        }
      else if (solver_family == "lethe-fluid-particles")
        {
          const Parameters::SizeOfSubsections size_of_subsections =
            Parameters::get_size_of_subsections(argv[1]);
          if (dim == 2)
            {
              ParameterHandler              prm;
              CFDDEMSimulationParameters<2> param;
              param.declare(prm, size_of_subsections);
              prm.parse_input(argv[1]);
            }

          else if (dim == 3)
            {
              ParameterHandler              prm;
              CFDDEMSimulationParameters<3> param;
              param.declare(prm, size_of_subsections);
              prm.parse_input(argv[1]);
            }
        }
      else
        {
          return 1;
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
  return 0;
}
