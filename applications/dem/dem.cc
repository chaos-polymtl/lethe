#include "dem/dem.h"

#include "core/dem_properties.h"

using namespace dealii;

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);


      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }

      const unsigned int dim = get_dimension(argv[1]);

      if (dim == 2)
        {
          ParameterHandler       prm;
          DEMSolverParameters<2> dem_parameters;
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(argv[1]);
          dem_parameters.parse(prm);

          DEMSolver<2> problem(dem_parameters);
          problem.solve();
        }

      else if (dim == 3)
        {
          ParameterHandler       prm;
          DEMSolverParameters<3> dem_parameters;
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(argv[1]);
          dem_parameters.parse(prm);

          DEMSolver<3> problem(dem_parameters);
          problem.solve();
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
