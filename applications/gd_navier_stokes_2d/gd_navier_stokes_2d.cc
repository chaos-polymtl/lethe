#include <deal.II/base/convergence_table.h>

#include <solvers/gd_navier_stokes.h>

int
main(int argc, char *argv[])
{
  try
    {
      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

      ParameterHandler                prm;
      NavierStokesSolverParameters<2> NSparam;
      NSparam.declare(prm);
      // Parsing of the file
      prm.parse_input(argv[1]);
      NSparam.parse(prm);

      GDNavierStokesSolver<2> problem_2d(NSparam,
                                         NSparam.femParameters.velocityOrder,
                                         NSparam.femParameters.pressureOrder);
      problem_2d.solve();
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
