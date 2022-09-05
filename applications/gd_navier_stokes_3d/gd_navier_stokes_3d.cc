#include "solvers/gd_navier_stokes.h"

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

      ParameterHandler        prm;
      SimulationParameters<3> NSparam;
      NSparam.declare(prm);
      // Parsing of the file
      prm.parse_input(argv[1]);
      NSparam.parse(prm);

      AssertThrow(NSparam.nitsche->number_solids == 0,
                  SolidWarning(NSparam.nitsche->number_solids,
                               "gd_navier_stokes_3d",
                               "gls_nitsche_navier_stokes_23 (2D solid) or "
                               "gls_nitsche_navier_stokes_33 (3D solid)"));

      GDNavierStokesSolver<3> problem_3d(NSparam);
      problem_3d.solve();
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
