#include <solvers/gd_navier_stokes.h>

#include <deal.II/base/convergence_table.h>

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
      SimulationParameters<2> NSparam;
      NSparam.declare(prm);
      // Parsing of the file
      prm.parse_input(argv[1]);
      NSparam.parse(prm);

      if (NSparam.nitsche->number_solids > 0)
        {
          std::cerr
            << "----------------------------------------------------"
            << std::endl
            << "Warning: you are using gd_navier_stokes_2d solver, but"
            << std::endl
            << "'number of solids' > 0 in 'subsection nitsche' (see documentation)."
            << std::endl
            << "To assemble the Nitsche restriction, use a gls_nitsche_navier_stokes solver instead."
            << std::endl
            << "----------------------------------------------------"
            << std::endl
            << std::endl;
        }

      GDNavierStokesSolver<2> problem_2d(NSparam);
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
