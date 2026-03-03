#include <solvers/solid_phase.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

#include <exception>
#include <fstream>
#include <iostream>

int
main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                                argv,
                                                numbers::invalid_unsigned_int);
      const MPI_Comm                   comm = MPI_COMM_WORLD;
      const unsigned int rank = Utilities::MPI::this_mpi_process(comm);


      ParameterHandler prm;
      SolidPhaseParameters::declare_parameters(prm);


      const bool print_prm =
        (argc > 1 && std::string(argv[1]) == "--print-prm");


      const std::string prm_file =
        (print_prm ? (argc > 2 ? std::string(argv[2]) : "solid_phase.prm") :
                     (argc > 1 ? std::string(argv[1]) : "solid_phase.prm"));


      if (print_prm)
        {
          if (rank == 0)
            {
              std::ofstream out(prm_file);
              prm.print_parameters(out, ParameterHandler::Text);
              std::cout << "Wrote default prm: " << prm_file << std::endl;
            }
          MPI_Barrier(comm);
          return 0;
        }


      bool prm_exists = false;
      if (rank == 0)
        {
          std::ifstream in(prm_file);
          prm_exists = in.good();

          if (!prm_exists)
            {
              std::ofstream out(prm_file);
              prm.print_parameters(out, ParameterHandler::Text);
              std::cout << "Parameter file not found. Wrote default: "
                        << prm_file << "\nEdit it and rerun.\n";
            }
        }

      prm_exists = Utilities::MPI::broadcast(comm, prm_exists, 0);

      MPI_Barrier(comm);
      if (!prm_exists)
        return 0;


      if (rank == 0)
        std::cout << "Reading parameters from: " << prm_file << std::endl;

      prm.parse_input(prm_file);

      SolidPhaseParameters parameters;
      parameters.parse_parameters(prm);


      SolidPhaseSolver<3> solver(parameters, comm);
      solver.run();
    }
  catch (std::exception &exc)
    {
      const unsigned int rank =
        dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

      std::cerr << "\n[rank " << rank << "] Exception:\n" << exc.what() << "\n";
      return 1;
    }
  catch (...)
    {
      const unsigned int rank =
        dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

      std::cerr << "\n[rank " << rank << "] Unknown exception!\n";
      return 1;
    }

  return 0;
}
