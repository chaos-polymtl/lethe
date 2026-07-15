#include <core/parameters.h>

#include <solvers/euler_euler_sf.h>
#include <solvers/euler_void_fraction.h>
#include <solvers/solid_phase.h>

#include <fem-dem/cfd_dem_simulation_parameters.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

#include <exception>
#include <fstream>
#include <iostream>

using namespace dealii;

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_init(argc,
                                                argv,
                                                numbers::invalid_unsigned_int);

      const MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
      const unsigned int rank =
        Utilities::MPI::this_mpi_process(mpi_communicator);

      constexpr int dim = 2;

      const bool print_prm =
        (argc > 1 && std::string(argv[1]) == "--print-prm");

      const std::string prm_file =
        print_prm ? (argc > 2 ? std::string(argv[2]) : "euler_euler.prm") :
                    (argc > 1 ? std::string(argv[1]) : "euler_euler.prm");

      ParameterHandler prm;

      SolidPhaseParameters<dim>::declare_parameters(prm);

      CFDDEMSimulationParameters<dim> fluid_parameters;

      Parameters::SizeOfSubsections size_of_subsections;
      size_of_subsections.boundary_conditions = 5;
      size_of_subsections.manifolds           = 1;

      fluid_parameters.declare(prm, size_of_subsections);

      if (print_prm)
        {
          if (rank == 0)
            {
              std::ofstream out(prm_file);
              prm.print_parameters(out, ParameterHandler::Text);
              std::cout << "Wrote default prm: " << prm_file << std::endl;
            }

          MPI_Barrier(mpi_communicator);
          return 0;
        }

      bool prm_exists = false;

      if (rank == 0)
        {
          std::ifstream in(prm_file);
          prm_exists = in.good();
        }

      prm_exists = Utilities::MPI::broadcast(mpi_communicator, prm_exists, 0);

      if (!prm_exists)
        {
          if (rank == 0)
            {
              std::ofstream out(prm_file);
              prm.print_parameters(out, ParameterHandler::Text);
              std::cout << "Parameter file not found. Wrote default: "
                        << prm_file << "\nEdit it and rerun.\n";
            }

          MPI_Barrier(mpi_communicator);
          return 0;
        }

      if (rank == 0)
        std::cout << "Reading parameters from: " << prm_file << std::endl;

      prm.parse_input(prm_file);

      SolidPhaseParameters<dim> solid_parameters;
      solid_parameters.parse_parameters(prm);

      fluid_parameters.parse(prm);

      EulerEulerOneWay<dim> problem(fluid_parameters, solid_parameters);

      problem.solve();
    }
  catch (std::exception &exc)
    {
      const unsigned int rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

      std::cerr << "\n[rank " << rank << "] Exception:\n" << exc.what() << "\n";
      return 1;
    }
  catch (...)
    {
      const unsigned int rank =
        Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

      std::cerr << "\n[rank " << rank << "] Unknown exception!\n";
      return 1;
    }

  return 0;
}