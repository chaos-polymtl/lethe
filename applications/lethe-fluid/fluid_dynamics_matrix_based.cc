// SPDX-FileCopyrightText: Copyright (c) 2019, 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "solvers/fluid_dynamics_matrix_based.h"

#include <core/utilities.h>

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      ConditionalOStream pcout(
        std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

      auto [options, args] = parse_args(argc, argv);

      // Print version information
      if (options["-V"])
        {
          pcout << "Running: " << concatenate_strings(argc, argv) << std::endl;

          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_version_info(pcout);

          return EXIT_SUCCESS;
        }

      if (args.empty())
        {
          pcout << "Usage: " << argv[0] << " input_file" << std::endl;
          return EXIT_FAILURE;
        }

      const std::string file_name(args[0]);

      const unsigned int                  dim = get_dimension(file_name);
      const Parameters::SizeOfSubsections size_of_subsections =
        Parameters::get_size_of_subsections(file_name);

      if (dim == 2)
        {
          ParameterHandler        prm;
          SimulationParameters<2> NSparam;
          NSparam.declare(prm, size_of_subsections);

          // Parsing of the file
          prm.parse_input(file_name);
          NSparam.parse(prm);

          // Remove old output files
          if (options["-R"])
            {
              std::string output_path =
                NSparam.simulation_control.output_folder;
              delete_vtu_and_pvd_files(output_path);
            }

          AssertThrow(NSparam.nitsche->number_solids == 0,
                      SolidWarning(NSparam.nitsche->number_solids,
                                   "lethe-fluid",
                                   "lethe-fluid-nitsche"));

          // Print parameters if needed
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_parameters_to_output_file(pcout, prm, file_name);

          FluidDynamicsMatrixBased<2> problem(NSparam);
          problem.solve();
        }

      else if (dim == 3)
        {
          ParameterHandler        prm;
          SimulationParameters<3> NSparam;
          NSparam.declare(prm, size_of_subsections);
          // Parsing of the file
          prm.parse_input(file_name);
          NSparam.parse(prm);

          // Remove old output files
          if (options["-R"])
            {
              std::string output_path =
                NSparam.simulation_control.output_folder;
              delete_vtu_and_pvd_files(output_path);
            }

          AssertThrow(NSparam.nitsche->number_solids == 0,
                      SolidWarning(NSparam.nitsche->number_solids,
                                   "lethe-fluid",
                                   "lethe-fluid-nitsche"));

          // Print parameters if needed
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_parameters_to_output_file(pcout, prm, file_name);

          FluidDynamicsMatrixBased<3> problem(NSparam);
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
