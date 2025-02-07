// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "solvers/fluid_dynamics_matrix_free.h"

#include <core/utilities.h>

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      if (argc == 1)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }

      const std::string file_name(argv[argc - 1]);
      const bool        print_parameters =
        (argc == 2) ? false : (std::string(argv[1]) == "--print-parameters");

      const unsigned int                  dim = get_dimension(file_name);
      const Parameters::SizeOfSubsections size_of_subsections =
        Parameters::get_size_of_subsections(file_name);

      ConditionalOStream pcout(std::cout,
                               (Utilities::MPI::this_mpi_process(
                                  MPI_COMM_WORLD) == 0) &&
                                 print_parameters);

      print_version_info(argc, argv, pcout);

      if (dim == 2)
        {
          ParameterHandler        prm;
          SimulationParameters<2> NSparam;
          NSparam.declare(prm, size_of_subsections);
          // Parsing of the file
          prm.parse_input(file_name);
          NSparam.parse(prm);

          print_parameters_to_output_file(pcout, prm);

          FluidDynamicsMatrixFree<2> problem(NSparam);
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

          print_parameters_to_output_file(pcout, prm);

          FluidDynamicsMatrixFree<3> problem(NSparam);
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
