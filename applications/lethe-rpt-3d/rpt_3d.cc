// SPDX-FileCopyrightText: Copyright (c) 2021 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <rpt/rpt.h>
#include <rpt/rpt_calculating_parameters.h>

#include <fstream>
#include <iostream>

using namespace dealii;

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      // Check the number of MPI processes
      int number_of_processes;
      MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
      AssertThrow(number_of_processes == 1,
                  ExcMessage(
                    "The rpt_3d application can only run with 1 MPI process."));

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

      ParameterHandler         prm;
      RPTCalculatingParameters rpt_parameters;
      rpt_parameters.declare(prm);

      // Parsing of the file
      prm.parse_input(file_name);
      rpt_parameters.parse(prm);

      // Print parameters if needed
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        print_parameters_to_output_file(pcout, prm, file_name);

      RPT<3> rpt(rpt_parameters);
      rpt.setup_and_calculate();
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
