// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "core/dem_properties.h"

#include "dem/ray_tracing.h"

using namespace dealii;

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

      const unsigned int dim = get_dimension(file_name);

      if (dim == 2)
        {
          ParameterHandler              prm;
          RayTracingSolverParameters<2> parameters;
          DEMSolverParameters<2>        dem_parameters;

          parameters.declare(prm);
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(file_name);
          parameters.parse(prm);
          dem_parameters.parse(prm);

          // Remove old output files
          if (options["-R"])
            {
              std::string output_path =
                parameters.simulation_control.output_folder;
              delete_vtu_and_pvd_files(output_path);
            }

          // Print parameters if needed
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_parameters_to_output_file(pcout, prm, file_name);

          const SolverType solver_type =
            parameters.model_parameters.solver_type;

          if (solver_type == dem)
            {
              RayTracingSolver<2> problem(parameters, dem_parameters);
              problem.solve();
            }
          else
            {
              AssertThrow(
                false,
                dealii::ExcMessage(
                  "While reading the solver type from the input file, "
                  "Lethe found a value different than \"dem\"."));
            }
        }

      else if (dim == 3)
        {
          ParameterHandler              prm;
          RayTracingSolverParameters<3> parameters;
          DEMSolverParameters<3>        dem_parameters;

          parameters.declare(prm);
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(file_name);
          parameters.parse(prm);
          dem_parameters.parse(prm);

          // Remove old output files
          if (options["-R"])
            {
              std::string output_path =
                parameters.simulation_control.output_folder;
              delete_vtu_and_pvd_files(output_path);
            }

          // Print parameters if needed
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_parameters_to_output_file(pcout, prm, file_name);

          RayTracingSolver<3> problem(parameters, dem_parameters);
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
