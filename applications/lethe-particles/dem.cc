// SPDX-FileCopyrightText: Copyright (c) 2020, 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "dem/dem.h"

#include "core/dem_properties.h"

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
          ParameterHandler       prm;
          DEMSolverParameters<2> dem_parameters;
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(file_name);
          dem_parameters.parse(prm);

          // Remove old output files
          if (options["-R"])
            {
              std::string output_path =
                dem_parameters.simulation_control.output_folder;
              delete_vtu_and_pvd_files(output_path);
            }

          // Print parameters if needed
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_parameters_to_output_file(pcout, prm, file_name);

          const DEM::SolverType solver_type =
            dem_parameters.model_parameters.solver_type;

          if (solver_type == DEM::SolverType::dem)
            {
              DEMSolver<2, DEM::DEMProperties::PropertiesIndex> problem(
                dem_parameters);
              problem.solve();
            }
          else if (solver_type == DEM::SolverType::dem_mp)
            {
              DEMSolver<2, DEM::DEMMPProperties::PropertiesIndex> problem(
                dem_parameters);
              problem.solve();
            }
          else
            {
              AssertThrow(
                false,
                dealii::ExcMessage(
                  "While reading the solver type from the input file, "
                  "Lethe found a value different than \"dem\" or \"dem_mp\"."
                  "lethe-particles does not support other solver type."));
            }
        }

      else if (dim == 3)
        {
          ParameterHandler       prm;
          DEMSolverParameters<3> dem_parameters;
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(file_name);
          dem_parameters.parse(prm);

          // Remove old output files
          if (options["-R"])
            {
              std::string output_path =
                dem_parameters.simulation_control.output_folder;
              delete_vtu_and_pvd_files(output_path);
            }

          // Print parameters if needed
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            print_parameters_to_output_file(pcout, prm, file_name);

          const DEM::SolverType solver_type =
            dem_parameters.model_parameters.solver_type;
          if (solver_type == DEM::SolverType::dem)
            {
              DEMSolver<3, DEM::DEMProperties::PropertiesIndex> problem(
                dem_parameters);
              problem.solve();
            }
          else if (solver_type == DEM::SolverType::dem_mp)
            {
              DEMSolver<3, DEM::DEMMPProperties::PropertiesIndex> problem(
                dem_parameters);
              problem.solve();
            }
          else
            {
              AssertThrow(
                false,
                dealii::ExcMessage(
                  "While reading the solver type from the input file, "
                  "Lethe found a value different than \"dem\" or \"dem_mp\"."));
            }
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
