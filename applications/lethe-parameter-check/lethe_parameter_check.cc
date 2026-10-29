// SPDX-FileCopyrightText: Copyright (c) 2023-2024, 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Executable used to test the validity of parameter files. It can also print
// the schema of the parameters of a solver family.


#include "solvers/simulation_parameters.h"

#include "dem/dem_solver_parameters.h"
#include "fem-dem/cfd_dem_simulation_parameters.h"

#include <core/parameters.h>
#include <core/utilities.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

namespace
{
  /**
   * @brief Print the usage of the application.
   *
   * @param[in] pcout Parallel output stream.
   * @param[in] application Name of the application.
   */
  void
  print_usage(const ConditionalOStream &pcout, const std::string &application)
  {
    pcout << "Usage:" << std::endl
          << "  " << application << " input_file solver_family" << std::endl
          << "      Check the validity of a parameter file." << std::endl
          << "  " << application << " --schema solver_family dimension"
          << std::endl
          << "      Print the schema of the parameters in the JSON format."
          << std::endl
          << "Solver family include: lethe-fluid, lethe-particles, "
             "lethe-fluid-particles"
          << std::endl;
  }

  /**
   * @brief Declare the parameters of a solver family and print their schema.
   *
   * A single subsection of each variable size subsection (bc 0 and
   * manifold 0) is declared. It is the template of the bc N and manifold N
   * subsections.
   *
   * @tparam ParametersType Parameters of the solver family.
   *
   * @param[in] solver_family Solver family.
   * @param[in] dim Dimension of the parameters.
   */
  template <typename ParametersType>
  void
  print_schema(const std::string &solver_family, const unsigned int dim)
  {
    Parameters::SizeOfSubsections size_of_subsections;
    size_of_subsections.boundary_conditions = 1;
    size_of_subsections.manifolds           = 1;

    ParameterHandler prm;
    ParametersType   parameters;
    parameters.declare(prm, size_of_subsections);

    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      print_parameter_schema(std::cout, prm, solver_family, dim);
  }

  /**
   * @brief Declare and parse the parameters of a parameter file. Parsing
   * throws if the parameter file is not valid.
   *
   * @tparam ParametersType Parameters of the solver family.
   *
   * @param[in] input_file Parameter file.
   * @param[in] size_of_subsections Size of the variable size subsections.
   */
  template <typename ParametersType>
  void
  check_parameter_file(const std::string                   &input_file,
                       const Parameters::SizeOfSubsections &size_of_subsections)
  {
    ParameterHandler prm;
    ParametersType   parameters;
    parameters.declare(prm, size_of_subsections);
    prm.parse_input(input_file);
    parameters.parse(prm);
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      const ConditionalOStream pcout(
        std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

      auto [options, args] = parse_args(argc, argv);

      const bool schema = options["--schema"];

      // A misspelled option must not be silently ignored
      const bool unknown_option =
        std::any_of(options.begin(), options.end(), [](const auto &option) {
          return option.first != "--schema";
        });

      if (unknown_option || args.size() != 2)
        {
          print_usage(pcout, argv[0]);
          return EXIT_FAILURE;
        }

      if (schema)
        {
          const std::string solver_family(args[0]);
          const int         dim = Utilities::string_to_int(args[1]);
          AssertThrow(dim == 2 || dim == 3,
                      ExcMessage("The dimension must be 2 or 3."));

          if (solver_family == "lethe-fluid")
            {
              if (dim == 2)
                print_schema<SimulationParameters<2>>(solver_family, dim);
              else
                print_schema<SimulationParameters<3>>(solver_family, dim);
            }
          else if (solver_family == "lethe-particles")
            {
              if (dim == 2)
                print_schema<DEMSolverParameters<2>>(solver_family, dim);
              else
                print_schema<DEMSolverParameters<3>>(solver_family, dim);
            }
          else if (solver_family == "lethe-fluid-particles")
            {
              if (dim == 2)
                print_schema<CFDDEMSimulationParameters<2>>(solver_family,
                                                            dim);
              else
                print_schema<CFDDEMSimulationParameters<3>>(solver_family,
                                                            dim);
            }
          else
            {
              print_usage(pcout, argv[0]);
              return EXIT_FAILURE;
            }

          return EXIT_SUCCESS;
        }

      const std::string  input_file(args[0]);
      const std::string  solver_family(args[1]);
      const unsigned int dim = get_dimension(input_file);
      AssertThrow(dim == 2 || dim == 3,
                  ExcMessage("The dimension must be 2 or 3."));

      if (solver_family == "lethe-fluid")
        {
          const Parameters::SizeOfSubsections size_of_subsections =
            Parameters::get_size_of_subsections(input_file);
          if (dim == 2)
            check_parameter_file<SimulationParameters<2>>(input_file,
                                                          size_of_subsections);
          else
            check_parameter_file<SimulationParameters<3>>(input_file,
                                                          size_of_subsections);
        }
      else if (solver_family == "lethe-particles")
        {
          // The DEM solver does not require any variable size subsection. Its
          // boundary conditions are declared in the "DEM boundary conditions"
          // subsection, and a file that declares no manifold leads to
          // subsections of size zero.
          const Parameters::SizeOfSubsections size_of_subsections =
            Parameters::get_size_of_subsections(input_file, false);
          if (dim == 2)
            check_parameter_file<DEMSolverParameters<2>>(input_file,
                                                         size_of_subsections);
          else
            check_parameter_file<DEMSolverParameters<3>>(input_file,
                                                         size_of_subsections);
        }
      else if (solver_family == "lethe-fluid-particles")
        {
          const Parameters::SizeOfSubsections size_of_subsections =
            Parameters::get_size_of_subsections(input_file);
          if (dim == 2)
            check_parameter_file<CFDDEMSimulationParameters<2>>(
              input_file, size_of_subsections);
          else
            check_parameter_file<CFDDEMSimulationParameters<3>>(
              input_file, size_of_subsections);
        }
      else
        {
          print_usage(pcout, argv[0]);
          return EXIT_FAILURE;
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
