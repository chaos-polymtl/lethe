// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "fem-dem/fluid_dynamics_vans_matrix_free.h"

#include <core/revision.h>
#include <core/utilities.h>

#include <deal.II/base/revision.h>

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

      if (print_parameters)
        {
          pcout << "Running: " << concatenate_strings(argc, argv) << std::endl;
          pcout << "  - deal.II (branch: " << DEAL_II_GIT_BRANCH
                << "; revision: " << DEAL_II_GIT_REVISION
                << "; short: " << DEAL_II_GIT_SHORTREV << ")" << std::endl;
          pcout << "  - Lethe (branch: " << LETHE_GIT_BRANCH
                << "; revision: " << LETHE_GIT_REVISION
                << "; short: " << LETHE_GIT_SHORTREV << ")" << std::endl;
          pcout << std::endl;
          pcout << std::endl;
        }

      if (dim == 2)
        {
          ParameterHandler              prm;
          CFDDEMSimulationParameters<2> NSparam;
          NSparam.declare(prm, size_of_subsections);
          // Parsing of the file
          prm.parse_input(file_name);
          NSparam.parse(prm);

          if (pcout.is_active())
            prm.print_parameters(pcout.get_stream(),
                                 ParameterHandler::OutputStyle::PRM |
                                   ParameterHandler::OutputStyle::Short |
                                   ParameterHandler::KeepDeclarationOrder
#if DEAL_II_VERSION_GTE(9, 7, 0)
                                   | ParameterHandler::KeepOnlyChanged
#endif
            );
          pcout << std::endl << std::endl;

          FluidDynamicsVANSMatrixFree<2> problem(NSparam);
          problem.solve();
        }

      else if (dim == 3)
        {
          ParameterHandler              prm;
          CFDDEMSimulationParameters<3> NSparam;
          NSparam.declare(prm, size_of_subsections);
          // Parsing of the file
          prm.parse_input(file_name);
          NSparam.parse(prm);

          if (pcout.is_active())
            prm.print_parameters(pcout.get_stream(),
                                 ParameterHandler::OutputStyle::PRM |
                                   ParameterHandler::OutputStyle::Short |
                                   ParameterHandler::KeepDeclarationOrder
#if DEAL_II_VERSION_GTE(9, 7, 0)
                                   | ParameterHandler::KeepOnlyChanged
#endif
            );
          pcout << std::endl << std::endl;

          FluidDynamicsVANSMatrixFree<3> problem(NSparam);
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
