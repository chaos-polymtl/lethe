/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 ---------------------------------------------------------------------*/

#include "solvers/fluid_dynamics_matrix_free.h"

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

      if (dim == 2)
        {
          ParameterHandler        prm;
          SimulationParameters<2> NSparam;
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
