// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/fluid_dynamics_block.h>

#include <deal.II/base/convergence_table.h>

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }

      const unsigned int                  dim = get_dimension(argv[1]);
      const Parameters::SizeOfSubsections size_of_subsections =
        Parameters::get_size_of_subsections(argv[1]);

      if (dim == 2)
        {
          ParameterHandler        prm;
          SimulationParameters<2> NSparam;
          NSparam.declare(prm, size_of_subsections);
          // Parsing of the file
          prm.parse_input(argv[1]);
          NSparam.parse(prm);

          AssertThrow(NSparam.nitsche->number_solids == 0,
                      SolidWarning(NSparam.nitsche->number_solids,
                                   "lethe-fluid-block",
                                   "lethe-fluid-nitsche"));

          FluidDynamicsBlock<2> problem(NSparam);
          problem.solve();
        }

      else if (dim == 3)
        {
          ParameterHandler        prm;
          SimulationParameters<3> NSparam;
          NSparam.declare(prm, size_of_subsections);
          // Parsing of the file
          prm.parse_input(argv[1]);
          NSparam.parse(prm);

          AssertThrow(NSparam.nitsche->number_solids == 0,
                      SolidWarning(NSparam.nitsche->number_solids,
                                   "lethe-fluid-block",
                                   "lethe-fluid-nitsche"));

          FluidDynamicsBlock<3> problem(NSparam);
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
