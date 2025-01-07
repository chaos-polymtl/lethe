// SPDX-FileCopyrightText: Copyright (c) 2020, 2023-2024 The Lethe Authors
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


      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }

      const unsigned int dim = get_dimension(argv[1]);

      if (dim == 2)
        {
          ParameterHandler       prm;
          DEMSolverParameters<2> dem_parameters;
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(argv[1]);
          dem_parameters.parse(prm);
          const DEM::SolverType solver_type =
            dem_parameters.model_parameters.solver_type;

          if (solver_type == DEM::SolverType::dem)
            {
              DEMSolver<2, DEM::DEMProperties::PropertiesIndex> problem(dem_parameters);
              problem.solve();
            }
          else
            {
              AssertThrow(
                false,
                dealii::ExcMessage(
                  "While reading the solver type from the input file, "
                  "Lethe found a value different than \"dem\". As of January 2025, "
                  "the lethe-particles application requires the uses of  "
                  "\"solver type = dem\", which is the default value."));
            }
        }

      else if (dim == 3)
        {
          ParameterHandler       prm;
          DEMSolverParameters<3> dem_parameters;
          dem_parameters.declare(prm);

          // Parsing of the file
          prm.parse_input(argv[1]);
          dem_parameters.parse(prm);
          // const DEM::SolverType solver_type =
          // dem_parameters.model_parameters.solver_type;
          const DEM::SolverType solver_type = DEM::SolverType::dem;
          if (solver_type == DEM::SolverType::dem)
            {
              DEMSolver<3, DEM::DEMProperties::PropertiesIndex> problem(dem_parameters);
              problem.solve();
            }
          else
            {
              AssertThrow(
                false,
                dealii::ExcMessage(
                  "While reading the solver type from the input file, "
                  "Lethe found a value different than \"dem\". As of January 2025, "
                  "the lethe-particles application requires the uses of  "
                  "\"solver type = dem\", which is the default value."));
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
