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
 * ---------------------------------------------------------------------

*
* Author: Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2019-
*/

#include "solvers/gls_navier_stokes.h"
#include "core/grids.h"
#include "solvers/solid_base.h"

int
main(int argc, char *argv[])
{
  try
    {
      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

      ParameterHandler                prm;
      NavierStokesSolverParameters<3> param;
      param.declare(prm);
      // Parsing of the file
      prm.parse_input(argv[1]);
      param.parse(prm);

      read_mesh_and_manifold(fluid_tria, param.mesh, param.manifolds_parameters, param.boundary_conditions);

      SolidBase<3,3> solid(param.nitsche,
                          fluid_tria,
                          param.fem_parameters.velocityOrder);
      

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
