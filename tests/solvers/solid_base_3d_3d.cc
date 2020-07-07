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
* Author: Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020-
*/

#include "deal.II/grid/grid_generator.h"
#include "solvers/solid_base.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

      MPI_Comm mpi_communicator = MPI_COMM_WOLRD;
      Parameters::Nitsche                                         param;
      std::shared_ptr<parallel::DistributedTriangulationBase<3>>  fluid_tria = std::make_shared<parallel::distributed::Triangulation<3>>(mpi_communicator,
            typename Triangulation<3>::MeshSmoothing(
              Triangulation<3>::smoothing_on_refinement |
              Triangulation<3>::smoothing_on_coarsening));

      // Mesh of the solid
      // param.solid_mesh.type = Parameters::Mesh::Type::gmsh;
      // param.solid_mesh.file_name = "50.msh";

      param.solid_mesh.type = Parameters::Mesh::Type::dealii;
      param.solid_mesh.grid_type = "hyper_shell";
      param.solid_mesh.grid_arguments = "0, 0 : 0 : 0.75 : 48 : true";

      // Mesh of the fluid
      GridGenerator::generate_from_name_and_arguments(*fluid_tria, "hyper_cube", "-1 : 1 : true");

      const unsigned int degree_velocity = 1;

      SolidBase<3,3> solid(param, fluid_tria, degree_velocity);
      solid.initial_setup();
      solid.setup_particles();
      solid.output_particles("output_particles");

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
