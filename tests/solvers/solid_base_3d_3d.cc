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

#include <mpi.h>

#include "../tests.h"
#include "solvers/solid_base.h"

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);

      MPI_Comm mpi_communicator(MPI_COMM_WORLD);

      initlog();

      Parameters::Nitsche                                        param;
      std::shared_ptr<parallel::DistributedTriangulationBase<3>> fluid_tria =
        std::make_shared<parallel::distributed::Triangulation<3>>(
          mpi_communicator,
          typename Triangulation<3>::MeshSmoothing(
            Triangulation<3>::smoothing_on_refinement |
            Triangulation<3>::smoothing_on_coarsening));

      // Mesh of the solid
      param.solid_mesh.type           = Parameters::Mesh::Type::dealii;
      param.solid_mesh.grid_type      = "hyper_ball";
      param.solid_mesh.grid_arguments = "0 , 0 , 0 :  0.75 : false";

      // Mesh of the fluid
      GridGenerator::hyper_cube(*fluid_tria, -1, 1);

      const unsigned int degree_velocity = 1;

      SolidBase<3, 3> solid(param, fluid_tria, degree_velocity);
      solid.initial_setup();
      solid.setup_particles();

      Particles::DataOut<3, 3>                       particles_out;
      std::shared_ptr<Particles::ParticleHandler<3>> solid_particle_handler =
        solid.get_solid_particle_handler();
      particles_out.build_patches(*solid_particle_handler);
      const std::string filename = ("particles.vtu");
      particles_out.write_vtu_in_parallel(filename, mpi_communicator);

      for (const auto &particle : (*solid_particle_handler))
        {
          deallog << "Particle index: " << particle.get_id() << std::endl;
          deallog << "Particle location: " << particle.get_location()
                  << std::endl;
          deallog << "Particle JxW: " << particle.get_properties()[0] << std::endl;
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
