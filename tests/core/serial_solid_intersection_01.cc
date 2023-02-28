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
 * This test generates a flat plane using a simplex mesh out of the deal.II grid
 * generator and calculates the cells that it cuts on an hyper_ball.
 * The list of cut cells/triangles combo is printed and a vtu file with the
 * intersections is generated.
 */

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/data_out.h>

// Lethe
#include <core/parameters.h>
#include <core/serial_solid.h>
#include <core/solid_objects_parameters.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // Generate a background fluid triangulation made of a sphere
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> fluid_tria =
    std::make_shared<parallel::distributed::Triangulation<3>>(
      mpi_communicator,
      typename Triangulation<3>::MeshSmoothing(
        Triangulation<3>::smoothing_on_refinement |
        Triangulation<3>::smoothing_on_coarsening));

  // Mesh of the fluid
  GridGenerator::hyper_ball(*fluid_tria, {0.2, 0, 0}, 2);

  // Generate the serial solid

  // Parameters for the Serial solid object
  auto param             = std::make_shared<Parameters::RigidSolidObject<3>>();
  param->solid_mesh.type = Parameters::Mesh::Type::dealii;
  param->solid_mesh.grid_type          = "hyper_rectangle";
  param->solid_mesh.grid_arguments     = "-1, -1 : 1, 1 : 1";
  param->solid_mesh.initial_refinement = 5;
  param->solid_mesh.simplex            = true;
  param->solid_mesh.translate          = false;
  param->solid_mesh.rotate             = false;

  std::shared_ptr<Mapping<3>> fluid_mapping =
    std::make_shared<MappingQGeneric<3>>(1);
  SerialSolid<2, 3> solid(param, 0);

  // Set-up the solid
  solid.initial_setup();


  // Calculate the intersections between the background triangulation and the
  // floating solid


  // Generate a VTU for debugging purposes which shows the intersection

  // Generate the particles
  //  Particles::DataOut<3, 3>                       particles_out;
  //  std::shared_ptr<Particles::ParticleHandler<3>> solid_particle_handler =
  //    solid.get_solid_particle_handler();
  //  particles_out.build_patches(*solid_particle_handler);
  //  const std::string filename = ("particles.vtu");
  //  particles_out.write_vtu_in_parallel(filename, mpi_communicator);

  //      deallog << "Particle location: " << particle.get_location() <<
  //      std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      test();
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
