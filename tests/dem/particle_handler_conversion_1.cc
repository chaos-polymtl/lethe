// SPDX-FileCopyrightText: Copyright (c) 2020-2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Inserting particles using volume insertion class.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/particle_handler_conversion.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Generate a dummy triangulation to store the particles
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr);
  int refinement_number = 1;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Define a particle handler with a particle and properties
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::DEMProperties::PropertiesIndex::n_properties);
  // Inserting one particle in contact with wall
  Point<dim> position1;
  for (int d = 0; d < dim; ++d)
    position1[d] = 0.5;
  int                      id = 0;
  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, particle_cell);

  // Set the properties of the particle manually to some values
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::type]    = 1;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::dp]      = 0.005;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::v_x]     = 1;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::v_y]     = 2;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::v_z]     = 3;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::omega_x] = 40;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::omega_y] = 50;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::omega_z] = 60;
  pit1->get_properties()[DEM::DEMProperties::PropertiesIndex::mass]    = 100;

  // Output the properties of particle
  deallog << "Output of the original particle handler" << std::endl;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      deallog << "Particle " << particle->get_id() << std::endl;
      deallog << "Position : " << particle->get_location() << std::endl;
      for (unsigned int i = 0;
           i < DEM::DEMProperties::PropertiesIndex::n_properties;
           ++i)
        {
          deallog << "Property " << i << " : " << particle->get_properties()[i]
                  << std::endl;
        }
    }


  // Create a second Particle Handler using a second set of properties
  Particles::ParticleHandler<dim> particle_handler_output(
    tr, mapping, DEM::CFDDEMProperties::PropertiesIndex::n_properties);

  // Fill the second particle handler using the first one
  convert_particle_handler<dim,
                           DEM::DEMProperties::PropertiesIndex,
                           DEM::CFDDEMProperties::PropertiesIndex>(
    tr, particle_handler, particle_handler_output);

  // Output the properties of particle
  deallog << "Output of the second particle handler" << std::endl;
  for (auto particle = particle_handler_output.begin();
       particle != particle_handler_output.end();
       ++particle)
    {
      deallog << "Particle " << particle->get_id() << std::endl;
      deallog << "Position : " << particle->get_location() << std::endl;
      for (unsigned int i = 0;
           i < DEM::CFDDEMProperties::PropertiesIndex::n_properties;
           ++i)
        {
          deallog << "Property " << i << " : " << particle->get_properties()[i]
                  << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<2>();
      test<3>();
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
