/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019-
 */

/**
 * @brief In this test, the performance of particle-particle fine search class
 * is evaluated.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <dem/dem_contact_manager.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_contact_info.h>
#include <dem/particle_particle_fine_search.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties());

  DEMContactManager<dim> container_manager;

  // Finding cell neighbors
  FindCellNeighbors<dim> cell_neighbor_object;
  cell_neighbor_object.find_cell_neighbors(
    triangulation,
    container_manager.cells_local_neighbor_list,
    container_manager.cells_ghost_neighbor_list);

  // Inserting two particles in contact
  Point<3> position1 = {0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.40499, 0, 0};
  int      id2       = 1;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());
  double       particle_diameter      = 0.005;
  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);

  pit1->get_properties()[DEM::PropertiesIndex::type]    = 0;
  pit1->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
  pit1->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass]    = 1;

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);
  pit2->get_properties()[DEM::PropertiesIndex::type]    = 0;
  pit2->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
  pit2->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::mass]    = 1;

  // Calling broad search
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      container_manager.particle_container[particle_iterator->get_id()] =
        particle_iterator;
    }

  container_manager.execute_particle_particle_broad_search(particle_handler);

  // Calling fine search
  container_manager.execute_particle_particle_fine_search(
    neighborhood_threshold);

  // Output
  for (auto adjacent_particles_iterator =
         container_manager.local_adjacent_particles.begin();
       adjacent_particles_iterator !=
       container_manager.local_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      auto information_iterator =
        (&adjacent_particles_iterator->second)->begin();
      while (information_iterator !=
             (&adjacent_particles_iterator->second)->end())
        {
          deallog << "The particle pair in contact are particles: "
                  << information_iterator->second.particle_one->get_id()
                  << " and "
                  << information_iterator->second.particle_two->get_id()
                  << std::endl;
          deallog << "Tangential overlap at the beginning of contact is: "
                  << information_iterator->second.tangential_overlap[0] << " "
                  << information_iterator->second.tangential_overlap[1] << " "
                  << information_iterator->second.tangential_overlap[2]
                  << std::endl;
          ++information_iterator;
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
