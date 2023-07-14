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
 * @brief In this test, the performance of particle-wall broad search is
 * investigated.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_wall_broad_search.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  tr.refine_global(refinement_number);

  MappingQ1<dim> mapping;


  Particles::ParticleHandler<dim> particle_handler(tr, mapping);

  // inserting three particles at x = -0.4 , x = 0.4 and x = 0.8
  // which means only particle 3 is located in a boundary cell
  Point<dim> position1 = {-0.4, 0, 0};
  int        id1       = 0;
  Point<dim> position2 = {0.4, 0, 0};
  int        id2       = 1;
  Point<dim> position3 = {0.8, 0, 0};
  int        id3       = 2;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(tr, particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);

  Particles::Particle<dim> particle3(position3, position3, id3);
  typename Triangulation<dim>::active_cell_iterator cell3 =
    GridTools::find_active_cell_around_point(tr, particle3.get_location());
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, cell3);

  // Calling find_boundary_cells_information function to find the information of
  // boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::vector<unsigned int>     outlet_boundaries;
  boundary_cells_object.build(
    tr,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Calling particle-wall broad search
  ParticleWallBroadSearch<dim> broad_search_object;
  // P-W broad search
  ParticleWallBroadSearch<dim> particle_wall_broad_search_object;
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    particle_wall_contact_list;
  particle_wall_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);
  broad_search_object.find_particle_wall_contact_pairs(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);

  // Output
  for (auto particle_wall_contact_list_iterator =
         particle_wall_contact_list.begin();
       particle_wall_contact_list_iterator != particle_wall_contact_list.end();
       ++particle_wall_contact_list_iterator)
    {
      auto particle_wall_candidate_content =
        &particle_wall_contact_list_iterator->second;
      for (auto particle_wall_candidate_content_iterator =
             particle_wall_candidate_content->begin();
           particle_wall_candidate_content_iterator !=
           particle_wall_candidate_content->end();
           ++particle_wall_candidate_content_iterator)
        {
          auto contact_information =
            particle_wall_candidate_content_iterator->second;
          auto particle_information = std::get<0>(contact_information);
          deallog << "Particle " << particle_information->get_id()
                  << " is located in a boundary cell" << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
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
