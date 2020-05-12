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

// Three particles are inserted manually in the x direction.
// We check if these particles appear correctly to each other
// as potential neighbhors

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/find_cell_neighbors.h>
#include <dem/pp_broad_search.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  // Generate a cube triangulation and refine it twice globally
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  Particles::ParticleHandler<dim> particle_handler;

  MappingQ1<dim>     mapping;
  PPBroadSearch<dim> broad_search_object;

  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
    cell_neighbor_list;

  // Finding cell neighbors list, it is required for finding the broad search
  // pairs
  FindCellNeighbors<dim> cell_neighbor_object;
  cell_neighbor_list = cell_neighbor_object.find_cell_neighbors(triangulation);

  // inserting three particles at x = -0.4 , x = 0.4 and x = 0.8
  // which means they are inserted in three adjacent cells in x direction
  Point<3> position1 = {-0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.4, 0, 0};
  int      id2       = 1;
  Point<3> position3 = {0.8, 0, 0};
  int      id3       = 2;

  // Manually insert the three particles
  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
                                   pt1_info = GridTools::find_active_cell_around_point(mapping,
                                                        triangulation,
                                                        position1);
  Particles::Particle<dim>         particle1(position1, pt1_info.second, id1);
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, pt1_info.first);

  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
                                   pt2_info = GridTools::find_active_cell_around_point(mapping,
                                                        triangulation,
                                                        position2);
  Particles::Particle<dim>         particle2(position2, pt2_info.second, id2);
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, pt2_info.first);

  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
                                   pt3_info = GridTools::find_active_cell_around_point(mapping,
                                                        triangulation,
                                                        position3);
  Particles::Particle<dim>         particle3(position3, pt3_info.second, id3);
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, pt3_info.first);

  // Calling broad search function
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
    broad_search_pairs;
  broad_search_object.find_PP_Contact_Pairs(particle_handler,
                                            cell_neighbor_list,
                                            broad_search_pairs);

  // Output
  for (auto pairs_iterator = broad_search_pairs.begin();
       pairs_iterator != broad_search_pairs.end();
       ++pairs_iterator)
    {
      deallog << "A pair is detected: particle "
              << pairs_iterator->first->get_id() << " and particle "
              << pairs_iterator->second->get_id() << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);
  test<3>();
}
