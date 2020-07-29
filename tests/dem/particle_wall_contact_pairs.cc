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

// In this test, the performance of particle-wall broad search is investigated

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/find_boundary_cells_information.h>
#include <dem/pw_broad_search.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim> void test() {
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr, -1 * hyper_cube_length, hyper_cube_length,
                            true);
  int refinement_number = 2;
  tr.refine_global(refinement_number);

  Particles::ParticleHandler<dim> particle_handler;

  // inserting three particles at x = -0.4 , x = 0.4 and x = 0.8
  // which means only particle 3 is located in a boundary cell
  Point<dim> position1 = {-0.4, 0, 0};
  int id1 = 0;
  Point<dim> position2 = {0.4, 0, 0};
  int id2 = 1;
  Point<dim> position3 = {0.8, 0, 0};
  int id3 = 2;

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
  std::vector<typename Triangulation<dim>::active_cell_iterator>
      boundary_cells_with_faces;
  std::map<int, boundary_cells_info_struct<dim>> boundary_cells_information;
  FindBoundaryCellsInformation<dim> boundary_cells_object;
  boundary_cells_information =
      boundary_cells_object.find_boundary_cells_information(
          boundary_cells_with_faces, tr);

  // Calling particle-wall broad search
  PWBroadSearch<dim> broad_search_object;
  std::map<std::pair<int, int>,
           std::tuple<typename Particles::ParticleIterator<dim>, Tensor<1, dim>,
                      Point<dim>>>
      pw_contact_list;
  broad_search_object.find_PW_Contact_Pairs(boundary_cells_information,
                                            particle_handler, pw_contact_list);

  // Output
  for (auto pw_contact_list_iterator = pw_contact_list.begin();
       pw_contact_list_iterator != pw_contact_list.end();
       ++pw_contact_list_iterator) {
    auto contact_information = pw_contact_list_iterator->second;
    auto particle_information = std::get<0>(contact_information);
    deallog << "Particle " << particle_information->get_id()
            << " is located in a boundary cell" << std::endl;
  }
}

int main(int argc, char **argv) {
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
      argc, argv, numbers::invalid_unsigned_int);
  test<3>();
}
