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

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/contact_search.h"
#include "dem/particle_insertion.h"

using namespace dealii;

template <int dim>
void
test()
{
  // Generate a cube triangulation and refine it twice globally
  Triangulation<dim, dim>              tr;
  Particles::ParticleHandler<dim, dim> particle_handler;
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  int cellNum = tr.n_active_cells();


  MappingQ1<dim> mapping;


  ContactSearch<dim, dim> cs1;



  std::pair<std::vector<
              std::set<typename Triangulation<dim, dim>::active_cell_iterator>>,
            std::vector<typename Triangulation<dim, dim>::active_cell_iterator>>
    cellNeighbor = cs1.findCellNeighbors(tr);

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
                           pt1_info = GridTools::find_active_cell_around_point(mapping, tr, position1);
  Particles::Particle<dim> particle1(position1, pt1_info.second, id1);
  Particles::ParticleIterator<dim, dim> pit1 =
    particle_handler.insert_particle(particle1, pt1_info.first);

  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
                           pt2_info = GridTools::find_active_cell_around_point(mapping, tr, position2);
  Particles::Particle<dim> particle2(position2, pt2_info.second, id2);
  Particles::ParticleIterator<dim, dim> pit2 =
    particle_handler.insert_particle(particle2, pt2_info.first);

  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
                           pt3_info = GridTools::find_active_cell_around_point(mapping, tr, position3);
  Particles::Particle<dim> particle3(position3, pt3_info.second, id3);
  Particles::ParticleIterator<dim, dim> pit3 =
    particle_handler.insert_particle(particle3, pt3_info.first);

  std::vector<std::pair<Particles::ParticleIterator<dim, dim>,
                        Particles::ParticleIterator<dim, dim>>>
    pairs;
  pairs = cs1.findContactPairs(particle_handler, tr, cellNeighbor.first);



  for (unsigned int i = 0; i != pairs.size(); i++)
    {
      deallog << "A pair is detected: particle " << pairs[i].first->get_id()
              << " and particle " << pairs[i].second->get_id() << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  test<3>();
}
