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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include "dem/boundary_cells_info_struct.h"
#include <deal.II/distributed/tria.h>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <iostream>
#include <vector>

using namespace dealii;

#ifndef PWBROADSEARCH_H_
#define PWBROADSEARCH_H_

/**
 * This class is used for broad particle-wall contact search. Broad search
 * is used to obtain all the particles located at boundary cells
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim, int spacedim = dim> class PWBroadSearch {
public:
  PWBroadSearch<dim, spacedim>();

  /**
   * Finds a vector of tuples (tuple of particles located in boundary cells,
   * normal vector of the boundary face and a point on the face) which shows the
   * candidate particle-wall collision pairs. These collision pairs will be
   * investigated in the fine search to check if they are in contact or not
   *
   * @param boundary_cells_information Information of the boundary cells and
   * faces. This is the output of the FindBoundaryCellsInformation class
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @return A vector of tuples. Each element of vector (tuple) contains a
   * particle located near boundaries, the normal vector of the corresponding
   * face boundary and a point on the boundary
   */

  std::vector<std::tuple<typename Particles::ParticleIterator<dim, spacedim>,
                         Tensor<1, dim>, Point<dim>>>
  find_PW_Contact_Pairs(
      std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
      Particles::ParticleHandler<dim, spacedim> &particle_handler);
};

#endif /* PWBROADSEARCH_H_ */
