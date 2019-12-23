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

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include <iostream>
#include <vector>


using namespace dealii;

#ifndef CONTACTSEARCH_H_
#  define CONTACTSEARCH_H_

class ContactSearch
{
public:
  ContactSearch();
  std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                        Particles::ParticleIterator<3, 3>>>
    findContactPairs(
      dealii::Particles::ParticleHandler<3, 3> &,
      const Triangulation<3, 3> &,
      std::vector<Triangulation<3>::active_cell_iterator>,
      std::vector<std::set<Triangulation<3>::active_cell_iterator>>);

  std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,
            std::vector<Triangulation<3>::active_cell_iterator>>
       findCellNeighbors(int, const Triangulation<3, 3> &);
  void fineSearch(
    std::vector<std::pair<Particles::ParticleIterator<3, 3>,
                          Particles::ParticleIterator<3, 3>>>,
    dealii::Particles::ParticleHandler<3, 3> &,
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                                     Particles::ParticleIterator<3, 3>>,
                           double,
                           Point<3>,
                           double,
                           Point<3>,
                           double,
                           double>> &,
    float);
  // std::directSearch();

  std::vector<Particles::ParticleIterator<3, 3>>
    pWSearch(dealii::Particles::ParticleHandler<3, 3> &);
};

#endif /* CONTACTSEARCH_H_ */
