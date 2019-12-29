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

template <int dim, int spacedim = dim>
class ContactSearch
{
public:
  ContactSearch<dim, spacedim>();

  std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                        Particles::ParticleIterator<dim, spacedim>>>
  findContactPairs(
    dealii::Particles::ParticleHandler<dim, spacedim> &,
    const Triangulation<dim, spacedim> &,
    std::vector<typename Triangulation<dim>::active_cell_iterator>,
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>);


  std::pair<
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>,
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
  findCellNeighbors(const Triangulation<dim, spacedim> &);


  void
  fineSearch(std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                                   Particles::ParticleIterator<dim, spacedim>>>,
             std::vector<
               std::tuple<std::pair<Particles::ParticleIterator<dim, spacedim>,
                                    Particles::ParticleIterator<dim, spacedim>>,
                          double,
                          Point<dim>,
                          double,
                          Point<dim>,
                          double,
                          double>> &,
             float);

  std::vector<Particles::ParticleIterator<dim, spacedim>>
  pWSearch(dealii::Particles::ParticleHandler<dim, spacedim> &);
};

#endif /* CONTACTSEARCH_H_ */
