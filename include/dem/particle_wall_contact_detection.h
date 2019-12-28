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

#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>


using namespace dealii;


#ifndef PWCONTACTDETECTION_H_
#  define PWCONTACTDETECTION_H_

template <int dim, int spacedim>
class ParticleWallContactDetection
{
public:
  ParticleWallContactDetection<dim, spacedim>();
  void
  boundaryCellsAndFaces(
    const Triangulation<dim, spacedim> &,
    std::vector<std::tuple<int,
                           typename Triangulation<dim>::active_cell_iterator,
                           int,
                           Point<dim>,
                           Point<dim>>> &);
  std::vector<std::tuple<
    std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
    Point<dim>,
    Point<dim>>>
  pwcontactlist(
    std::vector<std::tuple<int,
                           typename Triangulation<dim>::active_cell_iterator,
                           int,
                           Point<dim>,
                           Point<dim>>>,
    Particles::ParticleHandler<dim, spacedim> &);

  void
  pwFineSearch(
    std::vector<std::tuple<
      std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
      Point<dim>,
      Point<dim>>>,
    std::vector<std::tuple<
      std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
      Point<dim>,
      Point<dim>,
      double,
      double,
      double,
      Point<dim>,
      double>> &,
    float);

private:
  Point<dim> findProjection(Point<dim>, Point<dim>);
};

#endif /* PWCONTACTDETECTION_H_ */
