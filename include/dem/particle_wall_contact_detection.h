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

class ParticleWallContactDetection
{
public:
  ParticleWallContactDetection();
  void
  boundaryCellsAndFaces(
    const Triangulation<3, 3> &,
    std::vector<std::tuple<int,
                           Triangulation<3>::active_cell_iterator,
                           int,
                           Point<3>,
                           Point<3>>> &);
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                         Point<3>,
                         Point<3>>>
  pwcontactlist(std::vector<std::tuple<int,
                                       Triangulation<3>::active_cell_iterator,
                                       int,
                                       Point<3>,
                                       Point<3>>>,
                Particles::ParticleHandler<3, 3> &);

  void pwFineSearch(
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                           Point<3>,
                           Point<3>>>,
    Particles::ParticleHandler<3, 3> &,
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                           Point<3>,
                           Point<3>,
                           double,
                           double,
                           double,
                           Point<3>,
                           double>> &,
    float);

private:
  Point<3> findProjection(Point<3>, Point<3>);
  double   dotProduct(Point<3>, Point<3>);
  Point<3> crossProduct(Point<3>, Point<3>);
  double   vecValue(Point<3>);
};

#endif /* PWCONTACTDETECTION_H_ */
