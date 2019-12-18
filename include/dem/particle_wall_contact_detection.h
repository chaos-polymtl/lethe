/*
 * pwcontactdetection.h
 *
 *  Created on: Nov 26, 2019
 *      Author: shahab
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

class pwcontactdetection
{
public:
  pwcontactdetection();
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
