/*
 * DEMiterator.h
 *
 *  Created on: Sep 26, 2019
 *      Author: shahab
 */
#include "integration.h"
#include "read_input_script.h"

#ifndef DEMITERATOR_H_
#  define DEMITERATOR_H_

class DEM_iterator
{
public:
  DEM_iterator();
  void
  engine(
    int &,
    dealii::Particles::ParticleHandler<3, 3> &,
    const dealii::Triangulation<3, 3> &,
    int &,
    float &,
    ReadInputScript,
    std::pair<std::vector<std::set<Triangulation<3>::active_cell_iterator>>,
              std::vector<Triangulation<3>::active_cell_iterator>>,
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                                     Particles::ParticleIterator<3, 3>>,
                           std::vector<double>,
                           double,
                           std::vector<double>,
                           double,
                           std::vector<double>,
                           std::vector<double>,
                           double,
                           double>> &,
    std::vector<std::tuple<int,
                           Triangulation<3>::active_cell_iterator,
                           int,
                           Point<3>,
                           Point<3>>>,
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                           Point<3>,
                           Point<3>,
                           double,
                           double,
                           double,
                           Point<3>,
                           double>> &);

private:
  void forceReinit(Particles::ParticleHandler<3, 3> &);
  void checkSimBound(Particles::ParticleHandler<3, 3> &, ReadInputScript);
};

#endif /* DEMITERATOR_H_ */
