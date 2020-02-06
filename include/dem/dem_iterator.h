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

#include "dem/contact_force.h"
#include "dem/contact_info_struct.h"
#include "dem/contact_search.h"
#include "dem/dem_solver_parameters.h"
#include "dem/insertion_info_struct.h"
#include "dem/integrator.h"
#include "dem/particle_wall_contact_detection.h"
#include "dem/particle_wall_contact_force.h"
#include "dem/pp_broad_search.h"

#ifndef DEMITERATOR_H_
#define DEMITERATOR_H_

template <int dim, int spacedim> class DEM_iterator {
public:
  DEM_iterator<dim, spacedim>();
  void engine(
      dealii::Particles::ParticleHandler<dim, spacedim> &,
      const dealii::Triangulation<dim, spacedim> &, int &, float &,
      std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
          &,
      std::vector<std::map<int, Particles::ParticleIterator<dim, spacedim>>> &,
      std::vector<std::map<int, ContactInfoStruct<dim, spacedim>>> &,
      std::vector<
          std::tuple<int, typename Triangulation<dim>::active_cell_iterator,
                     int, Point<dim>, Point<dim>>>,
      std::vector<std::tuple<
          std::pair<Particles::ParticleIterator<dim, spacedim>, int>,
          Point<dim>, Point<dim>, double, double, double, Point<dim>, double>>
          &,
      std::vector<std::tuple<std::string, int>>,
      dealii::Particles::PropertyPool &, ContactSearch<dim, spacedim>,
      ParticleWallContactDetection<dim, spacedim>, ContactForce<dim, spacedim>,
      ParticleWallContactForce<dim, spacedim>, Integrator<dim, spacedim> *,
      double, int, int, Point<dim>, double, int, int, int, float, float, float,
      float, float, float, float, float, InsertionInfoStruct<dim, spacedim>,
      int, int, PPBroadSearch<dim, spacedim>);

private:
  void forceReinit(Particles::ParticleHandler<dim, spacedim> &);
  // void checkSimBound(Particles::ParticleHandler<3, 3> &, ReadInputScript);
};

#endif /* DEMITERATOR_H_ */
