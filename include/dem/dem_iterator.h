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
#include "dem/dem_solver_parameters.h"
#include "dem/insertion_info_struct.h"
#include "dem/integrator.h"
#include "dem/physical_info_struct.h"
#include "dem/pp_broad_search.h"
#include "dem/pp_contact_force.h"
#include "dem/pp_contact_info_struct.h"
#include "dem/pp_fine_search.h"
#include "dem/pw_broad_search.h"
#include "dem/pw_contact_force.h"
#include "dem/pw_contact_info_struct.h"
#include "dem/pw_fine_search.h"

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
      std::vector<std::map<int, pp_contact_info_struct<dim, spacedim>>> &,
      std::vector<boundary_cells_info_struct<dim>>,
      std::vector<std::map<int, pw_contact_info_struct<dim, spacedim>>> &,
      std::vector<std::pair<std::string, int>>,
      dealii::Particles::PropertyPool &, PPContactForce<dim, spacedim> *,
      PWContactForce<dim, spacedim> *, Integrator<dim, spacedim> *, double, int,
      int, physical_info_struct<dim>, insertion_info_struct<dim, spacedim>,
      Tensor<1, dim>, int, PPBroadSearch<dim, spacedim>,
      PPFineSearch<dim, spacedim>, PWBroadSearch<dim, spacedim>,
      PWFineSearch<dim, spacedim>);

private:
  void forceReinit(Particles::ParticleHandler<dim, spacedim> &);
  // void checkSimBound(Particles::ParticleHandler<3, 3> &, ReadInputScript);
};

#endif /* DEMITERATOR_H_ */
