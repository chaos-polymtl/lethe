/* ---------------------------------------------------------------------
  *
  * Copyright (C) 2019 - 2020 by the Lethe authors
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

#include <dem/particle_particle_contact_force.h>

using namespace DEM;

template <int dim>
inline void
ParticleParticleContactForce<dim>::find_effective_radius_and_mass(
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties)
{
  effective_mass = (particle_one_properties[DEM::PropertiesIndex::mass] *
                    particle_two_properties[DEM::PropertiesIndex::mass]) /
                   (particle_one_properties[DEM::PropertiesIndex::mass] +
                    particle_two_properties[DEM::PropertiesIndex::mass]);
  effective_radius = (particle_one_properties[DEM::PropertiesIndex::dp] *
                      particle_two_properties[DEM::PropertiesIndex::dp]) /
                     (2 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                           particle_two_properties[DEM::PropertiesIndex::dp]));
}

template class ParticleParticleContactForce<2>;
template class ParticleParticleContactForce<3>;
