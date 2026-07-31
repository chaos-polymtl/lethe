// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <dem/particle_handler_conversion.h>

// The DEM particle handlers are converted into CFD-DEM ones when a CFD-DEM
// simulation is started from the results of a DEM simulation. The multiphysics
// solvers are converted to one another so that the temperature of the
// particles is carried over.
template void
convert_particle_handler<2,
                         DEM::DEMProperties::PropertiesIndex,
                         DEM::CFDDEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2> &triangulation,
  const Particles::ParticleHandler<2>           &ph_in,
  Particles::ParticleHandler<2>                 &ph_out);

template void
convert_particle_handler<3,
                         DEM::DEMProperties::PropertiesIndex,
                         DEM::CFDDEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3> &triangulation,
  const Particles::ParticleHandler<3>           &ph_in,
  Particles::ParticleHandler<3>                 &ph_out);

template void
convert_particle_handler<2,
                         DEM::DEMMPProperties::PropertiesIndex,
                         DEM::CFDDEMMPProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2> &triangulation,
  const Particles::ParticleHandler<2>           &ph_in,
  Particles::ParticleHandler<2>                 &ph_out);

template void
convert_particle_handler<3,
                         DEM::DEMMPProperties::PropertiesIndex,
                         DEM::CFDDEMMPProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3> &triangulation,
  const Particles::ParticleHandler<3>           &ph_in,
  Particles::ParticleHandler<3>                 &ph_out);
