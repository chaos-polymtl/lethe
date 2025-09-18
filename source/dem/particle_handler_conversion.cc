// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <dem/particle_handler_conversion.h>

template <int,
          DEM::DEMProperties::PropertiesIndex,
          DEM::CFDDEMProperties::PropertiesIndex>
void
convert_particle_handler(
  const parallel::distributed::Triangulation<3> &triangulation,
  const Particles::ParticleHandler<3>           &ph_in,
  Particles::ParticleHandler<3>                 &ph_out);

template <int,
          DEM::DEMProperties::PropertiesIndex,
          DEM::CFDDEMProperties::PropertiesIndex>
void
convert_particle_handler(
  const parallel::distributed::Triangulation<2> &triangulation,
  const Particles::ParticleHandler<2>           &ph_in,
  Particles::ParticleHandler<2>                 &ph_out);
