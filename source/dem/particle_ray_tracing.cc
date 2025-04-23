// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later\

#include <dem/particle_ray_tracing.h>

template <int dim, typename PropertiesIndex>
ParticleRayTracing<dim,PropertiesIndex>::ParticleRayTracing(ParticleRayTracingParameters<dim> parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , photon_handler(triangulation, mapping, 1)
  , particle_handler(triangulation, mapping, PropertiesIndex::n_properties)
{}






template class ParticleRayTracing<2, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<2, DEM::DEMMPProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMMPProperties::PropertiesIndex>;