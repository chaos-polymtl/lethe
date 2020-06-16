/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out_faces.h>

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <core/grids.h>
#include <core/solutions_output.h>
#include <core/utilities.h>

#include "core/time_integration_utilities.h"

#include <deal.II/base/std_cxx14/memory.h>


/*
 * Constructor for the Navier-Stokes base class
 */
template <int dim, int spacedim>
SolidBase<dim, spacedim>::SolidBase(
  NavierStokesSolverParameters<dim> &p_nsparam):
    mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , solid_tria(dynamic_cast<parallel::DistributedTriangulationBase<dim> *>(
      new parallel::distributed::Triangulation<dim>(
        this->mpi_communicator,
        typename Triangulation<dim>::MeshSmoothing(
          Triangulation<dim>::smoothing_on_refinement |
          Triangulation<dim>::smoothing_on_coarsening))))
  , solid_dh(*this->solid_tria)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,solid_fe

}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::initial_setup()
{
  TimerOutput::Scope t(this->computing_timer, "initiale_setup");
  
  solid_fe = std_cxx14::make_unique<FE_Nothing<dim, spacedim>>();
  solid_dh.distribute_dofs(*solid_fe);
  solid_quadrature_formula =
    std_cxx14::make_unique<QGauss<dim>>(nsparam.fem_paremeters.velocityOrder + 1);
}


template <int dim, int scapedim>
void
SolidBase<dim, spacedim>::setup_particles(fluid_tria)
{
  
}

template <int dim, int spacedim>
void
SolidBase<dim, spacedim>::output_particles()
{
  
}

template <int dim, int spacedim>
Particles::ParticleHandler
SolidBase<pacedim>::generate_solid_particle_handler()
{
  initial_setup();
  setup_particles(fluid_tria);
  return solid_particle_handler;
}




// Pre-compile the 2D and 3D version with the types that can occur
template class NavierStokesBase<2, TrilinosWrappers::MPI::Vector, IndexSet>;
template class NavierStokesBase<3, TrilinosWrappers::MPI::Vector, IndexSet>;
template class NavierStokesBase<2,
                                TrilinosWrappers::MPI::BlockVector,
                                std::vector<IndexSet>>;
template class NavierStokesBase<3,
                                TrilinosWrappers::MPI::BlockVector,
                                std::vector<IndexSet>>;
