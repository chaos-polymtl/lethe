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
 * Authors: Carole-Anne Daunais, Valérie Bibeau, Polytechnique Montreal, 2020
 */

#ifndef lethe_solid_base_h
#define lethe_solid_base_h

// Dealii Includes

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h> 

// Fe
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 

// Distributed
#include <deal.II/distributed/tria_base.h> 
#include <deal.II/distributed/tria.h> 

// Particles
#include <deal.II/particles/data_out.h> 
#include <deal.II/particles/particle_handler.h> 
#include <deal.II/particles/particle_accessor.h>

// Lethe Includes
#include <core/parameters.h>
#include <solvers/navier_stokes_solver_parameters.h> 

// Std
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * A base class that generates a particle handler for the solid 
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @tparam VectorType  The Vector type used for the solvers
 *
 * @tparam DofsType the type of dof storage indices
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim, int spacedim = dim>
class SolidBase 
{
public:

  //Member functions
  SolidBase(NavierStokesSolverParameters<dim> &param,
  std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>> fluid_tria);
  void initial_setup();
  void setup_particles();
  void output_particles(std::string fprefix) const;
  std::shared_ptr<Particles::ParticleHandler<spacedim>>  generate_solid_particle_handler();

private:

  // Member variables
  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  std::shared_ptr<parallel::DistributedTriangulationBase<dim, spacedim>> solid_tria;
  std::shared_ptr<parallel::DistributedTriangulationBase<spacedim>>      fluid_tria;
  std::unique_ptr<DoFHandler<dim, spacedim>>                             solid_dh;
  std::unique_ptr<FESystem<dim, spacedim>>                               solid_fe;
  std::shared_ptr<Particles::ParticleHandler<spacedim>>                  solid_particle_handler;

  NavierStokesSolverParameters<dim>                                      param;
};

#endif
