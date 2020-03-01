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
 * Author: Bruno Blais, Shahab Golshan, Polytechnique Montreal, 2019-
 */

#include <deal.II/base/tensor.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/property_pool.h>

#include <dem/dem_iterator.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/pp_broad_search.h>
#include <dem/pp_contact_info_struct.h>
#include <dem/pp_fine_search.h>
#include <dem/pp_linear_force.h>
#include <dem/pp_nonlinear_force.h>
#include <dem/pw_broad_search.h>
#include <dem/pw_contact_info_struct.h>
#include <dem/pw_fine_search.h>
#include <dem/pw_linear_force.h>
#include <dem/pw_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>
#include <dem/visualization.h>

#ifndef LETHE_DEM_H
#define LETHE_DEM_H

/**
 * The DEM class which initializes all the required parameters and iterates over
 * the DEM iterator
 */

template <int dim> class DEMSolver {
public:
  DEMSolver(DEMSolverParameters<dim> dem_parameters);

  /**
   * Initialiazes all the required parameters and iterates over the DEM iterator
   * (DEM engine).
   */
  void solve();

private:
  MPI_Comm mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
  ConditionalOStream pcout;
  DEMSolverParameters<dim> parameters;
  parallel::distributed::Triangulation<dim> triangulation;
  Particles::PropertyPool property_pool;
  MappingQGeneric<dim> mapping;
  Particles::ParticleHandler<dim, dim> particle_handler;

  /**
   * Defines or reads the mesh based on the information provided by the user.
   * Gmsh files can also be read in this function.
   */
  void read_mesh();
};

#endif
