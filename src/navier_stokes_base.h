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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#ifndef LETHE_BASERNAVIERSTOKES_H
#define LETHE_BASERNAVIERSTOKES_H

// Dealii Includes

// Base
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// Lac
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

// Grid
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Numerics
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// Distributed
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

// Lethe Includes
#include "bdf.h"
#include "boundaryconditions.h"
#include "manifolds.h"
#include "navierstokessolverparameters.h"
#include "parameters.h"
#include "postprocessors.h"
#include "pvdhandler.h"
#include "simulationcontrol.h"

// Std
#include <fstream>
#include <iostream>

using namespace dealii;

/**
 * A base class for all the Navier-Stokes equation
 * This class regroups common facilities that are shared by all
 * the Navier-Stokes implementations to reduce code multiplicity
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 * @author Bruno Blais, 2019
 */

template <int dim>
class NavierStokesBase
{
public:
  NavierStokesBase(NavierStokesSolverParameters<dim> &nsparam,
                   const unsigned int                 degreeVelocity,
                   const unsigned int                 degreePressure);

protected:
  MPI_Comm           mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;


  parallel::distributed::Triangulation<dim> triangulation;
  DoFHandler<dim>                           dof_handler;
  FESystem<dim>                             fe;
  ConditionalOStream                        pcout;

  TimerOutput computing_timer;

  NavierStokesSolverParameters<dim> nsparam;

  // Finite element order used
  const unsigned int degreeVelocity_;
  const unsigned int degreePressure_;
  unsigned int       degreeQuadrature_;

  // Post-processing variables
  TableHandler enstrophy_table;
  TableHandler kinetic_energy_table;
};


template <int dim>
NavierStokesBase<dim>::NavierStokesBase(
  NavierStokesSolverParameters<dim> &p_nsparam,
  const unsigned int                 p_degreeVelocity,
  const unsigned int                 p_degreePressure)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , triangulation(this->mpi_communicator,
                  typename Triangulation<dim>::MeshSmoothing(
                    Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening))
  , dof_handler(this->triangulation)
  , fe(FE_Q<dim>(p_degreeVelocity), dim, FE_Q<dim>(p_degreePressure), 1)
  , pcout(std::cout,
          (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0))
  , computing_timer(this->mpi_communicator,
                    pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , nsparam(p_nsparam)
  , degreeVelocity_(p_degreeVelocity)
  , degreePressure_(p_degreePressure)
  , degreeQuadrature_(p_degreeVelocity + 1)
{
  //  pcout << "Base classe bitches " << std::endl;
  // Overide default value of quadrature point if they are specified
  if (nsparam.femParameters.quadraturePoints > 0)
    degreeQuadrature_ = nsparam.femParameters.quadraturePoints;

  // Change the behavior of the timer for situations when you don't want outputs
  if (nsparam.timer.type == Parameters::Timer::none)
    this->computing_timer.disable_output();

  //  pcout << "Base class is finished" << std::endl;
}


#endif
