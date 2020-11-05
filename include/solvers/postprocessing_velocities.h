/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_postprocessing_velocities_h
#define lethe_postprocessing_velocities_h

// Dealii Includes

// Base
#include <deal.II/base/utilities.h>


// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Numerics
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// Lethe Includes
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <solvers/flow_control.h>

#include "navier_stokes_solver_parameters.h"
#include "post_processors.h"

// Std
#include <fstream>
#include <iostream>
#include <type_traits>

using namespace dealii;

template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  void
  calculate_average_velocities(
    const VectorType &                       local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::PostProcessing &       post_processing,
    const DofsType &                         locally_owned_dofs,
    const MPI_Comm &                         mpi_communicator);

  const VectorType
  get_average_velocities();

  const VectorType
  get_average_velocities(const double bulk_velocity);


private:
  TrilinosScalar inv_range_time;
  TrilinosScalar dt;

  VectorType sum_velocity_dt;
  VectorType average_velocities;
};

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType &                       local_evaluation_point,
  const std::shared_ptr<SimulationControl> &simulation_control,
  const Parameters::PostProcessing &       post_processing,
  const DofsType &                         locally_owned_dofs,
  const MPI_Comm &                         mpi_communicator)
{
  double time = simulation_control->get_current_time();
  double total_time = time - post_processing.initial_time;
  dt = simulation_control->calculate_time_step();

  if (simulation_control->get_step_number() == 0)
  {
    // Reinitializing vectors with zeros and dt at t = 0
    sum_velocity_dt.reinit(locally_owned_dofs,
                           mpi_communicator);
    average_velocities.reinit(locally_owned_dofs,
                              mpi_communicator);
  }
  else if (abs(total_time) < 1e-6 || total_time > 0)
  {
    // Generating average velocities at each time from initial time
    inv_range_time = 1. / (total_time + dt);

    VectorType velocity_dt(locally_owned_dofs,
                           mpi_communicator);
    velocity_dt.equ(dt, local_evaluation_point);
    sum_velocity_dt += velocity_dt;

    if (simulation_control->is_output_iteration())
      average_velocities.equ(inv_range_time, sum_velocity_dt);
  }
}

template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::
get_average_velocities()
{
  return average_velocities;
}

// Function not tested yet
template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::
get_average_velocities(const double bulk_velocity)
{
  VectorType nondimensionalized_average_velocities(average_velocities);
  nondimensionalized_average_velocities /= bulk_velocity;
  return nondimensionalized_average_velocities;
}


#endif
