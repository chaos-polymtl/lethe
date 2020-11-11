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


// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

// Numerics
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// Lethe Includes
#include <core/parameters.h>
#include <core/simulation_control.h>
#include "navier_stokes_solver_parameters.h"
#include "post_processors.h"


using namespace dealii;


/**
 * @brief AverageVelocities. The AverageVelocities class calculates the
 * time-averaged velocities and pressure. The generated vector is output
 * with the solution and visualization is possible.
 */
template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  /**
   * @brief calculate_average_velocities. This function calculates time-averaged
   * velocities and pressure with dof vector with no ghost cell.
   *
   * @param local_evaluation_point. The vector solutions with no ghost cells
   *
   * @param simulation_control. The simulation information (time)
   *
   * @param post_processing. The parameters to start the processing
   *
   * @param locally_owned_dofs. The owned dofs
   *
   * @param mpi_communicator. The mpi communicator information
   */
  void
  calculate_average_velocities(
    const VectorType &                        local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::PostProcessing &        post_processing,
    const DofsType &                          locally_owned_dofs,
    const MPI_Comm &                          mpi_communicator);

  /**
   * @brief get_average_velocities. Gives the average of solutions.
   */
  const VectorType
  get_average_velocities();

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


#endif
