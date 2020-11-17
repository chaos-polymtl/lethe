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
 * time-averaged velocities and pressure (<u>, <v>, <w>, <p>). The generated
 * vector is output with the solution and visualization is possible.
 */
template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  AverageVelocities();
  /**
   * @brief calculate_average_velocities. This function calculates time-averaged
   * velocities and pressure with dof vector with no ghost cell.
   *
   * @param local_evaluation_point. The vector solutions with no ghost cells
   *
   * @param post_processing. The parameters to start the processing
   *
   * @param current_time. The current time in the simulation
   *
   * @param time_step. The current time step
   *
   * @param is_output_iteration. The information if output is required at this
   *                             iteration
   *
   * @param locally_owned_dofs. The owned dofs
   *
   * @param mpi_communicator. The mpi communicator information
   */
  void
  calculate_average_velocities(
    const VectorType &                local_evaluation_point,
    const Parameters::PostProcessing &post_processing,
    const double &                    current_time,
    const double &                    time_step,
    const bool &                      is_output_iteration,
    const DofsType &                  locally_owned_dofs,
    const MPI_Comm &                  mpi_communicator);


  /**
   * @brief calculate_reynolds_stress. This function calculates normal
   * time-averaged Reynold stresses and shear stress (<u'u'>, <v'v'>, <w'w'>
   * and <u'v'>.
   *
   * @param local_evaluation_point. The vector solutions with no ghost cells
   *
   * @param is_output_iteration. The information if output is required at this
   *                             iteration
   *
   */
  void
  calculate_reynolds_stress(const VectorType &local_evaluation_point,
                            const bool &      is_output_iteration);
  /**
   * @brief get_average_velocities. Gives the average of solutions.
   */
  const VectorType
  get_average_velocities();

  /**
   * @brief get_reynolds_stress. Gives the time-averaged normal Reynolds
   * stresses and shear stress
   */
  const VectorType
  get_reynolds_stress();

private:
  TrilinosScalar inv_range_time;
  TrilinosScalar dt_0;

  VectorType velocity_dt;
  VectorType sum_velocity_dt;
  VectorType average_velocities;

  VectorType reynolds_stress_dt;
  VectorType sum_reynolds_stress_dt;
  VectorType reynolds_stress;

  double dt;
  bool   average_calculation;
  double real_initial_time;
};

template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::get_average_velocities()
{
  return average_velocities;
}

template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::get_reynolds_stress()
{
  return reynolds_stress;
}

#endif
