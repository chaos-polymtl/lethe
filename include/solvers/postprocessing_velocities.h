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
  void
  calculate_average_velocities(
    const VectorType &                       local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::PostProcessing &       post_processing,
    const DofsType &                         locally_owned_dofs,
    const MPI_Comm &                         mpi_communicator);

  void
  calculate_reynolds_stress(
    const VectorType &                        local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const DofsType &                          locally_owned_dofs,
    const MPI_Comm &                          mpi_communicator);

  const VectorType
  get_average_velocities();

  const VectorType
  get_average_velocities(const double bulk_velocity);

  const VectorType
  get_reynolds_stress();

  const VectorType
  get_reynolds_stress(const double bulk_velocity);


private:
  double total_time;

  TrilinosScalar inv_range_time;
  TrilinosScalar dt;

  VectorType sum_velocity_dt;
  VectorType average_velocities;

  VectorType sum_reynolds_stress_dt;
  VectorType reynolds_stress;
};

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
  total_time = time - post_processing.initial_time;
  dt = simulation_control->calculate_time_step();

  if (total_time + dt >= -1e-6 && total_time < -1e-6)
  {
    // Reinitializing vectors before calculating average
    sum_velocity_dt.reinit(locally_owned_dofs,
                           mpi_communicator);
    average_velocities.reinit(locally_owned_dofs,
                              mpi_communicator);
  }
  else if (total_time >= -1e-6)
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


/**
 * @brief calculate_reynolds_stress. This function calculates normal
 * time-averaged Reynold stresses and shear stress (<u'u'>, <v'v'>, <w'w'>
 * and <u'v'>.
 *
 * @param local_evaluation_point. The vector solutions with no ghost cells
 *
 * @param simulation_control. The simulation information (time)
 *
 * @param locally_owned_dofs. The owned dofs
 *
 * @param mpi_communicator. The mpi communicator information
 */
template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stress(
  const VectorType &                        local_evaluation_point,
  const std::shared_ptr<SimulationControl> &simulation_control,
  const DofsType &                          locally_owned_dofs,
  const MPI_Comm &                          mpi_communicator)
{

  if (total_time + dt >= -1e-6 && total_time < -1e-6)
  {
    // Reinitializing vectors before calculating average
    sum_reynolds_stress_dt.reinit(locally_owned_dofs,
                                  mpi_communicator);
    reynolds_stress.reinit(locally_owned_dofs,
                           mpi_communicator);
  }
  else if (total_time >= -1e-6)
  {
    VectorType reynolds_stress_dt(locally_owned_dofs,
                                  mpi_communicator);

    if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      const unsigned int begin_index = local_evaluation_point.local_range().first;
      const unsigned int end_index = local_evaluation_point.local_range().second;

      for (unsigned int i = begin_index; i <= end_index; i++)
      {
        if ((i + 4) % 4 == 0)
        {
          // Calculating (u'u')*dt, (v'v')*dt (w'w')*dt and (u'v')*dt
          reynolds_stress_dt[i] =
            (local_evaluation_point[i] - average_velocities[i]) *
            (local_evaluation_point[i] - average_velocities[i]) * dt;
          reynolds_stress_dt[i + 1] =
            (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
            (local_evaluation_point[i + 1] - average_velocities[i + 1]) * dt;
          reynolds_stress_dt[i + 2] =
            (local_evaluation_point[i + 2] - average_velocities[i + 2]) *
            (local_evaluation_point[i + 2] - average_velocities[i + 2]) * dt;
          reynolds_stress_dt[i + 3] =
            (local_evaluation_point[i] - average_velocities[i]) *
            (local_evaluation_point[i + 1] - average_velocities[i + 1]) * dt;
        }
      }
      // Summation of all reynolds stress during simulation
      sum_reynolds_stress_dt += reynolds_stress_dt;
    }
    // Next condition not tested yet.
    else if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::BlockVector>)
    {
      unsigned int begin_index = local_evaluation_point.block(0).local_range().first;
      unsigned int end_index   = local_evaluation_point.block(0).local_range().second;

      reynolds_stress_dt.block(0) = local_evaluation_point.block(0);
      reynolds_stress_dt.block(0) -= average_velocities.block(0);
      reynolds_stress_dt.block(0).scale(reynolds_stress_dt.block(0));
      reynolds_stress_dt.block(0) *= dt;

      for (unsigned int i = begin_index; i <= end_index; i += 3)
        if ((i + 3) % 3 == 0)
          reynolds_stress_dt.block(1)[i/3] = reynolds_stress_dt.block(0)[i] *
                                             reynolds_stress_dt.block(0)[i + 1] / dt;

      sum_reynolds_stress_dt += reynolds_stress_dt;
    }
  }

  // Calculating time-averaged reynolds stress
  if (simulation_control->is_output_iteration())
    reynolds_stress.equ(inv_range_time, sum_reynolds_stress_dt);
}

/**
 * @brief get_average_velocities. Gives the average of solutions
 */
template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::get_average_velocities()
{
  return average_velocities;
}

/**
 * @brief get_average_velocities. Gives the nondimensionalized of solutions
 *
 * @param bulk_velocity. The bulk velocity calculated by the FlowControl class
 */
template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::
get_average_velocities(const double bulk_velocity)
{
  VectorType nondimensionalized_average_velocities(average_velocities);
  nondimensionalized_average_velocities /= bulk_velocity;

  return nondimensionalized_average_velocities;
}

// Function not tested yet
template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::get_reynolds_stress()
{
  return reynolds_stress;
}

// Function not tested yet
template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::get_reynolds_stress(
  const double bulk_velocity)
{
  VectorType nondimensionalized_reynolds_stress(reynolds_stress);
  nondimensionalized_reynolds_stress /= bulk_velocity;
  return nondimensionalized_reynolds_stress;
}


#endif
