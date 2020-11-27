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
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

// Dofs
#include <deal.II/dofs/dof_handler.h>

// Lethe Includes
#include <core/parameters.h>

using namespace dealii;

/**
 * @brief AverageVelocities. The AverageVelocities class calculates the
 * time-averaged velocities and pressure (<u>, <v>, <w>, <p>) and the
 * independent components of the Reynolds stresses tensor (<u'u'>, <v'v'>,
 * <w'w'>, <u'v'>, <v'w'>, <w'u'>). The generated vectors are displayable of
 * visualization software.
 *
 * Important : Time-averaging velocities and calculating reynolds stresses are
 * currently unavailable for mesh adaptation.
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
   * @param locally_owned_dofs. The owned dofs
   *
   * @param locally_owned_rs_components. The owned Reynolds stress components
   *
   * @param mpi_communicator. The mpi communicator information
   */
  void
  calculate_average_velocities(
    const VectorType &                local_evaluation_point,
    const Parameters::PostProcessing &post_processing,
    const double &                    current_time,
    const double &                    time_step);

  /**
   * @brief calculate_reynolds_stress. This function calculates normal and
   * other resolved time-averaged Reynold stresses (<u'u'>, <v'v'>, <w'w'>
   * and <u'v'>, <v'w'>, <w'u'>).
   *
   * @param local_evaluation_point. The vector solutions with no ghost cells
   */
  void
  calculate_reynolds_stresses(const VectorType &local_evaluation_point);

  /**
   * @brief get_average_velocities. Gives the average of solutions.
   */
  const VectorType &
  get_average_velocities()
  {
    return average_velocities;
  }

  /**
   * @brief get_reynolds_stress. Gives the time-averaged Reynolds
   * stresses.
   */
  const LinearAlgebra::distributed::Vector<double> &
  get_reynolds_stresses()
  {
    return reynolds_stresses;
  }

  void
  initialize_averaged_vectors(
    parallel::DistributedTriangulationBase<dim> &triangulation,
    const unsigned int &                         velocity_fem_degree,
    const DofsType &                             locally_owned_dofs,
    const DofsType &                             locally_relevant_dofs,
    const MPI_Comm &                             mpi_communicator);

  /**
   * @brief get_reynolds_stress. Gives the time-averaged Reynolds
   * stresses.
   */
  DoFHandler<dim> &
  get_reynolds_stress_handler()
  {
    return handler_rs;
  }

private:
  TrilinosScalar inv_range_time;
  TrilinosScalar dt_0;

  VectorType velocity_dt;
  VectorType sum_velocity_dt;
  VectorType average_velocities;

  LinearAlgebra::distributed::Vector<double> reynolds_stress_dt;
  LinearAlgebra::distributed::Vector<double> sum_reynolds_stress_dt;
  LinearAlgebra::distributed::Vector<double> reynolds_stresses;

  double dt;
  bool   average_calculation;
  double real_initial_time;

  DoFHandler<dim> handler_rs;
};

#endif
