// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_average_scalar_in_time_h
#define lethe_average_scalar_in_time_h

#include <core/parameters.h>

#include <deal.II/distributed/solution_transfer.h>

using namespace dealii;

/**
 * @brief Calculates the time-averaged scalar field.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved.
 *
 * @tparam VectorType The vector type used for the solvers (Trilinos (normal or block)
 * or deal.II vectors).
 *
 * @tparam DofsType The type of storage of indices which differs according to vector.
 *
 */

template <int dim>
class AverageScalar
{
public:
  /**
   * @brief Constructor that initializes the solution transfer objects for the scalar field.
   *
   * @param[in] dof_handler Used to initialize the solution transfer objects.
   * The solution transfer object is later used for mesh adaptation.
   */
  AverageScalar(const DoFHandler<dim> &dof_handler);


  /**
   * @brief Calculate time-averaged scalar field using vector with no ghost
   * cells.
   *
   * @param[in] local_evaluation_point The solution vector with no ghost cells.
   *
   * @param[in] post_processing The postprocessing parameters to get the initial
   * time used for averaging.
   *
   * @param[in] current_time The current time.
   *
   * @param[in] time_step The time step.
   */
  void
  calculate_average_scalar(const GlobalVectorType &local_evaluation_point,
                           const Parameters::PostProcessing &post_processing,
                           const double                      current_time,
                           const double                      time_step);

  /**
   * @brief Calculate time-averaged scalar field using vector with no ghost
   * cells.
   *
   * @param[in] locally_owned_dofs Degrees of freedom owned by a single
   * processor.
   *
   * @param[in] locally_relevant_dofs Degrees of freedom that live on locally
   * owned or ghost cells of a single processor.
   *
   * @param[in] mpi_communicator The communicator information.
   */
  void
  initialize_vectors(const IndexSet &locally_owned_dofs,
                     const IndexSet &locally_relevant_dofs,
                     const MPI_Comm &mpi_communicator);

  /**
   * @brief Getter for the average scalar.
   *
   */
  GlobalVectorType
  get_average_scalar()
  {
    return average_scalar;
  }

  /**
   * @brief Prepare average scalar object for dynamic mesh adaptation.
   */
  void
  update_average_scalar();

  /**
   * @brief Prepare average scalar object for dynamic mesh adaptation.
   */
  void
  prepare_for_mesh_adaptation();

  /**
   * @brief Reestablish solution vectors after dynamic mesh adaptation.
   */
  void
  post_mesh_adaptation();

  /**
   * @brief Save checkpoints to continue the averaging after the restart.
   *
   * @param[in] prefix Name for checkpointing files.
   *
   * @return Vector containing the average values.
   */
  std::vector<const GlobalVectorType *>
  save(const std::string &prefix);

  /**
   * @brief Read checkpoints to continuing averaging after restart.
   *
   * @param[in] prefix Name for checkpoint files.
   *
   * @return Vector containing the average values.
   */
  std::vector<GlobalVectorType *>
  read(const std::string &prefix);


  /**
   * @brief Sanitize the average scalar object after a checkpoint has been read by
   * resetting all of the locally_owned vectors using the locally_relevant
   * vectors. This is necessary because only the locally_relevant_vectors are
   * saved, but the calculation routines expect that the content of the
   * locally_owned vectors match that of the locally_relevant vectors.
   */
  void
  sanitize_after_restart()
  {
    sum_scalar_dt = sum_scalar_dt_with_ghost_cells;
  }

  /**
   * @brief Reinitialize the values of sum_scalar_dt_with_ghost_cells vector to 0 and the has_started_averaging flag to false. If the initial time for average temperature and heat flux parameter is greater than the simulation time, the checkpointed time-averaged heat flux is ignored after a restart. This allows users to restart a simulation with a different averaging start time.
   */
  void
  zero_average_after_restart();

private:
  /**
   * @brief Vector to store the scalar field multiplied by time step.
   *
   */
  GlobalVectorType scalar_dt;

  /**
   * @brief Vector to store the sum of all the scalar multiplied by time step.
   *
   */
  GlobalVectorType sum_scalar_dt;

  /**
   * @brief Ghosted vector to store the sum of all the scalar multiplied by time step.
   *
   */
  GlobalVectorType sum_scalar_dt_with_ghost_cells;

  /**
   * @brief Object to transfer the averaged scalar in case of mesh adaptation.
   *
   */
  SolutionTransfer<dim, GlobalVectorType> solution_transfer_sum_scalar_dt;

  /**
   * @brief Time averaged scalar.
   */
  GlobalVectorType average_scalar;

  /**
   * @brief Time averaged scalar containing the locally owned and locally relevant dofs for all processes.
   */
  GlobalVectorType average_scalar_with_ghost_cells;

  /**
   * @brief Difference between the simulation starting time and the initial averaging time.
   */
  double real_initial_time;

  /**
   * @brief Vector to store the sum of all the scalar field multiplied by time step.
   *
   */
  double total_time_for_average;

  /**
   * @brief Time step.
   *
   */
  double dt;

  /**
   * @brief First time step is stored in dt_0 and added to total_time_for_average to ensure the first average calculation does not result in a division by zero.
   *
   */
  double dt_0;

  /**
   * @brief Inverse of the total time for averaging.
   *
   */
  double inv_range_time;

  /**
   * @brief Track whether we are within the averaging time period.
   *
   */
  bool has_started_averaging;
};

#endif
