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

// Distributed
#include <deal.II/distributed/solution_transfer.h>


// Lethe Includes
#include <core/parameters.h>

using namespace dealii;

/**
 * @brief Calculates the time-averaged velocities and pressure \f$(\langle u
 * \rangle, \langle v \rangle, \langle w \rangle, \langle p \rangle)\f$ and the
 * independent components of the Reynolds stresses tensor \f$(\langle u'u'
 * \rangle, \langle v'v' \rangle, \langle w'w' \rangle, \langle u'v' \rangle,
 * \langle v'w' \rangle, \langle w'u' \rangle)\f$. 
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
template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  /**
   * @brief Constructor that initializes the solution transfer objects for the velocity
   * and the Reynolds stresses, the total time for average and average
   * calculation variables.
   *
   * @param[in] dof_handler Used to initialize the solution transfer objects.
   *
   */
  AverageVelocities(DoFHandler<dim> &dof_handler);

  /**
   * @brief Calculate time-averaged velocities and pressure using vector with no ghost
   * cells.
   *
   * @param[in] local_evaluation_point The solution vector with no ghost cells.
   *
   * @param[in] post_processing The postprocessing parameters to definie initial
   * time of averaging.
   *
   * @param[in] current_time The current time in the simulation.
   *
   * @param[in] time_step The current time step.
   *
   */
  void
  calculate_average_velocities(
    const VectorType                 &local_evaluation_point,
    const Parameters::PostProcessing &post_processing,
    const double                     &current_time,
    const double                     &time_step);


  void
  update_average_velocities();

  /**
   * @brief Calculate normal and other resolved time-averaged Reynold stresses
   * \f$(\langle u'u' \rangle, \langle v'v' \rangle, \langle w'w' \rangle)\f$
   * and \f$(\langle u'v' \rangle, \langle v'w' \rangle, \langle w'u'
   * \rangle)\f$.
   *
   * @param[in] local_evaluation_point The solution vector with no ghost cells.
   */
  void
  calculate_reynolds_stresses(const VectorType &local_evaluation_point);

  /**
   * @brief Give the average of solutions with ghost cells.
   *
   * @return The vector of average solutions.
   *
   */
  VectorType &
  get_average_velocities()
  {
    return get_av = average_velocities;
  }

  /**
   * @brief Give the time-averaged Reynolds normal stresses with ghost cells.
   *
   * @return Vector with average normal Reynolds stresses.
   *
   */
  const VectorType &
  get_reynolds_normal_stresses()
  {
    return get_rns = reynolds_normal_stresses;
  }

  /**
   * @brief Give the time-averaged Reynolds shear stresses with ghost cells.
   *
   * @return Vector with average shear Reynolds stresses.
   *
   */
  const VectorType &
  get_reynolds_shear_stresses()
  {
    return get_rss = reynolds_shear_stresses;
  }

  /**
   * @brief Initialize all the important vectors for the average velocities and
   * reynolds stresses calculation.
   *
   * @param[in] locally_owned_dofs The owned dofs.
   *
   * @param[in] locally_owned_rs_components The owned Reynolds stress
   * components.
   *
   * @param[in] dofs_per_vertex The number of dofs per vortex (dim for block
   * vectors).
   *
   * @param[in] mpi_communicator The communicator information.
   *
   */
  void
  initialize_vectors(const DofsType     &locally_owned_dofs,
                     const DofsType     &locally_relevant_dofs,
                     const unsigned int &dofs_per_vertex,
                     const MPI_Comm     &mpi_communicator);

  /**
   * @brief Initialize all the sum vectors for the average velocities and reynolds
   *  stresses storage.
   *
   * @param[in] locally_owned_dofs The owned dofs.
   *
   * @param[in] locally_owned_rs_components The owned Reynolds stress
   * components.
   *
   * @param[in] mpi_communicator The communicator information.
   */
  void
  initialize_checkpoint_vectors(const DofsType &locally_owned_dofs,
                                const DofsType &locally_relevant_dofs,
                                const MPI_Comm &mpi_communicator);


  /**
   * @brief Prepare average velocity object for dynamic mesh adaptation.
   */
  void
  prepare_for_mesh_adaptation();

  /**
   * @brief Reestablish solution vectors after dynamic mesh adaptation.
   */
  void
  post_mesh_adaptation();

  /**
   * @brief Save checkpoints to continuing averaging after restart.
   *
   * @param[in] prefix Name for checkpointing files.
   *
   * @return Vector with average values.
   */
  std::vector<const VectorType *>
  save(std::string prefix);

  /**
   * @brief Read checkpoints to continuing averaging after restart.
   *
   * @param[in] prefix Name for checkpoint files.
   *
   * @return Vector with average values.
   */
  std::vector<VectorType *>
  read(std::string prefix);

private:
  /**
   * @brief Inverse for total time for averaging.
   *
   */
  double inv_range_time;

  /**
   * @brief Store the first time step in case it varies.
   *
   */
  double dt_0;

  /**
   * @brief Vector to store the velocity multiplied by time step.
   *
   */
  VectorType velocity_dt;

  /**
   * @brief Vector to store the sum of all the velocities multiplied by time step.
   *
   */
  VectorType sum_velocity_dt;

  /**
   * @brief Vector to store average velocities.
   *
   */
  VectorType average_velocities;

  /**
   * @brief Getter vector for the average velocities.
   *
   */
  VectorType get_av;

  /**
   * @brief Vector to store the velocities multiplied by the time step to
   * calculate the normal stresses.
   *
   */
  VectorType reynolds_normal_stress_dt;

  /**
   * @brief Vector to store sum of the velocities multiplied by the time step to
   * calculate the normal stresses.
   *
   */
  VectorType sum_reynolds_normal_stress_dt;

  /**
   * @brief Vector with average Reynolds normal stresses.
   *
   */
  VectorType reynolds_normal_stresses;

  /**
   * @brief Getter vector for the average Reynolds normal stresses.
   *
   */
  VectorType get_rns;

  /**
   * @brief Vector to store the velocities multiplied by the time step to
   * calculate the shear stresses.
   *
   */
  VectorType reynolds_shear_stress_dt;

  /**
   * @brief Vector to store the sum the velocities multiplied by the time step to
   * calculate the shear stresses.
   *
   */
  VectorType sum_reynolds_shear_stress_dt;

  /**
   * @brief Vector with average average Reynolds shear stresses.
   *
   */
  VectorType reynolds_shear_stresses;

  /**
   * @brief Getter vector for the average Reynolds shear stresses.
   *
   */
  VectorType get_rss;

  /**
   * @brief Ghosted vector to store the sum of all the velocities multiplied by time step.
   *
   */
  VectorType sum_velocity_dt_with_ghost_cells;

  /**
   * @brief Ghosted vector to store sum of the velocities multiplied by the time step to
   * calculate the normal stresses.
   *
   */
  VectorType sum_rns_dt_with_ghost_cells;

  /**
   * @brief Ghosted vector to store the sum the velocities multiplied by the time step to
   * calculate the shear stresses.
   *
   */
  VectorType sum_rss_dt_with_ghost_cells;

  /**
   * @brief Object used to transfer the average velocities in case of mesh adaptation.
   *
   */
  parallel::distributed::SolutionTransfer<dim, VectorType>
    solution_transfer_sum_velocity_dt;

  /**
   * @brief Object used to transfer the average Reynolds normal stresses in case of
   * mesh adaptation.
   *
   */
  parallel::distributed::SolutionTransfer<dim, VectorType>
    solution_transfer_sum_reynolds_normal_stress_dt;

  /**
   * @brief Object used to transfer the average Reynolds shear stresses in case of
   * mesh adaptation.
   *
   */
  parallel::distributed::SolutionTransfer<dim, VectorType>
    solution_transfer_sum_reynolds_shear_stress_dt;

  /**
   * @brief Time step.
   *
   */
  double dt;

  /**
   * @brief Initial time of the simulation.
   *
   */
  double real_initial_time;

  /**
   * @brief Total period of time for averaging.
   *
   */
  double total_time_for_average;

  /**
   * @brief Track whether we are within the averaging time period.
   *
   */
  bool average_calculation;

  /**
   * @brief Number of degrees of freedom per vertex needed to calculate correctly
   * Reynolds stresses.
   *
   */
  unsigned int n_dofs_per_vertex;
};

#endif
