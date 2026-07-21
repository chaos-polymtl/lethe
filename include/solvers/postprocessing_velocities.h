// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_postprocessing_velocities_h
#define lethe_postprocessing_velocities_h

#include <core/parameters.h>
#include <core/vector.h>

#include <deal.II/distributed/solution_transfer.h>

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
   * @return Shared pointer to the vector of average solutions.
   *
   */
  std::shared_ptr<VectorType>
  get_average_velocities()
  {
    *get_av = average_velocities;
    return get_av;
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
   * reynolds stresses calculation. This includes the vector for the
   * checkpoints.
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
  save(const std::string &prefix);

  /**
   * @brief Read checkpoints to continuing averaging after restart.
   *
   * @param[in] prefix Name for checkpoint files.
   *
   * @return Vector with average values.
   */
  std::vector<VectorType *>
  read(const std::string &prefix);


  /**
   * @brief Sanitize the average velocity object after a checkpoint has been read by
   * resetting all of the locally_owned vectors using the locally_relevant
   * vectors. This is necessary because only the locally_relevant_vectors are
   * saved, but the calculation routines expect that the content of the
   * locally_owned vectors match that of the locally_relevant vectors.
   */
  void
  sanitize_after_restart()
  {
    sum_velocity_dt               = sum_velocity_dt_with_ghost_cells;
    sum_reynolds_normal_stress_dt = sum_rns_dt_with_ghost_cells;
    sum_reynolds_shear_stress_dt  = sum_rss_dt_with_ghost_cells;
  }

  /**
   * @brief Reinitialize the values of sum_velocities_dt_with_ghost_cells vector to 0 and the has_started_averaging flag to false. If the initial time for average velocities parameter is greater than the simulation time, the checkpointed time average velocities are ignored after a restart. This allows users to restart a simulation with a different averaging start time.
   */
  void
  zero_average_after_restart();

private:
  /**
   * @brief Type of the underlying non-block vector of VectorType.
   *
   * For block vectors, this is the type of an individual block. For regular
   * vectors, this is VectorType itself. The Reynolds stresses are always
   * accumulated block by block, hence the need for this type.
   */
  using UnderlyingVectorType =
    typename LetheVectorTraits<VectorType>::UnderlyingVectorType;

  /**
   * @brief Accumulate the Reynolds stresses multiplied by the time step over the locally owned degrees of freedom.
   *
   * This gathers the loop that is common to all the supported vector types.
   * The caller is responsible for extracting the relevant block of each vector
   * and for providing the mapping to the turbulence kinetic energy index,
   * which both differ between regular and block vectors.
   *
   * @param[in] local_solution Solution over which the fluctuations are computed.
   *
   * @param[in] local_average Time-averaged solution.
   *
   * @param[in] k_index Mapping from the index of the first velocity component
   * of a vertex to the index of the turbulence kinetic energy of that vertex.
   *
   * @param[in,out] rns_dt Normal Reynolds stresses multiplied by the time step.
   *
   * @param[in,out] rss_dt Shear Reynolds stresses multiplied by the time step.
   *
   * @param[in,out] k_dt Turbulence kinetic energy multiplied by the time step.
   */
  void
  accumulate_reynolds_stresses(const UnderlyingVectorType &local_solution,
                               const UnderlyingVectorType &local_average,
                               unsigned int (*k_index)(unsigned int),
                               UnderlyingVectorType &rns_dt,
                               UnderlyingVectorType &rss_dt,
                               UnderlyingVectorType &k_dt);

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
   * @brief Shared pointer to the getter vector for the average velocities.
   *
   */
  std::shared_ptr<VectorType> get_av;

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
  SolutionTransfer<dim, VectorType> solution_transfer_sum_velocity_dt;

  /**
   * @brief Object used to transfer the average Reynolds normal stresses in case of
   * mesh adaptation.
   *
   */
  SolutionTransfer<dim, VectorType>
    solution_transfer_sum_reynolds_normal_stress_dt;

  /**
   * @brief Object used to transfer the average Reynolds shear stresses in case of
   * mesh adaptation.
   *
   */
  SolutionTransfer<dim, VectorType>
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
  bool has_started_averaging;

  /**
   * @brief Number of degrees of freedom per vertex needed to calculate correctly
   * Reynolds stresses.
   *
   */
  unsigned int n_dofs_per_vertex;
};

#endif
