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
 * @brief AverageVelocities. The AverageVelocities class calculates the
 * time-averaged velocities and pressure (<u>, <v>, <w>, <p>) and the
 * independent components of the Reynolds stresses tensor (<u'u'>, <v'v'>,
 * <w'w'>, <u'v'>, <v'w'>, <w'u'>). The generated vectors are displayable of
 * visualization software.
 */
template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  AverageVelocities(DoFHandler<dim> &dof_handler);
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
   * @brief get_average_velocities. Gives the average of solutions with ghost
   * cells.
   */
  const VectorType &
  get_average_velocities()
  {
    return get_av = average_velocities;
  }

  /**
   * @brief get_reynolds_normal_stresses. Gives the time-averaged Reynolds
   * normal stresses with ghost cells.
   */
  const VectorType &
  get_reynolds_normal_stresses()
  {
    return get_rns = reynolds_normal_stresses;
  }

  /**
   * @brief get_reynolds_shear_stress. Gives the time-averaged Reynolds
   * shear stresses with ghost cells.
   */
  const VectorType &
  get_reynolds_shear_stresses()
  {
    return get_rss = reynolds_shear_stresses;
  }

  /**
   * @brief initialize_vectors. This function initializes all the important
   * vectors for the average velocities and reynolds stresses calculation.
   *
   * @param locally_owned_dofs. The owned dofs
   *
   * @param locally_owned_rs_components. The owned Reynolds stress components
   *
   * @param dofs_per_vertex. The number of dofs per vortex (dim for block
   *                         vectors)
   *
   * @param mpi_communicator. The mpi communicator information
   */
  void
  initialize_vectors(const DofsType &    locally_owned_dofs,
                     const DofsType &    locally_relevant_dofs,
                     const unsigned int &dofs_per_vertex,
                     const MPI_Comm &    mpi_communicator);

  /**
   * @brief initialize_checkpoint_vectors. This function initializes all the
   * sum vectors for the average velocities and reynolds stresses storage.
   *
   * @param locally_owned_dofs. The owned dofs
   *
   * @param locally_owned_rs_components. The owned Reynolds stress components
   *
   * @param mpi_communicator. The mpi communicator information
   */
  void
  initialize_checkpoint_vectors(const DofsType &locally_owned_dofs,
                                const DofsType &locally_relevant_dofs,
                                const MPI_Comm &mpi_communicator);


  /**
   * @brief Prepares average velocity object for dynamic mesh adaptation
   */
  void
  prepare_for_mesh_adaptation();

  /**
   * @brief Re-establish solution vectors after dynamic mesh adaptation
   */
  void
  post_mesh_adaptation();

  /**
   * @brief save & read. Save or read checkpoints to continuing averaging after
   * restart.
   */
  std::vector<const VectorType *>
  save(std::string prefix);

  std::vector<VectorType *>
  read(std::string prefix);

private:
  TrilinosScalar inv_range_time;
  TrilinosScalar dt_0;

  VectorType velocity_dt;
  VectorType sum_velocity_dt;
  VectorType average_velocities;
  VectorType get_av;

  VectorType reynolds_normal_stress_dt;
  VectorType sum_reynolds_normal_stress_dt;
  VectorType reynolds_normal_stresses;
  VectorType get_rns;

  VectorType reynolds_shear_stress_dt;
  VectorType sum_reynolds_shear_stress_dt;
  VectorType reynolds_shear_stresses;
  VectorType get_rss;

  VectorType sum_velocity_dt_with_ghost_cells;
  VectorType sum_rns_dt_with_ghost_cells;
  VectorType sum_rss_dt_with_ghost_cells;

  // Solution transfer for the three permanent velocity storage
  parallel::distributed::SolutionTransfer<dim, VectorType>
    solution_transfer_sum_velocity_dt;
  parallel::distributed::SolutionTransfer<dim, VectorType>
    solution_transfer_sum_reynolds_normal_stress_dt;
  parallel::distributed::SolutionTransfer<dim, VectorType>
    solution_transfer_sum_reynolds_shear_stress_dt;

  double       dt;
  double       real_initial_time;
  bool         average_calculation;
  unsigned int n_dofs_per_vertex;
};

#endif
