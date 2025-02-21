// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/vector.h>

#include <solvers/postprocessing_scalar.h>

#include <fstream>

template <int dim>
AverageScalar<dim>::AverageScalar(const DoFHandler<dim> &dof_handler)
  : solution_transfer_sum_scalar_dt(dof_handler)
  , total_time_for_average(0.0)
  , average_calculation(false)
{}

template <int dim>
void
AverageScalar<dim>::calculate_average_scalar(
  const GlobalVectorType           &local_evaluation_point,
  const Parameters::PostProcessing &post_processing,
  const double                      current_time,
  const double                      time_step)
{
  const double epsilon = 1e-6;
  const double initial_time =
    post_processing.initial_time_for_average_temp_and_hf;
  dt = time_step;

  // When averaging scalar begins
  if (current_time >= (initial_time - epsilon) && !average_calculation)
    {
      average_calculation = true;
      real_initial_time   = current_time;

      // Store the first dt value in case dt varies.
      dt_0 = dt;
    }

  // Calculate (scalar*dt) at each time step and accumulates the values
  scalar_dt.equ(dt, local_evaluation_point);

  sum_scalar_dt += scalar_dt;

  // Get the inverse of the time since the beginning of the time averaging
  total_time_for_average = (current_time - real_initial_time) + dt_0;
  inv_range_time         = 1. / total_time_for_average;

  // Calculate the average scalars.
  average_scalar.equ(inv_range_time, sum_scalar_dt);
}

template <int dim>
void
AverageScalar<dim>::initialize_vectors(const IndexSet &locally_owned_dofs,
                                       const IndexSet &locally_relevant_dofs,
                                       const MPI_Comm &mpi_communicator)
{
  // Reinitialisation of the average scalars to get the right length.
  scalar_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_scalar_dt.reinit(locally_owned_dofs, mpi_communicator);
  average_scalar.reinit(locally_owned_dofs, mpi_communicator);

  average_scalar_with_ghost_cells.reinit(locally_owned_dofs,
                                         locally_relevant_dofs,
                                         mpi_communicator);

  sum_scalar_dt_with_ghost_cells.reinit(locally_owned_dofs,
                                        locally_relevant_dofs,
                                        mpi_communicator);
}

template <int dim>
void
AverageScalar<dim>::update_average_scalar()
{
  if (total_time_for_average > 1e-16)
    {
      inv_range_time = 1.0 / total_time_for_average;
      // Calculate the average scalar.
      average_scalar.equ(inv_range_time, sum_scalar_dt);
    }
}

template <int dim>
void
AverageScalar<dim>::prepare_for_mesh_adaptation()
{
  sum_scalar_dt_with_ghost_cells = sum_scalar_dt;
  solution_transfer_sum_scalar_dt.prepare_for_coarsening_and_refinement(
    sum_scalar_dt_with_ghost_cells);
}

template <int dim>
void
AverageScalar<dim>::post_mesh_adaptation()
{
  solution_transfer_sum_scalar_dt.interpolate(sum_scalar_dt);

  sum_scalar_dt_with_ghost_cells = sum_scalar_dt;

  update_average_scalar();
}

template <int dim>
std::vector<const GlobalVectorType *>
AverageScalar<dim>::save(const std::string &prefix)
{
  sum_scalar_dt_with_ghost_cells = sum_scalar_dt;

  std::vector<const GlobalVectorType *> avg_scalar_set_transfer;
  avg_scalar_set_transfer.push_back(&sum_scalar_dt_with_ghost_cells);

  std::string   filename = prefix + ".averagescalar";
  std::ofstream output(filename.c_str());
  output << "Average scalar" << std::endl;
  output << "dt_0 " << dt_0 << std::endl;
  output << "Average_calculation_boolean " << average_calculation << std::endl;
  output << "Real_initial_time " << real_initial_time << std::endl;

  return avg_scalar_set_transfer;
}

template <int dim>
std::vector<GlobalVectorType *>
AverageScalar<dim>::read(const std::string &prefix)
{
  std::vector<GlobalVectorType *> sum_vectors;
  sum_vectors.push_back(&sum_scalar_dt_with_ghost_cells);

  std::string   filename = prefix + ".averagescalar";
  std::ifstream input(filename.c_str());
  AssertThrow(input, ExcFileNotOpen(filename));

  std::string buffer;
  std::getline(input, buffer);
  input >> buffer >> dt_0;
  input >> buffer >> average_calculation;
  input >> buffer >> real_initial_time;

  return sum_vectors;
}

template <int dim>
void
AverageScalar<dim>::reinit_average_after_restart(const IndexSet &locally_owned_dofs,
                                                 const IndexSet &locally_relevant_dofs,
                                                 const MPI_Comm &mpi_communicator)
{
  sum_scalar_dt_with_ghost_cells.reinit(locally_owned_dofs,
                                        locally_relevant_dofs,
                                        mpi_communicator);             
  average_calculation = false;
}

template class AverageScalar<2>;
template class AverageScalar<3>;
