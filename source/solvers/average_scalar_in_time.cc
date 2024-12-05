// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/vector.h>

#include <solvers/average_scalar_in_time.h>

#include <fstream>

template <int dim>
AverageScalarInTime<dim>::AverageScalarInTime(
  DoFHandler<dim> &dof_handler)
  : solution_transfer_sum_scalar_dt(dof_handler)
  , total_time_for_average(0.0)
  , average_calculation(false)
{}

template <int dim>
void
AverageScalarInTime<dim>::calculate_average_scalar(
  const GlobalVectorType           &local_evaluation_point,
  const Parameters::PostProcessing &post_processing,
  const double                      current_time,
  const double                      time_step)
{
  const double epsilon      = 1e-6;
  const double initial_time = post_processing.initial_time;
  dt                        = time_step;

  std::cout<<"allo 0"<<std::endl;
  // When averaging velocities begins
  if (current_time >= (initial_time - epsilon) && !average_calculation)
    {
      average_calculation = true;
      real_initial_time   = current_time;

      // Store the first dt value in case dt varies.
      dt_0 = dt;
    }
  std::cout<< dt <<std::endl;

  std::cout<<local_evaluation_point.size()<<std::endl;
  // Calculate (scalar*dt) at each time step and accumulates the values
  scalar_dt.equ(dt, local_evaluation_point);

  sum_scalar_dt += scalar_dt;

  // Get the inverse of the time since the beginning of the time averaging
  total_time_for_average = (current_time - real_initial_time) + dt_0;
  inv_range_time         = 1. / total_time_for_average;

  // Calculate the average velocities.
  average_scalar.equ(inv_range_time, sum_scalar_dt);

}

template <int dim>
void
AverageScalarInTime<dim>::initialize_vectors(
  const IndexSet  &locally_owned_dofs,
  const IndexSet  &locally_relevant_dofs,
  const MPI_Comm         &mpi_communicator)
{
  // Reinitialisation of the average velocity and reynolds stress vectors
  // to get the right length.
  scalar_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_scalar_dt.reinit(locally_owned_dofs, mpi_communicator);
  average_scalar.reinit(locally_owned_dofs, mpi_communicator);

}

template class AverageScalarInTime<2>;
template class AverageScalarInTime<3>;