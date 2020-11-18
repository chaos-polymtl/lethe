#include <solvers/postprocessing_velocities.h>

template <int dim, typename VectorType, typename DofsType>
AverageVelocities<dim, VectorType, DofsType>::AverageVelocities()
  : average_calculation(false)
{}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType &                local_evaluation_point,
  const Parameters::PostProcessing &post_processing,
  const double &                    current_time,
  const double &                    time_step,
  const bool &                      is_output_iteration,
  const DofsType &                  locally_owned_dofs,
  const MPI_Comm &                  mpi_communicator)
{
  const double epsilon      = 1e-6;
  const double initial_time = post_processing.initial_time;
  dt                        = time_step;

  // When averaging velocities begins
  if (current_time >= (initial_time - epsilon) && !average_calculation)
    {
      average_calculation = true;
      real_initial_time   = current_time;

      // Store the first dt value in case dt varies.
      dt_0 = dt;

      // Reinitialisation of the average velocity vectors to get the right
      // length
      velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
      sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
      average_velocities.reinit(locally_owned_dofs, mpi_communicator);

      // Reinitialisation of reynolds stress vectors if calculation is enable
      if (post_processing.calculate_reynolds_stress)
        {
          fluctuation.reinit(locally_owned_dofs, mpi_communicator);
          reynolds_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
          sum_reynolds_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
          reynolds_stress.reinit(locally_owned_dofs, mpi_communicator);
        }
    }

  // Calculate (u*dt) at each time step and accumulate the values
  velocity_dt.equ(dt, local_evaluation_point);
  sum_velocity_dt += velocity_dt;

  // Get the inverse of the total time with the first time step since
  // the weighted velocities are calculated with the first velocity when
  // total time = 0.
  inv_range_time = 1. / ((current_time - real_initial_time) + dt_0);


  // If calculate reynolds stress is enable, average velocities must be
  // calculated at every iteration since the average values are needed to get
  // reynolds stress. If it is disabled, average velocities only need to be
  // calculated when output is required.
  if (post_processing.calculate_reynolds_stress)
    {
      average_velocities.equ(inv_range_time, sum_velocity_dt);

      this->calculate_reynolds_stress(local_evaluation_point,
                                      is_output_iteration);
    }
  else if (is_output_iteration)
    {
      // Calculate average_velocities. (u*dt) / (total time + dt)
      // The sum of all weighted velocities in time / total time calculated
      average_velocities.equ(inv_range_time, sum_velocity_dt);
    }
}

// Since Trilinos vectors and block vectors data doesn't have the same
// structure, those vectors need to be processed in different ways.
// Same about the dimension, where 2D solutions are stored in (u,v,p) groups
// for vectors and [(u,v)][p] for block vectors and where 3D solutions
// are stored in (u,v,w,p) groups for vectors and [(u,v,w)][p] for
// block vectors.
template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stress(
  const VectorType &local_evaluation_point,
  const bool &      is_output_iteration)
{
  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      const unsigned int begin_index =
        local_evaluation_point.local_range().first;
      const unsigned int end_index =
        local_evaluation_point.local_range().second;

      for (unsigned int i = begin_index; i < end_index; i++)
        {
          if (dim == 2 && (i + 3) % 3 == 0)
            {
              reynolds_stress_dt[i] =
                (local_evaluation_point[i] - average_velocities[i]) *
                (local_evaluation_point[i] - average_velocities[i]) * dt;
              reynolds_stress_dt[i + 1] =
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                dt;
              reynolds_stress_dt[i + 2] =
                (local_evaluation_point[i] - average_velocities[i]) *
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                dt;
            }
          else if (dim == 3 && (i + 4) % 4 == 0)
            {
              // Calculating (u'u')*dt, (v'v')*dt (w'w')*dt and (u'v')*dt
              reynolds_stress_dt[i] =
                (local_evaluation_point[i] - average_velocities[i]) *
                (local_evaluation_point[i] - average_velocities[i]) * dt;
              reynolds_stress_dt[i + 1] =
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                dt;
              reynolds_stress_dt[i + 2] =
                (local_evaluation_point[i + 2] - average_velocities[i + 2]) *
                (local_evaluation_point[i + 2] - average_velocities[i + 2]) *
                dt;
              reynolds_stress_dt[i + 3] =
                (local_evaluation_point[i] - average_velocities[i]) *
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                dt;
            }
        }
      // Summation of all reynolds stress during simulation
      sum_reynolds_stress_dt += reynolds_stress_dt;
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::BlockVector>)
    {
      fluctuation.block(0) = local_evaluation_point.block(0);
      fluctuation.block(0) -= average_velocities.block(0);
      reynolds_stress_dt.block(0) = fluctuation.block(0);
      reynolds_stress_dt.block(0).scale(reynolds_stress_dt.block(0));
      reynolds_stress_dt.block(0) *= dt;

      unsigned int begin_index =
        local_evaluation_point.block(0).local_range().first;
      unsigned int end_index =
        local_evaluation_point.block(0).local_range().second;
      unsigned int size_velocity_block = reynolds_stress.block(0).size();

      for (unsigned int i = begin_index; i < end_index; i++)
        {
          if (dim == 2 && (i + 2) % 2 == 0)
            {
              reynolds_stress_dt.block(1)[i / 2 + size_velocity_block] =
                fluctuation.block(0)[i] * fluctuation.block(0)[i + 1] * dt;
            }
          else if (dim == 3 && (i + 3) % 3 == 0)
            {
              reynolds_stress_dt.block(1)[i / 3 + size_velocity_block] =
                fluctuation.block(0)[i] * fluctuation.block(0)[i + 1] * dt;
            }
        }

      sum_reynolds_stress_dt += reynolds_stress_dt;
    }
  // Calculating time-averaged reynolds stress
  if (is_output_iteration)
    reynolds_stress.equ(inv_range_time, sum_reynolds_stress_dt);
}


template class AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<2,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;
