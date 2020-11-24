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
  const DofsType &                  locally_owned_dofs,
  const MPI_Comm &                  mpi_communicator)
{
  const double epsilon      = 1e-6;
  const double initial_time = post_processing.initial_time;
  dt                        = time_step;
  //locally_owned_tensor_components = get_tensor_index_set(locally_owned_dofs);

  // When averaging velocities begins
  if (current_time >= (initial_time - epsilon) && !average_calculation)
    {
      average_calculation = true;
      real_initial_time   = current_time;

      // Store the first dt value in case dt varies.
      dt_0 = dt;

      // Reinitialisation of the average velocity and reynolds stress vectors
      // to get the right length.
      velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
      sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
      average_velocities.reinit(locally_owned_dofs, mpi_communicator);

      // Reinitialize local independant components for calculation
      // (local_average_solution as 4 components and there's 6 independant
      // components for the Reynolds tensor.

      reynolds_stress_dt.reinit(locally_owned_tensor_components,
                                mpi_communicator);
      sum_reynolds_stress_dt.reinit(locally_owned_tensor_components,
                                    mpi_communicator);
      reynolds_stresses.reinit(locally_owned_tensor_components,
                               mpi_communicator);
    }

  // Calculate (u*dt) at each time step and accumulate the values
  velocity_dt.equ(dt, local_evaluation_point);
  sum_velocity_dt += velocity_dt;

  // Get the inverse of the total time with the first time step since
  // the weighted velocities are calculated with the first velocity when
  // total time = 0.
  inv_range_time = 1. / ((current_time - real_initial_time) + dt_0);

  // Calculate the average velocities.
  average_velocities.equ(inv_range_time, sum_velocity_dt);

  this->calculate_reynolds_stresses(local_evaluation_point);
}

// Since Trilinos vectors and block vectors data doesn't have the same
// structure, those vectors need to be processed in different ways.
// Same about the dimension, where 2D solutions are stored in (u,v,p) groups
// for vectors and [(u,v)][p] for block vectors and where 3D solutions
// are stored in (u,v,w,p) groups for vectors and [(u,v,w)][p] for
// block vectors.
template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stresses(
  const VectorType &local_evaluation_point)
{
  unsigned int begin_index, end_index;

  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      begin_index = local_evaluation_point.local_range().first;
      end_index   = local_evaluation_point.local_range().second;

      for (unsigned int i = begin_index; i < end_index; i += dim + 1)
        {
          // Set related index from solution vector to reynolds tensor
          // which is stored in a data vector.
          unsigned int j = get_tensor_index(i);

          // u'u'*dt
          reynolds_stress_dt[j] =
            (local_evaluation_point[i] - average_velocities[i]) *
            (local_evaluation_point[i] - average_velocities[i]) * dt;

          // u'v'*dt
          reynolds_stress_dt[j + 1] =
            (local_evaluation_point[i] - average_velocities[i]) *
            (local_evaluation_point[i + 1] - average_velocities[i + 1]) * dt;

          // v'u'*dt
          reynolds_stress_dt[j + dim] = reynolds_stress_dt[j + dim];

          // v'v'*dt
          reynolds_stress_dt[j + dim + 1] =
            (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
            (local_evaluation_point[i + 1] - average_velocities[i + 1]) * dt;


          if (dim == 3)
            {
              // u'w'*dt
              reynolds_stress_dt[j + 2] =
                (local_evaluation_point[i] - average_velocities[i]) *
                (local_evaluation_point[i + 2] - average_velocities[i + 2]) * dt;

              // v'w'*dt
              reynolds_stress_dt[j + 5] =
                (local_evaluation_point[i + 1] - average_velocities[i + 1]) *
                (local_evaluation_point[i + 2] - average_velocities[i + 2]) *
                dt;

              // w'u'*dt
              reynolds_stress_dt[j + 6] = reynolds_stress_dt[j + 2];

              // w'v'*dt
              reynolds_stress_dt[j + 7] = reynolds_stress_dt[j + 5];

              // w'w'*dt
              reynolds_stress_dt[j + 8] =
                (local_evaluation_point[i + 2] - average_velocities[i + 2]) *
                (local_evaluation_point[i + 2] - average_velocities[i + 2]) *
                dt;
            }
        }
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::BlockVector>)
    {
      begin_index = local_evaluation_point.block(0).local_range().first;
      end_index   = local_evaluation_point.block(0).local_range().second;


    }

  // Sum of all reynolds stress during simulation.
  sum_reynolds_stress_dt += reynolds_stress_dt;

  // Calculate the reynolds stresses.
  reynolds_stresses.equ(inv_range_time, sum_reynolds_stress_dt);
}

template <int dim, typename VectorType, typename DofsType>
unsigned int
AverageVelocities<dim, VectorType, DofsType>::get_tensor_index(unsigned int i)
{
  unsigned int j;
  j = i / (dim + 1) * dim * dim;

  return j;
}

template <int dim, typename VectorType, typename DofsType>
IndexSet
AverageVelocities<dim, VectorType, DofsType>::get_tensor_index_set(
  const DofsType &locally_owned_dofs,
  const unsigned int &n_dofs)
{
  IndexSet locally_owned_tensor(9 * n_dofs / 4);

  unsigned int first_index, last_index;

  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      first_index = locally_owned_dofs.nth_index_in_set(0);
      last_index  = first_index + locally_owned_dofs.n_elements();
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::BlockVector>)
    {
      first_index = locally_owned_dofs[0].nth_index_in_set(0);
      last_index  = first_index + locally_owned_dofs[0].n_elements();
    }

  first_index = get_tensor_index(first_index);
  last_index  = get_tensor_index(last_index);

  std::cout << "----------------------" << std::endl;
  std::cout << "first_index" << first_index << "last_index" << last_index << std::endl;

  locally_owned_tensor.add_range(first_index, last_index);

  locally_owned_tensor_components = locally_owned_tensor;

  return locally_owned_tensor_components;
}


template class AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;


template class AverageVelocities<2,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;
