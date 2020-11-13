#include <solvers/postprocessing_velocities.h>


template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType &                        local_evaluation_point,
  const std::shared_ptr<SimulationControl> &simulation_control,
  const Parameters::PostProcessing &        post_processing,
  const DofsType &                          locally_owned_dofs,
  const MPI_Comm &                          mpi_communicator)
{
  double epsilon = 1e-6;
  double total_time =
    simulation_control->get_current_time() - post_processing.initial_time;
  double dt = simulation_control->get_time_step();

  // When averaging velocities begins (time = initial time)
  if (abs(total_time) < epsilon)
    {
      // Store the first dt value in case dt varies.
      dt_0 = dt;

      // Reinitialisation of the vectors to get the right length
      velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
      sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
      average_velocities.reinit(locally_owned_dofs, mpi_communicator);
    }

  // Calculate (u*dt) at each time step and accumulate the values
  velocity_dt.equ(dt, local_evaluation_point);
  sum_velocity_dt += velocity_dt;

  // Get the inverse of the total time with the first time step since
  // the weighted velocities are calculated with the first velocity when
  // total time = 0.
  inv_range_time = 1. / (total_time + dt_0);

  // Calculate average_velocities. (u*dt) / (total time + dt)
  // The sum of all weighted velocities in time / total time calculated
  average_velocities.equ(inv_range_time, sum_velocity_dt);
}

template <int dim, typename VectorType, typename DofsType>
const VectorType
AverageVelocities<dim, VectorType, DofsType>::get_average_velocities()
{
  return average_velocities;
}

template class AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<2,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;
