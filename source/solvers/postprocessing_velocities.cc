#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

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
  const double &                    time_step)
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
  unsigned int begin_index, end_index, n_dofs_per_node, ij_factor;

  const TrilinosWrappers::MPI::Vector *local_solution;
  const TrilinosWrappers::MPI::Vector *local_average;

  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
    {
      begin_index     = local_evaluation_point.local_range().first;
      end_index       = local_evaluation_point.local_range().second;
      n_dofs_per_node = dim + 1;
      ij_factor       = dim;
      local_solution  = &local_evaluation_point;
      local_average   = &average_velocities;
    }
  else if constexpr (std::is_same_v<VectorType,
                                    TrilinosWrappers::MPI::BlockVector>)
    {
      begin_index     = local_evaluation_point.block(0).local_range().first;
      end_index       = local_evaluation_point.block(0).local_range().second;
      n_dofs_per_node = dim;
      ij_factor       = dim + 1;
      local_solution  = &local_evaluation_point.block(0);
      local_average   = &average_velocities.block(0);
    }

  for (unsigned int i = begin_index; i < end_index; i += n_dofs_per_node)
    {
      // Set index from solution vector to reynolds stresses vector.
      // ij_factor is not the same for vector and block vector because
      // the number of solution in packages in the solution vector are not
      // the same, as explained in the comment above this function.
      unsigned int j = i * ij_factor / 2;

      // u'u'*dt
      reynolds_stress_dt[j] = ((*local_solution)[i] - (*local_average)[i]) *
                              ((*local_solution)[i] - (*local_average)[i]) * dt;

      // v'v'*dt
      reynolds_stress_dt[j + 1] =
        ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
        ((*local_solution)[i + 1] - (*local_average)[i + 1]) * dt;

      // u'v'*dt
      reynolds_stress_dt[j + dim] =
        ((*local_solution)[i] - (*local_average)[i]) *
        ((*local_solution)[i + 1] - (*local_average)[i + 1]) * dt;

      if (dim == 3)
        {
          // w'w'*dt
          reynolds_stress_dt[j + 2] =
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

          // v'w'*dt
          reynolds_stress_dt[j + 4] =
            ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

          // w'u'*dt
          reynolds_stress_dt[j + 5] =
            ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
            ((*local_solution)[i] - (*local_average)[i]) * dt;
        }
    }

  // Sum of all reynolds stresses during simulation.
  sum_reynolds_stress_dt += reynolds_stress_dt;

  // Calculate the reynolds stresses.
  reynolds_stresses.equ(inv_range_time, sum_reynolds_stress_dt);
}


template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::initialize_vectors(
  parallel::DistributedTriangulationBase<dim> &triangulation,
  const unsigned int &                         velocity_fem_degree,
  const DofsType &                             locally_owned_dofs,
  const DofsType &                             locally_relevant_dofs,
  const MPI_Comm &                             mpi_communicator)
{
  if (dim == 2)
    {
      FESystem<dim> fe_rs(FE_Q<dim>(velocity_fem_degree),
                          dim,
                          FE_Q<dim>(velocity_fem_degree),
                          1);

#if (DEAL_II_BOOST_VERSION_MINOR == 2)
      this->handler_rs.initialize(triangulation, fe_rs);
#else
      this->handler_rs.reinit(triangulation);
      this->handler_rs.distribute_dofs(fe_rs);
#endif
    }
  else if (dim == 3)
    {
      FESystem<dim> fe_rs(FE_Q<dim>(velocity_fem_degree),
                          dim,
                          FE_Q<dim>(velocity_fem_degree),
                          1,
                          FE_Q<dim>(velocity_fem_degree),
                          1,
                          FE_Q<dim>(velocity_fem_degree),
                          1);

#if (DEAL_II_BOOST_VERSION_MINOR == 2)
      this->handler_rs.initialize(triangulation, fe_rs);
#else
      this->handler_rs.reinit(triangulation);
      this->handler_rs.distribute_dofs(fe_rs);
#endif
    }

  IndexSet locally_owned_rs_components = handler_rs.locally_owned_dofs();
  IndexSet locally_relevant_rs_components;
  DoFTools::extract_locally_relevant_dofs(handler_rs,
                                          locally_relevant_rs_components);

  // Reinitialisation of the average velocity and reynolds stress vectors
  // to get the right length.
  velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  average_velocities.reinit(locally_owned_dofs,
                            locally_relevant_dofs,
                            mpi_communicator);

  // Reinitialize independent components of stress tensor vectors.
  reynolds_stress_dt.reinit(locally_owned_rs_components, mpi_communicator);
  sum_reynolds_stress_dt.reinit(locally_owned_rs_components, mpi_communicator);
  reynolds_stresses.reinit(locally_owned_rs_components,
                           locally_relevant_rs_components,
                           mpi_communicator);
}

template class AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<2,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;
