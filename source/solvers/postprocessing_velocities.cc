#include <solvers/postprocessing_velocities.h>

#include <fstream>

template <int dim, typename VectorType, typename DofsType>
AverageVelocities<dim, VectorType, DofsType>::AverageVelocities(
  DoFHandler<dim> &dof_handler)
  : solution_transfer_sum_velocity_dt(dof_handler)
  , solution_transfer_sum_reynolds_normal_stress_dt(dof_handler)
  , solution_transfer_sum_reynolds_shear_stress_dt(dof_handler)
  , total_time_for_average(0.0)
  , average_calculation(false)
{}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType                 &local_evaluation_point,
  const Parameters::PostProcessing &post_processing,
  const double                     &current_time,
  const double                     &time_step)
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

  // Get the inverse of the time since the beginning of the time averaging
  total_time_for_average = (current_time - real_initial_time) + dt_0;
  inv_range_time         = 1. / total_time_for_average;

  // Calculate the average velocities.
  average_velocities.equ(inv_range_time, sum_velocity_dt);

  this->calculate_reynolds_stresses(local_evaluation_point);
}


template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::update_average_velocities()
{
  // Use the inverse of the time since the beginning of the time averaging to
  // reevaluate the average velocity field and Reynolds stress.

  if (total_time_for_average > 1e-16)
    {
      inv_range_time = 1.0 / total_time_for_average;
      // Calculate the average velocities.
      average_velocities.equ(inv_range_time, sum_velocity_dt);
      reynolds_normal_stresses.equ(inv_range_time,
                                   sum_reynolds_normal_stress_dt);
      reynolds_shear_stresses.equ(inv_range_time, sum_reynolds_shear_stress_dt);
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
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stresses(
  const VectorType &local_evaluation_point)
{
  // TODO: generalize this for the use of
  // LinearAlegra::distributed::Vector<double>
  if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector> ||
                std::is_same_v<VectorType, TrilinosWrappers::MPI::BlockVector>)
    {
      unsigned int begin_index, end_index;

      // As explained in the comment above, about the structure of Trilinos
      // vectors, a function is required to get the index of the turbulence
      // kinetic energy. Trilinos vector : rss = (<u'u'>,<v'v'>,<w'w'>,k) =>
      // k_index = i + dim Trilinos blocks : rss = [(<u'u'>,<v'v'>,<w'w'>)][k])
      // => k_index = i / dim The k_index lambda is later changed for block
      // vectors.
      unsigned int (*k_index)(unsigned int) = [](unsigned int i) {
        return i + dim;
      };

      const TrilinosWrappers::MPI::Vector *local_solution, *local_average;
      TrilinosWrappers::MPI::Vector       *rns_dt, *rss_dt, *k_dt;

      if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
        {
          begin_index    = local_evaluation_point.local_range().first;
          end_index      = local_evaluation_point.local_range().second;
          local_solution = &local_evaluation_point;
          local_average  = &average_velocities;
          rns_dt         = &reynolds_normal_stress_dt;
          rss_dt         = &reynolds_shear_stress_dt;
          k_dt           = &reynolds_normal_stress_dt;
        }
      else if constexpr (std::is_same_v<VectorType,
                                        TrilinosWrappers::MPI::BlockVector>)
        {
          begin_index    = local_evaluation_point.block(0).local_range().first;
          end_index      = local_evaluation_point.block(0).local_range().second;
          local_solution = &local_evaluation_point.block(0);
          local_average  = &average_velocities.block(0);
          rns_dt         = &reynolds_normal_stress_dt.block(0);
          rss_dt         = &reynolds_shear_stress_dt.block(0);
          k_dt           = &reynolds_normal_stress_dt.block(1);
          k_index        = [](unsigned int i) { return i / dim; };
        }

      for (unsigned int i = begin_index; i < end_index; i += n_dofs_per_vertex)
        {
          // u'u'*dt
          (*rns_dt)[i] = ((*local_solution)[i] - (*local_average)[i]) *
                         ((*local_solution)[i] - (*local_average)[i]) * dt;

          // v'v'*dt
          (*rns_dt)[i + 1] =
            ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
            ((*local_solution)[i + 1] - (*local_average)[i + 1]) * dt;

          // u'v'*dt
          (*rss_dt)[i] = ((*local_solution)[i] - (*local_average)[i]) *
                         ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
                         dt;

          // k*dt = 1/2(u'u'+v'v')*dt (turbulence kinetic energy)
          // Note : k_dt and rns_dt are both pointers of
          // reynolds_normal_stress_dt for Trilinos vector (not block vectors)
          (*k_dt)[k_index(i)] = ((*rns_dt)[i] + (*rns_dt)[i + 1]) / 2;

          if (dim == 3)
            {
              // w'w'*dt
              (*rns_dt)[i + 2] =
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

              // v'w'*dt
              (*rss_dt)[i + 1] =
                ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

              // w'u'*dt
              (*rss_dt)[i + 2] =
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
                ((*local_solution)[i] - (*local_average)[i]) * dt;

              // k*dt = 1/2(u'u'+v'v'+w'w')*dt
              (*k_dt)[k_index(i)] = (*k_dt)[k_index(i)] + (*rns_dt)[i + 2] / 2;
            }
        }
    }

  if constexpr (std::is_same_v<VectorType,
                               LinearAlgebra::distributed::Vector<double>>)
    {
      unsigned int begin_index, end_index;
      unsigned int (*k_index)(unsigned int) = [](unsigned int i) {
        return i + dim;
      };
      const LinearAlgebra::distributed::Vector<double> *local_solution,
        *local_average;
      LinearAlgebra::distributed::Vector<double> *rns_dt, *rss_dt, *k_dt;

      local_evaluation_point.update_ghost_values();
      begin_index =
        local_evaluation_point.get_partitioner()->local_range().first;
      end_index =
        local_evaluation_point.get_partitioner()->local_range().second;
      local_solution = &local_evaluation_point;
      local_average  = &average_velocities;
      rns_dt         = &reynolds_normal_stress_dt;
      rss_dt         = &reynolds_shear_stress_dt;
      k_dt           = &reynolds_normal_stress_dt;

      for (unsigned int i = begin_index; i < end_index; i += n_dofs_per_vertex)
        {
          // u'u'*dt
          (*rns_dt)[i] = ((*local_solution)[i] - (*local_average)[i]) *
                         ((*local_solution)[i] - (*local_average)[i]) * dt;

          // v'v'*dt
          (*rns_dt)[i + 1] =
            ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
            ((*local_solution)[i + 1] - (*local_average)[i + 1]) * dt;

          // u'v'*dt
          (*rss_dt)[i] = ((*local_solution)[i] - (*local_average)[i]) *
                         ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
                         dt;

          // // k*dt = 1/2(u'u'+v'v')*dt (turbulence kinetic energy)
          // // Note : k_dt and rns_dt are both pointers of
          // // reynolds_normal_stress_dt for Trilinos vector (not block
          // vectors)
          (*k_dt)[k_index(i)] = ((*rns_dt)[i] + (*rns_dt)[i + 1]) / 2;

          if (dim == 3)
            {
              // w'w'*dt
              (*rns_dt)[i + 2] =
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

              // v'w'*dt
              (*rss_dt)[i + 1] =
                ((*local_solution)[i + 1] - (*local_average)[i + 1]) *
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) * dt;

              // w'u'*dt
              (*rss_dt)[i + 2] =
                ((*local_solution)[i + 2] - (*local_average)[i + 2]) *
                ((*local_solution)[i] - (*local_average)[i]) * dt;

              // k*dt = 1/2(u'u'+v'v'+w'w')*dt
              (*k_dt)[k_index(i)] = (*k_dt)[k_index(i)] + (*rns_dt)[i + 2] / 2;
            }
        }
    }

  // Sum of all reynolds stresses during simulation.
  sum_reynolds_normal_stress_dt += reynolds_normal_stress_dt;
  sum_reynolds_shear_stress_dt += reynolds_shear_stress_dt;

  // Calculate the reynolds stresses.
  reynolds_normal_stresses.equ(inv_range_time, sum_reynolds_normal_stress_dt);
  reynolds_shear_stresses.equ(inv_range_time, sum_reynolds_shear_stress_dt);
}


template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::initialize_vectors(
  const DofsType     &locally_owned_dofs,
  const DofsType     &locally_relevant_dofs,
  const unsigned int &dofs_per_vertex,
  const MPI_Comm     &mpi_communicator)
{
  // Save the number of dofs per vertex. If solution is in block vectors,
  // this is the number of dofs about velocity, dim.
  n_dofs_per_vertex = dofs_per_vertex;

  // Reinitialisation of the average velocity and reynolds stress vectors
  // to get the right length.
  velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  average_velocities.reinit(locally_owned_dofs, mpi_communicator);
  get_av.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  // Reinitialize independent components of stress tensor vectors.
  reynolds_normal_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_reynolds_normal_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  reynolds_normal_stresses.reinit(locally_owned_dofs, mpi_communicator);
  get_rns.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  reynolds_shear_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_reynolds_shear_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  reynolds_shear_stresses.reinit(locally_owned_dofs, mpi_communicator);
  get_rss.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::prepare_for_mesh_adaptation()
{
  get_av  = sum_velocity_dt;
  get_rns = sum_reynolds_normal_stress_dt;
  get_rss = sum_reynolds_shear_stress_dt;
  solution_transfer_sum_velocity_dt.prepare_for_coarsening_and_refinement(
    get_av);
  solution_transfer_sum_reynolds_normal_stress_dt
    .prepare_for_coarsening_and_refinement(get_rns);
  solution_transfer_sum_reynolds_shear_stress_dt
    .prepare_for_coarsening_and_refinement(get_rss);
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::post_mesh_adaptation()
{
  solution_transfer_sum_velocity_dt.interpolate(sum_velocity_dt);
  solution_transfer_sum_reynolds_normal_stress_dt.interpolate(
    sum_reynolds_normal_stress_dt);
  solution_transfer_sum_reynolds_shear_stress_dt.interpolate(
    sum_reynolds_shear_stress_dt);



  sum_velocity_dt_with_ghost_cells = sum_velocity_dt;
  sum_rns_dt_with_ghost_cells      = sum_reynolds_normal_stress_dt;
  sum_rss_dt_with_ghost_cells      = sum_reynolds_shear_stress_dt;

  update_average_velocities();
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::initialize_checkpoint_vectors(
  const DofsType &locally_owned_dofs,
  const DofsType &locally_relevant_dofs,
  const MPI_Comm &mpi_communicator)
{
  // Initializing vector with locally_relevant_dofs because writing checkpoint
  // needs to read ghost cells.
  sum_velocity_dt_with_ghost_cells.reinit(locally_owned_dofs,
                                          locally_relevant_dofs,
                                          mpi_communicator);
  sum_rns_dt_with_ghost_cells.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
  sum_rss_dt_with_ghost_cells.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
}

template <int dim, typename VectorType, typename DofsType>
std::vector<const VectorType *>
AverageVelocities<dim, VectorType, DofsType>::save(std::string prefix)
{
  sum_velocity_dt_with_ghost_cells = sum_velocity_dt;
  sum_rns_dt_with_ghost_cells      = sum_reynolds_normal_stress_dt;
  sum_rss_dt_with_ghost_cells      = sum_reynolds_shear_stress_dt;


  std::vector<const VectorType *> av_set_transfer;
  av_set_transfer.push_back(&sum_velocity_dt_with_ghost_cells);
  av_set_transfer.push_back(&sum_rns_dt_with_ghost_cells);
  av_set_transfer.push_back(&sum_rss_dt_with_ghost_cells);


  std::string   filename = prefix + ".averagevelocities";
  std::ofstream output(filename.c_str());
  output << "Average velocities" << std::endl;
  output << "dt_0 " << dt_0 << std::endl;
  output << "Average_calculation_boolean " << average_calculation << std::endl;
  output << "Real_initial_time " << real_initial_time << std::endl;

  return av_set_transfer;
}

template <int dim, typename VectorType, typename DofsType>
std::vector<VectorType *>
AverageVelocities<dim, VectorType, DofsType>::read(std::string prefix)
{
  std::vector<VectorType *> sum_vectors;
  sum_vectors.push_back(&sum_velocity_dt);
  sum_vectors.push_back(&sum_reynolds_normal_stress_dt);
  sum_vectors.push_back(&sum_reynolds_shear_stress_dt);


  std::string   filename = prefix + ".averagevelocities";
  std::ifstream input(filename.c_str());
  AssertThrow(input, ExcFileNotOpen(filename));

  std::string buffer;
  std::getline(input, buffer);
  input >> buffer >> dt_0;
  input >> buffer >> average_calculation;
  input >> buffer >> real_initial_time;

  return sum_vectors;
}


template class AverageVelocities<2, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>;

template class AverageVelocities<2,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 TrilinosWrappers::MPI::BlockVector,
                                 std::vector<IndexSet>>;

template class AverageVelocities<2,
                                 LinearAlgebra::distributed::Vector<double>,
                                 IndexSet>;

template class AverageVelocities<3,
                                 LinearAlgebra::distributed::Vector<double>,
                                 IndexSet>;
