// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>
#include <core/vector.h>

#include <solvers/postprocessing_velocities.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/fe.h>

#include <deal.II/lac/vector_operation.h>

#include <array>
#include <fstream>
#include <type_traits>
#include <vector>

template <int dim, typename VectorType, typename DofsType>
AverageVelocities<dim, VectorType, DofsType>::AverageVelocities(
  DoFHandler<dim> &dof_handler)
  : solution_transfer_sum_velocity_dt(dof_handler, true)
  , solution_transfer_sum_reynolds_normal_stress_dt(dof_handler, true)
  , solution_transfer_sum_reynolds_shear_stress_dt(dof_handler, true)
  , total_time_for_average(0.0)
  , has_started_averaging(false)
  , dof_handler(dof_handler)
{
  // Initialize average velocity solution vector shared pointer
  get_av = std::make_shared<VectorType>();
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType                      &local_evaluation_point,
  const Parameters::PostProcessing<dim> &post_processing,
  const double                          &current_time,
  const double                          &time_step)
{
  const double epsilon = 1e-6;
  const double initial_time =
    post_processing.initial_time_for_average_velocities;
  dt = time_step;

  // When averaging velocities begins
  if (current_time >= (initial_time - epsilon))
    {
      if (!has_started_averaging)
        {
          has_started_averaging = true;
          real_initial_time     = current_time;

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

// The Reynolds stresses are computed at every locally owned velocity support
// point from the degrees of freedom of its velocity components, which are
// gathered by gather_locally_owned_support_points(). The global indices of the
// degrees of freedom are also the global indices of the vectors, including the
// block vectors, whose blocks hold the velocity and the pressure degrees of
// freedom numbered component-wise. The stresses are stored in the slots of the
// solution: (<u'u'>, <v'v'>, <w'w'>, k) in the normal stress vector and
// (<u'v'>, <v'w'>, <w'u'>, 0) in the shear stress vector, where k is only
// stored at the pressure degrees of freedom that share a velocity support
// point. In 2D, the v slot of the shear stress vector is not used.
template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stresses(
  const VectorType &local_evaluation_point)
{
  for (const SupportPointDoFs &support_point : locally_owned_support_points)
    {
      const std::array<types::global_dof_index, dim> &velocity_dofs =
        support_point.velocity;

      // Velocity fluctuations u', v' and w'
      std::array<double, dim> fluctuation;
      for (int d = 0; d < dim; ++d)
        fluctuation[d] = local_evaluation_point(velocity_dofs[d]) -
                         average_velocities(velocity_dofs[d]);

      // u'u'*dt, v'v'*dt and w'w'*dt
      double sum_normal_stresses = 0.;
      for (int d = 0; d < dim; ++d)
        {
          const double normal_stress_dt = fluctuation[d] * fluctuation[d] * dt;
          reynolds_normal_stress_dt(velocity_dofs[d]) = normal_stress_dt;
          sum_normal_stresses += normal_stress_dt;
        }

      // u'v'*dt
      reynolds_shear_stress_dt(velocity_dofs[0]) =
        fluctuation[0] * fluctuation[1] * dt;

      if constexpr (dim == 3)
        {
          // v'w'*dt
          reynolds_shear_stress_dt(velocity_dofs[1]) =
            fluctuation[1] * fluctuation[2] * dt;

          // w'u'*dt
          reynolds_shear_stress_dt(velocity_dofs[2]) =
            fluctuation[2] * fluctuation[0] * dt;
        }

      // k*dt = 1/2(u'u'+v'v'+w'w')*dt (turbulent kinetic energy)
      if (support_point.pressure != numbers::invalid_dof_index)
        reynolds_normal_stress_dt(support_point.pressure) =
          0.5 * sum_normal_stresses;
    }

  reynolds_normal_stress_dt.compress(VectorOperation::insert);
  reynolds_shear_stress_dt.compress(VectorOperation::insert);

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
  const DofsType &locally_owned_dofs,
  const DofsType &locally_relevant_dofs,
  const MPI_Comm &mpi_communicator)
{
  // Reinitialisation of the average velocity and reynolds stress vectors
  // to get the right length.
  velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_velocity_dt.reinit(locally_owned_dofs, mpi_communicator);
  average_velocities.reinit(locally_owned_dofs, mpi_communicator);
  get_av->reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  // Reinitialize independent components of stress tensor vectors.
  reynolds_normal_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_reynolds_normal_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  reynolds_normal_stresses.reinit(locally_owned_dofs, mpi_communicator);
  get_rns.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  reynolds_shear_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  sum_reynolds_shear_stress_dt.reinit(locally_owned_dofs, mpi_communicator);
  reynolds_shear_stresses.reinit(locally_owned_dofs, mpi_communicator);
  get_rss.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

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

  // The Reynolds stresses are written using the global indices of the degrees
  // of freedom, which requires vectors that span all of them.
  AssertDimension(reynolds_normal_stress_dt.size(), dof_handler.n_dofs());

  gather_locally_owned_support_points();
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::
  gather_locally_owned_support_points()
{
  // Tolerance used to identify two unit support points of the reference cell.
  constexpr double support_point_tolerance = 1e-12;

  Assert(!dof_handler.has_hp_capabilities(), ExcNotImplemented());

  const FiniteElement<dim> &fe = dof_handler.get_fe();
  AssertDimension(fe.n_components(), dim + 1);
  Assert(fe.is_primitive(), ExcNotImplemented());

  // All the velocity components share one base element, so the degrees of
  // freedom of a velocity support point are the degrees of freedom of the dim
  // components with the same index within this base element.
  const unsigned int velocity_base = fe.component_to_base_index(0).first;
  for (int d = 1; d < dim; ++d)
    Assert(fe.component_to_base_index(d).first == velocity_base,
           ExcInternalError());

  const FiniteElement<dim> &velocity_fe = fe.base_element(velocity_base);
  const FiniteElement<dim> &pressure_fe =
    fe.base_element(fe.component_to_base_index(dim).first);
  AssertThrow(velocity_fe.has_support_points() &&
                pressure_fe.has_support_points(),
              ExcMessage(
                "The time-averaged Reynolds stresses require velocity and "
                "pressure finite elements with support points."));

  const unsigned int n_velocity_nodes = velocity_fe.n_dofs_per_cell();
  const unsigned int n_pressure_nodes = pressure_fe.n_dofs_per_cell();
  const std::vector<Point<dim>> &velocity_unit_support_points =
    velocity_fe.get_unit_support_points();
  const std::vector<Point<dim>> &pressure_unit_support_points =
    pressure_fe.get_unit_support_points();

  // Cell-local indices of the velocity components of every velocity node
  std::vector<std::array<unsigned int, dim>> local_velocity_dofs(
    n_velocity_nodes);
  for (unsigned int i = 0; i < n_velocity_nodes; ++i)
    for (int d = 0; d < dim; ++d)
      local_velocity_dofs[i][d] = fe.component_to_system_index(d, i);

  // Cell-local index of the pressure at every velocity node. Every pressure
  // node is paired with the first velocity node that has the same unit support
  // point and that is not paired yet. With equal-order elements, this pairs the
  // nodes with the same index. With bubble functions, which all have their
  // support point at the center of the cell, the Lagrange nodes are paired
  // first, since they are numbered before the bubble functions.
  std::vector<unsigned int> local_pressure_dofs(n_velocity_nodes,
                                                numbers::invalid_unsigned_int);
  for (unsigned int j = 0; j < n_pressure_nodes; ++j)
    for (unsigned int i = 0; i < n_velocity_nodes; ++i)
      if (local_pressure_dofs[i] == numbers::invalid_unsigned_int &&
          velocity_unit_support_points[i].distance(
            pressure_unit_support_points[j]) < support_point_tolerance)
        {
          local_pressure_dofs[i] = fe.component_to_system_index(dim, j);
          break;
        }

  // Gather the global indices of every locally owned velocity support point
  // once. A support point is identified by the degree of freedom of its first
  // velocity component. All the degrees of freedom of a support point belong to
  // the same mesh object, and therefore to the same process.
  const IndexSet   &locally_owned_dofs = dof_handler.locally_owned_dofs();
  std::vector<bool> is_gathered(locally_owned_dofs.n_elements(), false);
  std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());

  locally_owned_support_points.clear();
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      cell->get_dof_indices(dof_indices);
      for (unsigned int i = 0; i < n_velocity_nodes; ++i)
        {
          const types::global_dof_index first_velocity_dof =
            dof_indices[local_velocity_dofs[i][0]];
          if (!locally_owned_dofs.is_element(first_velocity_dof))
            continue;

          const types::global_dof_index index_within_owned_dofs =
            locally_owned_dofs.index_within_set(first_velocity_dof);
          if (is_gathered[index_within_owned_dofs])
            continue;
          is_gathered[index_within_owned_dofs] = true;

          SupportPointDoFs support_point;
          for (int d = 0; d < dim; ++d)
            {
              support_point.velocity[d] =
                dof_indices[local_velocity_dofs[i][d]];
              Assert(locally_owned_dofs.is_element(support_point.velocity[d]),
                     ExcInternalError());
            }

          support_point.pressure =
            (local_pressure_dofs[i] == numbers::invalid_unsigned_int) ?
              numbers::invalid_dof_index :
              dof_indices[local_pressure_dofs[i]];
          Assert(support_point.pressure == numbers::invalid_dof_index ||
                   locally_owned_dofs.is_element(support_point.pressure),
                 ExcInternalError());

          locally_owned_support_points.push_back(support_point);
        }
    }
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::prepare_for_mesh_adaptation()
{
  *get_av = sum_velocity_dt;
  get_rns = sum_reynolds_normal_stress_dt;
  get_rss = sum_reynolds_shear_stress_dt;
  solution_transfer_sum_velocity_dt.prepare_for_coarsening_and_refinement(
    *get_av);
  solution_transfer_sum_reynolds_normal_stress_dt
    .prepare_for_coarsening_and_refinement(get_rns);
  solution_transfer_sum_reynolds_shear_stress_dt
    .prepare_for_coarsening_and_refinement(get_rss);
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::post_mesh_adaptation()
{
  if constexpr (std::is_same_v<VectorType,
                               LinearAlgebra::distributed::Vector<double>>)
    {
      // When the vector are deal.II vectors, the solution transfer expects a
      // vector with ghost entries
      solution_transfer_sum_velocity_dt.interpolate(
        sum_velocity_dt_with_ghost_cells);
      solution_transfer_sum_reynolds_normal_stress_dt.interpolate(
        sum_rns_dt_with_ghost_cells);
      solution_transfer_sum_reynolds_shear_stress_dt.interpolate(
        sum_rss_dt_with_ghost_cells);

      sum_velocity_dt               = sum_velocity_dt_with_ghost_cells;
      sum_reynolds_normal_stress_dt = sum_rns_dt_with_ghost_cells;
      sum_reynolds_shear_stress_dt  = sum_rss_dt_with_ghost_cells;
    }
  else
    {
      // When the vectors are trilinos vectors, the solution transfer expects a
      // locally_owned vector
      solution_transfer_sum_velocity_dt.interpolate(sum_velocity_dt);
      solution_transfer_sum_reynolds_normal_stress_dt.interpolate(
        sum_reynolds_normal_stress_dt);
      solution_transfer_sum_reynolds_shear_stress_dt.interpolate(
        sum_reynolds_shear_stress_dt);



      sum_velocity_dt_with_ghost_cells = sum_velocity_dt;
      sum_rns_dt_with_ghost_cells      = sum_reynolds_normal_stress_dt;
      sum_rss_dt_with_ghost_cells      = sum_reynolds_shear_stress_dt;
    }



  update_average_velocities();
}

template <int dim, typename VectorType, typename DofsType>
std::vector<const VectorType *>
AverageVelocities<dim, VectorType, DofsType>::save(const std::string &prefix)
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
  output << "has_started_averaging_boolean " << has_started_averaging
         << std::endl;
  output << "Real_initial_time " << real_initial_time << std::endl;

  return av_set_transfer;
}

template <int dim, typename VectorType, typename DofsType>
std::vector<VectorType *>
AverageVelocities<dim, VectorType, DofsType>::read(const std::string &prefix)
{
  std::vector<VectorType *> sum_vectors;
  sum_vectors.push_back(&sum_velocity_dt_with_ghost_cells);
  sum_vectors.push_back(&sum_rns_dt_with_ghost_cells);
  sum_vectors.push_back(&sum_rss_dt_with_ghost_cells);


  std::string filename = prefix + ".averagevelocities";
  check_file_exists(filename,
                    "checkpoint file of the time-averaged velocities, given by "
                    "'subsection restart' - 'set filename'");
  std::ifstream input(filename.c_str());

  std::string buffer;
  std::getline(input, buffer);
  input >> buffer >> dt_0;
  input >> buffer >> has_started_averaging;
  input >> buffer >> real_initial_time;

  return sum_vectors;
}

template <int dim, typename VectorType, typename DofsType>
void
AverageVelocities<dim, VectorType, DofsType>::zero_average_after_restart()
{
  sum_velocity_dt_with_ghost_cells = 0.0;
  sum_rns_dt_with_ghost_cells      = 0.0;
  sum_rss_dt_with_ghost_cells      = 0.0;

  sum_velocity_dt               = 0.0;
  sum_reynolds_normal_stress_dt = 0.0;
  sum_reynolds_shear_stress_dt  = 0.0;

  has_started_averaging = false;
}

template class AverageVelocities<2, GlobalVectorType, IndexSet>;

template class AverageVelocities<3, GlobalVectorType, IndexSet>;

template class AverageVelocities<2,
                                 GlobalBlockVectorType,
                                 std::vector<IndexSet>>;

template class AverageVelocities<3,
                                 GlobalBlockVectorType,
                                 std::vector<IndexSet>>;
#ifndef LETHE_USE_LDV
template class AverageVelocities<2,
                                 LinearAlgebra::distributed::Vector<double>,
                                 IndexSet>;

template class AverageVelocities<3,
                                 LinearAlgebra::distributed::Vector<double>,
                                 IndexSet>;
#endif
