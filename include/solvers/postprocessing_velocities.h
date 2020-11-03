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

// Dealii Includes

// Base
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

// Lac
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

// Lac - Trilinos includes
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

// Grid
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Numerics
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

// Distributed
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>


// Lethe Includes
#include <core/bdf.h>
#include <core/boundary_conditions.h>
#include <core/manifolds.h>
#include <core/newton_non_linear_solver.h>
#include <core/parameters.h>
#include <core/physics_solver.h>
#include <core/pvd_handler.h>
#include <core/simulation_control.h>
#include <solvers/flow_control.h>

#include "navier_stokes_solver_parameters.h"
#include "post_processors.h"

// Std
#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim, typename VectorType, typename DofsType>
class AverageVelocities
{
public:
  VectorType
  calculate_average_velocities(
    const VectorType &                       local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const Parameters::PostProcessing &       post_processing,
    const DofsType &                         locally_owned_dofs,
    const MPI_Comm &                         mpi_communicator);

  VectorType
  nondimensionalize_average_velocities(const double bulk_velocity);

  VectorType
  calculate_reynolds_stresses(const VectorType & local_evaluation_point,
                              const DofsType &   locally_owned_dofs,
                              const MPI_Comm &   mpi_communicator);

  VectorType
  nondimensionalize_reynolds_stresses(const double bulk_velocity);

private:
  double time;
  double total_time;
  VectorType sum_velocity_dt;
  VectorType average_velocities;
  VectorType nondimensionalized_average_velocities;

  Vector<double> sum_normal_stresses_dt;
  Vector<double> sum_shear_stress_dt;
  VectorType normal_stresses;
  VectorType shear_stress;

  TrilinosScalar inv_range_time;
  TrilinosScalar dt;
};

template <int dim, typename VectorType, typename DofsType>
VectorType
AverageVelocities<dim, VectorType, DofsType>::calculate_average_velocities(
  const VectorType &                       local_evaluation_point,
  const std::shared_ptr<SimulationControl> &simulation_control,
  const Parameters::PostProcessing &       post_processing,
  const DofsType &                         locally_owned_dofs,
  const MPI_Comm &                         mpi_communicator)
{
  time = simulation_control->get_current_time();
  total_time = time - post_processing.initial_time;

  if (time <= 0)
  {
    // Reinitializing vectors with zeros at t = 0
    sum_velocity_dt.reinit(locally_owned_dofs,
                           mpi_communicator);
    average_velocities.reinit(locally_owned_dofs,
                              mpi_communicator);
  }
  else if (total_time >= 0)
  {
    // Generating average velocities at each time step at initial time
    inv_range_time = 1. / (total_time + simulation_control->calculate_time_step()) ;
    dt = simulation_control->calculate_time_step();

    VectorType velocity_dt(locally_owned_dofs,
                           mpi_communicator);
    velocity_dt.equ(dt, local_evaluation_point);

    sum_velocity_dt += velocity_dt;
    if (simulation_control->is_output_iteration())
      average_velocities.equ(inv_range_time, sum_velocity_dt);
  }
  return average_velocities;
}

// Function not tested yet
template <int dim, typename VectorType, typename DofsType>
VectorType
AverageVelocities<dim, VectorType, DofsType>::
nondimensionalize_average_velocities(const double bulk_velocity)
{
  const TrilinosScalar trilinos_inv_bulk_velocity = 1. / bulk_velocity;
  nondimensionalized_average_velocities = average_velocities;
  nondimensionalized_average_velocities.equ(trilinos_inv_bulk_velocity,
                                            average_velocities);
  return nondimensionalized_average_velocities;
}

//Reynolds stresses
template <int dim, typename VectorType, typename DofsType>
VectorType
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stresses(
  const VectorType & local_evaluation_point,
  const DofsType &   locally_owned_dofs,
  const MPI_Comm &   mpi_communicator)
{
  /*
  if (time - 0.0 < 1e-6)
  {
    // Reinitializing vectors with zeros at t = 0
    sum_normal_stresses_dt.reinit(3 * local_evaluation_point.size() / 4);
    normal_stresses = sum_normal_stresses_dt;
    sum_shear_stress_dt.reinit(local_evaluation_point.size() / 4);
    shear_stress = sum_shear_stress_dt;
  }
  else if (abs(total_time) < 1e-6 || total_time > 1e-6)
    {
      Vector<double> normal_stresses_dt(sum_normal_stresses_dt.size());
      Vector<double> shear_stress_dt(sum_shear_stress_dt.size());
      for (const auto &dof : locally_owned_dofs)
        {
          for (unsigned int q = 0; q < local_evaluation_point.size(); q += 4)
            {
              unsigned int i = 0;
              unsigned int j = 0;
                {
                  //normal_stresses_dt[i] = local_evaluation_point[q].add(-average_velocities[q]);
                  //normal_stresses_dt[i+1] = local_evaluation_point[q+1].add(-average_velocities[q+1]);
                  //normal_stresses_dt[i+2] = local_evaluation_point[q+2].add(-average_velocities[q+2]);

                  // Calculate u'*dt, v'*dt and w'*dt
                  normal_stresses_dt[i] = (local_evaluation_point[q] - average_velocities[q]) * trilinos_dt;
                  normal_stresses_dt[i+1] = (local_evaluation_point[q+1] - average_velocities[q+1]) * trilinos_dt;
                  normal_stresses_dt[i+2] = (local_evaluation_point[q+2] - average_velocities[q+2]) * trilinos_dt;
                  i += 3;

                  // Calulate u'v'*dt (u'*dt * u'*dt) / dt with normal_stresses_dt vectors)
                  shear_stress_dt[j] = normal_stresses_dt[i] * normal_stresses_dt[i+1] / trilinos_dt;
                  j++;
                }
            }
        }
      sum_normal_stresses_dt += normal_stresses_dt;
      normal_stresses.equ(trilinos_inv_range_time, sum_normal_stresses_dt);

      sum_shear_stress_dt += shear_stress_dt;
      shear_stress.equ(trilinos_inv_range_time, sum_shear_stress_dt);
    }
  return normal_stresses; */
}

template <int dim, typename VectorType, typename DofsType>
VectorType
AverageVelocities<dim, VectorType, DofsType>::nondimensionalize_reynolds_stresses(
  const double bulk_velocity)
{

}


#endif
