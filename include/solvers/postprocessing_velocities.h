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

  Vector<double>
  calculate_reynolds_stress(
    const VectorType &                        local_evaluation_point,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const DofsType &                          locally_owned_dofs,
    const MPI_Comm &                          mpi_communicator);

  VectorType
  nondimensionalize_reynolds_stress(const double bulk_velocity);

private:
  double time;
  double total_time;
  VectorType sum_velocity_dt;
  VectorType average_velocities;
  VectorType nondimensionalized_average_velocities;

  Vector<TrilinosScalar> sum_reynolds_stress_dt;
  Vector<TrilinosScalar> reynolds_stress;

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

  if (simulation_control->get_step_number() == 0)
  {
    // Reinitializing vectors with zeros and dt at t = 0
    sum_velocity_dt.reinit(locally_owned_dofs,
                           mpi_communicator);
    average_velocities.reinit(locally_owned_dofs,
                              mpi_communicator);
    dt = simulation_control->calculate_time_step();
  }
  else if (abs(total_time) < 1e-6 || total_time > 0)
  {
    // Generating average velocities at each time from initial time
    inv_range_time = 1. / (total_time + dt);

    std::cout << total_time + dt << std::endl;
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
  const TrilinosScalar inv_bulk_velocity = 1. / bulk_velocity;
  nondimensionalized_average_velocities = average_velocities;
  nondimensionalized_average_velocities.equ(inv_bulk_velocity, average_velocities);
  return nondimensionalized_average_velocities;
}

//Reynolds stress
template <int dim, typename VectorType, typename DofsType>
Vector<TrilinosScalar>
AverageVelocities<dim, VectorType, DofsType>::calculate_reynolds_stress(
  const VectorType &                        present_solution,
  const std::shared_ptr<SimulationControl> &simulation_control,
  const DofsType &                          locally_owned_dofs,
  const MPI_Comm &                          mpi_communicator)
{
  Vector<TrilinosScalar> solution = present_solution;

  if (simulation_control->get_step_number() == 0)
  {
    // Reinitializing vectors with zeros at t = 0
    sum_reynolds_stress_dt.reinit(solution.size());
    reynolds_stress.reinit(solution.size());
  }
  else if (abs(total_time) < 1e-6 || total_time > 1e-6)
    {
      Vector<TrilinosScalar> reynolds_stress_dt(reynolds_stress.size());

      // Won't work with blockVector
      for (unsigned int i = 0; i < solution.size(); i += 4)
        {
          // normal_stress_dt[i] = local_evaluation_point[q].add(-average_velocities[q]);
          // normal_stress_dt[i+1] = local_evaluation_point[q+1].add(-average_velocities[q+1]);
          // normal_stress_dt[i+2] = local_evaluation_point[q+2].add(-average_velocities[q+2]);

          // Calculate u'*dt, v'*dt and w'*dt
          reynolds_stress_dt[i] = (solution[i] - average_velocities[i]) *
                                  (solution[i] - average_velocities[i]) * dt;
          reynolds_stress_dt[i + 1] =
            (solution[i + 1] - average_velocities[i + 1]) *
            (solution[i + 1] - average_velocities[i + 1]) * dt;
          reynolds_stress_dt[i + 2] =
            (solution[i + 2] - average_velocities[i + 2]) *
            (solution[i + 2] - average_velocities[i + 2]) * dt;
          reynolds_stress_dt[i + 3] =
            (solution[i] - average_velocities[i]) *
            (solution[i + 1] - average_velocities[i + 1]) * dt;

          sum_reynolds_stress_dt += reynolds_stress_dt;

          if (simulation_control->is_output_iteration())
            reynolds_stress.equ(inv_range_time, sum_reynolds_stress_dt);
        }
    }

  return reynolds_stress;
}

template <int dim, typename VectorType, typename DofsType>
VectorType
AverageVelocities<dim, VectorType, DofsType>::nondimensionalize_reynolds_stress(
  const double bulk_velocity)
{

}


#endif
