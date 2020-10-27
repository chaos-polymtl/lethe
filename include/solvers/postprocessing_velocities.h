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
//#include <solvers/flow_control.h>

#include "navier_stokes_solver_parameters.h"
#include "post_processors.h"

// Std
#include <fstream>
#include <iostream>

using namespace dealii;


template <int dim, typename VectorType>
class PostprocessingVelocities
{
public:
  std::vector<Tensor<1, 3>>
  calculate_velocity_fluctuations(
    const DoFHandler<dim> &              dof_handler,
    const VectorType &                   present_solution,
    const Parameters::SimulationControl &simulation_control,
    const Parameters::FEM &              fem_parameters,
    const double &                       step_number,
    const MPI_Comm &                     mpi_communicator);

private:
  // double bulk_velocity = FlowControl<dim, VectorType>::bulk_velocity();
  std::vector<Tensor<1, 3>> velocity_fluctuations;
  Tensor<1, dim>            average_velocities;
};

template <int dim, typename VectorType>
std::vector<Tensor<1, 3>>
PostprocessingVelocities<dim, VectorType>::calculate_velocity_fluctuations(
  const DoFHandler<dim> &              dof_handler,
  const VectorType &                   evaluation_point,
  const Parameters::SimulationControl &simulation_control,
  const Parameters::FEM &              fem_parameters,
  const double &                       step_number,
  const MPI_Comm &                     mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const MappingQ<dim>       mapping(fe.degree, fem_parameters.qmapping_all);
  QGauss<dim>               quadrature_formula(fe.degree + 1);
  const unsigned int        n_q_points = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      present_velocity_values(n_q_points);

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);


  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          for (unsigned int q = 0; q < n_q_points; q++)
            {
              fe_values[velocities].get_function_gradients(
                evaluation_point, present_velocity_values);
              std::cout << "velocity values : " << present_velocity_values
                        << std::endl;
            }
        }
    }
  return present_velocity_values;
}



#endif
