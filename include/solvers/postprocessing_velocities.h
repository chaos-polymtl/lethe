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


template <int dim, typename VectorType>
class PostprocessingVelocities
{
public:
  void
  calculate_velocity_fluctuations(
    const DoFHandler<dim> &              dof_handler,
    const VectorType &                   evaluation_point,
    const Parameters::SimulationControl &simulation_control,
    const Parameters::FEM &              fem_parameters,
    const Parameters::PostProcessing &   post_processing,
    const double &                       current_time,
    const double &                       bulk_velocity,
    const MPI_Comm &                     mpi_communicator);
  std::vector<std::vector<double>>
  average_velocities();
  std::vector<std::vector<double>>
  velocity_fluctuations();
  void
  calculate_velocity_profiles(const double &bulk_velocity);


private:
  std::vector<Tensor<1, dim>> u_dt;
  // std::vector<Tensor<1, dim>> velocity_fluctuations;
  // std::vector<Tensor<1, dim>> average_velocities;
  std::vector<double>         velocity_fluctuations_u;
  std::vector<double>         velocity_fluctuations_v;
  std::vector<double>         velocity_fluctuations_w;
  std::vector<double>         average_velocities_u;
  std::vector<double>         average_velocities_v;
  std::vector<double>         average_velocities_w;
  std::vector<Tensor<1, dim>> velocity_profiles;
  std::vector<Tensor<1, dim>> reynolds_stress;
  std::vector<Tensor<1, dim>> shear_stress;
};

template <int dim, typename VectorType>
void
PostprocessingVelocities<dim, VectorType>::calculate_velocity_fluctuations(
  const DoFHandler<dim> &              dof_handler,
  const VectorType &                   evaluation_point,
  const Parameters::SimulationControl &simulation_control,
  const Parameters::FEM &              fem_parameters,
  const Parameters::PostProcessing &   post_processing,
  const double &                       current_time,
  const double &                       bulk_velocity,
  const MPI_Comm &                     mpi_communicator)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  const MappingQ<dim>       mapping(fe.degree, fem_parameters.qmapping_all);
  QGauss<dim>               quadrature_formula(fe.degree + 1);
  const unsigned int        n_q_points = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      present_velocity_values(n_q_points);
  double total_time = current_time - post_processing.initial_time;

  FEValues<dim> fe_values(mapping,
                          fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  int point = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          // std::cout << "size : " << dof_handler.end() << std::endl;
          for (unsigned int q = 0; q < n_q_points; q++)
            {
              if (abs(current_time - 0.0) < 1e-6)
                {
                  u_dt.push_back(present_velocity_values[q] *
                                 0); // Tensor of zeros
                  average_velocities_u.push_back(0.);
                  average_velocities_v.push_back(0.);
                  velocity_fluctuations_u.push_back(0.);
                  velocity_fluctuations_v.push_back(0.);
                  if (dim == 3)
                    {
                      average_velocities_w.push_back(0.);
                      velocity_fluctuations_w.push_back(0.);
                    }
                }
              else if (total_time > 1e-6)
                {
                  u_dt[point] +=
                    present_velocity_values[q] * simulation_control.dt;
                  average_velocities_u[point] = u_dt[point][0] / total_time;
                  average_velocities_v[point] = u_dt[point][1] / total_time;
                  velocity_fluctuations_u[point] =
                    present_velocity_values[q][0] - average_velocities_u[point];
                  velocity_fluctuations_v[point] =
                    present_velocity_values[q][1] - average_velocities_v[point];
                  if (dim == 3)
                    {
                      average_velocities_w[point] = u_dt[point][2] / total_time;
                      velocity_fluctuations_w[point] =
                        present_velocity_values[q][2] -
                        average_velocities_w[point];
                    }
                }

              /*
              u_dt                    = Utilities::MPI::sum(u_dt,
              mpi_communicator); average_velocities_u    =
              Utilities::MPI::sum(average_velocities_u, mpi_communicator);
              average_velocities_v    =
              Utilities::MPI::sum(average_velocities_v, mpi_communicator);
              average_velocities_w    =
              Utilities::MPI::sum(average_velocities_w, mpi_communicator);
              velocity_fluctuations_u =
              Utilities::MPI::sum(velocity_fluctuations_u, mpi_communicator);
              velocity_fluctuations_v =
              Utilities::MPI::sum(velocity_fluctuations_v, mpi_communicator);
              velocity_fluctuations_w =
              Utilities::MPI::sum(velocity_fluctuations_w, mpi_communicator);
              */

              // std::cout << "u_dt : " << u_dt[point] << std::endl;
              std::cout << "average_velocities : "
                        << average_velocities_u[point] << ","
                        << average_velocities_v[point] << ","
                        << average_velocities_w[point] << std::endl;
              std::cout << "velocity_fluctuations : "
                        << velocity_fluctuations_u[point] << ","
                        << velocity_fluctuations_v[point] << ","
                        << velocity_fluctuations_w[point] << std::endl;
              point++;
            }
        }
    }
}

template <int dim, typename VectorType>
std::vector<std::vector<double>>
PostprocessingVelocities<dim, VectorType>::average_velocities()
{
  std::vector<std::vector<double>> av;
  if (dim == 2)
    {
      av[0] = average_velocities_u;
      av[1] = average_velocities_v;
    }
  else if (dim == 3)
    {
      av[2] = average_velocities_w;
    }
  return av;
}

template <int dim, typename VectorType>
std::vector<std::vector<double>>
PostprocessingVelocities<dim, VectorType>::velocity_fluctuations()
{
  std::vector<std::vector<double>> vf;
  if (dim == 2)
    {
      vf[0] = velocity_fluctuations_u;
      vf[1] = average_velocities_v;
    }
  else if (dim == 3)
    {
      vf[2] = velocity_fluctuations_w;
    }
  return vf;
}


template <int dim, typename VectorType>
void
PostprocessingVelocities<dim, VectorType>::calculate_velocity_profiles(
  const double &bulk_velocity)
{}



#endif
