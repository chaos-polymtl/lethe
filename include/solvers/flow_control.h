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

// Dofs
#include <deal.II/dofs/dof_handler.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

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

// Lethe includes
#include <core/parameters.h>

#ifndef lethe_flow_control_h
#  define lethe_flow_control_h

using namespace dealii;

template <int dim, typename VectorType>
class FlowControl
{
public:
  Tensor<1, dim>
  calculate_beta(const DoFHandler<dim> &               dof_handler,
                 const VectorType &                    present_solution,
                 const Parameters::DynamicFlowControl &flow_control,
                 const Parameters::SimulationControl & simulation_control,
                 const Parameters::FEM &               fem_parameters,
                 const double &                        step_number,
                 const MPI_Comm &                      mpi_communicator);
  std::vector<double>
  flow_summary();
  double
  bulk_velocity();

private:
  // The coefficients are stored in the following fashion :
  // 0 - Flow control intended, n - n, 1n - n-1, n1 - n+1
  Tensor<1, dim> beta;
  double         beta_n;
  double         beta_n1     = 0;
  double         flow_rate_n = 0;
  double         flow_rate_1n;
  double         area;
};

template <int dim, typename VectorType>
Tensor<1, dim>
FlowControl<dim, VectorType>::calculate_beta(
  const DoFHandler<dim> &               dof_handler,
  const VectorType &                    present_solution,
  const Parameters::DynamicFlowControl &flow_control,
  const Parameters::SimulationControl & simulation_control,
  const Parameters::FEM &               fem_parameters,
  const double &                        step_number,
  const MPI_Comm &                      mpi_communicator)
{
  beta_n       = beta_n1;
  flow_rate_1n = flow_rate_n;

  const FiniteElement<dim> &fe = dof_handler.get_fe();

  const MappingQ<dim> mapping(fe.degree, fem_parameters.qmapping_all);
  QGauss<dim - 1>     face_quadrature_formula(fe.degree + 1);
  const unsigned int  n_q_points = face_quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  std::vector<Tensor<1, dim>>      velocity_values(n_q_points);
  Tensor<1, dim>                   normal_vector;

  FEFaceValues<dim> fe_face_values(mapping,
                                   fe,
                                   face_quadrature_formula,
                                   update_values | update_quadrature_points |
                                     update_JxW_values | update_normal_vectors);

  unsigned int boundary_id = flow_control.id_flow_control;

  // Resetting next flow rate and area at inlet prior new calculation
  flow_rate_n = 0;
  area        = 0;

  // Calculating area and volumetric flow rate at the inlet flow
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->at_boundary())
        {
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               face++)
            {
              if (cell->face(face)->at_boundary())
                {
                  fe_face_values.reinit(cell, face);
                  if (cell->face(face)->boundary_id() == boundary_id)
                    {
                      for (unsigned int q = 0; q < n_q_points; q++)
                        {
                          area += fe_face_values.JxW(q);
                          normal_vector = fe_face_values.normal_vector(q);
                          fe_face_values[velocities].get_function_values(
                            present_solution, velocity_values);
                          flow_rate_n += velocity_values[q] * normal_vector *
                                         fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }
    }

  area        = Utilities::MPI::sum(area, mpi_communicator);
  flow_rate_n = Utilities::MPI::sum(flow_rate_n, mpi_communicator);

  // Calculating the next beta coefficient
  const double dt          = simulation_control.dt;
  const double flow_rate_0 = flow_control.flow_rate;

  if (step_number == 1)
    beta_n1 = flow_control.beta_0;
  else if (step_number == 2)
    beta_n1 = beta_n - (flow_rate_0 - flow_rate_n) / (area * dt);
  else
    beta_n1 =
      beta_n - (flow_rate_0 - 2 * flow_rate_n + flow_rate_1n) / (area * dt);

  // Setting beta coefficient only to new time step
  if (flow_control.flow_direction == 0)
    beta[0] = beta_n1; // beta = f_x
  else if (flow_control.flow_direction == 1)
    beta[1] = beta_n1; // beta = f_y
  else if (flow_control.flow_direction == 2)
    beta[2] = beta_n1; // beta = f_z

  return beta;
}

template <int dim, typename VectorType>
std::vector<double>
FlowControl<dim, VectorType>::flow_summary()
{
  std::vector<double> summary{area, flow_rate_n, beta_n};
  return summary;
}

template <int dim, typename VectorType>
double
FlowControl<dim, VectorType>::bulk_velocity()
{
  double u_b = flow_rate_n / area;
  return u_b;
}

#endif
