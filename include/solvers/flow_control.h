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


/**
 * @brief FlowControl. The FlowControl class calculates the volumetric flow
 * rate at the inlet or at a selected boundary. It also allows to dynamically
 * control the flow with a beta coefficient calculated at each step time.
 */

template <int dim, typename VectorType>
class FlowControl
{
public:
  /**
   * @brief calculate_flow_rate. This function calculates the volumetric flow
   * rate at the selected boundary. It actually calculates the flow rate through
   * the summation of the value at each cell surface with the normal vector,
   * the velocity value and the area.
   *
   * @param dof_handler. The argument used to get finite elements.
   *
   * @param present_solution. The vector which contains all the values to
   *                          calculate the flow rate.
   *
   * @param flow_control. The flow control parameters
   *
   * @param fem_parameters. The FEM parameters
   *
   * @param mpi_communicator. The mpi communicator information
   */
  void
  calculate_flow_rate(const DoFHandler<dim> &               dof_handler,
                      const VectorType &                    present_solution,
                      const Parameters::DynamicFlowControl &flow_control,
                      const Parameters::FEM &               fem_parameters,
                      const MPI_Comm &                      mpi_communicator);

  /**
   * @brief calculate_beta. This function calculates a beta coefficient which
   * applies a force to the flow in order to adjust the flow rate to a desired
   * value through the previous flow rate. Once the flow rate is reached, the
   * algorithm calculates a new beta to keep a constant flow rate.
   *
   * @param flow_control. The flow control parameters
   *
   * @param simulation_control. The simulation control parameters
   *
   * @param step_number. The current step
   */
  void
  calculate_beta(const Parameters::DynamicFlowControl &flow_control,
                 const Parameters::SimulationControl & simulation_control,
                 const unsigned int &                  step_number);

  /**
   * @brief get_flow_rate. This function gives the flow rate of the
   * step time. Arguments are the same as calculate_flow_rate(...)
   */
  double
  get_flow_rate(const DoFHandler<dim> &               dof_handler,
                const VectorType &                    present_solution,
                const Parameters::DynamicFlowControl &flow_control,
                const Parameters::FEM &               fem_parameters,
                const MPI_Comm &                      mpi_communicator);

  /**
   * @brief get_beta. This function gives the beta coefficient of the
   * step time. Arguments are those of calculate_flow_rate (...) and
   * calculate_beta(...)
   */
  Tensor<1, dim>
  get_beta(const DoFHandler<dim> &               dof_handler,
           const VectorType &                    present_solution,
           const Parameters::DynamicFlowControl &flow_control,
           const Parameters::SimulationControl & simulation_control,
           const Parameters::FEM &               fem_parameters,
           const unsigned int &                  step_number,
           const MPI_Comm &                      mpi_communicator);

  /**
   * @brief flow_summary. This function gives information of the flow control
   */
  std::vector<double>
  flow_summary();

  /**
   * @brief bulk_velocity. This function calculates the bulk velocity which is the
   * velocity at the boundary.
   */
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

  // Variables used to improve convergence
  bool   adjusted         = false;
  double threshold_factor = 1.01; // 1%
};


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
