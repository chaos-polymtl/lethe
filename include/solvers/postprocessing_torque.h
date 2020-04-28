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
 * Author: Bruno Blais, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_postprocessing_torque_h
#define lethe_postprocessing_torque_h

// Base
#include <deal.II/base/quadrature_lib.h>

// Lac
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/vector.h>

// Dofs
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

// Fe
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

// Lethe includes
#include <core/boundary_conditions.h>
#include <core/parameters.h>


using namespace dealii;
/**
 * @brief Calculates the torques due to the fluid motion on every boundary conditions
 * @return std::vector of torques on each boundary condition
 * Post-processing function
 * This function calculates the torqueacting on each of the boundary conditions
 * within the domain. It generates a vector which size is the number of boundary
 * conditions.
 *
 * @param dof_handler The dof_handler used for the calculation.
 *
 * @param evaluation_point The solution at which the torque is calculated.
 *
 * @param physical_properties The parameters containing the required physical properties.
 *
 * @param fem_parameters The fem_parameters of the simulation.
 *
 * @param boundary_conditions The boundary conditions object.
 *
 * @param mpi_communicator The mpi communicator. It is used to reduce the torque calculation.
 */
template <int dim, typename VectorType>
std::vector<Tensor<1, 3>>
calculate_torques(
  const DoFHandler<dim> &                              dof_handler,
  const VectorType &                                   evaluation_point,
  const Parameters::PhysicalProperties &               physical_properties,
  const Parameters::FEM &                              fem_parameters,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const MPI_Comm &                                     mpi_communicator);

#endif
