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

#ifndef lethe_postprocessing_kinetic_energy_h
#define lethe_postprocessing_kinetic_energy_h

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
#include <core/parameters.h>
#include <solvers/boundary_conditions.h>


using namespace dealii;
/**
 * @brief Calculates the force due to the fluid motion on every boundary conditions
 * @return std::vector of forces on each boundary condition
 * Post-processing function
 * This function calculates the force acting on each of the boundary conditions
 * within the domain. It generates a vector which size is the number of boundary
 * conditions
 *
 * @param fe The finite element system object used for the calculation
 *
 * @param dof_handler The dof_handler used for the calculation
 *
 * @param evaluation_point The solution at which the force is calculated
 *
 * @param fem_parameters The fem_parameters of the simulation
 *
 * @param mpi_communicator The mpi communicator. It is used to reduce the force calculation
 */
template <int dim, typename VectorType>
double
calculate_kinetic_energy(const FESystem<dim> &  fe,
                         const DoFHandler<dim> &dof_handler,
                         const VectorType &     evaluation_point,
                         const Parameters::FEM &fem_parameters,
                         const MPI_Comm &       mpi_communicator);

#endif
