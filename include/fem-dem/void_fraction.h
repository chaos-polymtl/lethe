/* ---------------------------------------------------------------------
*
* Copyright (C) 2019 - 2024 by the Lethe authors
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
*/


#ifndef lethe_void_fraction_h
#  define lethe_void_fraction_h

#include <core/vector.h>

#include <fem-dem/parameters_cfd_dem.h>

#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

using namespace dealii;

/**
 * @brief Base for all of the void fraction calculators. The void fraction is used as an auxiliary fields for the solvers that solve the Volume-Averaged Navier-Stokes equations.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class VoidFractionBase
{
  VoidFractionBase(  parallel::DistributedTriangulationBase<dim> * triangulation,
                   const Parameters::VoidFractionParameters<dim> & input_parameters, const unsigned int fe_degree):
    dof_handler(*triangulation),
    parameters(input_parameters),
    fe(fe_degree)
  {;}

  /**
  * @brief Setup the degrees of freedom and establishes the constraints of the void fraction systems.
   *
   */
  void setup_dofs();

  /**
  * @brief Assembles the diagonal of the mass matrix to impose constraints on the system
  *
  *  @param diagonal_mass_matrix The matrix for which the diagonal entries will be filled with the diagonal of the mass matrix
   *
   *  @todo Establish if it is really necessary to keep this as a matrix and not as a vector since a diagonal is nothing more than a vector.
  */
  void assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &diagonal_mass_matrix);

  /**
   * @brief Calculates the void fraction
   *
   * @param time current time for which the void fraction is to be calculated.
   *
   */
  void
  calculate_void_fraction(const double time);



protected:

  /// Parameters for the calculation of the void fraction
  const Parameters::VoidFractionParameters<dim> parameters;

  /// DoFHandler that manages the void fraction
  DoFHandler<dim> dof_handler;

  /// FE for the void fraction. Currently only FE_Q (Lagrange polynomials) are supported.
  FE_Q<dim>       fe;

  /// Index set for the locally owned degree of freedoms
  IndexSet locally_owned_dofs;

  /// Index set for the locally relevant degree of freedoms
  IndexSet locally_relevant_dofs;

  // Solution of the void fraction at previous time steps
  std::vector<GlobalVectorType> previous_void_fraction;

  /// Fully distributed (including locally relevant) solution
  GlobalVectorType void_fraction_locally_relevant;

  /// Locally owned solution of the void fraction
  GlobalVectorType void_fraction_locally_owned;

  /// ??? Not fully sure of this yet
  TrilinosWrappers::SparseMatrix complete_system_matrix_void_fraction;
  GlobalVectorType               complete_system_rhs_void_fraction;

  /// Mass matrix used to constraint the value of the void fraction to be bounded
  TrilinosWrappers::SparseMatrix mass_matrix;

  /// Mass matrix diagonal used to constraint the value of the void fraction to be bounded
  GlobalVectorType               diagonal_of_mass_matrix;

  /// ??? BB Not fully sure of this yet
  IndexSet                       active_set;

  /// Preconditioner used for the solution of the smoothed L2 projection of the void fraction
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  /// Constraints used for the boundary conditions of the void fraction. Currently, this is only used to establish periodic void fractions.
  AffineConstraints<double>                          void_fraction_constraints;

  /// System matrix used to assembled the smoothed L2 projection of the void fraction
  TrilinosWrappers::SparseMatrix system_matrix_void_fraction;

  /// Right-hand side used to assemble the smoothed L2 projection of the void fraction
  GlobalVectorType               system_rhs_void_fraction;
};

/**
 * @brief Calculate the void fraction using an analytical function (time and space). Although the void fraction is known analytically, it is still calculated on the grid using a DofHandler to ensure consistency with the underlying finite element scheme.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class VoidFractionFunction : public VoidFractionBase<dim>
{

};

#endif