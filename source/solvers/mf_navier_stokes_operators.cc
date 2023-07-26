/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------*/

#include "solvers/mf_navier_stokes_operators.h"

/**
 * @brief Matrix free helper function
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam Number Abstract type for number across the class (i.e., double).
 * @param function Function to evaluate.
 * @param p_vectorized Batch of points to evaluate function at.
 * @return VectorizedArray<Number> Batch of evaluated values.
 */
template <int dim, typename Number>
VectorizedArray<Number>
evaluate_function(const Function<dim> &                      function,
                  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  VectorizedArray<Number> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];
      result[v] = function.value(p);
    }
  return result;
}

/**
 * @brief Matrix free helper function
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam Number Abstract type for number across the class (i.e., double).
 * @tparam components Number of solution components.
 * @param function Function to evaluate.
 * @param p_vectorized Batch of points to evaluate function at.
 * @return Tensor<1, components, VectorizedArray<Number>> Batch of evaluated values.
 */
template <int dim, typename Number, int components>
Tensor<1, components, VectorizedArray<Number>>
evaluate_function(const Function<dim> &                      function,
                  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  Tensor<1, components, VectorizedArray<Number>> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];
      for (unsigned int d = 0; d < components; ++d)
        result[d][v] = function.value(p, d);
    }
  return result;
}

template <int dim, typename number>
NavierStokesOperatorBase<dim, number>::NavierStokesOperatorBase()
{}

template <int dim, typename number>
NavierStokesOperatorBase<dim, number>::NavierStokesOperatorBase(
  const Mapping<dim> &             mapping,
  const DoFHandler<dim> &          dof_handler,
  const AffineConstraints<number> &constraints,
  const Quadrature<dim> &          quadrature,
  const SimulationParameters<dim> &nsparam,
  const unsigned int               mg_level)
{
  this->reinit(
    mapping, dof_handler, constraints, quadrature, nsparam, mg_level);
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::reinit(
  const Mapping<dim> &             mapping,
  const DoFHandler<dim> &          dof_handler,
  const AffineConstraints<number> &constraints,
  const Quadrature<dim> &          quadrature,
  const SimulationParameters<dim> &nsparam,
  const unsigned int               mg_level)
{
  this->system_matrix.clear();
  this->constraints.copy_from(constraints);

  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    (update_values | update_gradients | update_JxW_values |
     update_quadrature_points | update_hessians);
  additional_data.mg_level = mg_level;

  matrix_free.reinit(
    mapping, dof_handler, constraints, quadrature, additional_data);

  this->parameters = nsparam;

  this->compute_element_size();

  constrained_indices.clear();
  for (auto i : this->matrix_free.get_constrained_dofs())
    constrained_indices.push_back(i);
  constrained_values.resize(constrained_indices.size());

  if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
    {
      std::vector<types::global_dof_index> interface_indices;
      IndexSet                             refinement_edge_indices;
      refinement_edge_indices = get_refinement_edges(this->matrix_free);
      refinement_edge_indices.fill_index_vector(interface_indices);

      edge_constrained_indices.clear();
      edge_constrained_indices.reserve(interface_indices.size());
      edge_constrained_values.resize(interface_indices.size());
      const IndexSet &locally_owned =
        this->matrix_free.get_dof_handler().locally_owned_mg_dofs(
          this->matrix_free.get_mg_level());
      for (unsigned int i = 0; i < interface_indices.size(); ++i)
        if (locally_owned.is_element(interface_indices[i]))
          edge_constrained_indices.push_back(
            locally_owned.index_within_set(interface_indices[i]));

      this->has_edge_constrained_indices =
        Utilities::MPI::max(edge_constrained_indices.size(),
                            dof_handler.get_communicator()) > 0;

      if (this->has_edge_constrained_indices)
        {
          edge_constrained_cell.resize(matrix_free.n_cell_batches(), false);

          VectorType temp;
          matrix_free.initialize_dof_vector(temp);

          for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
            temp.local_element(edge_constrained_indices[i]) = 1.0;

          temp.update_ghost_values();

          FECellIntegrator integrator(matrix_free);

          for (unsigned int cell = 0; cell < matrix_free.n_cell_batches();
               ++cell)
            {
              integrator.reinit(cell);
              integrator.read_dof_values(temp);

              for (unsigned int i = 0; i < integrator.dofs_per_cell; ++i)
                if ((integrator.begin_dof_values()[i] ==
                     VectorizedArray<number>()) == false)
                  {
                    edge_constrained_cell[cell] = true;
                    break;
                  }
            }

          unsigned int count = 0;
          for (const auto i : edge_constrained_cell)
            if (i)
              count++;

          const unsigned int count_global =
            Utilities::MPI::sum(count, dof_handler.get_communicator());

          const unsigned int count_cells_global =
            Utilities::MPI::sum(matrix_free.n_cell_batches(),
                                dof_handler.get_communicator());

          if (Utilities::MPI::this_mpi_process(
                dof_handler.get_communicator()) == 0)
            std::cout << count_global << " " << count_cells_global << std::endl;
        }
    }
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::compute_element_size()
{
  FECellIntegrator phi(matrix_free);

  unsigned int n_cells =
    matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
  element_size.resize(n_cells);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (auto lane = 0u;
           lane < matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          double h_k = matrix_free.get_cell_iterator(cell, lane)->measure();

          if (dim == 2)
            element_size[cell][lane] = std::sqrt(4. * h_k / M_PI);
          else if (dim == 3)
            element_size[cell][lane] = std::pow(6 * h_k / M_PI, 1. / 3.);
        }
    }
}

template <int dim, typename number>
types::global_dof_index
NavierStokesOperatorBase<dim, number>::m() const
{
  if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
    return this->matrix_free.get_dof_handler().n_dofs(
      this->matrix_free.get_mg_level());
  else
    return this->matrix_free.get_dof_handler().n_dofs();
}

template <int dim, typename number>
number
NavierStokesOperatorBase<dim, number>::el(unsigned int, unsigned int) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::clear()
{}

// The matrix free object initialized the vector accordingly
template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::initialize_dof_vector(
  VectorType &vec) const
{
  matrix_free.initialize_dof_vector(vec);
}

// Performs an operator evaluation by looping over all cells
// and evluating the integrals in the matrix-free way
template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::vmult(VectorType &      dst,
                                             const VectorType &src) const
{
  // save values for edge constrained dofs and set them to 0 in src vector
  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    {
      edge_constrained_values[i] = std::pair<number, number>(
        src.local_element(edge_constrained_indices[i]),
        dst.local_element(edge_constrained_indices[i]));

      const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
        .local_element(edge_constrained_indices[i]) = 0.;
    }

  this->matrix_free.cell_loop(
    &NavierStokesOperatorBase::do_cell_integral_range, this, dst, src, true);

  // set constrained dofs as the sum of current dst value and src value
  for (unsigned int i = 0; i < constrained_indices.size(); ++i)
    dst.local_element(constrained_indices[i]) =
      src.local_element(constrained_indices[i]);

  // restoring edge constrained dofs in src and dst
  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    {
      const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
        .local_element(edge_constrained_indices[i]) =
        edge_constrained_values[i].first;
      dst.local_element(edge_constrained_indices[i]) =
        edge_constrained_values[i].first;
    }
}

// Performs the transposed operator evaluation. Since we have
// non-symmetric matrices this is different from the vmult call.
// TODO: implement this correctly (let's start with something symetric)
template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::Tvmult(VectorType &      dst,
                                              const VectorType &src) const
{
  this->vmult(dst, src);
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::vmult_interface_down(
  VectorType &      dst,
  VectorType const &src) const
{
  this->matrix_free.cell_loop(
    &NavierStokesOperatorBase::do_cell_integral_range, this, dst, src, true);

  // set constrained dofs as the sum of current dst value and src value
  for (unsigned int i = 0; i < constrained_indices.size(); ++i)
    dst.local_element(constrained_indices[i]) =
      src.local_element(constrained_indices[i]);
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::vmult_interface_up(
  VectorType &      dst,
  VectorType const &src) const
{
  if (has_edge_constrained_indices == false)
    {
      dst = number(0.);
      return;
    }

  dst = 0.0;

  // make a copy of src vector and set everything to 0 except edge
  // constrained dofs
  VectorType src_cpy;
  src_cpy.reinit(src, /*omit_zeroing_entries=*/false);

  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    src_cpy.local_element(edge_constrained_indices[i]) =
      src.local_element(edge_constrained_indices[i]);

  // do loop with copy of src
  this->matrix_free.cell_loop(&NavierStokesOperatorBase::do_cell_integral_range,
                              this,
                              dst,
                              src_cpy,
                              false);
}


// Computes the matrix efficiently using the optimized compute_matrix call.
// This function is usually only used to build the system matrix for the coarse
// grid level in the multigrid algorithm.
template <int dim, typename number>
const TrilinosWrappers::SparseMatrix &
NavierStokesOperatorBase<dim, number>::get_system_matrix() const
{
  if (system_matrix.m() == 0 && system_matrix.n() == 0)
    {
      const auto &dof_handler = this->matrix_free.get_dof_handler();

      TrilinosWrappers::SparsityPattern dsp(
        this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int ?
          dof_handler.locally_owned_mg_dofs(this->matrix_free.get_mg_level()) :
          dof_handler.locally_owned_dofs(),
        dof_handler.get_triangulation().get_communicator());

      if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
        MGTools::make_sparsity_pattern(dof_handler,
                                       dsp,
                                       this->matrix_free.get_mg_level(),
                                       this->constraints);
      else
        DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints);

      dsp.compress();
      system_matrix.reinit(dsp);

      MatrixFreeTools::compute_matrix(
        matrix_free,
        constraints,
        system_matrix,
        &NavierStokesOperatorBase::do_cell_integral_local,
        this);
    }
  return this->system_matrix;
}

template <int dim, typename number>
const MatrixFree<dim, number> &
NavierStokesOperatorBase<dim, number>::get_system_matrix_free() const
{
  return this->matrix_free;
}

template <int dim, typename number>
const AlignedVector<VectorizedArray<number>>
NavierStokesOperatorBase<dim, number>::get_element_size() const
{
  return this->element_size;
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::evaluate_non_linear_term(
  const VectorType &newton_step)
{
  const unsigned int n_cells = matrix_free.n_cell_batches();
  FECellIntegrator   phi(matrix_free);

  nonlinear_previous_values.reinit(n_cells, phi.n_q_points);
  nonlinear_previous_gradient.reinit(n_cells, phi.n_q_points);
  nonlinear_previous_hessian_diagonal.reinit(n_cells, phi.n_q_points);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values_plain(newton_step);
      phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                   EvaluationFlags::hessians);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          nonlinear_previous_values(cell, q)   = phi.get_value(q);
          nonlinear_previous_gradient(cell, q) = phi.get_gradient(q);
          nonlinear_previous_hessian_diagonal(cell, q) =
            phi.get_hessian_diagonal(q);
        }
    }
}

// The diagonal of the matrix is needed, e.g., for smoothers used in
// multigrid. Therefore, it is computed by a sequence of operator
// evaluations to unit basis vectors using the optimized compute_diagonal
// call.
template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::compute_inverse_diagonal(
  VectorType &diagonal) const
{
  this->matrix_free.initialize_dof_vector(diagonal);
  MatrixFreeTools::compute_diagonal(
    matrix_free,
    diagonal,
    &NavierStokesOperatorBase::do_cell_integral_local,
    this);

  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    diagonal.local_element(edge_constrained_indices[i]) = 0.0;

  for (auto &i : diagonal)
    i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0;
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::do_cell_integral_local(
  FECellIntegrator &integrator) const
{
  (void)integrator;

  AssertThrow(false,
              dealii::ExcMessage(
                "OperatorBase::do_face_integral() has not been implemented!"));
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::do_cell_integral_range(
  const MatrixFree<dim, number> &              matrix_free,
  VectorType &                                 dst,
  const VectorType &                           src,
  const std::pair<unsigned int, unsigned int> &range) const
{
  FECellIntegrator integrator(matrix_free, range);

  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      integrator.reinit(cell);

      integrator.read_dof_values(src);

      do_cell_integral_local(integrator);

      integrator.distribute_local_to_global(dst);
    }
}

template <int dim, typename number>
IndexSet
NavierStokesOperatorBase<dim, number>::get_refinement_edges(
  const MatrixFree<dim, number> &matrix_free)
{
  const unsigned int level = matrix_free.get_mg_level();

  std::vector<IndexSet> refinement_edge_indices;
  refinement_edge_indices.clear();
  const unsigned int nlevels =
    matrix_free.get_dof_handler().get_triangulation().n_global_levels();
  refinement_edge_indices.resize(nlevels);
  for (unsigned int l = 0; l < nlevels; l++)
    refinement_edge_indices[l] =
      IndexSet(matrix_free.get_dof_handler().n_dofs(l));

  MGTools::extract_inner_interface_dofs(matrix_free.get_dof_handler(),
                                        refinement_edge_indices);
  return refinement_edge_indices[level];
}

template class NavierStokesOperatorBase<2, double>;
template class NavierStokesOperatorBase<3, double>;

template <int dim, typename number>
NavierStokesSUPGPSPGOperator<dim, number>::NavierStokesSUPGPSPGOperator()
{}

template <int dim, typename number>
void
NavierStokesSUPGPSPGOperator<dim, number>::do_cell_integral_local(
  FECellIntegrator &integrator) const
{
  using FECellIntegratorType =
    FEEvaluation<dim, -1, 0, dim + 1, double, VectorizedArray<double>>;

  integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                      EvaluationFlags::hessians);

  for (unsigned int q = 0; q < integrator.n_q_points; ++q)
    {
      // TODO
    }

  integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

template class NavierStokesSUPGPSPGOperator<2, double>;
template class NavierStokesSUPGPSPGOperator<3, double>;
