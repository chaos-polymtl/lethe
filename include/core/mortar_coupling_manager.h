// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_core_mortar_coupling_manager_h
#define lethe_core_mortar_coupling_manager_h

#include <deal.II/base/config.h>

#if DEAL_II_VERSION_GTE(9, 7, 0)

#  include <core/parameters.h>
#  include <core/utilities.h>

#  include <deal.II/base/mpi_noncontiguous_partitioner.h>
#  include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>
#  include <deal.II/base/quadrature_lib.h>

#  include <deal.II/fe/fe_system.h>

#  include <deal.II/lac/trilinos_precondition.h>
#  include <deal.II/lac/trilinos_solver.h>
#  include <deal.II/lac/trilinos_sparse_matrix.h>
#  include <deal.II/lac/trilinos_sparsity_pattern.h>

#  include <deal.II/matrix_free/fe_evaluation.h>
#  include <deal.II/matrix_free/fe_point_evaluation.h>
#  include <deal.II/matrix_free/matrix_free.h>
#  include <deal.II/matrix_free/tools.h>

using namespace dealii;

/**
 * @brief Base class for the mortar manager
 * @param n_subdivisions Number of cells at the interface between inner and outer domains
 * @param n_quadrature_points Number of quadrature points per cell
 * @param radius Radius at the interface between inner and outer domains
 * @param rotate_pi Rotation angle for the inner domain
 */
template <int dim>
class MortarManager
{
public:
  MortarManager(const unsigned int n_subdivisions,
                const unsigned int n_quadrature_points,
                const double       radius,
                const double       rotate_pi);

  /**
   * @brief Verify if cells of the inner and outer domains are aligned
   */
  bool
  is_mesh_aligned() const;

  /**
   * @brief Returns the total number of quadrature points at the inner/outer boundary interface
   * // TODO keep only one version of this function
   */
  unsigned int
  get_n_points() const;

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the inerface
   * @param[in] angle_cell_center Angle between cell center and x-axis (in
   * radians)
   *
   * @return points Coordinate of quadrature points of the cell
   */
  unsigned int
  get_n_points(const double &angle_cell_center) const;

  /**
   * @brief Returns the indices of all quadrature points at both sides of the interface
   * @param[in] angle_cell_center Angular coordinate of cell center
   */
  std::vector<unsigned int>
  get_indices(const double &angle_cell_center) const;

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the inerface
   * @param[in] rad Angular coordinate of cell center
   * @param[out] points Coordinate of quadrature points of the cell
   */
  std::vector<Point<dim>>
  get_points(const double rad) const;

  /**
   * @brief Returns the coordinates of the quadrature points at the interface
   * @param[in] angle_cell_center Angle between cell center and x-axis (in
   * radians)
   *
   * @return points Coordinate of quadrature points of the cell
   */
  std::vector<Point<1>>
  get_points_ref(const double angle_cell_center) const;

  /**
   * @brief Returns the weights of the quadrature points at both sides of the interface
   * @param[in] angle_cell_center Angular coordinate of cell center
   *
   * @return points Angular weights of quadrature points of the cell
   */
  std::vector<double>
  get_weights(const double &angle_cell_center) const;

  /**
   * @brief Returns the normal vector for the quadrature points
   * @param[in] angle_cell_center Angular coordinate of cell center
   *
   * @return result Normal vectors of the cell quadrature points
   */
  std::vector<Tensor<1, dim, double>>
  get_normals(const double &angle_cell_center) const;

private:
  /**
   * @brief Returns the mesh alignement type and cell index
   * @param[in] rad Angular coordinate of cell center
   *
   * @return type Cell configuration type at the interface
   * type = 0: mesh aligned
   * type = 1: mesh not aligned, inner domain (allows rotation)
   * type = 2: mesh not aligned, outer domain (fixed)
   * @return id Index of the cell in which lies the rotated cell center
   */
  std::pair<unsigned int, unsigned int>
  get_config(const double &rad) const;

  const unsigned int n_subdivisions;
  const unsigned int n_quadrature_points;
  const double       radius;
  const double       rotate_pi;
  QGauss<1>          quadrature;
};


/**
 * @brief Compute inner product
 *
 * @param[in] grad Rank-1 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-0 tensor
 */
template <int dim, typename Number>
Number
contract(const Tensor<1, dim, Number> &grad,
         const Tensor<1, dim, Number> &normal)
{
  return grad * normal;
}

/**
 * @brief Compute inner product
 *
 * @param[in] grad Rank-2 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-1 tensor
 */
template <int dim, typename Number>
Tensor<1, dim, Number>
contract(const Tensor<2, dim, Number> &grad,
         const Tensor<1, dim, Number> &normal)
{
  return grad * normal;
}

/**
 * @brief Compute outer product for the multicomponent case
 *
 * @param[in] grad Rank-1 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-0 tensor
 *
 */
template <int n_components, int dim, typename Number>
Tensor<1, n_components, Number>
contract(const Tensor<1, n_components, Tensor<1, dim, Number>> &grad,
         const Tensor<1, dim, Number>                          &normal)
{
  Tensor<1, n_components, Number> result;

  for (int r = 0; r < n_components; ++r)
    result[r] = grad[r] * normal;

  return result;
}

/**
 * @brief Compute outer product
 *
 * @param[in] value Rank-0 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-1 tensor
 */
template <int dim, typename Number>
Tensor<1, dim, Number>
outer(const Number &value, const Tensor<1, dim, Number> &normal)
{
  return value * normal;
}


/**
 * @brief Compute outer product
 *
 * @param[in] value Rank-1 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-2 tensor
 */
template <int dim, typename Number>
Tensor<2, dim, Number>
outer(const Tensor<1, dim, Number> &value, const Tensor<1, dim, Number> &normal)
{
  Tensor<2, dim, Number> result;

  for (unsigned int c = 0; c < dim; ++c)
    result[c] = value[c] * normal;

  return result;
}

/**
 * @brief Compute outer product
 *
 * @param[in] value Rank-1 tensor for n_components
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-2 tensor for n_components
 */
template <int n_components, int dim, typename Number>
Tensor<1, n_components, Tensor<1, dim, Number>>
outer(const Tensor<1, n_components, Number> &value,
      const Tensor<1, dim, Number>          &normal)
{
  Tensor<1, n_components, Tensor<1, dim, Number>> result;

  for (unsigned int c = 0; c < n_components; ++c)
    result[c] = value[c] * normal;

  return result;
}

template <int dim, int dim_, typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<1, dim_, Tensor<1, dim, Number>> &v_gradient,
                        const Tensor<2, dim, Number>            &u_gradient,
                        const Number                            &factor)
{
  for (unsigned int d = 0; d < dim; ++d)
    v_gradient[d][d] += u_gradient[d][d] * factor;

  for (unsigned int e = 0; e < dim; ++e)
    for (unsigned int d = e + 1; d < dim; ++d)
      {
        const auto tmp = (u_gradient[d][e] + u_gradient[e][d]) * (factor * 0.5);
        v_gradient[d][e] += tmp;
        v_gradient[e][d] += tmp;
      }
}

template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<2, dim, Number>       &v_gradient,
                        const Tensor<2, dim, Number> &u_gradient,
                        const Number                 &factor)
{
  for (unsigned int d = 0; d < dim; ++d)
    v_gradient[d][d] += u_gradient[d][d] * factor;

  for (unsigned int e = 0; e < dim; ++e)
    for (unsigned int d = e + 1; d < dim; ++d)
      {
        const auto tmp = (u_gradient[d][e] + u_gradient[e][d]) * (factor * 0.5);
        v_gradient[d][e] += tmp;
        v_gradient[e][d] += tmp;
      }
}

template <int dim, typename Number>
class CouplingOperatorBase
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  CouplingOperatorBase(const Mapping<dim>              &mapping,
                       const DoFHandler<dim>           &dof_handler,
                       const AffineConstraints<Number> &constraints,
                       const Quadrature<dim>            quadrature,
                       const unsigned int               n_subdivisions,
                       const unsigned int               n_components,
                       const unsigned int               N,
                       const double                     radius,
                       const double                     rotate_pi,
                       const unsigned int               bid_0,
                       const unsigned int               bid_1,
                       const double                     sip_factor,
                       const std::vector<unsigned int>  relevant_dof_indices,
                       const double                     penalty_factor_grad);

  const AffineConstraints<Number> &
  get_affine_constraints() const;

  void
  vmult_add(VectorType &dst, const VectorType &src) const;

  void
  add_diagonal_entries(VectorType &diagonal) const;

  void
  add_sparsity_pattern_entries(SparsityPatternBase &dsp) const;

  void
  add_system_matrix_entries(
    TrilinosWrappers::SparseMatrix &system_matrix) const;

private:
  Number
  compute_penalty_factor(const unsigned int degree, const Number factor) const;

  Number
  compute_penalty_parameter(
    const typename Triangulation<dim>::cell_iterator &cell) const;

  double
  get_rad(const typename Triangulation<dim>::cell_iterator &cell,
          const typename Triangulation<dim>::face_iterator &face) const;

  std::vector<types::global_dof_index>
  get_dof_indices(
    const typename DoFHandler<dim>::active_cell_iterator &cell) const;


  const Mapping<dim>              &mapping;
  const DoFHandler<dim>           &dof_handler;
  const AffineConstraints<Number> &constraints;
  const Quadrature<dim>            quadrature;

  std::vector<std::tuple<std::vector<double>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         std::vector<Tensor<1, dim, Number>>,
                         Number>>
    all_intersections;

protected:
  Number penalty_factor;
  Number penalty_factor_grad;

  std::shared_ptr<MortarManager<dim>> mortar_manager_q;
  std::shared_ptr<MortarManager<dim>> mortar_manager_cell;

  const unsigned int n_components;
  const unsigned int N;

  const unsigned int bid_0;
  const unsigned int bid_1;

  const std::vector<unsigned int> relevant_dof_indices;
  const unsigned int              n_dofs_per_cell;

  Utilities::MPI::NoncontiguousPartitioner partitioner;
  Utilities::MPI::NoncontiguousPartitioner partitioner_cell;

  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> dof_indices_ghost;

  std::vector<Number>                 all_penalty_parameter;
  std::vector<Number>                 all_weights;
  std::vector<Point<dim, Number>>     all_points_ref;
  std::vector<Tensor<1, dim, Number>> all_normals;

  AffineConstraints<Number> constraints_extended;

  virtual void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const = 0;

  virtual void
  local_evaluate(const Vector<Number> &buffer,
                 const unsigned int    ptr_q,
                 const unsigned int    q_stride,
                 Number               *all_value_m) const = 0;

  virtual void
  local_integrate(Vector<Number>    &buffer,
                  const unsigned int ptr_q,
                  const unsigned int q_stride,
                  Number            *all_value_m,
                  Number            *all_value_p) const = 0;
};



template <typename T>
class BufferRW
{
public:
  BufferRW(T *ptr, const unsigned int offset)
    : ptr(ptr ? (ptr + offset) : nullptr)
  {}

  void
  write(const T &in)
  {
    ptr[0] = in;
    ptr += 1;
  }

  template <int dim>
  void
  write(const Tensor<1, dim, T> &in)
  {
    for (unsigned int i = 0; i < dim; ++i)
      ptr[i] = in[i];

    ptr += dim;
  }

  template <typename T0>
  T0
  read() const
  {
    T0 result = {};

    if (ptr)
      read(result);

    return result;
  }

private:
  mutable T *ptr;

  template <int dim>
  void
  read(Tensor<1, dim, T> &out) const
  {
    for (unsigned int i = 0; i < dim; ++i)
      out[i] = ptr[i];

    ptr += dim;
  }

  void
  read(T &out) const
  {
    out = ptr[0];
    ptr += 1;
  }
};



template <int dim, int n_components, typename Number>
class CouplingOperator : public CouplingOperatorBase<dim, Number>
{
public:
  using FEPointIntegrator = FEPointEvaluation<n_components, dim, dim, Number>;
  using value_type        = typename FEPointIntegrator::value_type;

  CouplingOperator(const Mapping<dim>              &mapping,
                   const DoFHandler<dim>           &dof_handler,
                   const AffineConstraints<Number> &constraints,
                   const Quadrature<dim>            quadrature,
                   const unsigned int               n_subdivisions,
                   const double                     radius,
                   const double                     rotate_pi,
                   const unsigned int               bid_0,
                   const unsigned int               bid_1,
                   const double                     sip_factor = 1.0,
                   const unsigned int first_selected_component = 0,
                   const double       penalty_factor_grad      = 1.0);

  /**
   * @brief Constructor of the class
   * @param[in] mapping Mapping of the domain
   * @param[in] dof_handler DoFHandler associated to the triangulation
   * @param[in] constraints Object with the constrains according to DoFs
   * @param[in] quadrature Required for local operations on cells
   * @param[in] mortar_parameters The information about the mortar method
   * control, including the rotor mesh parameters
   */
  CouplingOperator(const Mapping<dim>              &mapping,
                   const DoFHandler<dim>           &dof_handler,
                   const AffineConstraints<Number> &constraints,
                   const Quadrature<dim>           &quadrature,
                   const Parameters::Mortar<dim>   &mortar_parameters,
                   const unsigned int first_selected_component = 0,
                   const double       penalty_factor_grad      = 1.0);

  static std::vector<unsigned int>
  get_relevant_dof_indices(const FiniteElement<dim> &fe,
                           const unsigned int        first_selected_component);

  void
  local_reinit(
    const typename Triangulation<dim>::cell_iterator &cell,
    const ArrayView<const Point<dim, Number>>        &points) const override;

  void
  local_evaluate(const Vector<Number> &buffer,
                 const unsigned int    ptr_q,
                 const unsigned int    q_stride,
                 Number               *all_value_m) const override;

  void
  local_integrate(Vector<Number>    &buffer,
                  const unsigned int ptr_q,
                  const unsigned int q_stride,
                  Number            *all_value_m,
                  Number            *all_value_p) const override;

  const FESystem<dim>       fe_sub;
  mutable FEPointIntegrator phi_m;
};

#endif
#endif
