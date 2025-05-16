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
 */
template <int dim>
class MortarManager
{
public:
  MortarManager(const unsigned int n_subdivisions,
                const unsigned int n_quadrature_points,
                const double       radius,
                const double       rotation_angle);

  /**
   * @brief Verify if cells of the inner and outer domains are aligned
   */
  bool
  is_mesh_aligned() const;

  /**
   * @brief Returns the total number of quadrature points at the inner/outer boundary interface
   */
  unsigned int
  get_n_total_points() const;

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the interface
   *
   * @return points Coordinate of quadrature points of the cell
   */
  unsigned int
  get_n_points() const;

  /**
   * @brief Returns the indices of all quadrature points at both sides of the interface
   *
   * @param[in] face_center Face center
   */
  std::vector<unsigned int>
  get_indices(const Point<dim> &face_center) const;

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the interface
   *
   * @param[in] face_center Face center
   *
   * @return points Coordinate of quadrature points of the cell
   */
  std::vector<Point<dim>>
  get_points(const Point<dim> &face_center) const;

  /**
   * @brief Returns the coordinates of the quadrature points at the interface
   *
   * @param[in] face_center Face center
   *
   * @return points Coordinate of quadrature points of the cell
   */
  std::vector<Point<1>>
  get_points_ref(const Point<dim> &face_center) const;

  /**
   * @brief Returns the weights of the quadrature points at both sides of the interface
   *
   * @param[in] face_center Face center
   *
   * @return points Angular weights of quadrature points of the cell
   */
  std::vector<double>
  get_weights(const Point<dim> &face_center) const;

  /**
   * @brief Returns the normal vector for the quadrature points
   *
   * @param[in] face_center Face center
   *
   * @return result Normal vectors of the cell quadrature points
   */
  std::vector<Tensor<1, dim, double>>
  get_normals(const Point<dim> &face_center) const;

private:
  /**
   * @brief Returns the mesh alignment type and cell index
   *
   * @param[in] face_center Face center
   *
   * @return type Cell configuration type at the interface
   * type = 0: mesh aligned
   * type = 1: mesh not aligned, inner domain (allows rotation)
   * type = 2: mesh not aligned, outer domain (fixed)
   * @return id Index of the cell in which lies the rotated cell center
   */
  std::pair<unsigned int, unsigned int>
  get_config(const Point<dim> &face_center) const;

  /// Number of cells at the interface between inner and outer domains
  const unsigned int n_subdivisions;
  /// Number of quadrature points per cell
  const unsigned int n_quadrature_points;
  /// Radius at the interface between inner and outer domains
  const double radius;
  /// Rotation angle for the inner domain
  const double rotation_angle;
  /// Mortar quadrature
  QGauss<1> quadrature;
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

/**
 * @brief Compute scalar product
 *
 * @param[in, out] v_gradient Rank-1 tensor within rank-1 tensor where result is
 * stored
 * @param[in] u_gradient Rank-2 tensor
 * @param[in] factor Scalar factor
 */
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

/**
 * @brief Compute scalar product
 *
 * @param[in, out] v_gradient Rank-2 tensor where result is stored
 * @param[in] u_gradient Rank-2 tensor
 * @param[in] factor Scalar factor
 */
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
struct CouplingEvaluationData
{
  /// Penalty factor (akin to symmetric interior penalty factor in SIPG)
  Number                              penalty_factor;
  std::vector<Number>                 all_penalty_parameter;
  std::vector<Number>                 all_weights;
  std::vector<Tensor<1, dim, Number>> all_normals;
};



/**
 * @brief Base class for coupling evaluation routines.
 */
template <int dim, typename Number>
class CouplingEvaluationBase
{
public:
  virtual unsigned int
  data_size() const = 0;

  /**
   * @brief Return relevant dof indices.
   */
  virtual const std::vector<unsigned int> &
  get_relevant_dof_indices() const = 0;

  /**
   * @brief Set up mapping information
   *
   * @param[in] cell Cell iterator
   * @param[in] points List of points where FEPointIntegrator should be
   * evaluated
   */
  virtual void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const = 0;

  /**
   * @brief Evaluate values and gradients at the coupling entries
   *
   * @param[in] buffer Temporary vector where data is stored before being passes
   * to the system matrix
   * @param[in] ptr_q Pointer for the quadrature point index related to the
   * rotor-stator interface
   * @param[in] q_stride Pointer for the cell index in which the quadrature
   * point lies in
   * @param[in] all_value_m Number of values stored in the mortar side of the
   * interface
   */
  virtual void
  local_evaluate(const CouplingEvaluationData<dim, Number> &data,
                 const Vector<Number>                      &buffer,
                 const unsigned int                         ptr_q,
                 const unsigned int                         q_stride,
                 Number *all_value_m) const = 0;

  /**
   * @brief Perform integral of mortar elements at the rotor-stator interface
   *
   * @param[in] buffer Temporary vector where data is stored before being passes
   * to the system matrix
   * @param[in] ptr_q Pointer for the quadrature point index related to the
   * rotor-stator interface
   * @param[in] q_stride Pointer for the cell index in which the quadrature
   * point lies in
   * @param[in] all_value_m Number of values stored in the mortar side of the
   * interface
   * @param[in] all_value_p Number of values stored in the non-mortar side of
   * the interface
   */
  virtual void
  local_integrate(const CouplingEvaluationData<dim, Number> &data,
                  Vector<Number>                            &buffer,
                  const unsigned int                         ptr_q,
                  const unsigned int                         q_stride,
                  Number                                    *all_value_m,
                  Number *all_value_p) const = 0;
};



/**
 * @brief Base class for the Coupling Operator Base
 */
template <int dim, typename Number>
class CouplingOperator
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  CouplingOperator(
    const Mapping<dim>                                        &mapping,
    const DoFHandler<dim>                                     &dof_handler,
    const AffineConstraints<Number>                           &constraints,
    const Quadrature<dim>                                      quadrature,
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator,
    const unsigned int                                         n_subdivisions,
    const double                                               radius,
    const double                                               rotation_angle,
    const unsigned int                                         bid_rotor,
    const unsigned int                                         bid_stator,
    const double                                               sip_factor);

  /**
   * @brief Return object containing problem constraints
   *
   * @return AffineConstraints
   */
  const AffineConstraints<Number> &
  get_affine_constraints() const;

  /**
   * @brief Add matrix-vector multiplication
   *
   * @param[in, out] dst Destination vector holding the result
   * @param[in] src Input source vector
   */
  void
  vmult_add(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Add mortar coupling terms in diagonal entries
   *
   * @param[in, out] diagonal Matrix diagonal
   */
  void
  add_diagonal_entries(VectorType &diagonal) const;

  /**
   * @brief Add mortar coupling terms in the sparsity pattern
   *
   * @param[in, out] dsp Dynamic Sparsity Pattern object
   */
  void
  add_sparsity_pattern_entries(SparsityPatternBase &dsp) const;

  /**
   * @brief Add mortar coupling terms in the system matrix
   *
   * @param[in, out] system_matrix System matrix
   */
  void
  add_system_matrix_entries(
    TrilinosWrappers::SparseMatrix &system_matrix) const;

private:
  /**
   * @brief Compute penalty factor used in weak imposition of coupling at the rotor-stator interface
   *
   * @param[in] degree Polynomial degree of the FE approximation
   * @param[in] factor Penalty factor (akin to symmetric interior penalty factor
   * in SIPG)
   *
   * @return penalty factor value
   */
  Number
  compute_penalty_factor(const unsigned int degree, const Number factor) const;

  /**
   * @brief Compute penalty parameter in a cell
   *
   * @param[in] cell Cell iterator
   * @return Penalty parameter
   *
   * @return penalty parameter value
   */
  Number
  compute_penalty_parameter(
    const typename Triangulation<dim>::cell_iterator &cell) const;

  /**
   * @brief Returns angle of a point (cell center)
   *
   * @param[in] cell Cell iterator
   * @param[in] face Face iterator
   *
   * @return Angle in radians
   */
  Point<dim>
  get_face_center(const typename Triangulation<dim>::cell_iterator &cell,
                  const typename Triangulation<dim>::face_iterator &face) const;

  /**
   * @brief Returns dof indices
   *
   * @param[in] cell Cell iterator
   */
  std::vector<types::global_dof_index>
  get_dof_indices(
    const typename DoFHandler<dim>::active_cell_iterator &cell) const;


  /// Mapping of the domain
  const Mapping<dim> &mapping;
  /// DoFHandler associated to the triangulation
  const DoFHandler<dim> &dof_handler;
  /// Object with the constrains according to DoFs
  const AffineConstraints<Number> &constraints;
  /// Quadrature required for local operations on cells
  const Quadrature<dim> quadrature;

  std::vector<std::tuple<std::vector<double>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         std::vector<Tensor<1, dim, Number>>,
                         Number>>
    all_intersections;

protected:
  std::shared_ptr<MortarManager<dim>> mortar_manager_q;
  std::shared_ptr<MortarManager<dim>> mortar_manager_cell;

  /// Number of data points per quadrature point
  unsigned int N;

  /// Boundary ID of the inner domain (rotor)
  const unsigned int bid_rotor;
  /// Boundary ID of the outer domain (stator)
  const unsigned int bid_stator;

  /// List of relevant DoF indices per cell
  std::vector<unsigned int> relevant_dof_indices;
  /// Number of DoFs per cell
  unsigned int n_dofs_per_cell;

  Utilities::MPI::NoncontiguousPartitioner partitioner;
  Utilities::MPI::NoncontiguousPartitioner partitioner_cell;

  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> dof_indices_ghost;

  /// Vectors storing information at quadrature points for all cells at the
  /// rotor-stator interface
  std::vector<Point<dim, Number>>     all_points_ref;
  CouplingEvaluationData<dim, Number> data;

  /// Constraints extended according to mortar entries
  AffineConstraints<Number> constraints_extended;

  std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator;
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
class CouplingEvaluationSIPG : public CouplingEvaluationBase<dim, Number>
{
public:
  using FEPointIntegrator = FEPointEvaluation<n_components, dim, dim, Number>;
  using value_type        = typename FEPointIntegrator::value_type;

  CouplingEvaluationSIPG(const Mapping<dim>    &mapping,
                         const DoFHandler<dim> &dof_handler,
                         const unsigned int     first_selected_component = 0);

  unsigned int
  data_size() const override;

  const std::vector<unsigned int> &
  get_relevant_dof_indices() const override;

  void
  local_reinit(
    const typename Triangulation<dim>::cell_iterator &cell,
    const ArrayView<const Point<dim, Number>>        &points) const override;

  void
  local_evaluate(const CouplingEvaluationData<dim, Number> &data,
                 const Vector<Number>                      &buffer,
                 const unsigned int                         ptr_q,
                 const unsigned int                         q_stride,
                 Number *all_value_m) const override;

  void
  local_integrate(const CouplingEvaluationData<dim, Number> &data,
                  Vector<Number>                            &buffer,
                  const unsigned int                         ptr_q,
                  const unsigned int                         q_stride,
                  Number                                    *all_value_m,
                  Number *all_value_p) const override;

  const FESystem<dim> fe_sub;
  /// Interface to the evaluation of mortar coupling interpolated solution
  mutable FEPointIntegrator phi_m;

  std::vector<unsigned int> relevant_dof_indices;
};

#endif
#endif
