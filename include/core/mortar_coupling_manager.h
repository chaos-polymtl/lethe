// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_core_mortar_coupling_manager_h
#define lethe_core_mortar_coupling_manager_h

#include <deal.II/base/config.h>

#if DEAL_II_VERSION_GTE(9, 7, 0)

#  include <core/parameters.h>
#  include <core/utilities.h>

#  include <deal.II/base/mpi_noncontiguous_partitioner.h>
#  include <deal.II/base/quadrature_lib.h>

#  include <deal.II/fe/fe_system.h>

#  include <deal.II/lac/trilinos_sparse_matrix.h>

#  include <deal.II/matrix_free/fe_point_evaluation.h>

using namespace dealii;

/**
 * @brief Base class for the mortar manager
 */
template <int dim>
class MortarManagerBase
{
public:
  template <int dim2>
  MortarManagerBase(const unsigned int      n_subdivisions,
                    const Quadrature<dim2> &quadrature,
                    const double            radius,
                    const double            rotation_angle);

  /**
   * @brief Verify if cells of the inner and outer domains are aligned
   */
  bool
  is_mesh_aligned() const;

  /**
   * @brief Returns the total number of mortars
   */
  unsigned int
  get_n_total_mortars() const;

  /**
   * @brief Returns the number of mortars per face
   */
  unsigned int
  get_n_mortars() const;

  /**
   * @brief Returns the indices of all mortars at both sides of the interface
   *
   * @param[in] face_center Face center
   */
  std::vector<unsigned int>
  get_mortar_indices(const Point<dim> &face_center) const;

  /**
   * @brief Returns the total number of quadrature points at the inner/outer boundary interface
   */
  unsigned int
  get_n_total_points() const;

  /**
   * @brief Returns the coordinates of the quadrature points at both sides of the interface
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

protected:
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

  /**
   * @brief Returns the indices of all quadrature points at both sides of the interface
   *
   * @param[in] face_center Face center
   */
  std::vector<unsigned int>
  get_indices_internal(const Point<dim>  &face_center,
                       const unsigned int n_quadrature_points) const;

  /**
   * @brief Convert radiant to quadrature point in real space.
   */
  virtual Point<dim>
  from_1D(const double radiant) const = 0;

  /**
   * @brief Convert quadrature point in real space to radiant.
   */
  virtual double
  to_1D(const Point<dim> &point) const = 0;

  /**
   * @brief Return the normal for a given quadrature point.
   */
  virtual Tensor<1, dim, double>
  get_normal(const Point<dim> &point) const = 0;

  /// Number of cells at the interface between inner and outer domains
  const unsigned int n_subdivisions;
  /// Mortar quadrature
  Quadrature<1> quadrature;
  /// Number of quadrature points per cell
  const unsigned int n_quadrature_points;
  /// Radius at the interface between inner and outer domains
  const double radius;
  /// Rotation angle for the inner domain
  const double rotation_angle;
};

/**
 * @brief Compute the number of subdivisions at the rotor-stator interface and the rotor radius
 * @param[in] dof_handler DoFHandler associated to the triangulation
 * @param[in] mortar_parameters The information about the mortar method
 * control, including the rotor mesh parameters
 *
 * @return n_subdivisions Number of cells at the interface between inner
 * and outer domains
 * @return radius Radius at the interface between inner and outer domains
 */
template <int dim>
std::pair<unsigned int, double>
compute_n_subdivisions_and_radius(
  const Triangulation<dim>      &triangulation,
  const Parameters::Mortar<dim> &mortar_parameters);

/**
 * @brief Construct oversampled quadrature
 *
 * @param[in] quadrature Quadrature for local cell operations
 * @param[in] mortar_parameters The information about the mortar method
 * control, including the rotor mesh parameters
 *
 * @return Quadrature oversampled
 */
template <int dim>
Quadrature<dim>
construct_quadrature(const Quadrature<dim>         &quadrature,
                     const Parameters::Mortar<dim> &mortar_parameters);


template <int dim>
class MortarManagerCircle : public MortarManagerBase<dim>
{
public:
  template <int dim2>
  MortarManagerCircle(const unsigned int      n_subdivisions,
                      const Quadrature<dim2> &quadrature,
                      const double            radius,
                      const double            rotation_angle);

  template <int dim2>
  MortarManagerCircle(const Quadrature<dim2>        &quadrature,
                      const DoFHandler<dim2>        &dof_handler,
                      const Parameters::Mortar<dim> &mortar_parameters);

protected:
  Point<dim>
  from_1D(const double rad) const override;

  double
  to_1D(const Point<dim> &point) const override;

  Tensor<1, dim, double>
  get_normal(const Point<dim> &point) const override;
};


template <int dim>
template <int dim2>
MortarManagerBase<dim>::MortarManagerBase(const unsigned int n_subdivisions,
                                          const Quadrature<dim2> &quadrature_in,
                                          const double            radius,
                                          const double rotation_angle)
  : n_subdivisions(n_subdivisions)
  , quadrature(quadrature_in.get_tensor_basis()[0])
  , n_quadrature_points(quadrature.size())
  , radius(radius)
  , rotation_angle(rotation_angle)
{}


template <int dim>
template <int dim2>
MortarManagerCircle<dim>::MortarManagerCircle(
  const unsigned int      n_subdivisions,
  const Quadrature<dim2> &quadrature,
  const double            radius,
  const double            rotation_angle)
  : MortarManagerBase<dim>(n_subdivisions, quadrature, radius, rotation_angle)
{}


template <int dim>
template <int dim2>
MortarManagerCircle<dim>::MortarManagerCircle(
  const Quadrature<dim2>        &quadrature,
  const DoFHandler<dim2>        &dof_handler,
  const Parameters::Mortar<dim> &mortar_parameters)
  : MortarManagerCircle(
      compute_n_subdivisions_and_radius(dof_handler.get_triangulation(),
                                        mortar_parameters)
        .first,
      quadrature,
      compute_n_subdivisions_and_radius(dof_handler.get_triangulation(),
                                        mortar_parameters)
        .second,
      mortar_parameters.rotor_angular_velocity->value(Point<dim>()))
{}


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

/**
 * @brief Compute scalar product
 *
 * @param[in, out] v_gradient Rank-2 tensor where result is stored
 * @param[in] u_gradient Rank-2 tensor
 * @param[in] factor Scalar factor
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<1, 1, Number>       &v_gradient,
                        const Tensor<1, 1, Number> &u_gradient,
                        const Number               &factor)
{
  v_gradient[0] += u_gradient[0] * factor;
}



template <int dim, typename Number>
struct CouplingEvaluationData
{
  /// Penalty factor (akin to penalty factor in SIPG)
  Number penalty_factor;
  /// Penalty parameter in symmetric interior penalty Galerkin (SIPG) method
  std::vector<Number> all_penalty_parameter;
  /// Weights of quadrature points
  std::vector<Number> all_weights;
  // Normal vectors of quadrature points
  std::vector<Tensor<1, dim, Number>> all_normals;
};



/**
 * @brief Base class for coupling evaluation routines.
 */
template <int dim, typename Number>
class CouplingEvaluationBase
{
public:
  /**
   * Number of data points of type Number associated to a quadrature point.
   */
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
   * @param[in,out] buffer Temporary vector where data is stored before being
   * passes to the system matrix
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
   *
   * Notation referring to quantities on both mortar sides:
   * {{.}} = avg(.), and [.]   = jump(.)   *
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
  /**
   * @brief Constructor.
   *
   * @param[in] bid_m Boundary ID of the face whose outwards-pointing
   *   normal shows in the same direction as the normal provided by
   *   @p mortar_manager.
   * @param[in] bid_p Boundary ID of the face whose outwards-pointing
   *   normal shows in the opposite direction as the normal provided by
   *   @p mortar_manager.
   */
  CouplingOperator(
    const Mapping<dim>                                        &mapping,
    const DoFHandler<dim>                                     &dof_handler,
    const AffineConstraints<Number>                           &constraints,
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator,
    const std::shared_ptr<MortarManagerBase<dim>>              mortar_manager,
    const unsigned int                                         bid_m,
    const unsigned int                                         bid_p,
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
  template <typename VectorType>
  void
  vmult_add(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Add mortar coupling terms in diagonal entries
   *
   * @param[in, out] diagonal Matrix diagonal
   */
  template <typename VectorType>
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
   * @param[in] factor Penalty factor (akin to penalty factor in SIPG)
   *
   * @return penalty factor value
   * penalty_factor = (degree + 1)^2
   */
  Number
  compute_penalty_factor(const unsigned int degree, const Number factor) const;

  /**
   * @brief Compute penalty parameter in a cell
   *
   * @param[in] cell Cell iterator
   * @return Penalty parameter
   *
   * @return penalty parameter value from SIPG method
   * penalty_parameter = (A(∂Ω_e \ Γ_h)/2 + A(∂Ω_e ∩ Γ_h))/V(Ω_e)
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

  std::vector<std::tuple<std::vector<double>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         std::vector<Tensor<1, dim, Number>>,
                         Number>>
    all_intersections;

protected:
  /// Number of data points per quadrature point
  unsigned int q_data_size;

  /// Boundary ID of the inner domain (rotor)
  const unsigned int bid_m;
  /// Boundary ID of the outer domain (stator)
  const unsigned int bid_p;

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
  std::shared_ptr<MortarManagerBase<dim>>              mortar_manager;
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

  /// Finite element that matches the components `n_components` components
  /// starting at component with index `first_selected_component`
  const FESystem<dim> fe_sub;

  /// Interface to the evaluation of mortar coupling interpolated solution
  mutable FEPointIntegrator phi_m;

  /// Relevant dof indices
  std::vector<unsigned int> relevant_dof_indices;
};

#endif
#endif
