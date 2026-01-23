// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>

using namespace dealii;


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
  virtual ~CouplingEvaluationBase() = default;
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
   * @param[in] data Coupling evaluation data containing values of normals,
   * weights, and penalty parameters at each integration point
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
   * @param[in] data Coupling evaluation data containing values of normals,
   * weights, and penalty parameters at each integration point
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

/*-------------- Coupling evaluation Stokes -------------------------------*/
template <int dim, typename Number>
class CouplingEvaluationStokes : public CouplingEvaluationBase<dim, Number>
{
public:
  using FEPointIntegratorU = FEPointEvaluation<dim, dim, dim, Number>;
  using FEPointIntegratorP = FEPointEvaluation<1, dim, dim, Number>;

  using u_value_type = typename FEPointIntegratorU::value_type;

  CouplingEvaluationStokes(const Mapping<dim>    &mapping,
                           const DoFHandler<dim> &dof_handler,
                           const bool weak_pressure_gradient_term   = true,
                           const bool weak_velocity_divergence_term = true)
    : fe_sub_u(dof_handler.get_fe().base_element(
                 dof_handler.get_fe().component_to_base_index(0).first),
               dim)
    , fe_sub_p(dof_handler.get_fe().base_element(
                 dof_handler.get_fe().component_to_base_index(dim).first),
               1)
    , phi_u_m(mapping, fe_sub_u, update_values | update_gradients)
    , phi_p_m(mapping, fe_sub_p, update_values)
    , weak_pressure_gradient_term(weak_pressure_gradient_term)
    , weak_velocity_divergence_term(weak_velocity_divergence_term)
  {
    for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
      if (dof_handler.get_fe().system_to_component_index(i).first < dim)
        relevant_dof_indices.push_back(i);

    for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
      if (dof_handler.get_fe().system_to_component_index(i).first == dim)
        relevant_dof_indices.push_back(i);

    AssertDimension(dof_handler.get_fe().n_dofs_per_cell(),
                    relevant_dof_indices.size());
  }

  unsigned int
  data_size() const override
  {
    return 4 * dim;
  }

  const std::vector<unsigned int> &
  get_relevant_dof_indices() const override
  {
    return relevant_dof_indices;
  }

  void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const override
  {
    this->phi_u_m.reinit(cell, points);
    this->phi_p_m.reinit(cell, points);
  }

  void
  local_evaluate(const CouplingEvaluationData<dim, Number> &data,
                 const Vector<Number>                      &buffer,
                 const unsigned int                         ptr_q,
                 const unsigned int                         q_stride,
                 Number *all_value_m) const override
  {
    AssertDimension(buffer.size(),
                    fe_sub_u.n_dofs_per_cell() + fe_sub_p.n_dofs_per_cell());

    ArrayView<const Number> buffer_u(buffer.data() + 0,
                                     fe_sub_u.n_dofs_per_cell());
    ArrayView<const Number> buffer_p(buffer.data() + fe_sub_u.n_dofs_per_cell(),
                                     fe_sub_p.n_dofs_per_cell());

    this->phi_u_m.evaluate(buffer_u,
                           EvaluationFlags::values |
                             EvaluationFlags::gradients);
    this->phi_p_m.evaluate(buffer_p, EvaluationFlags::values);

    for (const auto q : this->phi_u_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        const auto normal = data.all_normals[q_index];

        const auto value_m    = this->phi_u_m.get_value(q);
        const auto gradient_m = contract(this->phi_u_m.get_gradient(q), normal);
        const auto p_value_m  = this->phi_p_m.get_value(q) * normal;

        BufferRW<Number> buffer_m(all_value_m, q * 4 * dim * q_stride);

        buffer_m.write(value_m);
        buffer_m.write(gradient_m);
        buffer_m.write(p_value_m);
      }
  }

  void
  local_integrate(const CouplingEvaluationData<dim, Number> &data,
                  Vector<Number>                            &buffer,
                  const unsigned int                         ptr_q,
                  const unsigned int                         q_stride,
                  Number                                    *all_value_m,
                  Number *all_value_p) const override
  {
    for (const auto q : this->phi_u_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        BufferRW<Number> buffer_m(all_value_m, q * 4 * dim * q_stride);
        BufferRW<Number> buffer_p(all_value_p, q * 4 * dim * q_stride);

        const auto value_m           = buffer_m.template read<u_value_type>();
        const auto value_p           = buffer_p.template read<u_value_type>();
        const auto normal_gradient_m = buffer_m.template read<u_value_type>();
        const auto normal_gradient_p = buffer_p.template read<u_value_type>();
        const auto normal_p_value_m  = buffer_m.template read<u_value_type>();
        const auto normal_p_value_p  = buffer_p.template read<u_value_type>();

        const auto JxW               = data.all_weights[q_index];
        const auto penalty_parameter = data.all_penalty_parameter[q_index];
        const auto normal            = data.all_normals[q_index];

        const auto u_value_avg  = (value_m + value_p) * 0.5;
        const auto u_value_jump = value_m - value_p;
        const auto u_gradient_avg =
          (normal_gradient_m - normal_gradient_p) * 0.5;
        const auto p_value_avg = (normal_p_value_m - normal_p_value_p) * 0.5;

        typename FEPointIntegratorU::value_type u_normal_gradient_avg_result =
          {};
        typename FEPointIntegratorU::value_type u_value_jump_result = {};
        typename FEPointIntegratorP::value_type p_value_jump_result = {};

        const double sigma = penalty_parameter * data.penalty_factor;

        if (true /*Laplace term*/)
          {
            // - (n avg(∇v), jump(u))
            u_normal_gradient_avg_result -= u_value_jump;

            // - (jump(v), avg(∇u) n)
            u_value_jump_result -= u_gradient_avg;

            // + (jump(v), σ jump(u))
            u_value_jump_result += sigma * u_value_jump;
          }

        if (weak_pressure_gradient_term)
          {
            // + (jump(v), avg(p) n)
            u_value_jump_result += p_value_avg;
          }
        else
          {
            // nothing to do
          }

        if (weak_velocity_divergence_term)
          {
            // + (jump(q), avg(u) n)
            if constexpr (dim == 1)
              p_value_jump_result += u_value_avg * normal[0];
            else
              p_value_jump_result += u_value_avg * normal;
          }
        else
          {
            // - (avg(q), jump(u) n)
            if constexpr (dim == 1)
              p_value_jump_result -= 0.5 * u_value_jump * normal[0];
            else
              p_value_jump_result -= 0.5 * u_value_jump * normal;
          }

        phi_u_m.submit_gradient(outer(u_normal_gradient_avg_result, normal) *
                                  0.5 * JxW,
                                q);
        phi_u_m.submit_value(u_value_jump_result * JxW, q);
        phi_p_m.submit_value(p_value_jump_result * JxW, q);
      }

    AssertDimension(buffer.size(),
                    fe_sub_u.n_dofs_per_cell() + fe_sub_p.n_dofs_per_cell());

    ArrayView<Number> buffer_u(buffer.data() + 0, fe_sub_u.n_dofs_per_cell());
    ArrayView<Number> buffer_p(buffer.data() + fe_sub_u.n_dofs_per_cell(),
                               fe_sub_p.n_dofs_per_cell());

    this->phi_u_m.test_and_sum(buffer_u,
                               EvaluationFlags::values |
                                 EvaluationFlags::gradients);
    this->phi_p_m.test_and_sum(buffer_p, EvaluationFlags::values);
  }

  const FESystem<dim>        fe_sub_u;
  const FESystem<dim>        fe_sub_p;
  mutable FEPointIntegratorU phi_u_m;
  mutable FEPointIntegratorP phi_p_m;

  const bool weak_pressure_gradient_term;
  const bool weak_velocity_divergence_term;

  std::vector<unsigned int> relevant_dof_indices;
};


/*-------------- Coupling evaluation SIPG -------------------------------*/
template <int dim, int n_components, typename Number>
class CouplingEvaluationSIPG : public CouplingEvaluationBase<dim, Number>
{
public:
  using FEPointIntegrator = FEPointEvaluation<n_components, dim, dim, Number>;
  using value_type        = typename FEPointIntegrator::value_type;

  CouplingEvaluationSIPG(const Mapping<dim>    &mapping,
                         const DoFHandler<dim> &dof_handler,
                         const unsigned int     first_selected_component = 0)
    : fe_sub(dof_handler.get_fe().base_element(
               dof_handler.get_fe()
                 .component_to_base_index(first_selected_component)
                 .first),
             n_components)
    , phi_m(mapping, fe_sub, update_values | update_gradients)
  {
    for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
      if ((first_selected_component <=
           dof_handler.get_fe().system_to_component_index(i).first) &&
          (dof_handler.get_fe().system_to_component_index(i).first <
           first_selected_component + n_components))
        relevant_dof_indices.push_back(i);
  }

  unsigned int
  data_size() const override
  {
    return n_components * 2;
  }

  const std::vector<unsigned int> &
  get_relevant_dof_indices() const override
  {
    return relevant_dof_indices;
  }

  void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const override
  {
    this->phi_m.reinit(cell, points);
  }

  void
  local_evaluate(const CouplingEvaluationData<dim, Number> &data,
                 const Vector<Number>                      &buffer,
                 const unsigned int                         ptr_q,
                 const unsigned int                         q_stride,
                 Number *all_values_m) const override
  {
    this->phi_m.evaluate(buffer,
                         EvaluationFlags::values | EvaluationFlags::gradients);

    for (const auto q : this->phi_m.quadrature_point_indices())
      {
        // Quadrature point index ('global' index within the rotor-stator
        // interface)
        const unsigned int q_index = ptr_q + q;

        // Normal, value, and gradient referring to the quadrature point
        const auto normal     = data.all_normals[q_index];
        const auto value_m    = this->phi_m.get_value(q);
        const auto gradient_m = contract(this->phi_m.get_gradient(q), normal);

        // Initialize buffer for 'negative' side of the interface (i.e. rotor),
        // where information is evaluated
        BufferRW<Number> buffer_m(all_values_m,
                                  q * 2 * n_components * q_stride);
        // Store values and gradients at the created buffer
        buffer_m.write(value_m);
        buffer_m.write(gradient_m);
      }
  }

  void
  local_integrate(const CouplingEvaluationData<dim, Number> &data,
                  Vector<Number>                            &buffer,
                  const unsigned int                         ptr_q,
                  const unsigned int                         q_stride,
                  Number                                    *all_values_m,
                  Number *all_values_p) const override
  {
    for (const auto q : this->phi_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;
        // Initialize buffer for both 'mortar' and 'non-mortar' sides of the
        // interface
        BufferRW<Number> buffer_m(all_values_m,
                                  q * 2 * n_components * q_stride);
        BufferRW<Number> buffer_p(all_values_p,
                                  q * 2 * n_components * q_stride);
        // Read shape functions values and gradients stored in the buffer
        const auto value_m           = buffer_m.template read<value_type>();
        const auto value_p           = buffer_p.template read<value_type>();
        const auto normal_gradient_m = buffer_m.template read<value_type>();
        const auto normal_gradient_p = buffer_p.template read<value_type>();

        const auto JxW               = data.all_weights[q_index];
        const auto penalty_parameter = data.all_penalty_parameter[q_index];
        const auto normal            = data.all_normals[q_index];

        // The expression for the jump on the mortar interface is
        // jump(u) = u_m * normal_m + u_p * normal_p. Since we are accessing
        // only the value of normal_m, we use a minus sign here because normal_p
        // = - normal_m
        const auto value_jump = outer((value_m - value_p), normal);

        // The expression for the average on the mortar interface is
        // avg(∇u).n = (∇u_m.normal_m + ∇u_p.normal_p) * 0.5. For the same
        // reason above, we include the negative sign here
        const auto gradient_normal_avg =
          (normal_gradient_m - normal_gradient_p) * 0.5;

        // SIPG penalty parameter
        const double sigma = penalty_parameter * data.penalty_factor;

        // - (n avg(∇v), jump(u))
        this->phi_m.submit_gradient(-value_jump * 0.5 * JxW, q);

        // + (jump(v), σ jump(u) - avg(∇u) n)
        this->phi_m.submit_value(
          (contract(value_jump, normal) * sigma - gradient_normal_avg) * JxW,
          q);
      }
    // Multiply previous terms by respective test functions values/gradients
    this->phi_m.test_and_sum(buffer,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);
  }

  /// Finite element that matches the components `n_components` components
  /// starting at component with index `first_selected_component`
  const FESystem<dim> fe_sub;

  /// Interface to the evaluation of mortar coupling interpolated solution
  mutable FEPointIntegrator phi_m;

  /// Relevant dof indices
  std::vector<unsigned int> relevant_dof_indices;
};
