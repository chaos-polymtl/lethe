// SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <fem-dem/fluid_dynamics_vans_matrix_free_operators.h>

template <int dim, typename number>
void
VANSOperator<dim, number>::evaluate_non_linear_term_and_calculate_tau(
  const VectorType &newton_step)
{
  // Assert that a correct stabilization method is used
  // Currently the VANS solver only supports pspg_supg
  AssertThrow(this->stabilization ==
                Parameters::Stabilization::NavierStokesStabilization::pspg_supg,
              ExcMessage(
                "PSPG-SUPG stabilization is the only stabilization method"
                " currently supported by the VANS matrix-free solver"));

  NavierStokesOperatorBase<dim, number>::
    evaluate_non_linear_term_and_calculate_tau(newton_step);

  // Evaluate the grad-div stabilization constant
  const unsigned int n_cells = this->matrix_free.n_cell_batches();
  FECellIntegrator   integrator(this->matrix_free);

  grad_div_gamma.reinit(n_cells, integrator.n_q_points);

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      integrator.reinit(cell);
      integrator.read_dof_values_plain(newton_step);

      // Integrator must update the values since the velocity
      // magnitude is used to calculate the grad-div stabilization constant.
      integrator.evaluate(EvaluationFlags::values);
      for (const auto q : integrator.quadrature_point_indices())
        {
          // Get the velocity magnitude to calculate the grad_div constant
          VectorizedArray<number> u_mag_squared = 0;
          for (unsigned int k = 0; k < dim; ++k)
            u_mag_squared +=
              Utilities::fixed_power<2>(integrator.get_value(q)[k]);
          VectorizedArray<number> u = std::sqrt(u_mag_squared);
          grad_div_gamma(cell, q) =
            kinematic_viscosity + cfd_dem_parameters.cstar * u;
        }
    }
}

template <int dim, typename number>
void
VANSOperator<dim, number>::compute_void_fraction(
  const LinearAlgebra::distributed::Vector<double> &void_fraction_solution,
  const DoFHandler<dim>                            &void_fraction_dof_handler)
{
  this->timer.enter_subsection("operator::compute_void_fraction");

  const unsigned int n_cells = this->matrix_free.n_cell_batches();
  FECellIntegrator   integrator(this->matrix_free);

  void_fraction.reinit(n_cells, integrator.n_q_points);
  void_fraction_gradient.reinit(n_cells, integrator.n_q_points);

  FEValues<dim> fe_values(*(this->matrix_free.get_mapping_info().mapping),
                          void_fraction_dof_handler.get_fe(),
                          this->matrix_free.get_quadrature(),
                          update_values | update_gradients);

  std::vector<double>         cell_void_fraction(fe_values.n_quadrature_points);
  std::vector<Tensor<1, dim>> cell_void_fraction_gradient(
    fe_values.n_quadrature_points);


  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (auto lane = 0u;
           lane < this->matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          fe_values.reinit(
            this->matrix_free.get_cell_iterator(cell, lane)
              ->as_dof_handler_iterator(void_fraction_dof_handler));

          fe_values.get_function_values(void_fraction_solution,
                                        cell_void_fraction);
          fe_values.get_function_gradients(void_fraction_solution,
                                           cell_void_fraction_gradient);

          for (const auto q : fe_values.quadrature_point_indices())
            {
              for (unsigned int c = 0; c < dim; ++c)
                void_fraction_gradient[cell][q][c][lane] =
                  cell_void_fraction_gradient[q][c];

              void_fraction[cell][q][lane] = cell_void_fraction[q];
            }
        }
    }

  this->timer.leave_subsection("operator::compute_void_fraction");
}

/**
 * The expressions calculated in this cell integral are:
 * (q,ε∇·δu) + (q,δu·∇ε)  (Continuity equation)
 * \+ (v, ε ∂t δu) + (v, ε (u·∇)δu) + (v, ε (δu·∇)u)
 * \- (∇·v, ε δp) - (v,p∇ε)
 * \+ ε*ν(∇v,∇δu)  +  ν(v,∇ε·∇δu)   (VANS equations),
 * plus three additional terms in the case of SUPG-PSPG stabilization:
 * \+ (ε ∂t δu + ε(u·∇)δu +  ε(δu·∇)u +  ε∇δp -  εν∆δu)τ·∇q (PSPG Jacobian)
 * \+ (ε ∂t δu + ε(u·∇)δu +  ε(δu·∇)u +  ε∇δp -  εν∆δu)τu·∇v (SUPG Jacobian P.1)
 * \+ (ε ∂t u  + ε(u·∇)u  +  ε∇p -  εν∆u -  εf )τδu·∇v (SUPG Jacobian P.2),
 * in the case of additional grad-div stabilization
 * \+ (∇·v,γ(ɛ∇·δu+δu·∇ɛ)) (grad-div term)
 */
template <int dim, typename number>
void
VANSOperator<dim, number>::do_cell_integral_local(
  FECellIntegrator &integrator) const
{
  if (this->enable_hessians_jacobian)
    integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                        EvaluationFlags::hessians);
  else
    integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

  const unsigned int cell = integrator.get_current_cell_index();

  // To identify whether the problem is transient or steady
  bool transient =
    (is_bdf(this->simulation_control->get_assembly_method())) ? true : false;

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  // Vector for BDF coefficients
  const Vector<double> *bdf_coefs;
  if (transient)
    bdf_coefs = &this->simulation_control->get_bdf_coefficients();

  for (const auto q : integrator.quadrature_point_indices())
    {
      Tensor<1, dim, VectorizedArray<number>> source_value;

      // Evaluate source term function if enabled
      if (this->forcing_function)
        source_value = this->forcing_terms(cell, q);

      // Add to source term the dynamic flow control force (zero if not
      // enabled)
      source_value += this->beta_force;

      // Gather the original value/gradient
      typename FECellIntegrator::value_type    value = integrator.get_value(q);
      typename FECellIntegrator::gradient_type gradient =
        integrator.get_gradient(q);
      typename FECellIntegrator::gradient_type hessian_diagonal;

      if (this->enable_hessians_jacobian)
        hessian_diagonal = integrator.get_hessian_diagonal(q);

      // Result value/gradient we will use
      typename FECellIntegrator::value_type    value_result;
      typename FECellIntegrator::gradient_type gradient_result;
      typename FECellIntegrator::hessian_type  hessian_result;

      // Gather previous values of the velocity and the pressure
      auto previous_values   = this->nonlinear_previous_values(cell, q);
      auto previous_gradient = this->nonlinear_previous_gradient(cell, q);
      auto previous_hessian_diagonal =
        this->nonlinear_previous_hessian_diagonal(cell, q);

      // Gather void fraction value and gradient
      auto vf_value    = this->void_fraction(cell, q);
      auto vf_gradient = this->void_fraction_gradient(cell, q);

      Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
      if (transient)
        previous_time_derivatives =
          this->time_derivatives_previous_solutions(cell, q);

      // Get stabilization parameter
      const auto tau = this->stabilization_parameter[cell][q];

      // Weak form Jacobian
      value_result[dim] = 0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          // ν(∇v,ɛ∇δu)
          gradient_result[i] = kinematic_viscosity * vf_value * gradient[i];
          // ν(v,∇ɛ∇δu)
          value_result[i] = kinematic_viscosity * vf_gradient * gradient[i];
          // -(∇·v,ɛδp)
          gradient_result[i][i] += -vf_value * value[dim];
          // -(v,δp∇ɛ)
          value_result[i] += -vf_gradient[i] * value[dim];
          // +(q,ɛ∇·δu)
          value_result[dim] += vf_value * gradient[i][i];
          // +(q,∇ɛ·δu)
          value_result[dim] += vf_gradient[i] * value[i];

          for (unsigned int k = 0; k < dim; ++k)
            {
              // +(v,ɛ(u·∇)δu + ɛ(δu·∇)u)
              value_result[i] +=
                vf_value * (gradient[i][k] * previous_values[k] +
                            previous_gradient[i][k] * value[k]);
            }
          // +(v,ɛ ∂t δu)
          if (transient)
            value_result[i] += vf_value * (*bdf_coefs)[0] * value[i];
        }

      // PSPG Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              // (-νɛ∆δu + ɛ(u·∇)δu + ɛ(δu·∇)u)·τ∇q
              gradient_result[dim][i] +=
                tau * vf_value *
                (-kinematic_viscosity * hessian_diagonal[i][k] +
                 gradient[i][k] * previous_values[k] +
                 previous_gradient[i][k] * value[k]);
            }
          // +ɛ(∂t δu)·τ∇q
          if (transient)
            gradient_result[dim][i] +=
              tau * vf_value * (*bdf_coefs)[0] * value[i];
        }
      // (ɛ∇δp)τ·∇q
      gradient_result[dim] += tau * vf_value * gradient[dim];

      // Grad-div stabilization Jacobian
      // (∇·v,γ(ɛ∇·δu+δu·∇ɛ))
      if (cfd_dem_parameters.grad_div == true)
        {
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  gradient_result[i][i] +=
                    grad_div_gamma(cell, q) *
                    (vf_value * gradient[k][k] + value[k] * vf_gradient[k]);
                }
            }
        }

      // SUPG Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              // Part 1
              for (unsigned int l = 0; l < dim; ++l)
                {
                  // +(ɛ(u·∇)δu + ɛ(δu·∇)u - νɛ∆δu)τ(u·∇)v
                  gradient_result[i][k] +=
                    tau * vf_value * previous_values[k] *
                    (gradient[i][l] * previous_values[l] +
                     vf_value * previous_gradient[i][l] * value[l] -
                     kinematic_viscosity * hessian_diagonal[i][l]);
                }
              // +(ɛ∇δp)τ(u·∇)v
              gradient_result[i][k] +=
                tau * vf_value * previous_values[k] * (gradient[dim][i]);

              // +(ɛ∂t δu)τ(u·∇)v
              if (transient)
                gradient_result[i][k] += tau * previous_values[k] * vf_value *
                                         ((*bdf_coefs)[0] * value[i]);

              // Part 2
              for (unsigned int l = 0; l < dim; ++l)
                {
                  // +(ɛ(u·∇)u - νɛ∆u)τ(δu·∇)v
                  gradient_result[i][k] +=
                    tau * value[k] * vf_value *
                    (previous_gradient[i][l] * previous_values[l] -
                     kinematic_viscosity * previous_hessian_diagonal[i][l]);
                }
              // +(ɛ∇p - ɛf)τ(δu·∇)v
              gradient_result[i][k] +=
                tau * value[k] * vf_value *
                (previous_gradient[dim][i] - source_value[i]);

              // +(ɛ∂t u)τ(δu·∇)v
              if (transient)
                gradient_result[i][k] += tau * value[k] * vf_value *
                                         ((*bdf_coefs)[0] * previous_values[i] +
                                          previous_time_derivatives[i]);
            }
        }

      integrator.submit_gradient(gradient_result, q);
      integrator.submit_value(value_result, q);

      if (this->enable_hessians_jacobian)
        integrator.submit_hessian(hessian_result, q);
    }

  if (this->enable_hessians_jacobian)
    integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients |
                         EvaluationFlags::hessians);
  else
    integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

/**
 * The expressions calculated in this cell integral are:
 * (q, ɛ ∇·u) + (q, u·∇ɛ) (Continuity equation)
 * \+(v,ɛ∂t u) + (v,ɛ(u·∇)u) - (∇·v,ɛp) - (v,p∇ɛ)
 * \+ ɛν(∇v,∇u) + ν(v,∇u∇ɛ) - (v,ɛf) (Weak form of VANS),
 * plus two additional terms in the case of SUPG-PSPG
 * stabilization:
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - ɛf)τ∇·q (PSPG term)
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - ɛf)τu·∇v (SUPG term),
 * With additional grad-div stabilization
 * (∇·v,γ(ɛ∇·u+u·∇ɛ)) (grad-div term)
 */
template <int dim, typename number>
void
VANSOperator<dim, number>::local_evaluate_residual(
  const MatrixFree<dim, number>               &matrix_free,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::pair<unsigned int, unsigned int> &range) const
{
  FECellIntegrator integrator(matrix_free);

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      integrator.reinit(cell);
      integrator.read_dof_values_plain(src);

      if (this->enable_hessians_residual)
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients |
                            EvaluationFlags::hessians);
      else
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients);

      // To identify whether the problem is transient or steady
      bool transient =
        (is_bdf(this->simulation_control->get_assembly_method())) ? true :
                                                                    false;

      // Vector for BDF coefficients
      const Vector<double> *bdf_coefs;
      if (transient)
        bdf_coefs = &this->simulation_control->get_bdf_coefficients();

      for (const auto q : integrator.quadrature_point_indices())
        {
          Tensor<1, dim, VectorizedArray<number>> source_value;

          // Evaluate source term function if enabled
          if (this->forcing_function)
            source_value = this->forcing_terms(cell, q);

          // Add to source term the dynamic flow control force (zero if not
          // enabled)
          source_value += this->beta_force;

          // Gather the original value/gradient
          typename FECellIntegrator::value_type value = integrator.get_value(q);
          typename FECellIntegrator::gradient_type gradient =
            integrator.get_gradient(q);
          typename FECellIntegrator::gradient_type hessian_diagonal;

          if (this->enable_hessians_residual)
            hessian_diagonal = integrator.get_hessian_diagonal(q);

          // Gather void fraction value and gradient
          auto vf_value    = this->void_fraction(cell, q);
          auto vf_gradient = this->void_fraction_gradient(cell, q);

          // Time derivatives of previous solutions
          Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
          if (transient)
            previous_time_derivatives =
              this->time_derivatives_previous_solutions(cell, q);

          // Get stabilization parameter
          const auto tau = this->stabilization_parameter[cell][q];

          // Result value/gradient we will use
          typename FECellIntegrator::value_type    value_result;
          typename FECellIntegrator::gradient_type gradient_result;
          typename FECellIntegrator::hessian_type  hessian_result;

          value_result[dim] = 0;
          // Weak form
          for (unsigned int i = 0; i < dim; ++i)
            {
              // ν(∇v,ɛ∇u)
              gradient_result[i] = kinematic_viscosity * vf_value * gradient[i];

              // -(∇·v,ɛp)
              gradient_result[i][i] += -vf_value * value[dim];

              // -(v,ɛf)
              value_result[i] = -vf_value * source_value[i];

              // -(v,p∇ɛ)
              value_result[i] += -vf_gradient[i] * value[dim];

              // +(v,ɛ∂t u)
              if (transient)
                value_result[i] += vf_value * (*bdf_coefs)[0] * value[i] +
                                   previous_time_derivatives[i];

              // +(q,ɛ∇·u)
              value_result[dim] += vf_value * gradient[i][i];
              // +(q,∇ɛ·u)
              value_result[dim] += value[i] * vf_gradient[i];

              for (unsigned int k = 0; k < dim; ++k)
                {
                  // ν(v,∇u∇ɛ)
                  value_result[i] +=
                    kinematic_viscosity * gradient[i][k] * vf_gradient[k];
                  // +(v,ɛ(u·∇)u)
                  value_result[i] += vf_value * gradient[i][k] * value[k];
                }
            }

          // PSPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  // (-νɛ∆u + ɛ(u·∇)u)·τ∇q
                  gradient_result[dim][i] +=
                    tau * vf_value *
                    (-kinematic_viscosity * hessian_diagonal[i][k] +
                     gradient[i][k] * value[k]);
                }
              // +(-ɛf)·τ∇q
              gradient_result[dim][i] += tau * (-vf_value * source_value[i]);

              // +(ɛ∂t u)·τ∇q
              if (transient)
                gradient_result[dim][i] +=
                  tau * (vf_value * (*bdf_coefs)[0] * value[i] +
                         previous_time_derivatives[i]);
            }
          // +ɛ(∇p)τ∇·q
          gradient_result[dim] += tau * vf_value * gradient[dim];

          // Grad-div stabilization
          // (∇·v,γ(ɛ∇·u+u·∇ɛ))
          if (cfd_dem_parameters.grad_div == true)
            {
              for (unsigned int i = 0; i < dim; ++i)
                {
                  for (unsigned int k = 0; k < dim; ++k)
                    {
                      gradient_result[i][i] +=
                        grad_div_gamma(cell, q) *
                        (vf_value * gradient[k][k] + value[k] * vf_gradient[k]);
                    }
                }
            }

          // SUPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  for (unsigned int l = 0; l < dim; ++l)
                    {
                      // (-νɛ∆u )τ(u·∇)v
                      gradient_result[i][k] += -tau * kinematic_viscosity *
                                               vf_value * value[k] *
                                               hessian_diagonal[i][l];

                      // + (ɛ(u·∇)u)τ(u·∇)v
                      gradient_result[i][k] +=
                        tau * vf_value * value[k] * gradient[i][l] * value[l];
                    }
                  // + (ɛ∇p - ɛf)τ(u·∇)v
                  gradient_result[i][k] += tau * value[k] * vf_value *
                                           (gradient[dim][i] - source_value[i]);

                  // + (ɛ∂t u)τ(u·∇)v
                  if (transient)
                    gradient_result[i][k] +=
                      tau * value[k] *
                      (vf_value * (*bdf_coefs)[0] * value[i] +
                       previous_time_derivatives[i]);
                }
            }

          integrator.submit_gradient(gradient_result, q);
          integrator.submit_value(value_result, q);
          if (this->enable_hessians_residual)
            integrator.submit_hessian(hessian_result, q);
        }

      if (this->enable_hessians_residual)
        integrator.integrate_scatter(EvaluationFlags::values |
                                       EvaluationFlags::gradients |
                                       EvaluationFlags::hessians,
                                     dst);
      else
        integrator.integrate_scatter(EvaluationFlags::values |
                                       EvaluationFlags::gradients,
                                     dst);
    }
}

template class VANSOperator<2, double>;
template class VANSOperator<3, double>;
template class VANSOperator<2, float>;
template class VANSOperator<3, float>;
