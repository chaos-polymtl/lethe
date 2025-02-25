// SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <fem-dem/fluid_dynamics_vans_matrix_free_operators.h>
#include <fem-dem/void_fraction.h>


template <int dim, typename number>
VANSOperator<dim, number>::VANSOperator() = default;

template <int dim, typename number>
void
VANSOperator<dim, number>::evaluate_non_linear_term_and_calculate_tau(
  const VectorType &newton_step)
{
  NavierStokesOperatorBase<dim, number>::
    evaluate_non_linear_term_and_calculate_tau(newton_step);
}

template <int dim, typename number>
void
VANSOperator<dim, number>::evaluate_void_fraction(
  const VoidFractionBase<dim> &void_fraction_manager)
{
  this->timer.enter_subsection("operator::evaluate_void_fraction");

  typename MatrixFree<dim, number>::AdditionalData additional_data;
  additional_data.mapping_update_flags = (update_values | update_gradients);

  // Matrix_free_void_fraction is reinit with the quadrature of the
  // Navier-Stokes equations and not it's own quadrature rule.
  /*
  matrix_free_void_fraction.reinit(
    (*void_fraction_manager.mapping),
    void_fraction_manager.dof_handler,
    void_fraction_manager.void_fraction_constraints,
    this->quadrature,
    additional_data);
    */

  const unsigned int n_cells = this->matrix_free.n_cell_batches();
  FECellIntegrator   integrator(this->matrix_free);

  void_fraction.reinit(n_cells, integrator.n_q_points);
  void_fraction_gradient.reinit(n_cells, integrator.n_q_points);


  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      integrator.reinit(cell);
      for (const auto q : integrator.quadrature_point_indices())
        {
          const Point<dim, VectorizedArray<number>> point_batch =
            integrator.quadrature_point(q);
          // Temporary work around, extract the void fraction function from the
          // void fraction manager
          void_fraction(cell, q) = evaluate_function<dim, number>(
            (void_fraction_manager.void_fraction_parameters->void_fraction),
            point_batch);
          void_fraction_gradient(cell, q) =
            evaluate_function_gradient<dim, number>(
              (void_fraction_manager.void_fraction_parameters->void_fraction),
              point_batch);
        }
    }

  this->timer.leave_subsection("operator::evaluate_void_fraction");
}

/**
 * The expressions calculated in this cell integral are:
 * (q,∇·δu) + (v, ε ∂t δu) + (v, ε (u·∇)δu) + (v, ε (δu·∇)u) - (∇·v, ε δp) +  ε
 * ν(∇v,∇δu) (Weak form Jacobian), plus three additional terms in the case of
 * SUPG-PSPG stabilization:
 * \+ (ε ∂t δu + ε(u·∇)δu +  ε(δu·∇)u +  ε∇δp -  εν∆δu)τ·∇q (PSPG Jacobian)
 * \+ (ε ∂t δu + ε(u·∇)δu +  ε(δu·∇)u +  ε∇δp -  εν∆δu)τu·∇v (SUPG Jacobian
 * Part 1)
 * \+ (ε ∂t u  + ε(u·∇)u  +  ε∇p -  εν∆u -  εf )τδu·∇v (SUPG Jacobian Part 2),
 * plus two additional terms in the case of full gls stabilization:
 * \+ (ɛ∂t δu +ɛ(u·∇)δu + ɛ(δu·∇)u + ɛ∇δp - νɛ∆δu)τ(−ν∆v) (GLS Jacobian)
 * \+ (∇·δu)τ'(∇·v) (LSIC Jacobian).
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
      const auto tau      = this->stabilization_parameter[cell][q];
      const auto tau_lsic = this->stabilization_parameter_lsic[cell][q];

      // Weak form Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          // ν(∇v,ɛ∇δu)
          gradient_result[i] =
            this->kinematic_viscosity * vf_value * gradient[i];
          // -(∇·v,ɛδp)
          gradient_result[i][i] += -vf_value * value[dim];
          // +(q,∇δu)
          value_result[dim] += gradient[i][i];

          for (unsigned int k = 0; k < dim; ++k)
            {
              // +(v,ɛ(u·∇)δu + ɛ(δu·∇)u)
              value_result[i] +=
                vf_value * gradient[i][k] * previous_values[k] +
                vf_value * previous_gradient[i][k] * value[k];
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
                tau * (-vf_value * this->kinematic_viscosity *
                         hessian_diagonal[i][k] +
                       vf_value * gradient[i][k] * previous_values[k] +
                       vf_value * previous_gradient[i][k] * value[k]);
            }
          // +ɛ(∂t δu)·τ∇q
          if (transient)
            gradient_result[dim][i] +=
              tau * vf_value * (*bdf_coefs)[0] * value[i];
        }
      // (ɛ∇δp)τ·∇q
      gradient_result[dim] += tau * vf_value * gradient[dim];

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
                     this->kinematic_viscosity * hessian_diagonal[i][l]);
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
                     this->kinematic_viscosity *
                       previous_hessian_diagonal[i][l]);
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

      if (this->stabilization ==
          Parameters::Stabilization::NavierStokesStabilization::gls)
        {
          // GLS Jacobian
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  if (this->enable_hessians_jacobian)
                    {
                      for (unsigned int l = 0; l < dim; ++l)
                        {
                          // +(ɛ(u·∇)δu + ɛ(δu·∇)u - νɛ∆δu)τ(−ν∆v)
                          hessian_result[i][k][k] +=
                            tau * -this->kinematic_viscosity * vf_value *
                            (gradient[i][l] * previous_values[l] +
                             previous_gradient[i][l] * value[l] -
                             this->kinematic_viscosity *
                               hessian_diagonal[i][l]);
                        }

                      // +(ɛ∇δp)τ(−ν∆v)
                      hessian_result[i][k][k] += tau *
                                                 -this->kinematic_viscosity *
                                                 (vf_value * gradient[dim][i]);

                      // +(ɛ∂t δu)τ(−ν∆v)
                      if (transient)
                        hessian_result[i][k][k] +=
                          tau * -this->kinematic_viscosity *
                          (vf_value * (*bdf_coefs)[0] * value[i]);
                    }

                  // LSIC term
                  // (∇·δu)τ'(∇·v)
                  gradient_result[i][i] += tau_lsic * gradient[k][k];
                }
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
 * (q, ∇·u) + (v,ɛ∂t u) + (v,ɛ(u·∇)u) - (∇·v,ɛp) + ɛν(∇v,∇u) - (v,ɛf) (Weak
 * form), plus two additional terms in the case of SUPG-PSPG stabilization:
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - ɛf)τ∇·q (PSPG term)
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - ɛf)τu·∇v (SUPG term),
 * plus two additional terms in the case of full gls stabilization:
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - ɛf)τ(−ν∆v) (GLS term)
 * \+ (∇·u)τ'(∇·v) (LSIC term).
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
          const auto tau      = this->stabilization_parameter[cell][q];
          const auto tau_lsic = this->stabilization_parameter_lsic[cell][q];

          // Result value/gradient we will use
          typename FECellIntegrator::value_type    value_result;
          typename FECellIntegrator::gradient_type gradient_result;
          typename FECellIntegrator::hessian_type  hessian_result;

          // Weak form
          for (unsigned int i = 0; i < dim; ++i)
            {
              // νɛ(∇v,∇u)
              gradient_result[i] =
                this->kinematic_viscosity * vf_value * gradient[i];
              // ν(v,∇u∇ɛ)
              value_result[i] =
                this->kinematic_viscosity * vf_gradient * gradient[i];
              // -(∇·v,ɛp)
              gradient_result[i][i] += -vf_value * value[dim];
              // -(v,p∇ɛ)
              value_result[i] += -vf_gradient[i] * value[dim];
              // +(v,-ɛf)
              value_result[i] = -vf_value * source_value[i];

              // +(v,ɛ∂t u)
              if (transient)
                value_result[i] += vf_value * (*bdf_coefs)[0] * value[i] +
                                   previous_time_derivatives[i];


              // +(q,∇·u)
              value_result[dim] += gradient[i][i];

              for (unsigned int k = 0; k < dim; ++k)
                {
                  // +(v,ɛ(u·∇)u)
                  value_result[i] += vf_value * gradient[i][k] * value[k];
                }
            }

          // PSPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  // (-νɛ∆u -ν∇u∇ɛ + ɛ(u·∇)u)·τ∇q
                  gradient_result[dim][i] +=
                    tau * (-vf_value * this->kinematic_viscosity *
                             hessian_diagonal[i][k] -
                           this->kinematic_viscosity * vf_gradient[k] *
                             gradient[i][k] +
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

          // SUPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  for (unsigned int l = 0; l < dim; ++l)
                    {
                      // (-νɛ∆u )τ(u·∇)v
                      gradient_result[i][k] +=
                        -tau * this->kinematic_viscosity * vf_value * value[k] *
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

          if (this->stabilization ==
              Parameters::Stabilization::NavierStokesStabilization::gls)
            {
              // GLS term
              for (unsigned int i = 0; i < dim; ++i)
                {
                  for (unsigned int k = 0; k < dim; ++k)
                    {
                      if (this->enable_hessians_residual)
                        {
                          for (unsigned int l = 0; l < dim; ++l)
                            {
                              // (-νɛ∆u + ɛ(u·∇)u)τ(−ν∆v)
                              hessian_result[i][k][k] +=
                                tau * -this->kinematic_viscosity *
                                (-vf_value * this->kinematic_viscosity *
                                   hessian_diagonal[i][l] +
                                 vf_value * gradient[i][l] * value[l]);
                            }
                          // + ɛ(∇p - f)τ(−ν∆v)
                          hessian_result[i][k][k] +=
                            tau * -this->kinematic_viscosity * vf_value *
                            (gradient[dim][i] - source_value[i]);

                          // + (ɛ∂t u)τ(−ν∆v)
                          if (transient)
                            hessian_result[i][k][k] +=
                              tau * -this->kinematic_viscosity *
                              (vf_value * (*bdf_coefs)[0] * value[i] +
                               previous_time_derivatives[i]);
                        }

                      // LSIC term
                      // (∇·u)τ'(∇·v)
                      gradient_result[i][i] += tau_lsic * gradient[k][k];
                    }
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
