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
          for (int k = 0; k < dim; ++k)
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
  const DoFHandler<dim>                            &void_fraction_dof_handler,
  const LinearAlgebra::distributed::Vector<double> &void_fraction_solution)
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
              for (int c = 0; c < dim; ++c)
                void_fraction_gradient[cell][q][c][lane] =
                  cell_void_fraction_gradient[q][c];

              void_fraction[cell][q][lane] = cell_void_fraction[q];
            }
        }
    }

  this->timer.leave_subsection("operator::compute_void_fraction");
}

template <int dim, typename number>
void
VANSOperator<dim, number>::compute_particle_fluid_force(
  const DoFHandler<dim>                            &fp_force_dof_handler,
  const LinearAlgebra::distributed::Vector<double> &fp_force_solution,
  const DoFHandler<dim>                            &fp_drag_dof_handler,
  const LinearAlgebra::distributed::Vector<double> &fp_drag_solution,
  const DoFHandler<dim> &particle_velocity_dof_handler,
  const LinearAlgebra::distributed::Vector<double> &particle_velocity_solution,
  const DoFHandler<dim> &momentum_transfer_coefficient_dof_handler,
  const LinearAlgebra::distributed::Vector<double>
    &momentum_transfer_coefficient_solution)
{
  this->timer.enter_subsection("operator::compute_particle_fluid_forces");

  const unsigned int n_cells = this->matrix_free.n_cell_batches();
  FECellIntegrator   integrator(this->matrix_free);

  particle_fluid_force.reinit(n_cells, integrator.n_q_points);
  particle_fluid_drag.reinit(n_cells, integrator.n_q_points);
  particle_velocity.reinit(n_cells, integrator.n_q_points);


  // We create one FE_values per field that we wish to interpolate. This comes
  // with a significant overhead, but that's life.
  FEValues<dim> fe_values_force(*(this->matrix_free.get_mapping_info().mapping),
                                fp_force_dof_handler.get_fe(),
                                this->matrix_free.get_quadrature(),
                                update_values);

  FEValues<dim> fe_values_drag(*(this->matrix_free.get_mapping_info().mapping),
                               fp_force_dof_handler.get_fe(),
                               this->matrix_free.get_quadrature(),
                               update_values);

  FEValues<dim> fe_values_particle_velocity(
    *(this->matrix_free.get_mapping_info().mapping),
    particle_velocity_dof_handler.get_fe(),
    this->matrix_free.get_quadrature(),
    update_values);

  FEValues<dim> fe_values_momentum_transfer_coefficient(
    *(this->matrix_free.get_mapping_info().mapping),
    momentum_transfer_coefficient_dof_handler.get_fe(),
    this->matrix_free.get_quadrature(),
    update_values);

  std::vector<Tensor<1, dim>> cell_fp_force(
    fe_values_force.n_quadrature_points);
  std::vector<Tensor<1, dim>> cell_fp_drag(fe_values_drag.n_quadrature_points);
  std::vector<Tensor<1, dim>> cell_particle_velocity(
    fe_values_particle_velocity.n_quadrature_points);
  std::vector<double> cell_momentum_transfer_coefficient(
    fe_values_force.n_quadrature_points);

  constexpr FEValuesExtractors::Vector vector_index(0);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (auto lane = 0u;
           lane < this->matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          // Reinit the particle-fluid force
          fe_values_force.reinit(
            this->matrix_free.get_cell_iterator(cell, lane)
              ->as_dof_handler_iterator(fp_force_dof_handler));

          fe_values_force[vector_index].get_function_values(fp_force_solution,
                                                            cell_fp_force);

          // Reinit the drag
          fe_values_drag.reinit(
            this->matrix_free.get_cell_iterator(cell, lane)
              ->as_dof_handler_iterator(fp_drag_dof_handler));

          fe_values_drag[vector_index].get_function_values(fp_drag_solution,
                                                           cell_fp_drag);

          // Reinit the particle velocity
          fe_values_particle_velocity.reinit(
            this->matrix_free.get_cell_iterator(cell, lane)
              ->as_dof_handler_iterator(particle_velocity_dof_handler));

          fe_values_particle_velocity[vector_index].get_function_values(
            particle_velocity_solution, cell_particle_velocity);

          // Reinit the momentum transfer coefficient
          fe_values_momentum_transfer_coefficient.reinit(
            this->matrix_free.get_cell_iterator(cell, lane)
              ->as_dof_handler_iterator(
                momentum_transfer_coefficient_dof_handler));

          fe_values_momentum_transfer_coefficient.get_function_values(
            momentum_transfer_coefficient_solution,
            cell_momentum_transfer_coefficient);

          for (const auto q : fe_values_force.quadrature_point_indices())
            {
              for (int c = 0; c < dim; ++c)
                {
                  // The force applied on the fluid from the particle is (-) the
                  // force applied on the particles by the fluid following
                  // Newton's third law.
                  particle_fluid_force[cell][q][c][lane] = -cell_fp_force[q][c];
                  particle_fluid_drag[cell][q][c][lane]  = -cell_fp_drag[q][c];
                  particle_velocity[cell][q][c][lane] =
                    cell_particle_velocity[q][c];
                }
              momentum_transfer_coefficient[cell][q][lane] =
                cell_momentum_transfer_coefficient[q];
            }
        }
    }

  this->timer.leave_subsection("operator::compute_particle_fluid_forces");
}

/**
 * Calculates the Jacobian of VANS equations. Here, we note  f_pf
 * the forces from the particles to the fluid.
 * The expressions calculated in this cell integral are:
 * (Continuity equation)
 * (q,ε∇·δu) + (q,δu·∇ε) +
 * (VANS equations)
 * \+ (v, ε ∂t δu)  -> Time derivative
 * \+ (v, ε (u·∇)δu) + (v, ε (δu·∇)u) -> Advection
 * \- (∇·v, ε δp) - (v,p∇ε) -> Pressure
 * \+ (v, βu) -> Drag force (this is positive because f_pf = - beta(u-v))
 * \+ ε*ν(∇v,∇δu)  +  ν(v,∇ε·∇δu)   -> Viscous term  ---> There are currently
 * two terms missing here which would be:
 *\+ ε*ν(∇v,∇δu^T)  +  ν(v,∇ε·∇δu^T)
 * plus three additional terms in the case of
 * SUPG-PSPG stabilization:
 * \+ (ε ∂t δu + ε(u·∇)δu +  ε(δu·∇)u +  ε∇δp -  εν∆δu +βδu)τ·∇q (PSPG
 * Jacobian)
 * \+ (ε ∂t δu + ε(u·∇)δu +  ε(δu·∇)u +  ε∇δp -  εν∆δu +βδu)τu·∇v (SUPG
 * Jacobian P.1)
 * \+ (ε ∂t u  + ε(u·∇)u  +  ε∇p -  εν∆u -  εf - f_pf )τδu·∇v (SUPG Jacobian
 * P.2), in the case of additional grad-div stabilization
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
      // Gather void fraction value and gradient
      // We gather the void fraction force since we will use it in the forcing
      // terms (beta and source terms).
      auto vf_value    = this->void_fraction(cell, q);
      auto vf_gradient = this->void_fraction_gradient(cell, q);

      Tensor<1, dim, VectorizedArray<number>> source_value;

      // Gather particle-fluid force, will be zero if they have not been
      // gathered.
      auto pf_force_value = this->particle_fluid_force(cell, q);
      auto pf_drag_value  = this->particle_fluid_drag(cell, q);
      auto p_velocity     = this->particle_velocity(cell, q);

      // Add to source term the particle-fluid force and the drag force
      source_value = pf_force_value + pf_drag_value;

      // Evaluate source term function if enabled
      if (this->forcing_function)
        source_value += vf_value * this->forcing_terms(cell, q);

      // Add to source term the dynamic flow control force (zero if not
      // enabled)
      source_value += vf_value * this->beta_force;

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

      // Calculate norm of the relative velocity and of the drag force and use
      // it to calculate the beta momentum exchange coefficient A tolerance is
      // added (1e-9) to prevent division by 0 and occurence of NaN.
      VectorizedArray<number> relative_velocity_norm_squared = 0.;
      VectorizedArray<number> drag_force_norm_squared        = 0.;

      for (unsigned int i = 0; i < dim; ++i)
        {
          drag_force_norm_squared +=
            Utilities::fixed_power<2>(pf_drag_value[i]);
          relative_velocity_norm_squared +=
            Utilities::fixed_power<2>(previous_values[i] - p_velocity[i]);
        }
      relative_velocity_norm_squared += 1e-20;
      // Since the drag force and the relative velocity are both first rank
      // tensor extracting beta is that directly define. We do it in a
      // projection fashion.
      VectorizedArray<number> beta_momentum_exchange =
        std::sqrt(drag_force_norm_squared / relative_velocity_norm_squared);


      Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
      if (transient)
        previous_time_derivatives =
          this->time_derivatives_previous_solutions(cell, q);

      // Get stabilization parameter
      const auto tau = this->stabilization_parameter[cell][q];

      // Weak form Jacobian
      value_result[dim] = 0;
      for (int i = 0; i < dim; ++i)
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
          // +(v, βu)
          value_result[i] += beta_momentum_exchange * value[i];

          for (int k = 0; k < dim; ++k)
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
      for (int i = 0; i < dim; ++i)
        {
          for (int k = 0; k < dim; ++k)
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

          // +βδu·τ∇q
          gradient_result[dim][i] += tau * beta_momentum_exchange * value[i];
        }
      // (ɛ∇δp)τ·∇q
      gradient_result[dim] += tau * vf_value * gradient[dim];

      // Grad-div stabilization Jacobian
      // (∇·v,γ(ɛ∇·δu+δu·∇ɛ))
      if (cfd_dem_parameters.grad_div == true)
        {
          for (int i = 0; i < dim; ++i)
            {
              for (int k = 0; k < dim; ++k)
                {
                  gradient_result[i][i] +=
                    grad_div_gamma(cell, q) *
                    (vf_value * gradient[k][k] + value[k] * vf_gradient[k]);
                }
            }
        }

      // SUPG Jacobian
      for (int i = 0; i < dim; ++i)
        {
          for (int k = 0; k < dim; ++k)
            {
              // Part 1
              for (int l = 0; l < dim; ++l)
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

              // +(βδu)τ(u·∇)v
              gradient_result[i][k] +=
                beta_momentum_exchange * value[i] * tau * previous_values[k];

              // Part 2
              for (int l = 0; l < dim; ++l)
                {
                  // +(ɛ(u·∇)u - νɛ∆u)τ(δu·∇)v
                  gradient_result[i][k] +=
                    tau * value[k] * vf_value *
                    (previous_gradient[i][l] * previous_values[l] -
                     kinematic_viscosity * previous_hessian_diagonal[i][l]);
                }
              // +(ɛ∇p - ɛf)τ(δu·∇)v
              gradient_result[i][k] +=
                tau * value[k] *
                (vf_value * previous_gradient[dim][i] - source_value[i]);

              // +(ɛ∂t u)τ(δu·∇)v
              if (transient)
                gradient_result[i][k] += tau * value[k] * vf_value *
                                         ((*bdf_coefs)[0] * previous_values[i] +
                                          previous_time_derivatives[i]);

              // +(βu)τ(δu·∇)v
              gradient_result[i][k] +=
                beta_momentum_exchange * previous_values[i] * tau * value[k];
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
 * \+ ɛν(∇v,∇u) + ν(v,∇u∇ɛ) - (v, f_pf) - (v,ɛf) (Weak form of VANS),
 * plus two additional terms in the case of SUPG-PSPG
 * stabilization:
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - f_pf - ɛf)τ∇·q (PSPG term)
 * \+ (ɛ∂t u +ɛ(u·∇)u + ɛ∇p - νɛ∆u - f_pf - ɛf)τu·∇v (SUPG term),
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
          // Gather particle-fluid force and drag force, will be zero if they
          // have not been gathered.
          auto pf_force_value = this->particle_fluid_force(cell, q);
          auto pf_drag_value  = this->particle_fluid_drag(cell, q);
          // Add to source term the particle-fluid force (zero if not enabled)
          // We divide this source by the void fraction value since it is
          // multiplied by the void fraction value within the assembler.
          Tensor<1, dim, VectorizedArray<number>> source_value =
            pf_force_value + pf_drag_value;

          // Gather void fraction value and gradient
          auto vf_value    = this->void_fraction(cell, q);
          auto vf_gradient = this->void_fraction_gradient(cell, q);

          // Evaluate source term function if enabled
          if (this->forcing_function)
            source_value += this->forcing_terms(cell, q) * vf_value;

          // Add to source term the dynamic flow control force (zero if
          // not enabled)
          source_value += this->beta_force * vf_value;

          // Gather the original value/gradient
          typename FECellIntegrator::value_type value = integrator.get_value(q);
          typename FECellIntegrator::gradient_type gradient =
            integrator.get_gradient(q);
          typename FECellIntegrator::gradient_type hessian_diagonal;

          if (this->enable_hessians_residual)
            hessian_diagonal = integrator.get_hessian_diagonal(q);

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
          for (int i = 0; i < dim; ++i)
            {
              // ν(∇v,ɛ∇u)
              gradient_result[i] = kinematic_viscosity * vf_value * gradient[i];

              // -(∇·v,ɛp)
              gradient_result[i][i] += -vf_value * value[dim];

              // -(v,ɛf)
              value_result[i] = -source_value[i];

              // -(v,p∇ɛ)
              value_result[i] += -vf_gradient[i] * value[dim];

              // +(v,ɛ∂t u)
              if (transient)
                value_result[i] += vf_value * ((*bdf_coefs)[0] * value[i] +
                                               previous_time_derivatives[i]);

              // +(q,ɛ∇·u)
              value_result[dim] += vf_value * gradient[i][i];
              // +(q,∇ɛ·u)
              value_result[dim] += value[i] * vf_gradient[i];

              for (int k = 0; k < dim; ++k)
                {
                  // ν(v,∇u∇ɛ)
                  value_result[i] +=
                    kinematic_viscosity * gradient[i][k] * vf_gradient[k];
                  // +(v,ɛ(u·∇)u)
                  value_result[i] += vf_value * gradient[i][k] * value[k];
                }
            }

          // PSPG term
          for (int i = 0; i < dim; ++i)
            {
              for (int k = 0; k < dim; ++k)
                {
                  // (-νɛ∆u + ɛ(u·∇)u)·τ∇q
                  gradient_result[dim][i] +=
                    tau * vf_value *
                    (-kinematic_viscosity * hessian_diagonal[i][k] +
                     gradient[i][k] * value[k]);
                }
              // +(-ɛf)·τ∇q
              gradient_result[dim][i] += tau * (-source_value[i]);

              // +(ɛ∂t u)·τ∇q
              if (transient)
                gradient_result[dim][i] +=
                  tau * vf_value *
                  ((*bdf_coefs)[0] * value[i] + previous_time_derivatives[i]);
            }
          // +ɛ(∇p)τ∇·q
          gradient_result[dim] += tau * vf_value * gradient[dim];

          // Grad-div stabilization
          // (∇·v,γ(ɛ∇·u+u·∇ɛ))
          if (cfd_dem_parameters.grad_div == true)
            {
              for (int i = 0; i < dim; ++i)
                {
                  for (int k = 0; k < dim; ++k)
                    {
                      gradient_result[i][i] +=
                        grad_div_gamma(cell, q) *
                        (vf_value * gradient[k][k] + value[k] * vf_gradient[k]);
                    }
                }
            }

          // SUPG term
          for (int i = 0; i < dim; ++i)
            {
              for (int k = 0; k < dim; ++k)
                {
                  for (int l = 0; l < dim; ++l)
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
                  gradient_result[i][k] +=
                    tau * value[k] *
                    (vf_value * gradient[dim][i] - source_value[i]);

                  // + (ɛ∂t u)τ(u·∇)v
                  if (transient)
                    gradient_result[i][k] += tau * value[k] * vf_value *
                                             ((*bdf_coefs)[0] * value[i] +
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
