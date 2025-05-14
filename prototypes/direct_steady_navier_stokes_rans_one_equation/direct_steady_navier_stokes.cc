// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "boundaryconditions.h"
#include "exactsolutions.h"
#include "forcingfunctions.h"

// Finally, this is as in previous programs:
using namespace dealii;

enum class SimulationCases
{
  MMS           = 0,
  TaylorCouette = 1,
  BFS           = 2,
};

template <int dim>
class DirectSteadyNavierStokes
{
public:
  DirectSteadyNavierStokes(const unsigned int degreeVelocity,
                           const unsigned int degreePressure,
                           const unsigned int degreeTurbulence);
  ~DirectSteadyNavierStokes();
  void
  runMMS();
  void
  runCouette();
  void
  runBFS();

  Function<dim> *exact_solution;
  Function<dim> *forcing_function;

  std::vector<Tensor<1, dim>>
  get_velocity()
  {
    std::vector<Tensor<1, dim>> velocities;
    QGauss<dim>                 quadrature_formula(degreeIntegration_ + 2);
    FEValues<dim>               fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points);
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        std::vector<Tensor<1, dim>> cell_velocities(
          fe_values.n_quadrature_points);
        fe_values[FEValuesExtractors::Vector(0)].get_function_values(
          present_solution, cell_velocities);
        velocities.insert(velocities.end(),
                          cell_velocities.begin(),
                          cell_velocities.end());
      }
    return velocities;
  }

  std::vector<Tensor<2, dim>>
  get_velocity_gradient()
  {
    std::vector<Tensor<2, dim>> velocity_gradients;
    QGauss<dim>                 quadrature_formula(degreeIntegration_ + 2);
    FEValues<dim>               fe_values(fe,
                            quadrature_formula,
                            update_gradients | update_quadrature_points);
    typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                   endc = dof_handler.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        std::vector<Tensor<2, dim>> cell_gradients(
          fe_values.n_quadrature_points);
        fe_values[FEValuesExtractors::Vector(0)].get_function_gradients(
          present_solution, cell_gradients);
        velocity_gradients.insert(velocity_gradients.end(),
                                  cell_gradients.begin(),
                                  cell_gradients.end());
      }
    return velocity_gradients;
  }

  std::vector<double>
  get_turbulent_k()
  {
    std::vector<double> k_values;
    QGauss<dim>         quadrature_formula(degreeIntegration_);
    FEValues<dim>       fe_values(fe_turbulence,
                            quadrature_formula,
                            update_values | update_quadrature_points);
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler_turbulence.begin_active(),
      endc = dof_handler_turbulence.end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit(cell);
        std::vector<double> cell_k_values(fe_values.n_quadrature_points);
        fe_values.get_function_values(present_solution_turbulence,
                                      cell_k_values);
        k_values.insert(k_values.end(),
                        cell_k_values.begin(),
                        cell_k_values.end());
      }
    return k_values;
  }
  bool add_turbulence_to_ns = false;

  inline double
  calculate_navier_stokes_gls_tau_steady(const double u_mag,
                                         const double kinematic_viscosity,
                                         const double h)
  {
    return 1. / std::sqrt(Utilities::fixed_power<2>(2. * u_mag / h) +
                          9 * Utilities::fixed_power<2>(
                                4 * kinematic_viscosity / (h * h)));
  }

private:
  void
  make_cube_grid(int refinementLevel);
  void
  refine_grid();
  void
  refine_mesh();
  void
  refine_mesh_uniform();
  void
  setup_dofs();
  void
  initialize_system();
  void
  initialize_system_turbulent_model();

  void
  assemble(const bool initial_step, const bool assemble_matrix);
  void
  assemble_system(const bool initial_step);
  void
  assemble_rhs(const bool initial_step);
  void
  solve(bool initial_step);
  void
  calculateL2Error();
  void
  output_results(const unsigned int cycle) const;
  void
  newton_iteration(const double       tolerance,
                   const unsigned int max_iteration,
                   const bool         is_initial_step,
                   const bool         output_result);

  void
  assemble_turbulent_model(const bool initial_step, const bool assemble_matrix);
  void
  assemble_turbulent_model_system(const bool initial_step);
  void
  assemble_turbulent_model_rhs(const bool initial_step);
  void
  solve_turbulent_model(const bool initial_step);
  void
  newton_iteration_turbulent_model(const double       tolerance,
                                   const unsigned int max_iteration,
                                   const bool         is_initial_step,
                                   const bool         output_result);


  std::vector<types::global_dof_index> dofs_per_block;

  double             viscosity_;
  const unsigned int degreeIntegration_;
  Triangulation<dim> triangulation;
  FESystem<dim>      fe;
  DoFHandler<dim>    dof_handler;
  FE_Q<dim>          fe_turbulence;
  DoFHandler<dim>    dof_handler_turbulence;

  AffineConstraints<double> zero_constraints;
  AffineConstraints<double> nonzero_constraints;

  BlockSparsityPattern      sparsity_pattern;
  SparsityPattern           sparsity_pattern_turbulence;
  BlockSparseMatrix<double> system_matrix;
  SparseMatrix<double>      system_matrix_turbulence;

  BlockVector<double> present_solution;
  Vector<double>      present_solution_turbulence;
  BlockVector<double> newton_update;
  Vector<double>      newton_update_turbulence;
  BlockVector<double> system_rhs;
  Vector<double>      system_rhs_turbulence;
  BlockVector<double> evaluation_point;
  Vector<double>      evaluation_point_turbulence;

  const SimulationCases simulationCase_ = SimulationCases::BFS;
  const bool            stabilized_     = false;
  const bool            iterative_      = false;
  std::vector<double>   L2ErrorU_;
  // Define Spalart-Allmaras constants
  // const double C_b1 = 0.1355;
  // const double C_b2 = 0.622;
  // const double k = 0.41;
  // const double sigma = 2./3.;
  // const double C_w1 = C_b1 * std::pow(k, -2.) + (1 + C_b2) / sigma;
  // const double C_w2 = 0.3;
  // const double C_w3 = 2.;
  // const double C_v1 = 7.1;
  const double C_viscosity     = 0.09;
  const double sigma_viscosity = 1.;


  const int initialSize_ = 5;
};


// Constructor
template <int dim>
DirectSteadyNavierStokes<dim>::DirectSteadyNavierStokes(
  const unsigned int degreeVelocity,
  const unsigned int degreePressure,
  const unsigned int degreeTurbulence)
  : viscosity_(1)
  , degreeIntegration_(degreeVelocity)
  , fe(FE_Q<dim>(degreeVelocity), dim, FE_Q<dim>(degreePressure), 1)
  , dof_handler(triangulation)
  , fe_turbulence(degreeTurbulence)
  , dof_handler_turbulence(triangulation)
{}


template <int dim>
DirectSteadyNavierStokes<dim>::~DirectSteadyNavierStokes()
{
  triangulation.clear();
}



template <int dim>
void
DirectSteadyNavierStokes<dim>::make_cube_grid(int refinementLevel)
{
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  triangulation.refine_global(refinementLevel);
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::refine_grid()
{
  triangulation.refine_global(1);
}


template <int dim>
void
DirectSteadyNavierStokes<dim>::setup_dofs()
{
  system_matrix.clear();

  dof_handler.distribute_dofs(fe);
  dof_handler_turbulence.distribute_dofs(fe_turbulence);

  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(dof_handler, block_component);
  dofs_per_block =
    DoFTools::count_dofs_per_fe_block(this->dof_handler, block_component);
  unsigned int dof_u = dofs_per_block[0];
  unsigned int dof_p = dofs_per_block[1];

  FEValuesExtractors::Vector velocities(0);
  {
    nonzero_constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
    if (simulationCase_ == SimulationCases::BFS)
      {
        VectorTools::interpolate_boundary_values(
          dof_handler,
          0,
          dealii::Functions::ZeroFunction<dim>(dim + 1),
          nonzero_constraints,
          fe.component_mask(velocities));

        std::vector<double> inlet_bc = {1., 0.0, 0.0};
        VectorTools::interpolate_boundary_values(
          dof_handler,
          1,
          dealii::Functions::ConstantFunction<dim>(inlet_bc),
          nonzero_constraints,
          fe.component_mask(velocities));
      }
    else
      {
        for (unsigned int id = 0; id < triangulation.get_boundary_ids().size();
             id++)
          {
            VectorTools::interpolate_boundary_values(
              dof_handler,
              id,
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              nonzero_constraints,
              fe.component_mask(velocities));

            if (simulationCase_ == SimulationCases::TaylorCouette)
              {
                VectorTools::interpolate_boundary_values(dof_handler,
                                                         1,
                                                         RotatingWall<dim>(),
                                                         nonzero_constraints,
                                                         fe.component_mask(
                                                           velocities));
              }
          }
      }
    nonzero_constraints.close();

    {
      zero_constraints.clear();
      DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
      if (simulationCase_ == SimulationCases::BFS)
        {
          VectorTools::interpolate_boundary_values(
            dof_handler,
            0,
            dealii::Functions::ZeroFunction<dim>(dim + 1),
            zero_constraints,
            fe.component_mask(velocities));

          std::vector<double> inlet_bc = {1.0, 0.0, 0.0};
          VectorTools::interpolate_boundary_values(
            dof_handler,
            1,
            dealii::Functions::ConstantFunction<dim>(inlet_bc),
            zero_constraints,
            fe.component_mask(velocities));
        }
      else
        {
          for (unsigned int id = 0;
               id < triangulation.get_boundary_ids().size();
               id++)
            {
              VectorTools::interpolate_boundary_values(
                dof_handler,
                id,
                dealii::Functions::ZeroFunction<dim>(dim + 1),
                zero_constraints,
                fe.component_mask(velocities));

              if (simulationCase_ == SimulationCases::TaylorCouette)
                {
                  VectorTools::interpolate_boundary_values(
                    dof_handler,
                    id,
                    dealii::Functions::ZeroFunction<dim>(dim + 1),
                    zero_constraints,
                    fe.component_mask(velocities));
                }
            }
        }
    }
    zero_constraints.close();
    std::cout << "   Number of active cells: " << triangulation.n_active_cells()
              << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << " (" << dof_u << '+' << dof_p << ')' << std::endl;
    std::cout << "   Number of degrees of freedom turbulence: "
              << dof_handler_turbulence.n_dofs() << std::endl;
  }
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::initialize_system()
{
  {
    BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
    sparsity_pattern.copy_from(dsp);
  }
  system_matrix.reinit(sparsity_pattern);
  present_solution.reinit(dofs_per_block);
  newton_update.reinit(dofs_per_block);
  system_rhs.reinit(dofs_per_block);
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::initialize_system_turbulent_model()
{
  {
    DynamicSparsityPattern dsp(dof_handler_turbulence.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler_turbulence, dsp);
    sparsity_pattern_turbulence.copy_from(dsp);
  }
  system_matrix_turbulence.reinit(sparsity_pattern_turbulence);
  present_solution_turbulence.reinit(dof_handler_turbulence.n_dofs());
  newton_update_turbulence.reinit(dof_handler_turbulence.n_dofs());
  system_rhs_turbulence.reinit(dof_handler_turbulence.n_dofs());
  evaluation_point_turbulence.reinit(dof_handler_turbulence.n_dofs());
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::assemble(const bool initial_step,
                                        const bool assemble_matrix)
{
  if (assemble_matrix)
    system_matrix = 0;
  system_rhs = 0;
  QGauss<dim>                      quadrature_formula(degreeIntegration_ + 2);
  FEValues<dim>                    fe_values(fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);
  const unsigned int               dofs_per_cell = fe.dofs_per_cell;
  const unsigned int               n_q_points    = quadrature_formula.size();
  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);
  FullMatrix<double>               local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>                   local_rhs(dofs_per_cell);
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<Tensor<1, dim>>          present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>>          present_velocity_gradients(n_q_points);
  std::vector<Tensor<1, dim>>          velocity_laplacians(n_q_points);
  std::vector<double>                  present_pressure_values(n_q_points);
  std::vector<Tensor<1, dim>>          pressure_gradient(n_q_points);
  std::vector<double>                  div_phi_u(dofs_per_cell);
  std::vector<Tensor<1, dim>>          phi_u(dofs_per_cell);
  std::vector<Tensor<2, dim>>          grad_phi_u(dofs_per_cell);
  std::vector<double>                  phi_p(dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();


  std::vector<double> turbulent_viscosity(n_q_points, 0.0);

  if (add_turbulence_to_ns)
    {
      const auto present_k = get_turbulent_k();
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Calculate turbulent viscosity from the turbulent k
          turbulent_viscosity[q] = C_viscosity * present_k[q] * present_k[q];
        }
    }

  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      local_matrix = 0;
      local_rhs    = 0;
      fe_values[velocities].get_function_values(evaluation_point,
                                                present_velocity_values);
      fe_values[velocities].get_function_gradients(evaluation_point,
                                                   present_velocity_gradients);
      fe_values[velocities].get_function_laplacians(evaluation_point,
                                                    velocity_laplacians);
      fe_values[pressure].get_function_values(evaluation_point,
                                              present_pressure_values);
      fe_values[pressure].get_function_gradients(evaluation_point,
                                                 pressure_gradient);
      forcing_function->vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_force);

      Tensor<1, dim> present_velocity_divergence;

      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const double u_mag =
            std::max(present_velocity_values[q].norm(), 1e-12);
          const double h = cell->measure();
          const double tau =
            calculate_navier_stokes_gls_tau_steady(u_mag, viscosity_, h);
          auto strong_residual =
            present_velocity_gradients[q] * present_velocity_values[q] +
            pressure_gradient[q] -
            viscosity_ * velocity_laplacians[q]; // TODO: Missing rhs force term
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              div_phi_u[k]  = fe_values[velocities].divergence(k, q);
              grad_phi_u[k] = fe_values[velocities].gradient(k, q);
              phi_u[k]      = fe_values[velocities].value(k, q);
              phi_p[k]      = fe_values[pressure].value(k, q);
            }
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              if (assemble_matrix)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (viscosity_ *
                           scalar_product(grad_phi_u[j], grad_phi_u[i]) +
                         present_velocity_gradients[q] * phi_u[j] * phi_u[i] +
                         grad_phi_u[j] * present_velocity_values[q] * phi_u[i] -
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                        fe_values.JxW(q);
                      // Add the contribution of the turbulent viscosity

                      local_matrix(i, j) +=
                        turbulent_viscosity[q] *
                        scalar_product(grad_phi_u[j], grad_phi_u[i]) *
                        fe_values.JxW(q);
                    }
                }
              const unsigned int component_i =
                fe.system_to_component_index(i).first;
              double present_velocity_divergence =
                trace(present_velocity_gradients[q]);
              local_rhs(i) +=
                (-viscosity_ * scalar_product(present_velocity_gradients[q],
                                              grad_phi_u[i]) -
                 present_velocity_gradients[q] * present_velocity_values[q] *
                   phi_u[i] +
                 present_pressure_values[q] * div_phi_u[i] +
                 present_velocity_divergence * phi_p[i]) *
                fe_values.JxW(q);

              // Add the contribution of the turbulent viscosity
              // local_rhs(i) -=
              //   turbulent_viscosity[q] *
              //   scalar_product(present_velocity_gradients[q], grad_phi_u[i])
              //   * fe_values.JxW(q);


              // SUPG term
              // local_rhs(i) += -tau *
              //                 (strong_residual *
              //                  (grad_phi_u[i] * present_velocity_values[q]))
              //                  *
              //                 fe_values.JxW(q);

              local_rhs(i) += fe_values.shape_value(i, q) *
                              rhs_force[q](component_i) * fe_values.JxW(q);
            }
        }


      cell->get_dof_indices(local_dof_indices);
      const AffineConstraints<double> &constraints_used =
        initial_step ? nonzero_constraints : zero_constraints;
      if (assemble_matrix)
        {
          constraints_used.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      system_rhs);
        }
      else
        {
          constraints_used.distribute_local_to_global(local_rhs,
                                                      local_dof_indices,
                                                      system_rhs);
        }
    }
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::assemble_turbulent_model(
  const bool initial_step,
  const bool assemble_matrix)
{
  const double C_mu    = this->C_viscosity;
  const double sigma_k = this->sigma_viscosity;
  if (assemble_matrix)
    system_matrix_turbulence = 0;
  system_rhs_turbulence = 0.;
  QGauss<dim>         quadrature_formula(degreeIntegration_ + 1);
  FEValues<dim>       fe_values(fe_turbulence,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients);
  FEFaceValues<dim>   fe_face_values(fe_turbulence,
                                   QGauss<dim - 1>(degreeIntegration_ + 1),
                                   update_values | update_quadrature_points |
                                     update_JxW_values);
  const unsigned int  n_face_q_points = fe_face_values.n_quadrature_points;
  const unsigned int  dofs_per_cell   = fe_turbulence.dofs_per_cell;
  const unsigned int  n_q_points      = quadrature_formula.size();
  FullMatrix<double>  local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>      local_rhs(dofs_per_cell);
  std::vector<double> local_turbulence_values(n_q_points);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  std::vector<double>                  present_k_values(n_q_points, 1.);
  std::vector<Tensor<1, dim>>          present_k_gradients(n_q_points);
  std::vector<double>                  phi_k(dofs_per_cell);
  std::vector<Tensor<1, dim>>          grad_phi_k(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_turbulence
                                                          .begin_active(),
                                                 endc =
                                                   dof_handler_turbulence.end();

  const std::vector<Tensor<1, dim>> present_velocity =
    get_velocity(); // Get the velocity from the Navier-Stokes solution
  const std::vector<Tensor<2, dim>> present_velocity_gradient =
    get_velocity_gradient(); // Get the velocity gradient from the
                             // Navier-Stokes solution
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      local_matrix = 0.;
      local_rhs    = 0.;
      fe_values.get_function_values(evaluation_point_turbulence,
                                    present_k_values);
      fe_values.get_function_gradients(evaluation_point_turbulence,
                                       present_k_gradients);
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          const auto shear_rate =
            0.5 * (present_velocity_gradient[q] +
                   transpose(present_velocity_gradient[q]));
          const double shear_rate_squared =
            scalar_product(shear_rate, shear_rate);
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              grad_phi_k[k] = fe_values.shape_grad(k, q);
              phi_k[k]      = fe_values.shape_value(k, q);
            }
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              if (assemble_matrix)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (scalar_product(present_velocity[q], grad_phi_k[j]) *
                           phi_k[i] +
                         2. * C_mu / sigma_k * phi_k[j] *
                           present_k_gradients[q] * grad_phi_k[i] +
                         C_mu / sigma_k * present_k_values[q] *
                           present_k_values[q] * grad_phi_k[j] * grad_phi_k[i] -
                         2. * C_mu * phi_k[j] * shear_rate_squared * phi_k[i]) *
                        fe_values.JxW(q);
                    }
                }

              local_rhs(i) +=
                (present_velocity[q] * present_k_gradients[q] * phi_k[i] +
                 C_mu / sigma_k * present_k_values[q] * present_k_values[q] *
                   present_k_gradients[q] * grad_phi_k[i] -
                 C_mu / sigma_k * present_k_values[q] * present_k_values[q] *
                   phi_k[i] * shear_rate_squared) *
                fe_values.JxW(q);
            }
        }
      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix_turbulence.add(local_dof_indices[i],
                                       local_dof_indices[j],
                                       local_matrix(i, j));

      for (const unsigned int i : fe_values.dof_indices())
        system_rhs_turbulence(local_dof_indices[i]) += local_rhs(i);
    }
  std::map<types::global_dof_index, double> boundary_values_turbulence;
  double                                    contour_value = 0.;

  VectorTools::interpolate_boundary_values(
    dof_handler_turbulence,
    0,
    dealii::Functions::ConstantFunction<dim>(0., 1),
    boundary_values_turbulence);

  VectorTools::interpolate_boundary_values(
    dof_handler_turbulence,
    1,
    dealii::Functions::ConstantFunction<dim>(10, 1),
    boundary_values_turbulence);

  //         VectorTools::interpolate_boundary_values(
  // dof_handler_turbulence,
  // 2,
  // dealii::Functions::ConstantFunction<dim>(0., 1),
  // boundary_values_turbulence);


  MatrixTools::apply_boundary_values(boundary_values_turbulence,
                                     system_matrix_turbulence,
                                     present_solution_turbulence,
                                     system_rhs_turbulence);
}



template <int dim>
void
DirectSteadyNavierStokes<dim>::assemble_system(const bool initial_step)
{
  assemble(initial_step, true);
}
template <int dim>
void
DirectSteadyNavierStokes<dim>::assemble_rhs(const bool initial_step)
{
  assemble(initial_step, false);
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::assemble_turbulent_model_system(
  const bool initial_step)
{
  assemble_turbulent_model(initial_step, true);
}
template <int dim>
void
DirectSteadyNavierStokes<dim>::assemble_turbulent_model_rhs(
  const bool initial_step)
{
  assemble_turbulent_model(initial_step, false);
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::solve(const bool initial_step)
{
  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : zero_constraints;
  SparseDirectUMFPACK direct;
  direct.initialize(system_matrix);
  direct.vmult(newton_update, system_rhs);
  constraints_used.distribute(newton_update);
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::solve_turbulent_model(const bool initial_step)
{
  // const AffineConstraints<double> &constraints_used =
  //   initial_step ? nonzero_constraints : zero_constraints;
  SparseDirectUMFPACK direct;
  direct.initialize(system_matrix_turbulence);
  direct.vmult(newton_update_turbulence, system_rhs_turbulence);
  // constraints_used.distribute(newton_update_turbulence);
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::refine_mesh()
{
  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  FEValuesExtractors::Vector velocity(0);
  KellyErrorEstimator<dim>::estimate(
    dof_handler,
    QGauss<dim - 1>(degreeIntegration_ + 1),
    typename std::map<types::boundary_id, const Function<dim, double> *>(),
    present_solution,
    estimated_error_per_cell,
    fe.component_mask(velocity));
  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                  estimated_error_per_cell,
                                                  0.15,
                                                  0.0);
  triangulation.prepare_coarsening_and_refinement();
  SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  triangulation.execute_coarsening_and_refinement();
  setup_dofs();
  BlockVector<double> tmp(dofs_per_block);
  Vector<double>      tmp_turbulence(dof_handler_turbulence.n_dofs());

#if DEAL_II_VERSION_GTE(9, 7, 0)
  solution_transfer.interpolate(tmp);
#else
  solution_transfer.interpolate(tmp, present_solution);
#endif

  nonzero_constraints.distribute(tmp);
  initialize_system();
  initialize_system_turbulent_model();
  present_solution            = tmp;
  present_solution_turbulence = tmp_turbulence;
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::refine_mesh_uniform()
{
  SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler);
  solution_transfer.prepare_for_coarsening_and_refinement(present_solution);
  triangulation.refine_global(1);
  setup_dofs();
  BlockVector<double> tmp(dofs_per_block);
  Vector<double>      tmp_turbulence(dof_handler_turbulence.n_dofs());

#if DEAL_II_VERSION_GTE(9, 7, 0)
  solution_transfer.interpolate(tmp);
#else
  solution_transfer.interpolate(tmp, present_solution);
#endif

  nonzero_constraints.distribute(tmp);
  initialize_system();
  initialize_system_turbulent_model();
  present_solution            = tmp;
  present_solution_turbulence = tmp_turbulence;
}


template <int dim>
void
DirectSteadyNavierStokes<dim>::newton_iteration(
  const double       tolerance,
  const unsigned int max_iteration,
  const bool         is_initial_step,
  const bool /*output_result*/)
{
  double current_res;
  double last_res;
  bool   first_step = is_initial_step;
  {
    unsigned int outer_iteration = 0;
    last_res                     = 1.0;
    current_res                  = 1.0;
    while ((first_step || (current_res > tolerance)) &&
           outer_iteration < max_iteration)
      {
        if (first_step)
          {
            initialize_system();
            evaluation_point = present_solution;
            assemble_system(first_step);
            current_res = system_rhs.l2_norm();
            std::cout << "Newton iteration: " << outer_iteration
                      << "  - Residual:  " << current_res << std::endl;
            solve(first_step);
            present_solution = newton_update;
            nonzero_constraints.distribute(present_solution);
            first_step       = false;
            evaluation_point = present_solution;
            assemble_rhs(first_step);
            current_res = system_rhs.l2_norm();
            last_res    = current_res;
          }
        else
          {
            std::cout << "Newton iteration: " << outer_iteration
                      << "  - Residual:  " << current_res << std::endl;
            evaluation_point = present_solution;
            assemble_system(first_step);
            solve(first_step);
            present_solution = newton_update;
            for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
              {
                evaluation_point = present_solution;
                evaluation_point.add(alpha, newton_update);
                nonzero_constraints.distribute(evaluation_point);
                assemble_rhs(first_step);
                current_res = system_rhs.l2_norm();
                std::cout << "\t\talpha = " << std::setw(6) << alpha
                          << std::setw(0) << " res = " << current_res
                          << std::endl;
                if (current_res < last_res)
                  break;
              }
            {
              present_solution = evaluation_point;
              last_res         = current_res;
            }
          }
        ++outer_iteration;
      }
  }
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::newton_iteration_turbulent_model(
  const double       tolerance,
  const unsigned int max_iteration,
  const bool         is_initial_step,
  const bool /*output_result*/)
{
  double current_res;
  double last_res;
  bool   first_step = is_initial_step;
  {
    unsigned int outer_iteration = 0;
    last_res                     = 1.0;
    current_res                  = 1.0;
    std::cout << "\n\nTURBULENCE NEWTON SOLVER:\n" << std::endl;
    while ((first_step || (current_res > tolerance)) &&
           outer_iteration < max_iteration)
      {
        if (first_step)
          {
            initialize_system_turbulent_model();
            evaluation_point_turbulence = present_solution_turbulence;
            assemble_turbulent_model_system(first_step);
            current_res = system_rhs_turbulence.l2_norm();
            std::cout << "Newton iteration: " << outer_iteration
                      << "  - Residual:  " << current_res << std::endl;
            solve_turbulent_model(first_step);
            first_step                  = false;
            present_solution_turbulence = newton_update_turbulence;
            evaluation_point_turbulence = present_solution_turbulence;
            current_res                 = system_rhs_turbulence.l2_norm();
            last_res                    = current_res;
          }
        else
          {
            evaluation_point_turbulence = present_solution_turbulence;
            assemble_turbulent_model_system(first_step);
            solve_turbulent_model(first_step);
            for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
              {
                evaluation_point_turbulence = present_solution_turbulence;
                evaluation_point_turbulence.add(alpha,
                                                newton_update_turbulence);
                // nonzero_constraints.distribute(evaluation_point_turbulence);
                assemble_turbulent_model_rhs(first_step);
                current_res = system_rhs_turbulence.l2_norm();
                std::cout << "\t\talpha = " << std::setw(6) << alpha
                          << std::setw(0) << " res = " << current_res
                          << std::endl;
                // Check if the residual is decreasing
                // If it is not, we need to reduce alpha
                // If it is, we can break the loop
                if (current_res < last_res)
                  break;
              }
            {
              present_solution_turbulence = evaluation_point_turbulence;
              last_res                    = current_res;
            }
          }
        ++outer_iteration;
        // Print system matrix and rhs
        system_matrix_turbulence.print(std::cout);
        system_rhs_turbulence.print(std::cout);
      }
  }
  std::cout << "Final residual: " << current_res << std::endl;
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::output_results(const unsigned int cycle) const
{
  std::vector<std::string> solution_names(dim, "velocity");
  solution_names.push_back("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(present_solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);
  data_out.build_patches(1);

  std::string filenamesolution = "solution-";
  if (add_turbulence_to_ns)
    filenamesolution += "turbulence-";
  filenamesolution += ('0' + cycle);
  filenamesolution += ".vtk";

  std::cout << "Writing file : " << filenamesolution << std::endl;
  std::ofstream outputSolution(filenamesolution.c_str());

  data_out.write_vtk(outputSolution);
  if (add_turbulence_to_ns)
    {
      DataOut<dim>             data_out_turbulence;
      std::vector<std::string> solution_names_turbulence(1, "turbulence");
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation_turbulence(
          1, DataComponentInterpretation::component_is_scalar);
      data_out_turbulence.attach_dof_handler(dof_handler_turbulence);
      data_out_turbulence.add_data_vector(
        present_solution_turbulence,
        solution_names_turbulence,
        DataOut<dim>::type_dof_data,
        data_component_interpretation_turbulence);

      data_out_turbulence.build_patches(1);

      std::string filenameturbulence = "turbulence-result-";
      filenameturbulence += ('0' + cycle);
      filenameturbulence += ".vtk";
      std::ofstream outputTurbulence(filenameturbulence.c_str());
      data_out_turbulence.write_vtk(outputTurbulence);
      std::cout << "Writing file : " << filenameturbulence << std::endl;
    }
}

// Find the l2 norm of the error between the finite element sol'n and the
// exact sol'n
template <int dim>
void
DirectSteadyNavierStokes<dim>::calculateL2Error()
{
  QGauss<dim>   quadrature_formula(fe.degree + 2);
  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);
  const FEValuesExtractors::Scalar pressure(dim);


  const unsigned int dofs_per_cell =
    fe.dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();
  double             l2errorU   = 0.;

  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim + 1));


  std::vector<Tensor<1, dim>> local_velocity_values(n_q_points);
  std::vector<double>         local_pressure_values(n_q_points);

  double maxPressure = -DBL_MAX;
  // Get the maximal value of the pressure
  for (auto icell = dof_handler.begin_active(); icell != dof_handler.end();
       ++icell)
    {
      fe_values.reinit(icell);
      fe_values[pressure].get_function_values(present_solution,
                                              local_pressure_values);

      for (unsigned int i = 0; i < local_pressure_values.size(); ++i)
        {
          maxPressure = std::max(local_pressure_values[i], maxPressure);
        }
    }

  // loop over elements
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values[velocities].get_function_values(present_solution,
                                                local_velocity_values);
      fe_values[pressure].get_function_values(present_solution,
                                              local_pressure_values);

      // Retrieve the effective "connectivity matrix" for this element
      cell->get_dof_indices(local_dof_indices);

      // Get the exact solution at all gauss points
      exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                        q_exactSol);

      for (unsigned int q = 0; q < n_q_points; q++)
        {
          // Find the values of x and u_h (the finite element solution) at the
          // quadrature points
          double ux_sim   = local_velocity_values[q][0];
          double ux_exact = q_exactSol[q][0];

          double uy_sim   = local_velocity_values[q][1];
          double uy_exact = q_exactSol[q][1];

          l2errorU +=
            (ux_sim - ux_exact) * (ux_sim - ux_exact) * fe_values.JxW(q);
          l2errorU +=
            (uy_sim - uy_exact) * (uy_sim - uy_exact) * fe_values.JxW(q);
        }
    }
  std::cout << "L2Error is : " << std::sqrt(l2errorU) << std::endl;
  L2ErrorU_.push_back(std::sqrt(l2errorU));
}


template <int dim>
void
DirectSteadyNavierStokes<dim>::runMMS()
{
  make_cube_grid(initialSize_);
  exact_solution   = new ExactSolutionMMS<dim>;
  forcing_function = new MMSSineForcingFunction<dim>;
  viscosity_       = 1.;
  setup_dofs();

  //    compute_initial_guess();
  for (unsigned int cycle = 0; cycle < 1; cycle++)
    {
      if (cycle != 0)
        refine_mesh_uniform();
      newton_iteration(1.e-6, 5, true, true);
      output_results(cycle);

      newton_iteration_turbulent_model(1.e-3, 5, true, true);
      add_turbulence_to_ns = true;
      newton_iteration(1.e-6, 5, true, true);
      output_results(cycle);
      calculateL2Error();
    }
  std::ofstream output_file("./L2Error.dat");
  for (unsigned int i = 0; i < L2ErrorU_.size(); ++i)
    {
      output_file << i + initialSize_ << " " << L2ErrorU_[i] << std::endl;
    }
  output_file.close();
}

template <int dim>
void
DirectSteadyNavierStokes<dim>::runBFS()
{
  viscosity_ = 1;
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file("backward-facing-step.msh");

  grid_in.read_msh(input_file);

  forcing_function = new NoForce<dim>;
  exact_solution   = new ExactSolutionTaylorCouette<dim>;
  setup_dofs();

  for (int cycle = 0; cycle < 1; cycle++)
    {
      if (cycle != 0)
        refine_mesh();
      newton_iteration(1.e-3, 50, true, true);
      output_results(cycle);
      newton_iteration_turbulent_model(1.e-3, 5, true, true);
      add_turbulence_to_ns = true;
      newton_iteration(1.e-3, 5, true, true);
      output_results(cycle);
      calculateL2Error();
    }

  std::ofstream output_file("./L2Error.dat");
  for (unsigned int i = 0; i < L2ErrorU_.size(); ++i)
    {
      output_file << i + initialSize_ << " " << L2ErrorU_[i] << std::endl;
    }
  output_file.close();
}


template <int dim>
void
DirectSteadyNavierStokes<dim>::runCouette()
{
  viscosity_ = 10;
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file("taylorcouette.msh");

  grid_in.read_msh(input_file);


  static const SphericalManifold<dim> boundary;

  triangulation.set_all_manifold_ids_on_boundary(0);
  triangulation.set_manifold(0, boundary);

  forcing_function = new NoForce<dim>;
  exact_solution   = new ExactSolutionTaylorCouette<dim>;
  setup_dofs();



  for (int cycle = 0; cycle < 4; cycle++)
    {
      if (cycle != 0)
        refine_mesh();
      newton_iteration(1.e-10, 50, true, true);
      output_results(cycle);
      calculateL2Error();
    }

  std::ofstream output_file("./L2Error.dat");
  for (unsigned int i = 0; i < L2ErrorU_.size(); ++i)
    {
      output_file << i + initialSize_ << " " << L2ErrorU_[i] << std::endl;
    }
  output_file.close();
}

int
main()
{
  try
    {
      DirectSteadyNavierStokes<2> problem_2d(2, 1, 1);
      //        problem_2d.runCouette();
      // problem_2d.runMMS();
      problem_2d.runBFS();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
