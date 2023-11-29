/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */


// @sect3{Include files}

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

//#include "boundaryconditions.h"
//#include "exactsolutions.h"
//#include "forcingfunctions.h"

#include <cmath>
#include <fstream>
#include <iostream>

// Finally, this is as in previous programs:
using namespace dealii;

enum SimulationCases
{
  MMS        = 0,
  Couette    = 1,
  Poiseuille = 2,
  Cavity     = 3
};

template <int dim>
class CouetteTopVelocity : public Function<dim>
{
public:
  CouetteTopVelocity()
    : Function<dim>(dim)
  {}
  virtual void
  vector_value(const Point<dim> &point, Vector<double> &values) const override;
};

template <int dim>
void
CouetteTopVelocity<dim>::vector_value(const Point<dim> & /*point*/,
                                      Vector<double> &values) const
{
  values[0] = 1;
  values[1] = 0;
  if (dim == 3)
    {
      values[2] = 0;
    }
}

template <int dim>
class MMSSineForcingFunction : public Function<dim>
{
public:
  MMSSineForcingFunction()
    : Function<dim>(2){};
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;
};
template <int dim>
void
MMSSineForcingFunction<dim>::vector_value(const Point<dim> &p,
                                          Vector<double> &  values) const
{
  assert(dim == 2);
  const double a = M_PI;

  double x = p[0];
  double y = p[1];
  values(0) =
    (2 * a * a * (-sin(a * x) * sin(a * x) + cos(a * x) * (cos(a * x))) *
       sin(a * y) * cos(a * y) -
     4 * a * a * sin(a * x) * sin(a * x) * sin(a * y) * cos(a * y) - 2.0 * x) *
      (-1.) +
    a * std::pow(sin(a * x), 3.) * std::pow(sin(a * y), 2.) * std::cos(a * x);
  values(1) =
    (2 * a * a * (sin(a * y) * (sin(a * y)) - cos(a * y) * cos(a * y)) *
       sin(a * x) * cos(a * x) +
     4 * a * a * sin(a * x) * sin(a * y) * sin(a * y) * cos(a * x) - 2.0 * y) *
      (-1) +
    a * std::pow(sin(a * x), 2.) * std::pow(sin(a * y), 3.) * std::cos(a * y);
}

template <int dim>
class NoForce : public Function<dim>
{
public:
  NoForce()
    : Function<dim>(2){};
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;
};
template <int dim>
void
NoForce<dim>::vector_value(const Point<dim> & /*p*/,
                           Vector<double> &values) const
{
  assert(dim == 2);
  values(0) = 0.;
  values(1) = 0.;
}

template <int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
  ExactSolutionMMS()
    : Function<dim>(2)
  {}
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;
};
template <int dim>
void
ExactSolutionMMS<dim>::vector_value(const Point<dim> &p,
                                    Vector<double> &  values) const
{
  const double a = M_PI;
  double       x = p[0];
  double       y = p[1];
  values(0)      = sin(a * x) * sin(a * x) * cos(a * y) * sin(a * y);
  values(1)      = -cos(a * x) * sin(a * x) * sin(a * y) * sin(a * y);
  // values(2)      = -2 + x * x + y * y;
}


template <int dim>
class ChorinNavierStokes
{
public:
  ChorinNavierStokes(const unsigned int    degreeVelocity,
                     const unsigned int    degreePressure,
                     const SimulationCases sim_case);
  ~ChorinNavierStokes();
  void
  run();

  Function<dim> *exact_solution;
  Function<dim> *forcing_function;

private:
  void
  make_cube_grid(int refinementLevel);
  void
  refine_mesh();
  void
  refine_mesh_uniform();
  void
  setup_dofs();
  void
  initialize_system();
  void
  set_pressure_reference();
  void
  solve_cfd(const double tolerance);
  void
  assemble_predictor();
  void
  solve_init_velocity_eq();
  void
  assemble_pressure_eq();
  void
  solve_pressure_eq();
  void
  assemble_corrector();
  void
  solve_corrector();
  void
  calculateL2Error();
  void
  output_results(const unsigned int cycle) const;

  SimulationCases    simulation_case;
  double             viscosity;
  Triangulation<dim> triangulation;
  double             k_l;
  double             simulation_time;

  // FE system for velocity
  FESystem<dim>             fe_velocity;
  DoFHandler<dim>           dof_handler_velocity;
  SparsityPattern           sparsity_pattern_velocity;
  AffineConstraints<double> velocity_constraints;


  // FE system for pressure
  FESystem<dim>             fe_pressure;
  DoFHandler<dim>           dof_handler_pressure;
  SparsityPattern           sparsity_pattern_pressure;
  AffineConstraints<double> pressure_constraints;
  AffineConstraints<double> delta_pressure_constraints;



  // Components required for equation 1
  SparseMatrix<double> predictor_eq_system_matrix;
  Vector<double>       predicted_velocity;
  Vector<double>       predictor_eq_system_rhs;

  // Components required for equation 2
  SparseMatrix<double> pressure_eq_system_matrix;
  Vector<double>       pressure;
  Vector<double>       delta_pressure;

  Vector<double> pressure_eq_system_rhs;

  // Components required for equation 3
  SparseMatrix<double> corrector_eq_system_matrix;
  Vector<double>       corrected_velocity;
  Vector<double>       corrector_eq_system_rhs;

  Vector<double> velocity_solution;
  Vector<double> velocity_change;


  bool pressure_initialized;
};


// Constructor
template <int dim>
ChorinNavierStokes<dim>::ChorinNavierStokes(const unsigned int degreeVelocity,
                                            const unsigned int degreePressure,
                                            const SimulationCases sim_case)
  : simulation_case(sim_case)
  , viscosity(1.0)
  // Initialise FE system for velocity
  , fe_velocity(FE_Q<dim>(degreeVelocity), dim)
  , dof_handler_velocity(triangulation)
  // Initialise FE system for pressure
  , fe_pressure(FE_Q<dim>(degreePressure), 1)
  , dof_handler_pressure(triangulation)
{}

template <int dim>
ChorinNavierStokes<dim>::~ChorinNavierStokes()
{
  triangulation.clear();
}

template <int dim>
void
ChorinNavierStokes<dim>::make_cube_grid(int refinementLevel)
{
  // Colorize BCs from the get go.
  GridGenerator::hyper_cube(triangulation, -1, 1, true);

  // Global refinement
  triangulation.refine_global(refinementLevel);
}

template <int dim>
void
ChorinNavierStokes<dim>::setup_dofs()
{
  // Distribute DOFs
  dof_handler_velocity.distribute_dofs(fe_velocity);
  dof_handler_pressure.distribute_dofs(fe_pressure);

  // Output information
  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Number of velocity degrees of freedom: "
            << dof_handler_velocity.n_dofs() << std::endl
            << "   Number of pressure degrees of freedom: "
            << dof_handler_pressure.n_dofs() << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::initialize_system()
{
  // Setup boundary conditions
  if (simulation_case == SimulationCases::Couette)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               CouetteTopVelocity<dim>(),
                                               velocity_constraints);

      set_pressure_reference();
    }

  else if (simulation_case == SimulationCases::Poiseuille)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);

      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               0,
                                               Functions::ConstantFunction<dim>(
                                                 1.),
                                               pressure_constraints);

      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               pressure_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               delta_pressure_constraints);

      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               delta_pressure_constraints);
    }

  else if (simulation_case == SimulationCases::Cavity)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               0,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               1,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               CouetteTopVelocity<dim>(),
                                               velocity_constraints);

      set_pressure_reference();
    }

  else if (simulation_case == SimulationCases::MMS)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               0,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               1,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
      set_pressure_reference();
    }


  else
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               0,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               velocity_constraints);
    }


  velocity_constraints.close();
  pressure_constraints.close();
  delta_pressure_constraints.close();

  std::cout << "Number of pressure constraints       : "
            << pressure_constraints.n_constraints() << std::endl;
  std::cout << "Number of delta pressure constraints : "
            << pressure_constraints.n_constraints() << std::endl;



  // Build sparsity patterns for velocity and pressure
  DynamicSparsityPattern dsp_velocity(dof_handler_velocity.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_velocity,
                                  dsp_velocity,
                                  velocity_constraints,
                                  false);
  sparsity_pattern_velocity.copy_from(dsp_velocity);

  DynamicSparsityPattern dsp_pressure(dof_handler_pressure.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_pressure,
                                  dsp_pressure,
                                  pressure_constraints,
                                  false);
  sparsity_pattern_pressure.copy_from(dsp_pressure);

  // Initialise matrices/vectors for each equation
  predictor_eq_system_matrix.reinit(sparsity_pattern_velocity);
  predicted_velocity.reinit(dof_handler_velocity.n_dofs());
  predictor_eq_system_rhs.reinit(dof_handler_velocity.n_dofs());

  pressure_eq_system_matrix.reinit(sparsity_pattern_pressure);
  pressure.reinit(dof_handler_pressure.n_dofs());
  delta_pressure.reinit(dof_handler_pressure.n_dofs());
  pressure_eq_system_rhs.reinit(dof_handler_pressure.n_dofs());

  corrector_eq_system_matrix.reinit(sparsity_pattern_velocity);
  corrected_velocity.reinit(dof_handler_velocity.n_dofs());
  corrector_eq_system_rhs.reinit(dof_handler_velocity.n_dofs());

  velocity_solution.reinit(dof_handler_velocity.n_dofs());
  velocity_change.reinit(dof_handler_velocity.n_dofs());
}

template <int dim>
void
ChorinNavierStokes<dim>::set_pressure_reference()
{
  // First find candidates for DoF indices to constrain for each velocity
  // component.
  types::global_dof_index pressure_id[1];
  {
    pressure_id[0] = numbers::invalid_dof_index;

    unsigned int n_left_to_find = 1;

    std::vector<types::global_dof_index> local_dof_indices(
      fe_pressure.dofs_per_cell);
    typename DoFHandler<dim>::active_cell_iterator cell;
    for (const auto &cell : dof_handler_pressure.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);

          for (unsigned int i = 0; i < fe_pressure.dofs_per_cell; ++i)
            {
              const types::global_dof_index idx = local_dof_indices[i];

              if (pressure_constraints.can_store_line(idx) &&
                  !pressure_constraints.is_constrained(idx))
                {
                  pressure_id[0] = idx;
                  --n_left_to_find;
                }

              // are we done searching?
              if (n_left_to_find == 0)
                break; // exit inner loop
            }

          if (n_left_to_find == 0)
            break; // exit outer loop
        }
  }

  std::cout << "Adding constraint to the following dof " << pressure_id[0]
            << std::endl;
  pressure_constraints.add_line(pressure_id[0]);
  delta_pressure_constraints.add_line(pressure_id[0]);
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_predictor()
{
  predictor_eq_system_matrix = 0;
  predictor_eq_system_rhs    = 0;

  QGauss<dim>   quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values_velocity(fe_velocity,
                                   quadrature_velocity,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values);

  FEValues<dim> fe_values_pressure(fe_pressure,
                                   quadrature_velocity,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values);

  const unsigned int velocity_dofs_per_cell = fe_velocity.n_dofs_per_cell();
  const unsigned int n_q_points             = quadrature_velocity.size();

  const FEValuesExtractors::Vector velocities(0);

  // Initialise cell contribution matrices
  FullMatrix<double> cell_matrix(velocity_dofs_per_cell,
                                 velocity_dofs_per_cell);
  Vector<double>     cell_rhs(velocity_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(
    velocity_dofs_per_cell);

  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
  std::vector<Tensor<1, dim>> present_pressure_gradients(n_q_points);


  std::vector<Tensor<1, dim>> phi_u(velocity_dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(velocity_dofs_per_cell);

  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim));
  Tensor<1, dim>              force;


  // Set up iterators over both DoF_handlers
  typename DoFHandler<dim>::active_cell_iterator cell_u =
    dof_handler_velocity.begin_active();
  const auto endc = dof_handler_velocity.end();

  typename DoFHandler<dim>::active_cell_iterator cell_p =
    dof_handler_pressure.begin_active();

  // Iterate over each cell
  for (; cell_u != endc; ++cell_u, ++cell_p)
    {
      // Reset each cell contributions
      fe_values_velocity.reinit(cell_u);
      fe_values_pressure.reinit(cell_p);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Get values, shapes, gradients, etc.
      fe_values_velocity[velocities].get_function_values(
        velocity_solution, present_velocity_values);
      fe_values_pressure.get_function_gradients(pressure,
                                                present_pressure_gradients);

      // Force
      forcing_function->vector_value_list(
        fe_values_velocity.get_quadrature_points(), rhs_force);

      // Iterate through quadrature points
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
            {
              phi_u[i]      = fe_values_velocity[velocities].value(i, q);
              grad_phi_u[i] = fe_values_velocity[velocities].gradient(i, q);
            }
          // Establish the force vector
          for (int i = 0; i < dim; ++i)
            {
              force[i] = rhs_force[q](i);
            }

          for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
                {
                  // Build matrix A1 by cell components
                  cell_matrix(i, j) +=
                    ((phi_u[i] * phi_u[j]) / k_l +
                     +viscosity * scalar_product(grad_phi_u[i], grad_phi_u[j]) +
                     grad_phi_u[j] * present_velocity_values[q] * phi_u[i]) *
                    fe_values_velocity.JxW(q);
                }

              // Build vector F1 by cell components
              cell_rhs(i) +=
                (1. / k_l * phi_u[i] * present_velocity_values[q] -
                 phi_u[i] * present_pressure_gradients[q] + phi_u[i] * force) *
                fe_values_velocity.JxW(q); // k_l * phi_i * f_l * dx
            }
        }

      // Transfer cell components to global matrix/vector
      cell_u->get_dof_indices(local_dof_indices);

      velocity_constraints.distribute_local_to_global(
        cell_matrix,
        cell_rhs,
        local_dof_indices,
        predictor_eq_system_matrix,
        predictor_eq_system_rhs);
    }
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_init_velocity_eq()
{
  SparseDirectUMFPACK init_velocity_eq_direct;
  init_velocity_eq_direct.initialize(predictor_eq_system_matrix);
  init_velocity_eq_direct.vmult(predicted_velocity, predictor_eq_system_rhs);
  velocity_constraints.distribute(predicted_velocity);

  std::cout << "    Predictor Equation (1) solved using direct solver."
            << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_pressure_eq()
{
  pressure_eq_system_matrix = 0;
  pressure_eq_system_rhs    = 0;

  // Pressure FE system
  QGauss<dim>   quadrature_pressure(fe_pressure.degree + 1);
  FEValues<dim> fe_values_pressure(fe_pressure,
                                   quadrature_pressure,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values);

  // Velocity FE system
  FEValues<dim> fe_values_velocity(fe_velocity,
                                   quadrature_pressure,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values);

  const FEValuesExtractors::Vector velocities(0);

  const unsigned int pressure_dofs_per_cell = fe_pressure.n_dofs_per_cell();
  const unsigned int n_q_pressure           = quadrature_pressure.size();

  // Initialise cell contribution matrices
  FullMatrix<double> cell_matrix(pressure_dofs_per_cell,
                                 pressure_dofs_per_cell);
  Vector<double>     cell_rhs(pressure_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(
    pressure_dofs_per_cell);

  std::vector<double> div_tentative_velocity(n_q_pressure);

  // Set up iterators over both DoF_handlers
  typename DoFHandler<dim>::active_cell_iterator cell_u =
    dof_handler_velocity.begin_active();

  typename DoFHandler<dim>::active_cell_iterator cell_p =
    dof_handler_pressure.begin_active();
  const auto endc = dof_handler_pressure.end();

  // Iterate over each cell
  for (; cell_p != endc; ++cell_p, ++cell_u)
    {
      // Reset each cell contributions
      fe_values_pressure.reinit(cell_p);
      fe_values_velocity.reinit(cell_u);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Get function values
      fe_values_velocity[velocities].get_function_divergences(
        predicted_velocity, div_tentative_velocity);

      // Iterate over quadrature points
      for (unsigned int q = 0; q < n_q_pressure; ++q)
        {
          for (const unsigned int i : fe_values_pressure.dof_indices())
            for (const unsigned int j : fe_values_pressure.dof_indices())
              cell_matrix(i, j) +=
                (fe_values_pressure.shape_grad(i,
                                               q) * // grad phi_i(x_q)
                 fe_values_pressure.shape_grad(j,
                                               q) * // grad phi_j(x_q)
                 fe_values_pressure.JxW(q));        // dx

          for (const unsigned int i : fe_values_pressure.dof_indices())
            cell_rhs(i) -= (fe_values_pressure.shape_value(i, q) *
                            // phi_i(x_q)
                            div_tentative_velocity[q] * // div u*
                            fe_values_pressure.JxW(q)) /
                           k_l; // dx / k_l
        }

      // Transfer cell components to global matrix/vector
      cell_p->get_dof_indices(local_dof_indices);

      if (pressure_initialized)
        {
          delta_pressure_constraints.distribute_local_to_global(
            cell_matrix,
            cell_rhs,
            local_dof_indices,
            pressure_eq_system_matrix,
            pressure_eq_system_rhs);
        }
      else
        {
          pressure_constraints.distribute_local_to_global(
            cell_matrix,
            cell_rhs,
            local_dof_indices,
            pressure_eq_system_matrix,
            pressure_eq_system_rhs);
        }
    }
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_pressure_eq()
{
  SparseDirectUMFPACK pressure_eq_direct;
  pressure_eq_direct.initialize(pressure_eq_system_matrix);
  pressure_eq_direct.vmult(delta_pressure, pressure_eq_system_rhs);

  if (pressure_initialized)
    {
      delta_pressure_constraints.distribute(delta_pressure);
    }
  else
    {
      pressure_constraints.distribute(delta_pressure);
      pressure_initialized = true;
    }
  // delta_pressure *= 0.3;
  pressure += delta_pressure;


  std::cout << "    Pressure  Equation (2) solved using direct solver."
            << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_corrector()
{
  corrector_eq_system_matrix = 0;
  corrector_eq_system_rhs    = 0;

  // Velocity FE system
  QGauss<dim>   quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values_velocity(fe_velocity,
                                   quadrature_velocity,
                                   update_values | update_quadrature_points |
                                     update_gradients | update_JxW_values);

  // Pressure FE system
  FEValues<dim> fe_values_pressure(fe_pressure,
                                   quadrature_velocity,
                                   update_values | update_gradients |
                                     update_JxW_values);

  const unsigned int velocity_dofs_per_cell = fe_velocity.n_dofs_per_cell();
  const unsigned int n_q_points_velocity    = quadrature_velocity.size();

  const FEValuesExtractors::Vector velocities(0);

  // Initialise cell contribution matrices
  FullMatrix<double> cell_matrix(velocity_dofs_per_cell,
                                 velocity_dofs_per_cell);
  Vector<double>     cell_rhs(velocity_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(
    velocity_dofs_per_cell);

  std::vector<Tensor<1, dim>> tentative_velocity_values(n_q_points_velocity);
  std::vector<Tensor<1, dim>> pressure_gradients(n_q_points_velocity);

  std::vector<Tensor<1, dim>> phi_u(velocity_dofs_per_cell);

  // Set up iterators over both DoF_handlers
  typename DoFHandler<dim>::active_cell_iterator cell_u =
    dof_handler_velocity.begin_active();

  typename DoFHandler<dim>::active_cell_iterator cell_p =
    dof_handler_pressure.begin_active();
  const auto endc = dof_handler_pressure.end();

  // Iterate over each cell
  for (; cell_p != endc; ++cell_p, ++cell_u)
    {
      // Reset each cell contributions
      fe_values_velocity.reinit(cell_u);
      fe_values_pressure.reinit(cell_p);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Get function values
      fe_values_velocity[velocities].get_function_values(
        predicted_velocity, tentative_velocity_values);
      fe_values_pressure.get_function_gradients(delta_pressure,
                                                pressure_gradients);

      for (unsigned int q_velocity = 0; q_velocity < n_q_points_velocity;
           ++q_velocity)
        {
          for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
            phi_u[i] = fe_values_velocity[velocities].value(i, q_velocity);

          for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
                {
                  // Build matrix A3 by cell components
                  cell_matrix(i, j) +=
                    (phi_u[i] *                           // phi_i(x_q)
                     phi_u[j] *                           // phi_j(x_q)
                     fe_values_velocity.JxW(q_velocity)); // dx
                }

              // Build vector F3 by cell components
              cell_rhs(i) +=
                (phi_u[i] * tentative_velocity_values[q_velocity] *
                 fe_values_velocity.JxW(q_velocity)) -
                (k_l * phi_u[i] * pressure_gradients[q_velocity] *
                 fe_values_velocity.JxW(q_velocity)); // phi_i(x_q) * dx
            }
        }

      // Transfer cell components to global matrix/vector
      cell_u->get_dof_indices(local_dof_indices);

      velocity_constraints.distribute_local_to_global(
        cell_matrix,
        cell_rhs,
        local_dof_indices,
        corrector_eq_system_matrix,
        corrector_eq_system_rhs);
    }
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_corrector()
{
  // Solve AU=F directly
  SparseDirectUMFPACK new_velocity_eq_direct;
  new_velocity_eq_direct.initialize(corrector_eq_system_matrix);
  new_velocity_eq_direct.vmult(corrected_velocity, corrector_eq_system_rhs);
  velocity_constraints.distribute(corrected_velocity);


  std::cout << "    Velocity Correction Equation (3) directly solved."
            << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::refine_mesh()
{}

template <int dim>
void
ChorinNavierStokes<dim>::refine_mesh_uniform()
{}

template <int dim>
void
ChorinNavierStokes<dim>::output_results(const unsigned int cycle) const
{
  // Output
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler_pressure, pressure, "pressure");
  data_out.add_data_vector(dof_handler_pressure,
                           delta_pressure,
                           "delta_pressure");


  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  {
    std::vector<std::string> velocity_solution_names(dim, "velocity");
    data_out.add_data_vector(dof_handler_velocity,
                             velocity_solution,
                             velocity_solution_names,
                             data_component_interpretation);
  }

  {
    std::vector<std::string> velocity_solution_names(dim, "predicted_velocity");
    data_out.add_data_vector(dof_handler_velocity,
                             predicted_velocity,
                             velocity_solution_names,
                             data_component_interpretation);
  }

  {
    std::vector<std::string> velocity_solution_names(dim, "corrected_velocity");
    data_out.add_data_vector(dof_handler_velocity,
                             corrected_velocity,
                             velocity_solution_names,
                             data_component_interpretation);
  }



  data_out.build_patches();
  std::ofstream output("output-" + Utilities::int_to_string(cycle, 4) + ".vtu");
  data_out.write_vtu(output);
}

// Find the l2 norm of the error between the finite element solution and the
// exact solution
template <int dim>
void
ChorinNavierStokes<dim>::calculateL2Error()
{
  QGauss<dim>   quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values(fe_velocity,
                          quadrature_velocity,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  double l2errorU = 0;

  const FEValuesExtractors::Vector velocities(0);
  const unsigned int velocity_dofs_per_cell = fe_velocity.n_dofs_per_cell();
  const unsigned int n_q_points             = quadrature_velocity.size();
  std::vector<types::global_dof_index> local_dof_indices(
    velocity_dofs_per_cell);

  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
  std::vector<Vector<double>> q_exactSol(n_q_points, Vector<double>(dim));



  // loop over elements
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_velocity
                                                          .begin_active(),
                                                 endc =
                                                   dof_handler_velocity.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      fe_values[velocities].get_function_values(velocity_solution,
                                                present_velocity_values);

      // Retrieve the effective "connectivity matrix" for this element
      cell->get_dof_indices(local_dof_indices);

      // Get the exact solution at all gauss points
      exact_solution->vector_value_list(fe_values.get_quadrature_points(),
                                        q_exactSol);

      for (unsigned int q = 0; q < n_q_points; q++)
        {
          // Find the values of x and u_h (the finite element solution) at the
          // quadrature points
          double ux_sim   = present_velocity_values[q][0];
          double ux_exact = q_exactSol[q][0];

          double uy_sim   = present_velocity_values[q][1];
          double uy_exact = q_exactSol[q][1];

          l2errorU +=
            (ux_sim - ux_exact) * (ux_sim - ux_exact) * fe_values.JxW(q);
          l2errorU +=
            (uy_sim - uy_exact) * (uy_sim - uy_exact) * fe_values.JxW(q);
        }
    }
  std::cout << "L2Error is : " << std::sqrt(l2errorU) << std::endl;
}



template <int dim>
void
ChorinNavierStokes<dim>::solve_cfd(const double tolerance)
{
  double       residual = 1;
  unsigned int cycle    = 0;
  while (residual > tolerance)
    {
      // At each time step
      std::cout << " Time: " << simulation_time << " seconds" << std::endl;

      assemble_predictor();
      solve_init_velocity_eq();
      assemble_pressure_eq();
      solve_pressure_eq();
      assemble_corrector();
      solve_corrector();
      output_results(cycle);

      velocity_change = corrected_velocity;
      velocity_change -= velocity_solution;
      double residual_velocity = velocity_change.l2_norm();
      double residual_pressure = delta_pressure.l2_norm();

      std::cout << "    Residual velocity : " << residual_velocity << std::endl;
      std::cout << "    Residual pressure : " << residual_pressure << std::endl;
      residual = std::max(residual_pressure, residual_velocity);

      velocity_solution = corrected_velocity;
      simulation_time   = simulation_time + k_l;
      cycle++;
    }
}

template <int dim>
void
ChorinNavierStokes<dim>::run()
{
  if (simulation_case == SimulationCases::MMS)
    {
      forcing_function = new MMSSineForcingFunction<dim>;
      exact_solution   = new ExactSolutionMMS<dim>;
    }
  else
    forcing_function = new NoForce<dim>;


  viscosity = 1;
  // Set time parameters for simulation
  k_l                  = 0.01;
  double tolerance     = 1e-3;
  pressure_initialized = false;

  simulation_time = 0;

  make_cube_grid(5);
  setup_dofs();
  initialize_system();
  solve_cfd(tolerance);
  calculateL2Error();
}

int
main()
{
  try
    {
      ChorinNavierStokes<2> problem_2d(
        2, 1, SimulationCases::MMS); // degreeVelocity, degreePressure

      problem_2d.run();
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
