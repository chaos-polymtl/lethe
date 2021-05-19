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
  MMS     = 0,
  Couette = 1,
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
CouetteTopVelocity<dim>::vector_value(const Point<dim> &point,
                                      Vector<double> &  values) const
{
  values[0] = 1;
  values[1] = 0;
  if (dim == 3)
    {
      values[2] = 0;
    }
}

template <int dim>
class ChorinNavierStokes
{
public:
  ChorinNavierStokes(const unsigned int degreeVelocity,
                     const unsigned int degreePressure);
  ~ChorinNavierStokes();
  void
  runMMS();
  void
  runCouette();
  void
  runPoiseulle();
  void
  runCavity();

  Function<dim> *exact_solution;
  Function<dim> *forcing_function;

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
  assemble_init_velocity_eq();
  void
  solve_init_velocity_eq();
  void
  assemble_pressure_eq();
  void
  solve_pressure_eq();
  void
  assemble_new_velocity_eq();
  void
  solve_new_velocity_eq();
  void
  calculateL2Error();
  void
  output_results(const unsigned int cycle) const;

  double             viscosity;
  Triangulation<dim> triangulation;
  double             timestep;
  double             simulation_time;

  bool couette_case;
  bool poiseuille_case;
  bool cavity_case;
  bool MMS_case;

  // FE system for velocity
  FESystem<dim>   fe_velocity;
  DoFHandler<dim> dof_handler_velocity;
  SparsityPattern sparsity_pattern_velocity;

  // FE system for pressure
  FESystem<dim>   fe_pressure;
  DoFHandler<dim> dof_handler_pressure;
  SparsityPattern sparsity_pattern_pressure;

  // Components required for equation 1
  SparseMatrix<double> init_velocity_eq_system_matrix;
  Vector<double>       tentative_velocity;
  Vector<double>       init_velocity_eq_system_rhs;

  // Components required for equation 2
  SparseMatrix<double> pressure_eq_system_matrix;
  Vector<double>       pressure_solution;
  Vector<double>       pressure_eq_system_rhs;

  // Components required for equation 3
  SparseMatrix<double> new_velocity_eq_system_matrix;
  Vector<double>       next_velocity;
  Vector<double>       new_velocity_eq_system_rhs;

  Vector<double> velocity_solution;
};


// Constructor
template <int dim>
ChorinNavierStokes<dim>::ChorinNavierStokes(const unsigned int degreeVelocity,
                                            const unsigned int degreePressure)
  : viscosity(10.0)

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

  //  // Define boundary ids according to test case
  //  if (couette_case)
  //    {
  //      for (auto &face : triangulation.active_face_iterators())
  //        {
  //          // if face is on top of cube, set id to 1
  //          if (std::fabs(face->center()(1) - (1.0)) < 1e-12)
  //            face->set_boundary_id(1);
  //        }
  //    }
  //  else if (poiseulle_case)
  //    {
  //      for (auto &face : triangulation.active_face_iterators())
  //        {
  //          // bottom face id = 0 (automatically)
  //          // right face id = 1
  //          if (std::fabs(face->center()(0) - (1.0)) < 1e-12)
  //            face->set_boundary_id(1);
  //          // top face id = 2
  //          if (std::fabs(face->center()(1) - (1.0)) < 1e-12)
  //            face->set_boundary_id(2);
  //          // left face id = 3
  //          if (std::fabs(face->center()(0) - (-1.0)) < 1e-12)
  //            face->set_boundary_id(3);
  //        }
  //    }

  // Global refinement
  triangulation.refine_global(refinementLevel);
}

template <int dim>
void
ChorinNavierStokes<dim>::refine_grid()
{
  triangulation.refine_global(1);
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
  // Build sparsity patterns
  DynamicSparsityPattern dsp_velocity(dof_handler_velocity.n_dofs(),
                                      dof_handler_velocity.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp_velocity);
  sparsity_pattern_velocity.copy_from(dsp_velocity);

  DynamicSparsityPattern dsp_pressure(dof_handler_pressure.n_dofs(),
                                      dof_handler_pressure.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp_pressure);
  sparsity_pattern_pressure.copy_from(dsp_pressure);

  // Initialise matrices/vectors for each equation
  init_velocity_eq_system_matrix.reinit(sparsity_pattern_velocity);
  tentative_velocity.reinit(dof_handler_velocity.n_dofs());
  init_velocity_eq_system_rhs.reinit(dof_handler_velocity.n_dofs());

  pressure_eq_system_matrix.reinit(sparsity_pattern_pressure);
  pressure_solution.reinit(dof_handler_pressure.n_dofs());
  pressure_eq_system_rhs.reinit(dof_handler_pressure.n_dofs());

  new_velocity_eq_system_matrix.reinit(sparsity_pattern_velocity);
  next_velocity.reinit(dof_handler_velocity.n_dofs());
  new_velocity_eq_system_rhs.reinit(dof_handler_velocity.n_dofs());

  velocity_solution.reinit(dof_handler_velocity.n_dofs());
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_init_velocity_eq()
{
  QGauss<dim>   quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values_velocity(fe_velocity,
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

  std::vector<Tensor<1, dim>> phi_u(velocity_dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi_u(velocity_dofs_per_cell);

  std::vector<Tensor<1, dim>> force(velocity_dofs_per_cell);

  // Iterate over each cell
  for (const auto &cell : dof_handler_velocity.active_cell_iterators())
    {
      // Reset each cell contributions
      fe_values_velocity.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Get function values
      fe_values_velocity[velocities].get_function_values(
        velocity_solution, present_velocity_values);

      // Iterate through quadrature points
      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
        {
          for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
            {
              phi_u[i] = fe_values_velocity[velocities].value(i, q_index);
              grad_phi_u[i] =
                fe_values_velocity[velocities].gradient(i, q_index);
            }

          for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
                {
                  // Build matrix A1 by cell components
                  cell_matrix(i, j) +=
                    ((phi_u[i] * phi_u[j]) +
                     (timestep *
                      scalar_product(phi_u[i],
                                     present_velocity_values[q_index] *
                                       grad_phi_u[j])) +
                     (viscosity * timestep *
                      scalar_product(grad_phi_u[i], grad_phi_u[j]))) *
                    fe_values_velocity.JxW(q_index);
                }

              // Build vector F1 by cell components
              cell_rhs(i) +=
                (phi_u[i] * present_velocity_values[q_index] *
                 fe_values_velocity.JxW(q_index)) // phi_i * u* * dx
                + (timestep * phi_u[i] * force[i] *
                   fe_values_velocity.JxW(q_index)); // k_l * phi_i * f_l * dx
            }
        }

      // Transfer cell components to global matrix/vector
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
        for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
          {
            init_velocity_eq_system_matrix.add(local_dof_indices[i],
                                               local_dof_indices[j],
                                               cell_matrix(i, j));
          }

      for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
        init_velocity_eq_system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  // Add boundary conditions
  std::map<types::global_dof_index, double> boundary_values;
  if (couette_case)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               CouetteTopVelocity<dim>(),
                                               boundary_values);
    }

  else if (poiseuille_case)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
    }

  else if (cavity_case)
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               0,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               1,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               3,
                                               CouetteTopVelocity<dim>(),
                                               boundary_values);
    }

  else
    {
      VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                               0,
                                               Functions::ZeroFunction<dim>(
                                                 dim),
                                               boundary_values);
    }
  MatrixTools::apply_boundary_values(boundary_values,
                                     init_velocity_eq_system_matrix,
                                     tentative_velocity,
                                     init_velocity_eq_system_rhs);
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_init_velocity_eq()
{
  // Solve AU=F directly
  SparseDirectUMFPACK init_velocity_eq_direct;
  init_velocity_eq_direct.initialize(init_velocity_eq_system_matrix);
  init_velocity_eq_direct.vmult(tentative_velocity,
                                init_velocity_eq_system_rhs);

  std::cout << "    Initial Velocity Equation (1) directly solved."
            << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_pressure_eq()
{
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
        tentative_velocity, div_tentative_velocity);

      // Iterate over quadrature points
      for (unsigned int q_pressure = 0; q_pressure < n_q_pressure; ++q_pressure)
        {
          for (const unsigned int i : fe_values_pressure.dof_indices())
            for (const unsigned int j : fe_values_pressure.dof_indices())
              cell_matrix(i, j) +=
                (fe_values_pressure.shape_grad(i,
                                               q_pressure) * // grad phi_i(x_q)
                 fe_values_pressure.shape_grad(j,
                                               q_pressure) * // grad phi_j(x_q)
                 fe_values_pressure.JxW(q_pressure));        // dx

          for (const unsigned int i : fe_values_pressure.dof_indices())
            cell_rhs(i) -=
              (fe_values_pressure.shape_value(i, q_pressure) * // phi_i(x_q)
               div_tentative_velocity[q_pressure] *            // div u*
               fe_values_pressure.JxW(q_pressure)) /
              timestep; // dx / k_l
        }

      // Transfer cell components to global matrix/vector
      cell_p->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values_pressure.dof_indices())
        for (const unsigned int j : fe_values_pressure.dof_indices())
          pressure_eq_system_matrix.add(local_dof_indices[i],
                                        local_dof_indices[j],
                                        cell_matrix(i, j));

      for (const unsigned int i : fe_values_pressure.dof_indices())
        pressure_eq_system_rhs(local_dof_indices[i]) += cell_rhs(i);

      // Add boundary conditions (couette)
      if (couette_case || cavity_case)
        {
          // Target bottom left corner at (-1,-1)
          Point<dim> target_point;
          target_point(0) = -1.0;
          target_point(1) = -1.0;

          std::map<unsigned int, double> boundary_condition;

          for (unsigned int vertex = 0;
               vertex < GeometryInfo<dim>::vertices_per_cell;
               ++vertex)
            {
              if (target_point.distance(cell_p->vertex(vertex)) <
                  1e-3 * cell_p->diameter())
                boundary_condition[cell_p->vertex_dof_index(vertex, 0)] = 0;
            }
        }
    }

  // Add boundary conditions (poiseulle)
  std::map<types::global_dof_index, double> boundary_values;
  if (poiseuille_case)
    {
      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                               0,
                                               Functions::ConstantFunction<dim>(
                                                 1.),
                                               boundary_values);
      MatrixTools::apply_boundary_values(boundary_values,
                                         pressure_eq_system_matrix,
                                         pressure_solution,
                                         pressure_eq_system_rhs);
    }

  //  if (couette_case)
  //    {
  //      VectorTools::interpolate_boundary_values(dof_handler_pressure,
  //                                               0,
  //                                               Functions::ConstantFunction<dim>(
  //                                                 0.),
  //                                               boundary_values);
  //      MatrixTools::apply_boundary_values(boundary_values,
  //                                         pressure_eq_system_matrix,
  //                                         pressure_solution,
  //                                         pressure_eq_system_rhs);
  //    }

  //  else
  //    {
  //      VectorTools::interpolate_boundary_values(dof_handler_pressure,
  //                                               0,
  //                                               Functions::ZeroFunction<dim>(),
  //                                               boundary_values);
  //      MatrixTools::apply_boundary_values(boundary_values,
  //                                         pressure_eq_system_matrix,
  //                                         pressure_solution,
  //                                         pressure_eq_system_rhs);
  //    }
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_pressure_eq()
{
  // Solve AU=F directly
  SparseDirectUMFPACK pressure_eq_direct;
  pressure_eq_direct.initialize(pressure_eq_system_matrix);
  pressure_eq_direct.vmult(pressure_solution, pressure_eq_system_rhs);

  std::cout << "    Pressure Equation (2) directly solved." << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_new_velocity_eq()
{
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
        tentative_velocity, tentative_velocity_values);
      fe_values_pressure.get_function_gradients(pressure_solution,
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
                (timestep * phi_u[i] * pressure_gradients[q_velocity] *
                 fe_values_velocity.JxW(q_velocity)); // phi_i(x_q) * dx
            }
        }

      // Transfer cell components to global matrix/vector
      cell_u->get_dof_indices(local_dof_indices);

      for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
        for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
          {
            new_velocity_eq_system_matrix.add(local_dof_indices[i],
                                              local_dof_indices[j],
                                              cell_matrix(i, j));
          }

      for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
        new_velocity_eq_system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  //  // Add boundary conditions
  //  std::map<types::global_dof_index, double> boundary_values;
  //  if (couette_case)
  //    {
  //      VectorTools::interpolate_boundary_values(dof_handler_velocity,
  //                                               0,
  //                                               Functions::ZeroFunction<dim>(
  //                                                 dim),
  //                                               boundary_values);
  //      VectorTools::interpolate_boundary_values(
  //        dof_handler_velocity,
  //        1,
  //        Functions::ConstantFunction<dim>(1., dim),
  //        boundary_values);
  //    }

  //  else if (poiseulle_case)
  //    {
  //      VectorTools::interpolate_boundary_values(dof_handler_velocity,
  //                                               2,
  //                                               Functions::ZeroFunction<dim>(
  //                                                 dim),
  //                                               boundary_values);
  //      VectorTools::interpolate_boundary_values(dof_handler_velocity,
  //                                               3,
  //                                               Functions::ZeroFunction<dim>(
  //                                                 dim),
  //                                               boundary_values);
  //    }

  //  else
  //    {
  //      VectorTools::interpolate_boundary_values(dof_handler_velocity,
  //                                               0,
  //                                               Functions::ZeroFunction<dim>(
  //                                                 dim),
  //                                               boundary_values);
  //    }

  //  MatrixTools::apply_boundary_values(boundary_values,
  //                                     new_velocity_eq_system_matrix,
  //                                     next_velocity,
  //                                     new_velocity_eq_system_rhs);
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_new_velocity_eq()
{
  // Solve AU=F directly
  SparseDirectUMFPACK new_velocity_eq_direct;
  new_velocity_eq_direct.initialize(new_velocity_eq_system_matrix);
  new_velocity_eq_direct.vmult(next_velocity, new_velocity_eq_system_rhs);

  std::cout << "    New Velocity Equation (3) directly solved." << std::endl;
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

  data_out.add_data_vector(dof_handler_pressure, pressure_solution, "pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  std::vector<std::string> velocity_solution_names(dim, "velocity");
  data_out.add_data_vector(dof_handler_velocity,
                           velocity_solution,
                           velocity_solution_names,
                           data_component_interpretation);

  data_out.build_patches();
  std::ofstream output("output-" + Utilities::int_to_string(cycle, 4) + ".vtu");
  data_out.write_vtu(output);
}

// Find the l2 norm of the error between the finite element solution and the
// exact solution
template <int dim>
void
ChorinNavierStokes<dim>::calculateL2Error()
{}

template <int dim>
void
ChorinNavierStokes<dim>::runMMS()
{
  couette_case    = false;
  poiseuille_case = false;
  MMS_case        = true;
  cavity_case     = false;
}

template <int dim>
void
ChorinNavierStokes<dim>::runCouette()
{
  // Set time parameters for simulation
  timestep              = 0.1;
  const double end_time = 10;

  simulation_time              = 0;
  const unsigned int end_cycle = end_time / timestep;

  couette_case    = true;
  poiseuille_case = false;
  MMS_case        = false;
  cavity_case     = false;

  std::cout << " Test case: Couette" << std::endl;

  // Increase time
  for (unsigned int cycle = 0; cycle < end_cycle; cycle++)
    {
      std::cout << " Time: " << simulation_time << " seconds" << std::endl;

      // For first iteration in time
      if (cycle == 0.0)
        {
          make_cube_grid(4);
          refine_grid();
          setup_dofs();
          initialize_system();
        }

      // At each time step
      assemble_init_velocity_eq();
      solve_init_velocity_eq();
      assemble_pressure_eq();
      solve_pressure_eq();
      assemble_new_velocity_eq();
      solve_new_velocity_eq();
      output_results(cycle);

      velocity_solution = next_velocity;
      simulation_time   = simulation_time + timestep;
    }
}

template <int dim>
void
ChorinNavierStokes<dim>::runPoiseulle()
{
  // Set time parameters for simulation
  timestep              = 0.1;
  const double end_time = 10;

  simulation_time              = 0;
  const unsigned int end_cycle = end_time / timestep;

  couette_case    = false;
  poiseuille_case = true;
  MMS_case        = false;
  cavity_case     = false;

  std::cout << " Test case: Poiseulle" << std::endl;

  // Increase time
  for (unsigned int cycle = 0; cycle < end_cycle; cycle++)
    {
      std::cout << " Time: " << simulation_time << " seconds" << std::endl;
      // For first iteration in time
      if (cycle == 0)
        {
          make_cube_grid(4);
          refine_grid();
          setup_dofs();
          initialize_system();
        }

      // At each time step
      assemble_init_velocity_eq();
      solve_init_velocity_eq();
      assemble_pressure_eq();
      solve_pressure_eq();
      assemble_new_velocity_eq();
      solve_new_velocity_eq();
      output_results(cycle);

      velocity_solution = next_velocity;
      simulation_time   = simulation_time + timestep;
    }
}

template <int dim>
void
ChorinNavierStokes<dim>::runCavity()
{
  // Set time parameters for simulation
  timestep              = 0.1;
  const double end_time = 10;

  simulation_time              = 0;
  const unsigned int end_cycle = end_time / timestep;

  couette_case    = false;
  poiseuille_case = false;
  MMS_case        = false;
  cavity_case     = true;
  std::cout << " Test case: Cavity" << std::endl;

  // Increase time
  for (unsigned int cycle = 0; cycle < end_cycle; cycle++)
    {
      std::cout << " Time: " << simulation_time << " seconds" << std::endl;

      // For first iteration in time
      if (cycle == 0.0)
        {
          make_cube_grid(4);
          refine_grid();
          setup_dofs();
          initialize_system();
        }

      // At each time step
      assemble_init_velocity_eq();
      solve_init_velocity_eq();
      assemble_pressure_eq();
      solve_pressure_eq();
      assemble_new_velocity_eq();
      solve_new_velocity_eq();
      output_results(cycle);

      velocity_solution = next_velocity;
      simulation_time   = simulation_time + timestep;
    }
}

int
main()
{
  try
    {
      ChorinNavierStokes<2> problem_2d(1, 1); // degreeVelocity, degreePressure

      problem_2d.runCouette();
      //      problem_2d.runPoiseulle();
      //            problem_2d.runCavity();
      // problem_2d.runMMS();
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
