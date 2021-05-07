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
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
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

//#include "boundaryconditions.h"
//#include "exactsolutions.h"
//#include "forcingfunctions.h"

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include <deal.II/lac/sparse_direct.h>

// Finally, this is as in previous programs:
using namespace dealii;

enum SimulationCases
{
  MMS           = 0,
  TaylorCouette = 1,
};

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
  run();

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
  assemble_eq1();
  void
  solve_eq1();
  void
  assemble_eq2();
  void
  solve_eq2();
  void
  assemble_eq3();
  void
  solve_eq3();
  void
  calculateL2Error();
  void
  output_results(const unsigned int cycle) const;

  double             viscosity;
  Triangulation<dim> triangulation;
  double             timestep;
  std::vector<types::global_dof_index> dofs_per_block;
  const FEValuesExtractors::Vector velocities;

  // FE system for velocity
  FESystem<dim>      fe_velocity;
  DoFHandler<dim>    dof_handler_velocity;
  BlockSparsityPattern    sparsity_pattern_velocity;

  // FE system for pressure
  FESystem<dim>      fe_pressure;
  DoFHandler<dim>    dof_handler_pressure;
  SparsityPattern    sparsity_pattern_pressure;

  // Components required for equation 1
  BlockSparseMatrix<double> eq1_system_matrix;
  BlockVector<double> tentative_velocity;
  BlockVector<double> eq1_system_rhs;

  // Components required for equation 2
  SparseMatrix<double> eq2_system_matrix;
  Vector<double> pressure_solution;
  Vector<double> eq2_system_rhs;

  // Components required for equation 3
  BlockSparseMatrix<double> eq3_system_matrix;
  BlockVector<double> next_velocity;
  BlockVector<double> eq3_system_rhs;

  BlockVector<double> velocity_solution;
};


// Constructor
template <int dim>
ChorinNavierStokes<dim>::ChorinNavierStokes(const unsigned int degreeVelocity,
                                            const unsigned int degreePressure)
  : viscosity(1)

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
  GridGenerator::hyper_cube(triangulation, -1, 1);
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

  // Establish dofs_per_block
  std::vector<unsigned int> block_component(dim, 0);
  DoFRenumbering::component_wise(dof_handler_velocity, block_component);
  dofs_per_block = DoFTools::count_dofs_per_fe_block(dof_handler_velocity, block_component);

  // Output information
  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Number of velocity degrees of freedom: " << dof_handler_velocity.n_dofs()
            << std::endl
            << "   Number of pressure degrees of freedom: " << dof_handler_pressure.n_dofs()
            << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::initialize_system()
{
  // Build sparsity patterns
  BlockDynamicSparsityPattern dsp_velocity(dofs_per_block, dofs_per_block);
  DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp_velocity);
  sparsity_pattern_velocity.copy_from(dsp_velocity);

  DynamicSparsityPattern dsp_pressure(dof_handler_pressure.n_dofs(), dof_handler_pressure.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp_pressure);
  sparsity_pattern_pressure.copy_from(dsp_pressure);

  // Initialise matrices/vectors for each equation
  eq1_system_matrix.reinit(sparsity_pattern_velocity);
  tentative_velocity.reinit(dofs_per_block);
  eq1_system_rhs.reinit(dofs_per_block);

  eq2_system_matrix.reinit(sparsity_pattern_pressure);
  pressure_solution.reinit(dof_handler_pressure.n_dofs());
  eq2_system_rhs.reinit(dof_handler_pressure.n_dofs());

  eq3_system_matrix.reinit(sparsity_pattern_velocity);
  next_velocity.reinit(dofs_per_block);
  eq3_system_rhs.reinit(dofs_per_block);

  velocity_solution.reinit(dofs_per_block);
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_eq1()
{
  QGauss<dim> quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values_velocity(fe_velocity,
                      quadrature_velocity,
                      update_values | update_quadrature_points |
                       update_gradients | update_JxW_values);
  
  const unsigned int velocity_dofs_per_cell = fe_velocity.n_dofs_per_cell();
  const unsigned int n_q_points             = quadrature_velocity.size();

  const FEValuesExtractors::Vector velocities(0);

  // Initialise cell contribution matrices
  FullMatrix<double> cell_matrix(velocity_dofs_per_cell, velocity_dofs_per_cell);
  Vector<double>     cell_rhs(velocity_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(velocity_dofs_per_cell);

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
      fe_values_velocity[velocities].get_function_values(velocity_solution, present_velocity_values);

      // Iterate through quadrature points
      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
      {
        for (unsigned int i=0; i < velocity_dofs_per_cell; ++i)
            {
              phi_u[i]      = fe_values_velocity[velocities].value(i, q_index);
              grad_phi_u[i] = fe_values_velocity[velocities].gradient(i, q_index);
            }

        for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
              {
                // Build matrix A1 by cell components
                cell_matrix(i, j) += ((phi_u[i] * phi_u[j])
                                      + (timestep * scalar_product(phi_u[i], present_velocity_values[q_index] * grad_phi_u[j]))
                                      + (viscosity * timestep * scalar_product(grad_phi_u[i], grad_phi_u[j])))
                                      * fe_values_velocity.JxW(q_index); 
              }

            // Build vector F1 by cell components
            cell_rhs(i) += (phi_u[i] * present_velocity_values[q_index] * fe_values_velocity.JxW(q_index))       // phi_i * u* * dx
                            + (timestep * phi_u[i] * force[i] * fe_values_velocity.JxW(q_index));       // k_l * phi_i * f_l * dx
          }
      }

      // Transfer cell components to global matrix/vector
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i=0; i < velocity_dofs_per_cell; ++i)
        for (unsigned int j=0; j < velocity_dofs_per_cell; ++j)
          {
            eq1_system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));
          }

      for (unsigned int i=0; i < velocity_dofs_per_cell; ++i)
        eq1_system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  // Add boundary conditions
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                           0,
                                           Functions::ZeroFunction<dim>(dim),
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     eq1_system_matrix,
                                     tentative_velocity,
                                     eq1_system_rhs);
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_eq1()
{
  // Solve AU=F directly
  SparseDirectUMFPACK eq1_direct;
  eq1_direct.initialize(eq1_system_matrix);
  eq1_direct.vmult(tentative_velocity, eq1_system_rhs);

  std::cout << "    Equation 1 directly solved." << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_eq2()
{
  // Velocity FE system
  QGauss<dim> quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values_velocity(fe_velocity,
                      quadrature_velocity,
                      update_values | update_quadrature_points |
                       update_gradients | update_JxW_values);

  // Pressure FE system                     
  QGauss<dim> quadrature_pressure(fe_pressure.degree + 1);
  FEValues<dim> fe_pressure_values(fe_pressure,
                        quadrature_pressure,
                        update_values | update_quadrature_points |
                         update_gradients | update_JxW_values);
  
  const unsigned int dofs_per_cell = fe_pressure.n_dofs_per_cell();
  const unsigned int n_q_pressure  = quadrature_pressure.size();

  // Initialise cell contribution matrices
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<double> div_tentative_velocity(n_q_pressure);  

 // Iterate over each cell
 for (const auto &cell : dof_handler_pressure.active_cell_iterators())
    {
      // Reset each cell contributions
      fe_pressure_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Get function values
      fe_values_velocity[velocities].get_function_divergences(tentative_velocity,
                                                              div_tentative_velocity);

      // Iterate over quadrature points
      for (const unsigned int q_pressure : fe_pressure_values.quadrature_point_indices())
        {
          // Get equivalant velocity quadrature point                                                 >>>>>> LINE OF CODE NOT CORRECT!   
          //const auto q_velocity = fe_values_velocity.quadrature_point(q_pressure);

          for (const unsigned int i : fe_pressure_values.dof_indices())
            for (const unsigned int j : fe_pressure_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_pressure_values.shape_grad(i, q_pressure) * // grad phi_i(x_q)
                 fe_pressure_values.shape_grad(j, q_pressure) * // grad phi_j(x_q)
                 fe_pressure_values.JxW(q_pressure));           // dx

          for (const unsigned int i : fe_pressure_values.dof_indices())
            cell_rhs(i) += (fe_pressure_values.shape_value(i, q_pressure) * // phi_i(x_q)
                            div_tentative_velocity[q_pressure] *            // div u*                  >>>>>> Needs to be q_velocity
                            fe_pressure_values.JxW(q_pressure)) / timestep; // dx / k_l
        }

      // Transfer cell components to global matrix/vector
      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_pressure_values.dof_indices())
        for (const unsigned int j : fe_pressure_values.dof_indices())
          eq2_system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));
      for (const unsigned int i : fe_pressure_values.dof_indices())
        eq2_system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  // Add boundary conditions
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     eq2_system_matrix,
                                     pressure_solution,
                                     eq2_system_rhs);
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_eq2()
{
  // Solve AU=F using a preconditioned Conjugate Gradient algorithm
  SolverControl            eq2_solver_control(1000, 1e-12);
  SolverCG<Vector<double>> eq2_solver(eq2_solver_control);

  // Define Jacobi preconditioner
  PreconditionJacobi<SparseMatrix<double>> preconditioner;
  preconditioner.initialize(eq2_system_matrix, 1.2);

  eq2_solver.solve(eq2_system_matrix, pressure_solution, eq2_system_rhs, preconditioner);
  std::cout << "   " << eq2_solver_control.last_step()
            << " CG iterations needed to obtain convergence in Equation 2." << std::endl;
}

template <int dim>
void
ChorinNavierStokes<dim>::assemble_eq3()
{
  // Velocity FE system
  QGauss<dim> quadrature_velocity(fe_velocity.degree + 1);
  FEValues<dim> fe_values_velocity(fe_velocity,
                      quadrature_velocity,
                      update_values | update_quadrature_points |
                       update_gradients | update_JxW_values);

  // Pressure FE system
  QGauss<dim> quadrature_pressure(fe_pressure.degree + 1);
  FEValues<dim> fe_values_pressure(fe_pressure,
                        quadrature_pressure,
                        update_values | update_gradients | update_JxW_values);
  
  const unsigned int velocity_dofs_per_cell = fe_velocity.n_dofs_per_cell();
  const unsigned int n_q_points_velocity    = quadrature_velocity.size();
  const unsigned int n_q_points_pressure    = quadrature_pressure.size();

  const FEValuesExtractors::Vector velocities(0);

  // Initialise cell contribution matrices
  FullMatrix<double> cell_matrix(velocity_dofs_per_cell, velocity_dofs_per_cell);
  Vector<double>     cell_rhs(velocity_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(velocity_dofs_per_cell);

  std::vector<Tensor<1, dim>> tentative_velocity_values(n_q_points_velocity);
  std::vector<Tensor<1, dim>> pressure_gradients(n_q_points_pressure);

  std::vector<Tensor<1, dim>> phi_u(velocity_dofs_per_cell);

  // Iterate over each cell
  for (const auto &cell : dof_handler_velocity.active_cell_iterators())
    {
      // Reset each cell contributions
      fe_values_velocity.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      // Get function values
      fe_values_velocity[velocities].get_function_values(tentative_velocity, tentative_velocity_values);
      fe_values_pressure.get_function_gradients(pressure_solution, pressure_gradients);

      for (unsigned int q_index = 0; q_index < n_q_points_velocity; ++q_index)
      {
        for (unsigned int i=0; i < velocity_dofs_per_cell; ++i)
            phi_u[i]      = fe_values_velocity[velocities].value(i, q_index);

        for (unsigned int i = 0; i < velocity_dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < velocity_dofs_per_cell; ++j)
                {
                  // Build matrix A3 by cell components
                  cell_matrix(i, j) +=  (phi_u[i] *                          // phi_i(x_q)
                                         phi_u[j] *                          // phi_j(x_q)
                                          fe_values_velocity.JxW(q_index));  // dx
                }

            // Build vector F3 by cell components
            cell_rhs(i) += (phi_u[i] * tentative_velocity_values[q_index] * fe_values_velocity.JxW(q_index))
                            - (timestep * phi_u[i] * pressure_gradients[q_index] * fe_values_velocity.JxW(q_index));       // phi_i(x_q) * dx
          }
      }

      // Transfer cell components to global matrix/vector
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int i=0; i < velocity_dofs_per_cell; ++i)
        for (unsigned int j=0; j < velocity_dofs_per_cell; ++j)
          {
            eq3_system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i, j));
          }

      for (unsigned int i=0; i < velocity_dofs_per_cell; ++i)
        eq3_system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  // Add boundary conditions
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler_velocity,
                                           0,
                                           Functions::ZeroFunction<dim>(dim),
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     eq3_system_matrix,
                                     next_velocity,
                                     eq3_system_rhs);
}

template <int dim>
void
ChorinNavierStokes<dim>::solve_eq3()
{
  // Solve AU=F directly
  SparseDirectUMFPACK eq3_direct;
  eq3_direct.initialize(eq3_system_matrix);
  eq3_direct.vmult(next_velocity, eq3_system_rhs);

  std::cout << "    Equation 3 directly solved." << std::endl;
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
  // Output pressure
  DataOut<dim> data_out_pressure;
  data_out_pressure.attach_dof_handler(dof_handler_pressure);
  data_out_pressure.add_data_vector(pressure_solution, "pressure");
  data_out_pressure.build_patches();
  std::ofstream output_pressure("output_pressure-" + Utilities::int_to_string(cycle, 4) + ".vtu");
  data_out_pressure.write_vtu(output_pressure);
  
  // Output velocity
  DataOut<dim> data_out_velocity;
  data_out_velocity.attach_dof_handler(dof_handler_velocity);
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
  std::vector<std::string> solution_names(dim, "velocity");
  data_out_velocity.add_data_vector(velocity_solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           data_component_interpretation);
  data_out_velocity.build_patches();
  std::ofstream output_velocity("output_velocity-" + Utilities::int_to_string(cycle, 4) + ".vtu");
  data_out_velocity.write_vtu(output_velocity);
}

// Find the l2 norm of the error between the finite element solution and the exact solution
template <int dim>
void
ChorinNavierStokes<dim>::calculateL2Error()
{}

// Run function while building to test
template <int dim>
void
ChorinNavierStokes<dim>::run()
{
  timestep = 0.1;

  make_cube_grid(3);
  refine_grid();
  setup_dofs();
  initialize_system();

  assemble_eq1();
  solve_eq1();
  assemble_eq2();
  solve_eq2();
  assemble_eq3();
  solve_eq3();
  output_results(0);

  velocity_solution = next_velocity;
}

template <int dim>
void
ChorinNavierStokes<dim>::runMMS()
{}


template <int dim>
void
ChorinNavierStokes<dim>::runCouette()
{}

int
main()
{
  try
    {
      ChorinNavierStokes<2> problem_2d(1, 1); // degreeVelocity, degreePressure
      problem_2d.run();
      //problem_2d.runCouette();
      //problem_2d.runMMS();
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