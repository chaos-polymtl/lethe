/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2016 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Heidelberg, 1999
 */


// @sect3{Include files}

// The first few (many?) include files have already been used in the previous
// example, so we will not explain their meaning here again.
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;
enum simCase
{
  TaylorCouette,
  MMS
};

template <int dim>
class DGHeat
{
public:
  DGHeat(simCase scase, int refinementLevel);
  void
  run();
  double
  getL2Error()
  {
    return L2Error_;
  }

private:
  void
  make_grid();
  void
  make_cube_grid();
  void
  make_ring_grid();
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results() const;
  void
  calculateL2Error();

  Triangulation<dim> triangulation;
  FE_Q<dim>          fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
  Point<dim>     center;

  simCase simulationCase_;
  int     refinementLevel_;
  double  L2Error_;
};



template <int dim>
class RightHandSideMMS : public Function<dim>
{
public:
  RightHandSideMMS()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
};

template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  BoundaryValues()
    : Function<dim>()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;
};

template <int dim>
double
RightHandSideMMS<dim>::value(const Point<dim> &p,
                             const unsigned int /*component*/) const
{
  double return_value = 0.0;
  double x            = p(0);
  double y            = p(1);

  return_value = -2. * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);

  return return_value;
}

// As boundary values, we choose $x^2+y^2$ in 2D, and $x^2+y^2+z^2$ in 3D. This
// happens to be equal to the square of the vector from the origin to the
// point at which we would like to evaluate the function, irrespective of the
// dimension. So that is what we return:
template <int dim>
double
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  return p.square();
}

template <int dim>
DGHeat<dim>::DGHeat(simCase scase, int refinementLevel)
  : fe(1)
  , dof_handler(triangulation)
  , simulationCase_(scase)
  , refinementLevel_(refinementLevel)
{}

template <int dim>
void
DGHeat<dim>::make_grid()
{
  if (simulationCase_ == MMS)
    make_cube_grid();
  else if (simulationCase_ == TaylorCouette)
    make_ring_grid();
}

template <int dim>
void
DGHeat<dim>::make_cube_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(refinementLevel_);

  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}

template <int dim>
void
DGHeat<dim>::make_ring_grid()
{
  const double inner_radius = 0.25, outer_radius = 1.0;
  if (dim == 2)
    center = Point<dim>(0, 0);
  GridGenerator::hyper_shell(
    triangulation, center, inner_radius, outer_radius, 10, true);

  triangulation.refine_global(refinementLevel_);

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;

  std::cout << "Number of total cells: " << triangulation.n_cells()
            << std::endl;
}


template <int dim>
void
DGHeat<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
DGHeat<dim>::assemble_system()
{
  QGauss<dim> quadrature_formula(5);

  RightHandSideMMS<dim> right_hand_side;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;

      for (unsigned int q_index = 0; q_index < n_q_points; ++q_index)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              // Stiffness Matrix
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) *
                 fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));

            if (simulationCase_ == MMS)
              {
                // Right Hand Side
                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_index) *
                   right_hand_side.value(fe_values.quadrature_point(q_index)) *
                   fe_values.JxW(q_index));
              }
          }


      // Assemble global matrix
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            system_matrix.add(local_dof_indices[i],
                              local_dof_indices[j],
                              cell_matrix(i, j));

          system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }


  std::map<types::global_dof_index, double> boundary_values;
  if (simulationCase_ == TaylorCouette)
    {
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values);
      VectorTools::interpolate_boundary_values(
        dof_handler, 1, Functions::ConstantFunction<dim>(1.), boundary_values);
    }

  if (simulationCase_ == MMS)
    {
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values);
    }

  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}

template <int dim>
void
DGHeat<dim>::solve()
{
  SolverControl solver_control(10000, 1e-12);
  SolverCG<>    solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  // We have made one addition, though: since we suppress output from the
  // linear solvers, we have to print the number of iterations by hand.
  std::cout << "   " << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}

template <int dim>
void
DGHeat<dim>::output_results() const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  std::string dimension(dim == 2 ? "solution-2d-case-" : "solution-3d-case-");

  std::string fname = dimension + Utilities::int_to_string(simulationCase_) +
                      "-" + Utilities::int_to_string(refinementLevel_) + ".vtk";

  std::ofstream output(fname.c_str());
  data_out.write_vtk(output);
}

// Find the l2 norm of the error between the finite element sol'n and the exact
// sol'n
template <int dim>
void
DGHeat<dim>::calculateL2Error()
{
  QGauss<dim>   quadrature_formula(5);
  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell =
    fe.dofs_per_cell; // This gives you dofs per cell
  std::vector<types::global_dof_index> local_dof_indices(
    dofs_per_cell); //  Local connectivity

  const unsigned int n_q_points = quadrature_formula.size();

  double l2error = 0.;

  // loop over elements
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (; cell != endc; ++cell)
    {
      fe_values.reinit(cell);

      // Retrieve the effective "connectivity matrix" for this element
      cell->get_dof_indices(local_dof_indices);


      for (unsigned int q = 0; q < n_q_points; q++)
        {
          const double x = fe_values.quadrature_point(q)[0];
          const double y = fe_values.quadrature_point(q)[1];
          if (dim > 2)
            const double z = fe_values.quadrature_point(q)[2];

          const double r       = std::sqrt(x * x + y * y);
          const double lnratio = std::log(1. / 0.25);
          double       u_exact = 0.;
          if (simulationCase_ == TaylorCouette)
            u_exact = 1. / (lnratio)*std::log(r / 0.25);
          if (simulationCase_ == MMS)
            u_exact = -sin(M_PI * x) * std::sin(M_PI * y);
          double u_sim = 0;

          // Find the values of x and u_h (the finite element solution) at the
          // quadrature points
          for (unsigned int i = 0; i < dofs_per_cell; i++)
            {
              u_sim +=
                fe_values.shape_value(i, q) * solution[local_dof_indices[i]];
            }
          l2error += (u_sim - u_exact) * (u_sim - u_exact) * fe_values.JxW(q);
          //       std::cout << " x = " << x << " y = " << y <<  " r = " << r <<
          //       "   u_exact = " << u_exact << "   u_sim=" << u_sim <<
          //       std::endl;
        }
    }
  std::cout << "L2Error is : " << std::sqrt(l2error) << std::endl;
  L2Error_ = std::sqrt(l2error);
}



template <int dim>
void
DGHeat<dim>::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
  calculateL2Error();
}

int
main()
{
  deallog.depth_console(0);

  // Taylor couette
  {
    int                 nsim = 5;
    std::vector<double> l2error;
    std::vector<int>    size;
    for (int m = 1; m < 1 + nsim; ++m)
      {
        std::cout << "Solving Taylor-Couette problem 2D - "
                  << " with mesh - " << m << std::endl;
        DGHeat<2> taylorCouette_problem_2d(TaylorCouette, m);
        taylorCouette_problem_2d.run();
        l2error.push_back(taylorCouette_problem_2d.getL2Error());
        size.push_back(m);
      }
    std::ofstream output_file("./L2Error-TaylorCouette.dat");
    for (unsigned int i = 0; i < size.size(); ++i)
      {
        output_file << size[i] << " " << l2error[i] << std::endl;
      }
  }

  // MMS
  {
    int                 nsim = 8;
    std::vector<double> l2error;
    std::vector<int>    size;
    for (int m = 2; m < 1 + nsim; ++m)
      {
        std::cout << "Solving MMS problem 2D - "
                  << " with mesh - " << m << std::endl;
        DGHeat<2> mms_problem_2d(MMS, m);
        mms_problem_2d.run();
        l2error.push_back(mms_problem_2d.getL2Error());
        size.push_back(m);
      }
    std::ofstream output_file("./L2Error-MMS.dat");
    for (unsigned int i = 0; i < size.size(); ++i)
      {
        output_file << size[i] << " " << l2error[i] << std::endl;
      }
  }

  return 0;
}
