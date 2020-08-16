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
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

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

#include <deal.II/meshworker/mesh_loop.h>

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
struct ScratchData
{
  ScratchData(const Mapping<dim> &      mapping,
              const FiniteElement<dim> &fe,
              const unsigned int        quadrature_degree,
              const UpdateFlags         update_flags = update_values |
                                               update_gradients |
                                               update_quadrature_points |
                                               update_JxW_values,
              const UpdateFlags interface_update_flags =
                update_values | update_gradients | update_quadrature_points |
                update_JxW_values | update_normal_vectors)
    : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
    , fe_interface_values(mapping,
                          fe,
                          QGauss<dim - 1>(quadrature_degree),
                          interface_update_flags)
  {}


  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(scratch_data.fe_values.get_mapping(),
                          scratch_data.fe_values.get_fe(),
                          scratch_data.fe_interface_values.get_quadrature(),
                          scratch_data.fe_interface_values.get_update_flags())
  {}

  FEValues<dim>          fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};



struct CopyDataFace
{
  FullMatrix<double>                   cell_matrix;
  std::vector<types::global_dof_index> joint_dof_indices;
};



struct CopyData
{
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
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
  if (p.square() > 0.9)
    return 1;
  else
    return 0.;
}

template <int dim>
class DGHeat
{
public:
  DGHeat(simCase scase, unsigned int initial_level, unsigned int final_level);
  void
  run();

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
  output_results(unsigned int it) const;
  void
  calculateL2Error();

  Triangulation<dim>  triangulation;
  FE_DGQ<dim>         fe;
  const MappingQ<dim> mapping;
  DoFHandler<dim>     dof_handler;


  RightHandSideMMS<dim> right_hand_side;


  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
  Point<dim>     center;

  ConvergenceTable error_table;

  simCase simulation_case;

  unsigned int initial_refinement_level;
  unsigned int number_refinement;
};



template <int dim>
DGHeat<dim>::DGHeat(simCase      scase,
                    unsigned int initial_refinement,
                    unsigned int number_refinement)
  : fe(1)
  , mapping(1)
  , dof_handler(triangulation)
  , simulation_case(scase)
  , initial_refinement_level(initial_refinement)
  , number_refinement(number_refinement)
{}

template <int dim>
void
DGHeat<dim>::make_grid()
{
  if (simulation_case == MMS)
    make_cube_grid();
  else if (simulation_case == TaylorCouette)
    make_ring_grid();
}

template <int dim>
void
DGHeat<dim>::make_cube_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(initial_refinement_level);

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

  triangulation.refine_global(initial_refinement_level);

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
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


template <int dim>
void
DGHeat<dim>::assemble_system()
{
  using Iterator = typename DoFHandler<dim>::active_cell_iterator;
  const BoundaryValues<dim> boundary_function;

  auto cell_worker = [&](const Iterator &  cell,
                         ScratchData<dim> &scratch_data,
                         CopyData &        copy_data) {
    const unsigned int n_dofs = scratch_data.fe_values.get_fe().dofs_per_cell;
    copy_data.reinit(cell, n_dofs);
    scratch_data.fe_values.reinit(cell);

    const auto &q_points = scratch_data.fe_values.get_quadrature_points();

    const FEValues<dim> &      fe_v = scratch_data.fe_values;
    const std::vector<double> &JxW  = fe_v.get_JxW_values();

    std::vector<double> f(q_points.size());
    right_hand_side.value_list(q_points, f);

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
      {
        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data.cell_matrix(i, j) +=
                  fe_v.shape_grad(i, point)   // \nabla \phi_i
                  * fe_v.shape_grad(j, point) // \nabla \phi_j
                  * JxW[point];               // dx
              }

            if (simulation_case == MMS)
              // Right Hand Side
              copy_data.cell_rhs(i) +=
                (fe_v.shape_value(i, point) * f[point] * JxW[point]);
          }
      }
  };

  auto boundary_worker = [&](const Iterator &    cell,
                             const unsigned int &face_no,
                             ScratchData<dim> &  scratch_data,
                             CopyData &          copy_data) {
    scratch_data.fe_interface_values.reinit(cell, face_no);

    const FEFaceValuesBase<dim> &fe_face =
      scratch_data.fe_interface_values.get_fe_face_values(0);

    const auto &       q_points     = fe_face.get_quadrature_points();
    const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
    const std::vector<double> &JxW  = fe_face.get_JxW_values();

    const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors();
    std::vector<double>                g(q_points.size());

    boundary_function.value_list(q_points, g);

    double h;
    if (dim == 2)
      h = std::sqrt(4. * cell->measure() / M_PI);
    else if (dim == 3)
      h = pow(6 * cell->measure() / M_PI, 1. / 3.);


    const double beta = 10.;

    for (unsigned int point = 0; point < q_points.size(); ++point)
      {
        for (unsigned int i = 0; i < n_facet_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_facet_dofs; ++j)
              {
                copy_data.cell_matrix(i, j) +=
                  -normals[point] *
                  fe_face.shape_grad(i, point)    // n*\nabla \phi_i
                  * fe_face.shape_value(j, point) // \phi_j
                  * JxW[point];                   // dx

                copy_data.cell_matrix(i, j) +=
                  -fe_face.shape_value(i, point) // \phi_i
                  * fe_face.shape_grad(j, point) *
                  normals[point] // n*\nabla \phi_j
                  * JxW[point];  // dx

                copy_data.cell_matrix(i, j) +=
                  beta * 1. / h * fe_face.shape_value(i, point) // \phi_i
                  * fe_face.shape_value(j, point) * JxW[point]; // dx
              }
          }

        if (simulation_case == TaylorCouette)
          for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              copy_data.cell_rhs(i) += beta * 1. / h *
                                       fe_face.shape_value(i, point) // \phi_i
                                       * g[point]                    // g
                                       * JxW[point];                 // dx
              copy_data.cell_rhs(i) +=
                -normals[point] *
                fe_face.shape_grad(i, point) // n*\nabla \phi_i
                * g[point]                   // g
                * JxW[point];                // dx
            }
      }
  };

  auto face_worker = [&](const Iterator &    cell,
                         const unsigned int &f,
                         const unsigned int &sf,
                         const Iterator &    ncell,
                         const unsigned int &nf,
                         const unsigned int &nsf,
                         ScratchData<dim> &  scratch_data,
                         CopyData &          copy_data) {
    FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;

    fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

    const auto &q_points = fe_iv.get_quadrature_points();

    copy_data.face_data.emplace_back();
    CopyDataFace &copy_data_face = copy_data.face_data.back();

    const unsigned int n_dofs        = fe_iv.n_current_interface_dofs();
    copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices();

    copy_data_face.cell_matrix.reinit(n_dofs, n_dofs);

    const std::vector<double> &        JxW     = fe_iv.get_JxW_values();
    const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();


    double h;
    if (dim == 2)
      h = std::sqrt(4. * cell->measure() / M_PI);
    else if (dim == 3)
      h = pow(6 * cell->measure() / M_PI, 1. / 3.);

    const double beta = 10.;

    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      {
        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data_face.cell_matrix(i, j) +=
                  -normals[qpoint] * fe_iv.average_gradient(i, qpoint) *
                  fe_iv.jump(j, qpoint) * JxW[qpoint];

                copy_data_face.cell_matrix(i, j) +=
                  -fe_iv.jump(i, qpoint) // \phi_i
                  * fe_iv.average_gradient(j, qpoint) *
                  normals[qpoint] // n*\nabla \phi_j
                  * JxW[qpoint];  // dx

                copy_data_face.cell_matrix(i, j) +=
                  beta * 1. / h * fe_iv.jump(i, qpoint) *
                  fe_iv.jump(j, qpoint) * JxW[qpoint];
              }
          }
      }
  };

  AffineConstraints<double> constraints;

  auto copier = [&](const CopyData &c) {
    constraints.distribute_local_to_global(c.cell_matrix,
                                           c.cell_rhs,
                                           c.local_dof_indices,
                                           system_matrix,
                                           system_rhs);

    for (auto &cdf : c.face_data)
      {
        constraints.distribute_local_to_global(cdf.cell_matrix,
                                               cdf.joint_dof_indices,
                                               system_matrix);
      }
  };

  const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1;

  ScratchData<dim> scratch_data(mapping, fe, n_gauss_points);
  CopyData         copy_data;

  MeshWorker::mesh_loop(dof_handler.begin_active(),
                        dof_handler.end(),
                        cell_worker,
                        copier,
                        scratch_data,
                        copy_data,
                        MeshWorker::assemble_own_cells |
                          MeshWorker::assemble_boundary_faces |
                          MeshWorker::assemble_own_interior_faces_once,
                        boundary_worker,
                        face_worker);
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
DGHeat<dim>::output_results(unsigned int it) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  data_out.build_patches();

  std::string dimension(dim == 2 ? "solution-2d-case-" : "solution-3d-case-");

  std::string fname = dimension + Utilities::int_to_string(simulation_case) +
                      "-" + Utilities::int_to_string(it) + ".vtk";

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

          const double r       = std::sqrt(x * x + y * y);
          const double lnratio = std::log(1. / 0.25);
          double       u_exact = 0.;
          if (simulation_case == TaylorCouette)
            u_exact = 1. / (lnratio)*std::log(r / 0.25);
          if (simulation_case == MMS)
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
        }
    }


  std::cout << "L2Error is : " << std::sqrt(l2error) << std::endl;
  error_table.add_value("error", std::sqrt(l2error));
  error_table.add_value("cells", triangulation.n_global_active_cells());
}



template <int dim>
void
DGHeat<dim>::run()
{
  make_grid();
  for (unsigned int it = 0; it < number_refinement; ++it)
    {
      if (it > 0)
        triangulation.refine_global(1);
      setup_system();
      assemble_system();
      solve();
      output_results(it);
      calculateL2Error();
    }

  error_table.omit_column_from_convergence_rate_evaluation("cells");
  error_table.evaluate_all_convergence_rates(
    ConvergenceTable::reduction_rate_log2);

  error_table.set_scientific("error", true);

  error_table.write_text(std::cout);
}

int
main()
{
  deallog.depth_console(0);

  // Taylor couette
  {
    std::cout << "Solving Taylor-Couette problem 2D  " << std::endl;
    DGHeat<2> taylorCouette_problem_2d(TaylorCouette, 1, 6);
    taylorCouette_problem_2d.run();
  }

  // MMS
  {
    std::cout << "Solving MMS problem 2D " << std::endl;
    DGHeat<2> mms_problem_2d(MMS, 2, 8);
    mms_problem_2d.run();
  }
  return 0;
}
