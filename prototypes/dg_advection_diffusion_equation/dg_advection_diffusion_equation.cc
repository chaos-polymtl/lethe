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
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
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
  SmithHutton,
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
  RightHandSideMMS(const double Pe_diff)
    : Function<dim>()
    , Pe_diff(Pe_diff)

  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const;

  const double Pe_diff;
};

template <int dim>
double
RightHandSideMMS<dim>::value(const Point<dim> &p,
                             const unsigned int /*component*/) const
{
  double return_value = 0.0;
  double x            = p(0);
  double y            = p(1);

  return_value =
    1. / Pe_diff * 2. * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);


  return return_value;
}

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
BoundaryValues<dim>::value(const Point<dim> &p,
                           const unsigned int /*component*/) const
{
  double alpha = 10.;
  double x     = p(0);
  double y     = p(1);

  if (y < 1e-16)
    return 1 + std::tanh(alpha * (2 * x + 1));
  else
    return 1 + std::tanh(alpha * (-1.));
}

template <int dim>
class VelocityFieldMMS : public Function<dim>
{
public:
  VelocityFieldMMS(double Pe)
    : Function<dim>(dim)
    , Pe(Pe)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const override;

  // virtual void
  // vector_value(const Point<dim> &p, Vector<double> &value) const override;
private:
  const double Pe;
};


template <int dim>
double
VelocityFieldMMS<dim>::value(const Point<dim> & p,
                             const unsigned int component) const
{
  double x = p(0);
  double y = p(1);

  if (component == 0)
    return Pe * std::sin(M_PI * x) * std::sin(M_PI * x) * std::sin(M_PI * y) *
           std::cos(M_PI * y);
  else // (component==1)
    return -Pe * std::cos(M_PI * x) * std::sin(M_PI * x) * std::sin(M_PI * y) *
           std::sin(M_PI * y);
}

template <int dim>
class VelocityFieldSmithHutton : public Function<dim>
{
public:
  VelocityFieldSmithHutton(double Pe)
    : Function<dim>(dim)
    , Pe(Pe)
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component) const override;

  // virtual void
  // vector_value(const Point<dim> &p, Vector<double> &value) const override;
private:
  const double Pe;
};


template <int dim>
double
VelocityFieldSmithHutton<dim>::value(const Point<dim> & p,
                                     const unsigned int component) const
{
  double x = p(0);
  double y = p(1);

  if (component == 0)
    return Pe * 2 * y * (1. - x * x);
  else // (component==1)
    return -Pe * 2 * x * (1. - y * y);
}

template <int dim>
class DGAdvectionDiffusion
{
public:
  DGAdvectionDiffusion(simCase      scase,
                       unsigned int initial_level,
                       unsigned int final_level,
                       const double Pe);
  void
  run();

private:
  void
  make_grid();
  void
  make_cube_grid();
  void
  make_smith_hutton_grid();
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

  const double Pe_diff;

  RightHandSideMMS<dim> right_hand_side_function;

  std::unique_ptr<Function<dim>> velocity_function;
  // VelocityFieldMMS<dim> velocity_function;



  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
  Point<dim>     center;

  ConvergenceTable error_table;

  simCase simulation_case;

  unsigned int initial_refinement_level;
  unsigned int number_refinement;

  const double beta;
};



template <int dim>
DGAdvectionDiffusion<dim>::DGAdvectionDiffusion(simCase      scase,
                                                unsigned int initial_refinement,
                                                unsigned int number_refinement,
                                                const double Pe)
  : fe(2)
  , mapping(2)
  , dof_handler(triangulation)
  , Pe_diff(Pe)
  , right_hand_side_function(Pe_diff)
  , simulation_case(scase)
  , initial_refinement_level(initial_refinement)
  , number_refinement(number_refinement)
  , beta(5.)
{
  if (simulation_case == MMS)
    velocity_function = std::make_unique<VelocityFieldMMS<dim>>(1);

  else if (simulation_case == SmithHutton)
    velocity_function = std::make_unique<VelocityFieldSmithHutton<dim>>(1);
}

template <int dim>
void
DGAdvectionDiffusion<dim>::make_grid()
{
  if (simulation_case == MMS)
    make_cube_grid();

  else if (simulation_case == SmithHutton)
    make_smith_hutton_grid();
}

template <int dim>
void
DGAdvectionDiffusion<dim>::make_cube_grid()
{
  GridGenerator::hyper_cube(triangulation, -1, 1);
  triangulation.refine_global(initial_refinement_level);

  GridTools::distort_random(0.1, triangulation, true);

  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}

template <int dim>
void
DGAdvectionDiffusion<dim>::make_smith_hutton_grid()
{
  GridIn<dim> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file("smith_hutton.msh");
  grid_in.read_msh(input_file);

  triangulation.refine_global(initial_refinement_level);

  GridTools::distort_random(0.1, triangulation, true);


  std::cout << "   Number of active cells: " << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: " << triangulation.n_cells()
            << std::endl;
}

template <int dim>
void
DGAdvectionDiffusion<dim>::setup_system()
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
DGAdvectionDiffusion<dim>::assemble_system()
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
    right_hand_side_function.value_list(q_points, f);

    std::vector<Vector<double>> velocity_list(q_points.size(),
                                              Vector<double>(dim));

    velocity_function->vector_value_list(q_points, velocity_list);

    Tensor<1, dim> velocity_q;


    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
      {
        // Establish the velocity
        for (int i = 0; i < dim; ++i)
          velocity_q[i] = velocity_function->value(q_points[point], i);

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data.cell_matrix(i, j) +=
                  1. / Pe_diff * fe_v.shape_grad(i, point) // \nabla \phi_i
                  * fe_v.shape_grad(j, point)              // \nabla \phi_j
                  * JxW[point];                            // dx

                copy_data.cell_matrix(i, j) +=
                  -velocity_q * fe_v.shape_grad(i, point) *
                  fe_v.shape_value(j, point) * JxW[point];
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

    std::vector<Vector<double>> velocity_list(q_points.size(),
                                              Vector<double>(dim));

    velocity_function->vector_value_list(q_points, velocity_list);

    Tensor<1, dim> velocity_q;

    double h;
    if (dim == 2)
      h = std::sqrt(4. * cell->measure() / M_PI);
    else if (dim == 3)
      h = pow(6 * cell->measure() / M_PI, 1. / 3.);


    for (unsigned int point = 0; point < q_points.size(); ++point)
      {
        //        if (cell->face(face_no)->boundary_id() == 0)
        {
          // Establish the velocity
          for (int i = 0; i < dim; ++i)
            velocity_q[i] = velocity_function->value(q_points[point], i);

          const double velocity_dot_n = velocity_q * normals[point];

          for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              for (unsigned int j = 0; j < n_facet_dofs; ++j)
                {
                  if (cell->face(face_no)->boundary_id() == 0)
                    {
                      copy_data.cell_matrix(i, j) +=
                        -1. / Pe_diff * normals[point] *
                        fe_face.shape_grad(i, point)    // n*\nabla \phi_i
                        * fe_face.shape_value(j, point) // \phi_j
                        * JxW[point];                   // dx

                      copy_data.cell_matrix(i, j) +=
                        -1. / Pe_diff * fe_face.shape_value(i, point) // \phi_i
                        * normals[point] *
                        fe_face.shape_grad(j, point)
                        // n*\nabla \phi_j
                        * JxW[point]; // dx

                      copy_data.cell_matrix(i, j) +=
                        beta * 1. / h / Pe_diff *
                        fe_face.shape_value(i, point)                 // \phi_i
                        * fe_face.shape_value(j, point) * JxW[point]; // dx
                    }

                  //                  if (velocity_dot_n < 0)
                  //                    copy_data.cell_matrix(i, j) +=
                  //                      -velocity_dot_n *
                  //                      fe_face.shape_value(i, point) *
                  //                      fe_face.shape_value(j, point) *
                  //                      JxW[point];

                  if (velocity_dot_n > 0)
                    {
                      copy_data.cell_matrix(i, j) +=
                        fe_face.shape_value(i, point)   // \phi_i
                        * fe_face.shape_value(j, point) // \phi_j
                        * velocity_dot_n                // \beta . n
                        * JxW[point];                   // dx
                    }
                }
            }

          if (simulation_case == SmithHutton)
            for (unsigned int i = 0; i < n_facet_dofs; ++i)
              {
                if (cell->face(face_no)->boundary_id() == 0)
                  {
                    copy_data.cell_rhs(i) +=
                      beta * 1. / h / Pe_diff *
                      fe_face.shape_value(i, point) // \phi_i
                      * g[point]                    // g
                      * JxW[point];                 // dx
                    copy_data.cell_rhs(i) +=
                      -1. / Pe_diff * normals[point] *
                      fe_face.shape_grad(i, point) // n*\nabla \phi_i
                      * g[point]                   // g
                      * JxW[point];                // dx
                  }

                if (velocity_dot_n < 0)
                  copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) *
                                           g[point] * velocity_dot_n *
                                           JxW[point];
              }
        }
        //        else
        {
          //            // Establish the velocity
          //            for (int i = 0; i < dim; ++i)
          //              velocity_q[i] =
          //              velocity_function->value(q_points[point], i);

          //            const double velocity_dot_n = velocity_q *
          //            normals[point]; for (unsigned int i = 0; i <
          //            n_facet_dofs; ++i)
          //              {
          //                for (unsigned int j = 0; j < n_facet_dofs; ++j)
          //                  {
          //                    if (velocity_dot_n > 0)
          //                      copy_data.cell_matrix(i, j) +=
          //                        -velocity_dot_n * fe_face.shape_value(i,
          //                        point) * fe_face.shape_value(j, point) *
          //                        JxW[point];
          //                  }

          //                if (velocity_dot_n < 0)
          //                  copy_data.cell_rhs(i) += -fe_face.shape_value(i,
          //                  point) *
          //                                           g[point] *
          //                                           velocity_dot_n *
          //                                           JxW[point];
          //              }
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

    std::vector<Vector<double>> velocity_list(q_points.size(),
                                              Vector<double>(dim));

    velocity_function->vector_value_list(q_points, velocity_list);

    Tensor<1, dim> velocity_q;

    double h;
    if (dim == 2)
      h = std::sqrt(4. * cell->measure() / M_PI);
    else if (dim == 3)
      h = pow(6 * cell->measure() / M_PI, 1. / 3.);

    for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint)
      {
        // Establish the velocity
        for (int i = 0; i < dim; ++i)
          velocity_q[i] = velocity_function->value(q_points[qpoint], i);

        const double velocity_dot_n = velocity_q * normals[qpoint];

        const double gamma = velocity_dot_n > 0 ? 0.5 : -0.5;

        for (unsigned int i = 0; i < n_dofs; ++i)
          {
            for (unsigned int j = 0; j < n_dofs; ++j)
              {
                copy_data_face.cell_matrix(i, j) +=
                  -1. / Pe_diff * normals[qpoint] *
                  fe_iv.average_gradient(i, qpoint) * fe_iv.jump(j, qpoint) *
                  JxW[qpoint];

                copy_data_face.cell_matrix(i, j) +=
                  -1. / Pe_diff * fe_iv.jump(i, qpoint) // \phi_i
                  * fe_iv.average_gradient(j, qpoint) *
                  normals[qpoint] // n*\nabla \phi_j
                  * JxW[qpoint];  // dx

                copy_data_face.cell_matrix(i, j) +=
                  beta * 1. / h / Pe_diff * fe_iv.jump(i, qpoint) *
                  fe_iv.jump(j, qpoint) * JxW[qpoint];

                copy_data_face.cell_matrix(i, j) +=
                  fe_iv.jump(i, qpoint) // [\phi_i]
                  * fe_iv.shape_value((velocity_dot_n > 0), j, qpoint) *
                  velocity_dot_n * JxW[qpoint];

                //                copy_data_face.cell_matrix(i, j) +=
                //                  fe_iv.average(i, qpoint)
                //                  * velocity_dot_n *
                //                  fe_iv.jump(j,
                //                             qpoint)
                //                  * JxW[qpoint];

                //                copy_data_face.cell_matrix(i, j) +=
                //                  gamma * fe_iv.jump(i, qpoint) // [\phi_i]
                //                  * velocity_dot_n *
                //                  fe_iv.jump(j,
                //                             qpoint)
                //                  * JxW[qpoint];
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
DGAdvectionDiffusion<dim>::solve()
{
  // SolverControl solver_control(50000, 1e-12);
  // SolverGMRES<> solver(solver_control);
  // solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
  //
  //// We have made one addition, though: since we suppress output from the
  //// linear solvers, we have to print the number of iterations by hand.
  // std::cout << "   " << solver_control.last_step()
  //          << " GMRES iterations needed to obtain convergence." << std::endl;

  SparseDirectUMFPACK direct;
  direct.initialize(system_matrix);
  direct.vmult(solution, system_rhs);
}

template <int dim>
void
DGAdvectionDiffusion<dim>::output_results(unsigned int it) const
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
DGAdvectionDiffusion<dim>::calculateL2Error()
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

          double u_exact = 0.;
          if (simulation_case == MMS)
            u_exact = sin(M_PI * x) * std::sin(M_PI * y);
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
  error_table.add_value("error", std::sqrt(l2error));
  error_table.add_value("cells", triangulation.n_global_active_cells());
}



template <int dim>
void
DGAdvectionDiffusion<dim>::run()
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
      if (simulation_case == MMS)
        calculateL2Error();
    }

  if (simulation_case == MMS)
    {
      error_table.omit_column_from_convergence_rate_evaluation("cells");
      error_table.evaluate_all_convergence_rates(
        ConvergenceTable::reduction_rate_log2);

      error_table.set_scientific("error", true);

      error_table.write_text(std::cout);
    }
}

int
main()
{
  deallog.depth_console(0);

  // MMS
  {
    std::cout << "Solving MMS problem 2D " << std::endl;
    DGAdvectionDiffusion<2> mms_problem_2d(MMS, 2, 6, 10);
    mms_problem_2d.run();
  }

  // Smith Hutton
  {
    std::cout << "Solving Smith-Hutton problem 2D " << std::endl;
    DGAdvectionDiffusion<2> problem_2d(SmithHutton, 2, 6, 1e12);
    problem_2d.run();
  }
  return 0;
}
