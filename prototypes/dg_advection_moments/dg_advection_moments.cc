// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/time_stepping.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <numbers>

using namespace dealii;

// ============================================================================
// 1. EQUATION DATA
// ============================================================================

// Velocity: v = 2*pi in 1D
template <int dim>
class VelocityField : public Function<dim>
{
public:
  VelocityField()
    : Function<dim>(dim)
  {}
  void
  vector_value(const Point<dim> & /*p*/, Vector<double> &values) const override
  {
    if constexpr (dim == 1)
      values[0] = 2.0 * std::numbers::pi;
    else
      values[0] = 1.0;
  }
};

// Exact Solution: u = sin(x - 2*pi*t)
template <int dim>
class AnalyticalSolution : public Function<dim>
{
public:
  AnalyticalSolution()
    : Function<dim>(1)
  {}
  double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    return std::sin(p(0) - 2.0 * std::numbers::pi * this->get_time());
  }
};

// Initial Condition: t=0
template <int dim>
class InitialCondition : public Function<dim>
{
public:
  InitialCondition()
    : Function<dim>(1)
  {}
  double
  value(const Point<dim> &p, const unsigned int = 0) const override
  {
    return std::sin(p(0));
  }
};

// ============================================================================
// 2. MESH WORKER DATA STRUCTURES
// ============================================================================

template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim>        &mapping,
              const FiniteElement<dim>  &fe,
              const Quadrature<dim>     &quadrature,
              const Quadrature<dim - 1> &quadrature_face)
    : fe_values(mapping,
                fe,
                quadrature,
                update_values | update_gradients | update_JxW_values |
                  update_quadrature_points)
    , fe_interface_values(mapping,
                          fe,
                          quadrature_face,
                          update_values | update_JxW_values |
                            update_normal_vectors | update_quadrature_points)
  {}

  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
    , fe_interface_values(scratch_data.fe_interface_values.get_mapping(),
                          scratch_data.fe_interface_values.get_fe(),
                          scratch_data.fe_interface_values.get_quadrature(),
                          scratch_data.fe_interface_values.get_update_flags())
  {}

  FEValues<dim>          fe_values;
  FEInterfaceValues<dim> fe_interface_values;
};

struct CopyDataFace
{
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> joint_dof_indices;
};

struct CopyData
{
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;
  std::vector<CopyDataFace>            face_data;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    face_data.clear();
    cell_rhs.reinit(dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }
};

// ============================================================================
// 3. MAIN SOLVER CLASS
// ============================================================================

template <int dim>
class DGAdvection
{
public:
  DGAdvection(const unsigned int degree, const unsigned int refinements)
    : fe_degree(degree)
    , n_refinements(refinements)
    , dof_handler(triangulation)
    , fe(degree)
    , mapping(degree)
  {}

  void
  run();

private:
  void
  make_grid();
  void
  setup_system();
  void
  assemble_mass_matrix();
  void
  assemble_rhs(const Vector<double> &current_solution,
               double                time,
               Vector<double>       &rhs);
  void
  apply_inverse_mass(Vector<double> &dst, const Vector<double> &src);
  void
  output_results(unsigned int step, double time);

  const unsigned int fe_degree;
  const unsigned int n_refinements;

  Triangulation<dim>              triangulation;
  DoFHandler<dim>                 dof_handler;
  FE_DGQ<dim>                     fe;
  MappingQ<dim>                   mapping;
  AffineConstraints<double>       constraints;
  SparseMatrix<double>            mass_matrix;
  SparsityPattern                 sparsity_pattern;
  std::vector<FullMatrix<double>> inverse_local_mass;
  Vector<double>                  solution;
  VelocityField<dim>              velocity;
};

// ----------------------------------------------------------------------------
// Grid & Setup
// ----------------------------------------------------------------------------

template <int dim>
void
DGAdvection<dim>::make_grid()
{
  // 1. Create the grid
  // CRITICAL FIX: The 'true' argument colorizes the boundaries.
  // Left = 0, Right = 1. Without this, all IDs are 0.
  GridGenerator::hyper_cube(triangulation, 0, 2.0 * std::numbers::pi, true);

  // 2. Apply Periodicity
  // Connect ID 0 (Left) to ID 1 (Right)
  std::vector<
    GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
    periodicity_vector;

  GridTools::collect_periodic_faces(triangulation,
                                      0,
                                      1,
                                      0,
                                      periodicity_vector);

  std::cout << "Periodic faces found: " << periodicity_vector.size() << std::endl;

  triangulation.add_periodicity(periodicity_vector);

  // 3. Refine
  triangulation.refine_global(n_refinements);
}

template <int dim>
void
DGAdvection<dim>::setup_system()
{
  dof_handler.distribute_dofs(fe);

  constraints.clear();
  constraints.close();

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  mass_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());

  assemble_mass_matrix(); // Build the inverse mass blocks immediately
}

template <int dim>
void
DGAdvection<dim>::assemble_mass_matrix()
{
  // Quadrature for Mass Matrix (exact integration)
  QGauss<dim> q_mass(fe_degree + 1);
  FEValues<dim> fe_values(mapping, fe, q_mass, update_values | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q           = q_mass.size();

  inverse_local_mass.resize(triangulation.n_active_cells());
  FullMatrix<double> Mlocal(dofs_per_cell, dofs_per_cell);

  unsigned int cell_index = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      Mlocal = 0;
      fe_values.reinit(cell);
      for (unsigned int q = 0; q < n_q; ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            Mlocal(i, j) += fe_values.shape_value(i, q) *
                            fe_values.shape_value(j, q) * fe_values.JxW(q);

      Mlocal.gauss_jordan(); // Invert in place
      inverse_local_mass[cell_index] = Mlocal;
      cell_index++;
    }
}

// ----------------------------------------------------------------------------
// RHS Assembly
// ----------------------------------------------------------------------------

template <int dim>
void DGAdvection<dim>::assemble_rhs(const Vector<double> &current_solution,
                                    double /*time*/,
                                    Vector<double> &rhs_vector)
{
  rhs_vector = 0;

  QGauss<dim>     q_cell(fe_degree + 1);
  QGauss<dim - 1> q_face(fe_degree + 1);

  FEValues<dim> fe_values(mapping, fe, q_cell, update_values | update_gradients | update_JxW_values | update_quadrature_points);
  FEInterfaceValues<dim> fe_iv(mapping, fe, q_face, update_values | update_JxW_values | update_normal_vectors | update_quadrature_points);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  // 1. VOLUME INTEGRALS (Loop over all cells)
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    cell_rhs = 0;
    fe_values.reinit(cell);
    cell->get_dof_indices(local_dof_indices);

    Vector<double> sol_loc(dofs_per_cell);
    cell->get_dof_values(current_solution, sol_loc);

    Vector<double> vel(dim);
    Tensor<1, dim> vel_tensor;

    for (unsigned int q = 0; q < q_cell.size(); ++q)
    {
      velocity.vector_value(fe_values.quadrature_point(q), vel);
      for(int d=0; d<dim; ++d) vel_tensor[d] = vel[d];

      double u_val = 0;
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        u_val += fe_values.shape_value(i, q) * sol_loc[i];

      // Volume Term: (grad v, a*u)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        cell_rhs(i) += (fe_values.shape_grad(i, q) * vel_tensor * u_val) * fe_values.JxW(q);
    }
    rhs_vector.add(local_dof_indices, cell_rhs);
  }

  // 2. FLUX INTEGRALS (Manual Loop over Faces)
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
    {
      auto neighbor = cell->neighbor(f);

      // CRITICAL: This checks if we should process this interface.
      // For periodic faces, 'neighbor' exists, and we only process if index > cell index
      // to avoid double counting.
      if (neighbor.state() == IteratorState::valid && neighbor->index() > cell->index())
      {
         const unsigned int neighbor_face_no = cell->neighbor_of_neighbor(f);

         // Initialize FEInterfaceValues for the pair (cell, neighbor)
         fe_iv.reinit(cell, f, 0, neighbor, neighbor_face_no, 0);

         const unsigned int n_face_q = q_face.size();

         // Prepare local RHS vectors for both sides
         Vector<double> rhs_in(dofs_per_cell);
         Vector<double> rhs_out(dofs_per_cell);
         rhs_in = 0; rhs_out = 0;

         // Get Solution values on both sides
         std::vector<double> u_in(n_face_q);
         std::vector<double> u_out(n_face_q);
         fe_iv.get_fe_face_values(0).get_function_values(current_solution, u_in);
         fe_iv.get_fe_face_values(1).get_function_values(current_solution, u_out);

         Vector<double> vel(dim);
         Tensor<1, dim> vel_tensor;

         for (unsigned int q = 0; q < n_face_q; ++q)
         {
           velocity.vector_value(fe_iv.quadrature_point(q), vel);
           for(int d=0; d<dim; ++d) vel_tensor[d] = vel[d];

           const double v_dot_n = vel_tensor * fe_iv.normal_vector(q);

           // === LAX-FRIEDRICHS FLUX ===
           const double u_avg  = 0.5 * (u_in[q] + u_out[q]);
           const double u_jump = u_in[q] - u_out[q];

           double flux = u_avg * v_dot_n + 0.5 * std::abs(v_dot_n) * u_jump;

           // Update Residuals
           // Cell (In): - phi * Flux
           for (unsigned int i = 0; i < dofs_per_cell; ++i)
             rhs_in(i) -= fe_iv.shape_value(0, i, q) * flux * fe_iv.JxW(q);

           // Neighbor (Out): + phi * Flux
           for (unsigned int i = 0; i < dofs_per_cell; ++i)
             rhs_out(i) += fe_iv.shape_value(1, i, q) * flux * fe_iv.JxW(q);
         }

         // Add to global vector
         cell->get_dof_indices(local_dof_indices);
         rhs_vector.add(local_dof_indices, rhs_in);

         neighbor->get_dof_indices(local_dof_indices);
         rhs_vector.add(local_dof_indices, rhs_out);
      }
    }
  }
}

template <int dim>
void
DGAdvection<dim>::apply_inverse_mass(Vector<double>       &dst,
                                     const Vector<double> &src)
{
  dst = 0;
  Vector<double>                       src_loc(fe.dofs_per_cell);
  Vector<double>                       dst_loc(fe.dofs_per_cell);
  std::vector<types::global_dof_index> indices(fe.dofs_per_cell);

  unsigned int cell_index = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      cell->get_dof_indices(indices);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        src_loc[i] = src[indices[i]];

      inverse_local_mass[cell_index].vmult(dst_loc, src_loc);

      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
        dst[indices[i]] = dst_loc[i];

      cell_index++;
    }
}

template <int dim>
void
DGAdvection<dim>::output_results(unsigned int step, double time)
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "u");

  // Add exact solution
  Vector<double> exact(dof_handler.n_dofs());
  AnalyticalSolution<dim> exact_func;
  exact_func.set_time(time);
  VectorTools::interpolate(dof_handler, exact_func, exact);
  data_out.add_data_vector(exact, "u_exact");

  data_out.build_patches(mapping, fe_degree);
  std::ofstream out("output/solution-" + std::to_string(step) + ".vtk");
  data_out.write_vtk(out);
}

template <int dim>
void
DGAdvection<dim>::run()
{
  make_grid();
  setup_system();

  InitialCondition<dim> initial_condition;
  VectorTools::project(dof_handler,
                       constraints,
                       QGauss<dim>(fe_degree + 1),
                       initial_condition,
                       solution);

  // Time Stepping Setup
  const double final_time = 2.0;
  const double h          = GridTools::minimal_cell_diameter(triangulation);
  // Max speed v = 2*pi approx 6.28
  const double max_v = 2.0 * std::numbers::pi;
  // CFL 0.1
  const double dt      = 0.05 * h / (max_v * (2 * fe_degree + 1));
  const int    n_steps = static_cast<int>(final_time / dt);

  std::cout << "h = " << h << " dt = " << dt << " steps = " << n_steps
            << std::endl;

  output_results(0, 0.0);

  TimeStepping::ExplicitRungeKutta<Vector<double>> rk(
    TimeStepping::RK_CLASSIC_FOURTH_ORDER);

  double         time = 0.0;
  Vector<double> rhs(dof_handler.n_dofs());
  Vector<double> k(dof_handler.n_dofs());

  for (int step = 1; step <= n_steps; ++step)
    {
      rk.evolve_one_time_step(
        [&](const double t, const Vector<double> &y) {
          assemble_rhs(y, t, rhs);
          apply_inverse_mass(k, rhs);
          return k;
        },
        time,
        dt,
        solution);

      time += dt;

      if (step % 500 == 0 || step == n_steps)
        {
          AnalyticalSolution<dim> exact_func;
          exact_func.set_time(time);
          Vector<double> diff(triangulation.n_active_cells());
          VectorTools::integrate_difference(mapping,
                                            dof_handler,
                                            solution,
                                            exact_func,
                                            diff,
                                            QGauss<dim>(fe_degree + 1),
                                            VectorTools::L2_norm);
          const double err = VectorTools::compute_global_error(
            triangulation, diff, VectorTools::L2_norm);

          std::cout << "Step " << step << " t=" << time << " L2=" << err
                    << std::endl;
        }

      if (step % 1000 == 0)
        output_results(step, time);
    }
}

int
main()
{

  try
    {
      DGAdvection<1> test(1, 6);
      test.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << exc.what() << std::endl;
      return 1;
    }
  return 0;
}