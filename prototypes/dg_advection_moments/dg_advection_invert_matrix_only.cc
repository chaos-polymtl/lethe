// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/* *
 * 1D DG Advection Solver (Explicit Time Integration)
 * Equation: du/dt + 2*pi * du/dx = 0
 * Domain: x in [0, 2*pi]
 * Method: Explicit RK3 with precomputed Inverse Mass Matrix
 */

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// --- Boundary and Initial Conditions ---
template <int dim>
class InflowCondition : public Function<dim>
{
public:
  InflowCondition()
    : Function<dim>()
  {}
  virtual double
  value(const Point<dim> & /*p*/, const unsigned int /*c*/ = 0) const override
  {
    return -std::sin(2 * numbers::PI * this->get_time());
  }
};

template <int dim>
class InitialValues : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int /*c*/ = 0) const override
  {
    return std::sin(p[0]);
  }
};

// --- MeshWorker Data Structures ---
template <int dim>
struct ScratchData
{
  ScratchData(const FiniteElement<dim>  &fe,
              const Quadrature<dim>     &quadrature,
              const Quadrature<dim - 1> &face_quadrature)
    : fe_values(fe,
                quadrature,
                update_values | update_gradients | update_JxW_values |
                  update_quadrature_points)
    , fe_face_values(fe,
                     face_quadrature,
                     update_values | update_normal_vectors | update_JxW_values |
                       update_quadrature_points)
    , fe_face_values_neighbor(fe, face_quadrature, update_values)
  {}

  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                update_values | update_gradients | update_JxW_values |
                  update_quadrature_points)
    , fe_face_values(scratch_data.fe_face_values.get_fe(),
                     scratch_data.fe_face_values.get_quadrature(),
                     update_values | update_normal_vectors | update_JxW_values |
                       update_quadrature_points)
    , fe_face_values_neighbor(
        scratch_data.fe_face_values_neighbor.get_fe(),
        scratch_data.fe_face_values_neighbor.get_quadrature(),
        update_values)
  {}

  FEValues<dim>     fe_values;
  FEFaceValues<dim> fe_face_values;
  FEFaceValues<dim> fe_face_values_neighbor;
};

struct CopyData
{
  struct CellData
  {
    std::vector<types::global_dof_index> local_dof_indices;
    Vector<double>                       cell_rhs;
  };
  CellData cell_data[2];
};

// --- Main Class ---
template <int dim>
class AdvectionDG
{
public:
  AdvectionDG();
  void
  run();

private:
  void
  setup_system();
  void
  assemble_inverse_mass_matrix(); // Changed from assemble_mass_matrix
  void
  assemble_rhs();
  void
  output_results(const unsigned int step_number);

  Triangulation<dim> triangulation;
  FE_DGQ<dim>        fe;
  DoFHandler<dim>    dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> inverse_mass_matrix; // Stores M^-1 directly

  Vector<double> solution;
  Vector<double> system_rhs;
  Vector<double> solution_update;

  InflowCondition<dim> inflow_function;
  const double         velocity_val = 2.0 * numbers::PI;
};

template <int dim>
AdvectionDG<dim>::AdvectionDG()
  : fe(1)
  , dof_handler(triangulation)
{}

template <int dim>
void
AdvectionDG<dim>::setup_system()
{
  GridGenerator::hyper_cube(triangulation, 0, 2.0 * numbers::PI);
  triangulation.refine_global(6);
  dof_handler.distribute_dofs(fe);

  std::cout << "Domain: [0, 2pi] | Cells: " << triangulation.n_active_cells()
            << " | DoFs: " << dof_handler.n_dofs() << std::endl;

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  inverse_mass_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
  solution_update.reinit(dof_handler.n_dofs());
}

template <int dim>
void
AdvectionDG<dim>::assemble_inverse_mass_matrix()
{
  std::cout << "Assembling and inverting Mass Matrix..." << std::endl;

  inverse_mass_matrix = 0;
  QGauss<dim>   quadrature(fe.degree + 1);
  FEValues<dim> fe_values(fe, quadrature, update_values | update_JxW_values);

  const unsigned int                   dofs_per_cell = fe.dofs_per_cell;
  FullMatrix<double>                   cell_mass(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_mass = 0;

      // 1. Compute Local Mass Matrix
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            cell_mass(i, j) += fe_values.shape_value(i, q) *
                               fe_values.shape_value(j, q) * fe_values.JxW(q);

      // 2. Invert Local Mass Matrix explicitly
      cell_mass.gauss_jordan();

      // 3. Store M^-1 in the global sparse matrix
      // Since this is DG, DoFs are not shared, so we can just add/set.
      cell->get_dof_indices(indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          inverse_mass_matrix.add(indices[i], indices[j], cell_mass(i, j));
    }
}

template <int dim>
void
AdvectionDG<dim>::assemble_rhs()
{
  system_rhs = 0;

  const QGauss<dim>     quadrature(fe.degree + 1);
  const QGauss<dim - 1> face_quadrature(fe.degree + 1);

  auto cell_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        ScratchData<dim>                                     &scratch,
        CopyData                                             &copy) {
      scratch.fe_values.reinit(cell);
      const unsigned int n_q  = scratch.fe_values.get_quadrature().size();
      const unsigned int dofs = scratch.fe_values.get_fe().dofs_per_cell;

      copy.cell_data[0].cell_rhs.reinit(dofs);
      copy.cell_data[0].local_dof_indices.resize(dofs);
      cell->get_dof_indices(copy.cell_data[0].local_dof_indices);

      std::vector<double> current_sol(n_q);
      scratch.fe_values.get_function_values(this->solution, current_sol);

      Tensor<1, dim> velocity;
      velocity[0] = this->velocity_val;

      for (unsigned int q = 0; q < n_q; ++q)
        for (unsigned int i = 0; i < dofs; ++i)
          copy.cell_data[0].cell_rhs(i) +=
            (velocity * scratch.fe_values.shape_grad(i, q)) * current_sol[q] *
            scratch.fe_values.JxW(q);
    };

  auto boundary_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        ScratchData<dim>                                     &scratch,
        CopyData                                             &copy) {
      scratch.fe_face_values.reinit(cell, face_no);
      const unsigned int n_q  = scratch.fe_face_values.get_quadrature().size();
      const unsigned int dofs = scratch.fe_face_values.get_fe().dofs_per_cell;

      copy.cell_data[0].local_dof_indices.resize(dofs);
      cell->get_dof_indices(copy.cell_data[0].local_dof_indices);

      std::vector<double> u_inner(n_q);
      scratch.fe_face_values.get_function_values(this->solution, u_inner);

      Tensor<1, dim> velocity;
      velocity[0] = this->velocity_val;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double beta_dot_n =
            velocity * scratch.fe_face_values.normal_vector(q);
          double u_flux = (beta_dot_n < 0) ?
                            this->inflow_function.value(
                              scratch.fe_face_values.quadrature_point(q)) :
                            u_inner[q];

          for (unsigned int i = 0; i < dofs; ++i)
            copy.cell_data[0].cell_rhs(i) -=
              u_flux * beta_dot_n * scratch.fe_face_values.shape_value(i, q) *
              scratch.fe_face_values.JxW(q);
        }
    };

  auto face_worker =
    [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int                                    face_no,
        const unsigned int                                    subface_no,
        const typename DoFHandler<dim>::active_cell_iterator &neighbor,
        const unsigned int                                    neighbor_face_no,
        const unsigned int neighbor_subface_no,
        ScratchData<dim>  &scratch,
        CopyData          &copy) {
      if constexpr (dim == 1)
        {
          scratch.fe_face_values.reinit(cell, face_no);
          scratch.fe_face_values_neighbor.reinit(neighbor, neighbor_face_no);
        }
      else
        {
          scratch.fe_face_values.reinit(cell, face_no, subface_no);
          scratch.fe_face_values_neighbor.reinit(neighbor,
                                                 neighbor_face_no,
                                                 neighbor_subface_no);
        }

      const unsigned int n_q  = scratch.fe_face_values.get_quadrature().size();
      const unsigned int dofs = scratch.fe_face_values.get_fe().dofs_per_cell;

      copy.cell_data[1].cell_rhs.reinit(dofs);
      copy.cell_data[1].local_dof_indices.resize(dofs);
      neighbor->get_dof_indices(copy.cell_data[1].local_dof_indices);

      std::vector<double> u_inner(n_q), u_outer(n_q);
      scratch.fe_face_values.get_function_values(this->solution, u_inner);
      scratch.fe_face_values_neighbor.get_function_values(this->solution,
                                                          u_outer);

      Tensor<1, dim> velocity;
      velocity[0] = this->velocity_val;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          const double beta_dot_n =
            velocity * scratch.fe_face_values.normal_vector(q);
          double u_flux = (beta_dot_n > 0) ? u_inner[q] : u_outer[q];

          for (unsigned int i = 0; i < dofs; ++i)
            copy.cell_data[0].cell_rhs(i) -=
              u_flux * beta_dot_n * scratch.fe_face_values.shape_value(i, q) *
              scratch.fe_face_values.JxW(q);

          for (unsigned int i = 0; i < dofs; ++i)
            copy.cell_data[1].cell_rhs(i) +=
              u_flux * beta_dot_n *
              scratch.fe_face_values_neighbor.shape_value(i, q) *
              scratch.fe_face_values.JxW(q);
        }
    };

  auto copier = [&](const CopyData &copy) {
    for (unsigned int i = 0; i < copy.cell_data[0].local_dof_indices.size();
         ++i)
      system_rhs(copy.cell_data[0].local_dof_indices[i]) +=
        copy.cell_data[0].cell_rhs(i);

    if (copy.cell_data[1].local_dof_indices.size() > 0)
      for (unsigned int i = 0; i < copy.cell_data[1].local_dof_indices.size();
           ++i)
        system_rhs(copy.cell_data[1].local_dof_indices[i]) +=
          copy.cell_data[1].cell_rhs(i);
  };

  ScratchData<dim> scratch_data(fe, quadrature, face_quadrature);
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
AdvectionDG<dim>::run()
{
  setup_system();
  assemble_inverse_mass_matrix(); // Compute M^-1 once

  VectorTools::interpolate(dof_handler, InitialValues<dim>(), solution);
  output_results(0);

  double time     = 0;
  double end_time = 1.0;
  double dt       = 0.001;

  Vector<double> old_solution, rk_step_1, rk_step_2;
  old_solution.reinit(solution.size());
  rk_step_1.reinit(solution.size());
  rk_step_2.reinit(solution.size());

  unsigned int timestep_number = 0;
  std::cout << "Starting Explicit RK3 loop..." << std::endl;

  while (time < end_time)
    {
      timestep_number++;
      old_solution = solution;

      // RK Stage 1
      inflow_function.set_time(time);
      assemble_rhs();
      // EXPLICIT UPDATE: du/dt = M^-1 * RHS
      inverse_mass_matrix.vmult(solution_update, system_rhs);
      rk_step_1 = old_solution;
      rk_step_1.add(dt, solution_update);

      // RK Stage 2
      solution = rk_step_1;
      inflow_function.set_time(time + dt);
      assemble_rhs();
      inverse_mass_matrix.vmult(solution_update, system_rhs);
      rk_step_2 = old_solution;
      rk_step_2 *= 0.75;
      rk_step_2.add(0.25, rk_step_1);
      rk_step_2.add(0.25 * dt, solution_update);

      // RK Stage 3
      solution = rk_step_2;
      inflow_function.set_time(time + 0.5 * dt);
      assemble_rhs();
      inverse_mass_matrix.vmult(solution_update, system_rhs);
      solution = old_solution;
      solution *= (1.0 / 3.0);
      solution.add(2.0 / 3.0, rk_step_2);
      solution.add(2.0 / 3.0 * dt, solution_update);

      time += dt;

      if (timestep_number % 100 == 0)
        {
          std::cout << "Time: " << time << std::endl;
          output_results(timestep_number);
        }
    }
}

template <int dim>
void
AdvectionDG<dim>::output_results(const unsigned int step_number)
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "u");
  data_out.build_patches();
  std::string filename =
    "solution-" + Utilities::int_to_string(step_number, 5) + ".vtu";
  std::ofstream output(filename);
  data_out.write_vtu(output);
}

int
main()
{
  try
    {
      AdvectionDG<1> problem;
      problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << exc.what() << std::endl;
      return 1;
    }
  return 0;
}
