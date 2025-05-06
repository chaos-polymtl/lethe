// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief PoissonOperator: test on rotated hypershell-hypershell geometry.
 */

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

#include "./tests.h"

using namespace dealii;


template <int dim>
class PoissonMatrixBased
{
public:
  PoissonMatrixBased(const unsigned int fe_degree,
                     const unsigned int mapping_degree,
                     const unsigned int n_global_refinements,
                     const double       radius)
    : fe_degree(fe_degree)
    , mapping_degree(mapping_degree)
    , n_global_refinements(n_global_refinements)
    , radius(radius)
    , comm(MPI_COMM_WORLD)
    , tria(comm)
    , dof_handler(tria)
    , fe(fe_degree)
    , mapping(mapping_degree)
    , quadrature(fe_degree + 1)
  {}

  void
  generate_grid()
  {
    // generate merged grid
    hyper_shell_with_hyper_shell(radius, tria);
    tria.refine_global(n_global_refinements);
  }

  void
  setup_dofs()
  {
    dof_handler.distribute_dofs(fe);

    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    constraints.reinit(dof_handler.locally_owned_dofs(), locally_relevant_dofs);
    DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
    DoFTools::make_zero_boundary_constraints(dof_handler, 3, constraints);
    // DoFTools::make_zero_boundary_constraints(dof_handler, 1, constraints);
    // DoFTools::make_zero_boundary_constraints(dof_handler, 2, constraints);



    coupling = std::make_shared<CouplingOperator<dim, 1, double>>(
      mapping,
      dof_handler,
      constraints,
      quadrature,
      4 * Utilities::pow(2, n_global_refinements + 1),
      0.5 * radius,
      0,
      1,
      2);



    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);

    coupling->add_sparsity_pattern_entries(dsp);

    constraints.close();


    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }

  void
  assemble_matrix_and_rhs()
  {
    FEValues<dim> fe_values(fe,
                            quadrature,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        cell_matrix = 0;
        cell_rhs    = 0;

        for (const unsigned int q_index : fe_values.quadrature_point_indices())
          {
            for (const unsigned int i : fe_values.dof_indices())
              {
                for (const unsigned int j : fe_values.dof_indices())
                  cell_matrix(i, j) += (               // a(x_q)
                    fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                    fe_values.JxW(q_index));           // dx

                cell_rhs(i) +=
                  (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                   1.0 *                               // f(x)
                   fe_values.JxW(q_index));            // dx
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }
  }

  void
  solve()
  {
    SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
    SolverCG<Vector<double>> solver(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);
  }

  void
  output_results()
  {
    DataOut<dim> data_out;

    DataOutBase::VtkFlags flags;
    flags.write_higher_order_cells = true;
    data_out.set_flags(flags);

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");

    Vector<double> ranks(tria.n_active_cells());
    // ranks = Utilities::MPI::this_mpi_process();
    // data_out.add_data_vector(ranks, "ranks");
    data_out.build_patches(mapping,
                           fe_degree + 1,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);
    data_out.write_vtu_in_parallel("poisson_dg.vtu", MPI_COMM_WORLD);
  }


private:
  const unsigned int fe_degree;
  const unsigned int mapping_degree;
  const unsigned int n_global_refinements;
  const double       radius;
  const MPI_Comm     comm;

  parallel::distributed::Triangulation<dim> tria;
  DoFHandler<dim>                           dof_handler;

  FE_Q<dim>     fe;
  MappingQ<dim> mapping;
  QGauss<dim>   quadrature;

  AffineConstraints<double> constraints;

  SparseMatrix<double>                              system_matrix;
  SparsityPattern                                   sparsity_pattern;
  Vector<double>                                    solution;
  Vector<double>                                    system_rhs;
  std::shared_ptr<CouplingOperator<dim, 1, double>> coupling;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  PoissonMatrixBased<2> problem(1, 1, 5, 1);
  problem.generate_grid();
  problem.setup_dofs();
  problem.assemble_matrix_and_rhs();
  problem.solve();
  problem.output_results();
}
