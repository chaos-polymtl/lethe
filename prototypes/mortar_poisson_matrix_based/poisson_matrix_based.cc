// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief PoissonOperator: test matrix-based assembly.
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

#include "./generate_grids.h"
#include "./mortar_operator.h"

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
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(comm))
    , this_mpi_process(Utilities::MPI::this_mpi_process(comm))
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
    // global mesh refinement
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
    constraints.close();

    // create mortar manager
    const auto mortar_manager = std::make_shared<MortarManagerCircle<dim>>(
      6 * Utilities::pow(2, n_global_refinements),
      0.5 * radius,
      quadrature,
      0.0);

    // create coupling evaluator
    const std::shared_ptr<CouplingEvaluationBase<dim, double>>
      mortar_coupling_evaluator =
        std::make_shared<CouplingEvaluationSIPG<dim, 1, double>>(mapping,
                                                                 dof_handler);

    // create coupling operator
    mortar_coupling_operator =
      std::make_shared<CouplingOperator<dim, double>>(mapping,
                                                      dof_handler,
                                                      constraints,
                                                      mortar_coupling_evaluator,
                                                      mortar_manager,
                                                      1,
                                                      2,
                                                      1.0);

    // create sparsity pattern
    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);

    // add coupling entries in sparsity pattern
    mortar_coupling_operator->add_sparsity_pattern_entries(dsp);
    constraints.close();

    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               comm,
                                               locally_relevant_dofs);

    system_matrix.reinit(dof_handler.locally_owned_dofs(),
                         dof_handler.locally_owned_dofs(),
                         dsp,
                         comm);
    solution.reinit(dof_handler.locally_owned_dofs(), comm);
    system_rhs.reinit(dof_handler.locally_owned_dofs(), comm);
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
        if (cell->subdomain_id() == this_mpi_process)
          {
            fe_values.reinit(cell);

            cell_matrix = 0;
            cell_rhs    = 0;

            for (const unsigned int q_index :
                 fe_values.quadrature_point_indices())
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
            constraints.distribute_local_to_global(cell_matrix,
                                                   cell_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
          }
      }
    // add coupling entries in stiffness matrix
    mortar_coupling_operator->add_system_matrix_entries(system_matrix);

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  void
  solve()
  {
    ConditionalOStream pcout(std::cout,
                             Utilities::MPI::this_mpi_process(comm) == 0);

    ReductionControl reduction_control(10000, 1e-20, 1e-6);
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(reduction_control);

    TrilinosWrappers::PreconditionILU preconditioner;
    preconditioner.initialize(system_matrix);

    solver.solve(system_matrix, solution, system_rhs, preconditioner);

    constraints.distribute(solution);

    pcout << reduction_control.last_step() << std::endl;
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
    ranks = Utilities::MPI::this_mpi_process(comm);
    data_out.add_data_vector(ranks, "ranks");
    data_out.build_patches(mapping,
                           fe_degree + 1,
                           DataOut<dim>::CurvedCellRegion::curved_inner_cells);
    data_out.write_vtu_in_parallel("poisson_MB.vtu", MPI_COMM_WORLD);
  }


private:
  const unsigned int fe_degree;
  const unsigned int mapping_degree;
  const unsigned int n_global_refinements;
  const double       radius;
  const MPI_Comm     comm;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;

  parallel::distributed::Triangulation<dim> tria;
  DoFHandler<dim>                           dof_handler;

  FE_Q<dim>     fe;
  MappingQ<dim> mapping;
  QGauss<dim>   quadrature;

  AffineConstraints<double> constraints;

  TrilinosWrappers::SparseMatrix                 system_matrix;
  TrilinosWrappers::MPI::Vector                  solution;
  TrilinosWrappers::MPI::Vector                  system_rhs;
  std::shared_ptr<CouplingOperator<dim, double>> mortar_coupling_operator;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);


  try
    {
      PoissonMatrixBased<2> problem(1, 1, 5, 1);
      problem.generate_grid();
      problem.setup_dofs();
      problem.assemble_matrix_and_rhs();
      problem.solve();
      problem.output_results();
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
