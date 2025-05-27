// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
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
  setup_system()
  {
    // distribute dofs
    dof_handler.distribute_dofs(fe);

    // get dofs sets
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

    // setup constraints
    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
    DoFTools::make_zero_boundary_constraints(dof_handler, 3, constraints);
    constraints.close();

    // create mortar manager
    const auto mortar_manager = std::make_shared<MortarManagerCircle<dim>>(
      6 * Utilities::pow(2, n_global_refinements),
      quadrature,
      0.5 * radius,
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
                                                      1000.0);

    // create sparsity pattern
    DynamicSparsityPattern dsp(locally_relevant_dofs);

    DoFTools::make_sparsity_pattern(dof_handler,
                                    dsp,
                                    constraints,
                                    /*keep_constrained_dofs = */ true);

    // add coupling entries in sparsity pattern
    mortar_coupling_operator->add_sparsity_pattern_entries(dsp);
    constraints.close();
    sparsity_pattern.copy_from(dsp);

    // distribute sparsity pattern in MPI
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               locally_owned_dofs,
                                               comm,
                                               locally_relevant_dofs);

    // initialze matrices
    system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, comm);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, comm);
    delta_solution.reinit(locally_owned_dofs, locally_relevant_dofs, comm);
    system_rhs.reinit(locally_owned_dofs, comm);
    previous_solution = solution;

    // apply BCs to solution vector for first iteration
    constraints.distribute(solution);
  }

  void
  assemble_system()
  {
    /* Matrix assembly for the nonlinear Poisson equation of the form -∇.(exp(u)
     * ∇(u)) = f */

    system_matrix = 0;
    system_rhs    = 0;

    FEValues<dim> fe_values(fe,
                            quadrature,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    // dofs per cell
    const unsigned int                   dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // quadrature points
    const unsigned int n_q_points = fe_values.n_quadrature_points;

    // initialize cell matrix and RHS
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    // initialize vectors for previous solution
    std::vector<double>         previous_values(n_q_points);
    std::vector<Tensor<1, dim>> previous_gradients(n_q_points);

    // loop over locally owned cells
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->subdomain_id() == this_mpi_process)
          {
            cell_matrix = 0.;
            cell_rhs    = 0.;

            fe_values.reinit(cell);

            // get previous values and gradients
            fe_values.get_function_values(solution, previous_values);
            fe_values.get_function_gradients(solution, previous_gradients);

            // loop over quadrature points
            for (unsigned int q = 0; q < n_q_points; q++)
              {
                // exp(u)
                const double nonlinearity = std::exp(previous_values[q]);
                // Jacobian determinant
                const double dx = fe_values.JxW(q);

                for (const unsigned int i : fe_values.dof_indices())
                  {
                    // test function values and gradients
                    const double         phi_i = fe_values.shape_value(i, q);
                    const Tensor<1, dim> grad_phi_i =
                      fe_values.shape_grad(i, q);

                    for (const unsigned int j : fe_values.dof_indices())
                      {
                        // shape function values and gradients
                        const double phi_j = fe_values.shape_value(j, q);
                        const Tensor<1, dim> grad_phi_j =
                          fe_values.shape_grad(j, q);

                        // (∇v, ∇δu) - (v, exp(u)*δu)
                        cell_matrix(i, j) += (grad_phi_i * grad_phi_j -
                                              phi_i * nonlinearity * phi_j) *
                                             dx;
                      }

                    // -(∇v, ∇u) + (v, exp(u))
                    cell_rhs(i) += (-grad_phi_i * previous_gradients[q] +
                                    phi_i * nonlinearity) *
                                   dx;
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
    mortar_coupling_operator->add_system_rhs_entries(system_rhs);

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  double
  compute_residual()
  {
    TrilinosWrappers::MPI::Vector residual;
    TrilinosWrappers::MPI::Vector evaluation_point(system_rhs);
    TrilinosWrappers::MPI::Vector local_newton_update(system_rhs);
    local_newton_update = delta_solution;
    TrilinosWrappers::MPI::Vector local_evaluation_point;

    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();

    residual.reinit(locally_owned_dofs, comm);
    local_evaluation_point.reinit(locally_owned_dofs,
                                  locally_relevant_dofs,
                                  comm);

    evaluation_point = solution;

    local_evaluation_point = solution;

    FEValues<dim> fe_values(fe,
                            quadrature,
                            update_values | update_gradients |
                              update_JxW_values | update_quadrature_points);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
    const unsigned int n_q_points    = fe_values.n_quadrature_points;

    Vector<double> cell_residual(dofs_per_cell);

    std::vector<double>         values(n_q_points);
    std::vector<Tensor<1, dim>> gradients(n_q_points);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            cell_residual = 0.0;
            fe_values.reinit(cell);

            fe_values.get_function_values(local_evaluation_point, values);
            fe_values.get_function_gradients(local_evaluation_point, gradients);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double nonlinearity = std::exp(values[q]);
                const double dx           = fe_values.JxW(q);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const double         phi_i = fe_values.shape_value(i, q);
                    const Tensor<1, dim> grad_phi_i =
                      fe_values.shape_grad(i, q);

                    cell_residual(i) +=
                      (grad_phi_i * gradients[q] - phi_i * nonlinearity) * dx;
                  }
              }

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_residual,
                                                   local_dof_indices,
                                                   residual);
          }
      }

    residual.compress(VectorOperation::add);
    residual.update_ghost_values();

    return residual.l2_norm();
  }

  void
  solve_linear()
  {
    // initialize GMRES solver
    ReductionControl reduction_control(10000, 1e-10, 1e-6);
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(reduction_control);

    // initialize ILU preconditioner
    TrilinosWrappers::PreconditionILU preconditioner;
    preconditioner.initialize(system_matrix);

    // solve linear system
    solver.solve(system_matrix, delta_solution, system_rhs, preconditioner);
  }

  void
  solve_non_linear()
  {
    // conditional output stream for MPI run
    ConditionalOStream pcout(std::cout,
                             Utilities::MPI::this_mpi_process(comm) == 0);

    // iteration parameters
    double       error = 1e6;
    double       tol   = 1e-6;
    unsigned int iter  = 0;

    // output results on first iteration
    output_results(iter);

    // iteration loop
    while (error > tol)
      {
        // update iteration
        iter++;

        // assemble and solve system
        assemble_system();
        solve_linear();

        // update solution
        solution += delta_solution;
        constraints.distribute(solution);

        // calculate error
        pcout << "   Iter " << iter << " - Delta solution norm, Linfty norm: "
              << delta_solution.linfty_norm()
              << " L2 norm: " << delta_solution.l2_norm() << std::endl;

        error = delta_solution.linfty_norm();

        // output iteration results
        output_results(iter);

        // update solution
        previous_solution = solution;
      }
  }

  void
  output_results(unsigned int iter)
  {
    DataOut<dim> data_out;
    std::string  filename =
      "poisson_nonlinear_MB-" + Utilities::int_to_string(iter) + ".vtu";

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
    data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
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
  SparsityPattern                                sparsity_pattern;
  TrilinosWrappers::MPI::Vector                  solution;
  TrilinosWrappers::MPI::Vector                  delta_solution;
  TrilinosWrappers::MPI::Vector                  previous_solution;
  TrilinosWrappers::MPI::Vector                  system_rhs;
  std::shared_ptr<CouplingOperator<dim, double>> mortar_coupling_operator;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  PoissonMatrixBased<2> problem(1, 1, 4, 1);

  // generate grid and setup dofs
  problem.generate_grid();
  problem.setup_system();

  // solve nonlinear problem
  problem.solve_non_linear();
}
