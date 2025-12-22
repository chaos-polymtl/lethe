// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief PoissonOperator: test matrix-based assembly for a nonlinear Poisson problem of the form
 * -∇.(∇u) = exp(u). Based on prototypes/matrix_based_non_linear_poisson */

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
  }

  void
  assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;

    FEValues<dim> fe_values(fe,
                            quadrature,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    // dofs per cell
    const unsigned int                   dofs_per_cell = fe.n_dofs_per_cell();
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // number of quadrature points
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

                    // +(∇v, ∇u) - (v, exp(u))
                    cell_rhs(i) += (grad_phi_i * previous_gradients[q] -
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

    // add coupling entries in system matrix and RHS
    mortar_coupling_operator->add_system_matrix_entries(system_matrix);
    mortar_coupling_operator->vmult_add(system_rhs, solution);

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    system_rhs *= -1.0;
  }

  void
  compute_update()
  {
    // initialize GMRES solver
    ReductionControl reduction_control(10000, 1e-10, 1e-6);
    SolverGMRES<TrilinosWrappers::MPI::Vector> solver(reduction_control);

    // initialize ILU preconditioner
    TrilinosWrappers::PreconditionILU preconditioner;
    preconditioner.initialize(system_matrix);

    TrilinosWrappers::MPI::Vector completely_distributed_solution(
      dof_handler.locally_owned_dofs(), comm);

    // solve linear system
    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 preconditioner);

    // update solution
    constraints.distribute(completely_distributed_solution);
    delta_solution = completely_distributed_solution;

    TrilinosWrappers::MPI::Vector local_solution(system_rhs);
    local_solution = solution;
    local_solution += completely_distributed_solution;
    solution = local_solution;
  }

  void
  solve_non_linear()
  {
    // conditional output stream for MPI run
    ConditionalOStream pcout(std::cout,
                             Utilities::MPI::this_mpi_process(comm) == 0);

    // iteration parameters
    double       error = 1e10;
    double       tol   = 1e-10;
    unsigned int it    = 0;

    // output results on first iteration
    output_results(it);

    // iteration loop
    while (error > tol)
      {
        // update iteration
        it++;

        // assemble and solve linear system
        assemble_system();
        compute_update();

        // get vectors without any ghosts
        TrilinosWrappers::MPI::Vector local_delta_solution(system_rhs);
        local_delta_solution = delta_solution;

        // calculate error
        error = local_delta_solution.l2_norm();

        // print norms
        pcout << "   Iter " << it << " - Delta solution norm, Linfty norm: "
              << local_delta_solution.linfty_norm()
              << " L2 norm: " << local_delta_solution.l2_norm() << std::endl;

        // output iteration results
        output_results(it);
      }
  }

  void
  output_results(unsigned int it)
  {
    DataOut<dim> data_out;
    std::string  filename =
      "poisson_nonlinear_MB-" + Utilities::int_to_string(it) + ".vtu";

    solution.update_ghost_values();

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
  TrilinosWrappers::MPI::Vector                  solution;
  TrilinosWrappers::MPI::Vector                  delta_solution;
  TrilinosWrappers::MPI::Vector                  system_rhs;
  std::shared_ptr<CouplingOperator<dim, double>> mortar_coupling_operator;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  PoissonMatrixBased<2> problem(1, 1, 3, 1);

  // generate grid and setup dofs
  problem.generate_grid();
  problem.setup_system();

  // solve nonlinear problem
  problem.solve_non_linear();
}
