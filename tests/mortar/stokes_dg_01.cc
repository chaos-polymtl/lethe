// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief StokesOperator: test on rotated hypercube-cylinder geometry.
 */

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/physics/transformations.h>

#include <fstream>

#include "./tests.h"

using namespace dealii;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int fe_degree            = 3;
  const unsigned int mapping_degree       = 3;
  const unsigned int dim                  = 2;
  const unsigned int n_global_refinements = 3;
  const double       radius               = 0.75;
  const double       outer_radius         = 1.0;
  const double       rotate               = 3.0;
  const double       rotate_pi            = 2 * numbers::PI * rotate / 360.0;
  const bool         rotate_triangulation = true;
  const MPI_Comm     comm                 = MPI_COMM_WORLD;
  const std::string  grid                 = "hyper_cube";
  const double       sip_factor           = 1.0;

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(comm) == 0);

  FESystem<dim> fe(FE_DGQ<dim>(fe_degree), dim, FE_DGQ<dim>(fe_degree - 1), 1);
  MappingQ<dim> mapping_q(mapping_degree);
  QGauss<dim>   quadrature(fe_degree + 1);

  parallel::distributed::Triangulation<dim> tria(comm);
  if (grid == "hyper_cube")
    GridGenerator::hyper_cube(tria, -outer_radius, +outer_radius);
  else
    AssertThrow(false, ExcNotImplemented());

  tria.refine_global(n_global_refinements);

  MappingQCache<dim> mapping(mapping_degree);

  mapping.initialize(mapping_q, tria);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  const IndexSet            locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);
  constraints.reinit(dof_handler.locally_owned_dofs(), locally_relevant_dofs);

  if (true /*TODO: better solution!*/)
    {
      unsigned int min_index = numbers::invalid_unsigned_int;

      std::vector<types::global_dof_index> dof_indices;

      // Loop over the cells to identify the min index
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              const auto &fe = cell->get_fe();

              dof_indices.resize(fe.n_dofs_per_cell());
              cell->get_dof_indices(dof_indices);

              for (unsigned int i = 0; i < dof_indices.size(); ++i)
                if (fe.system_to_component_index(i).first == dim)
                  min_index = std::min(min_index, dof_indices[i]);
            }
        }

      // Necessary to find the min across all cores.
      min_index =
        Utilities::MPI::min(min_index, dof_handler.get_communicator());

      if (locally_relevant_dofs.is_element(min_index))
        constraints.add_line(min_index);

      std::cout << min_index << std::endl;
    }


  constraints.close();

  GeneralStokesOperatorDG<dim, double> op(mapping,
                                          dof_handler,
                                          constraints,
                                          quadrature);

  LinearAlgebra::distributed::Vector<double> rhs, solution;
  op.initialize_dof_vector(rhs);
  op.initialize_dof_vector(solution);

  std::shared_ptr<Function<dim, Number>> exact_solution;

  exact_solution = std::make_shared<FunctionFromFunctionObjects<dim>>(
    [&](const auto &p, const auto c) {
      const double a = numbers::PI;
      const double x = p[0];
      const double y = p[1];

      if (c == 0)
        return std::sin(a * x) * std::sin(a * x) * std::cos(a * y) *
               std::sin(a * y);
      else if (c == 1)
        return -std::cos(a * x) * std::sin(a * x) * std::sin(a * y) *
               std::sin(a * y);
      else if (c == 2)
        return std::sin(a * x) * std::sin(a * y);

      AssertThrow(false, ExcNotImplemented());

      return 0.0;
    },
    dim + 1);

  if (true)
    {
      std::shared_ptr<Function<dim, Number>> rhs_func;

      rhs_func = std::make_shared<FunctionFromFunctionObjects<dim>>(
        [&](const auto &p, const auto c) {
          const double a = numbers::PI;
          const double x = p[0];
          const double y = p[1];

          if (c == 0)
            return 2 * a * a *
                     (std::sin(a * x) * std::sin(a * x) -
                      std::cos(a * x) * std::cos(a * x)) *
                     std::sin(a * y) * std::cos(a * y) +
                   4 * a * a * std::sin(a * x) * std::sin(a * x) *
                     std::sin(a * y) * std::cos(a * y) +
                   a * std::sin(a * y) * std::cos(a * x);
          else if (c == 1)
            return -2 * a * a *
                     (std::sin(a * y) * std::sin(a * y) -
                      std::cos(a * y) * std::cos(a * y)) *
                     std::sin(a * x) * std::cos(a * x) -
                   4 * a * a * std::sin(a * x) * std::sin(a * y) *
                     std::sin(a * y) * std::cos(a * x) +
                   a * std::sin(a * x) * std::cos(a * y);
          else if (c == 2)
            return 0.0;

          AssertThrow(false, ExcNotImplemented());

          return 0.0;
        },
        dim + 1);

      VectorTools::create_right_hand_side(
        mapping, dof_handler, quadrature, *rhs_func, rhs, constraints);
    }

  rhs.compress(VectorOperation::add);

  ReductionControl        reduction_control(10000, 1e-20, 1e-6);
  SolverGMRES<VectorType> solver(reduction_control);

  if (false)
    {
      // 0) no preconditioner
      PreconditionIdentity preconditioner;
      solution = 0.0;
      solver.solve(op, solution, rhs, preconditioner);
      pcout << reduction_control.last_step();
    }
  if (false)
    {
      // 1) with preconditioner: inverse diagonal
      DiagonalMatrix<LinearAlgebra::distributed::Vector<double>> preconditioner;
      op.compute_inverse_diagonal(preconditioner.get_vector());
      solution = 0.0;
      solver.solve(op, solution, rhs, preconditioner);
      pcout << reduction_control.last_step();
    }
  if (false)
    {
      // 2) with preconditioner: algebraic multigrid
      TrilinosWrappers::PreconditionILU preconditioner;
      preconditioner.initialize(op.get_system_matrix());
      solution = 0.0;
      solver.solve(op, solution, rhs, preconditioner);
      pcout << reduction_control.last_step();
    }
  if (true)
    {
      // 3) with preconditioner: direct solver
      TrilinosWrappers::SolverDirect preconditioner;
      preconditioner.initialize(op.get_system_matrix());
      solution = 0.0;
      solver.solve(op, solution, rhs, preconditioner);
      pcout << reduction_control.last_step();
    }

  if (true)
    {
      DataOut<dim> data_out;

      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);

      std::vector<std::string> labels(dim + 1, "u");
      labels[dim] = "p";

      std::vector<std::string> labels_ana(dim + 1, "ana_u");
      labels_ana[dim] = "ana_p";

      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretation(
          dim + 1, DataComponentInterpretation::component_is_part_of_vector);
      data_component_interpretation[dim] =
        DataComponentInterpretation::component_is_scalar;

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(dof_handler,
                               solution,
                               labels,
                               data_component_interpretation);

      Vector<double> ranks(tria.n_active_cells());
      ranks = Utilities::MPI::this_mpi_process(comm);
      data_out.add_data_vector(ranks, "ranks");

      LinearAlgebra::distributed::Vector<double> analytical_solution;
      op.initialize_dof_vector(analytical_solution);
      VectorTools::interpolate(mapping,
                               dof_handler,
                               *exact_solution,
                               analytical_solution);
      data_out.add_data_vector(dof_handler,
                               analytical_solution,
                               labels_ana,
                               data_component_interpretation);

      data_out.build_patches(
        mapping,
        fe_degree + 1,
        DataOut<dim>::CurvedCellRegion::curved_inner_cells);
      data_out.write_vtu_in_parallel("poisson_dg.vtu", MPI_COMM_WORLD);
    }

  if (true)
    {
      solution.update_ghost_values();

      const ComponentSelectFunction<dim> u_mask(std::make_pair(0, dim),
                                                dim + 1);
      const ComponentSelectFunction<dim> p_mask(dim, dim + 1);

      Vector<float> norm_per_cell(tria.n_active_cells());
      VectorTools::integrate_difference(dof_handler,
                                        solution,
                                        *exact_solution,
                                        norm_per_cell,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::L2_norm,
                                        &u_mask);
      const double error_L2_norm_u =
        VectorTools::compute_global_error(tria,
                                          norm_per_cell,
                                          VectorTools::L2_norm);

      VectorTools::integrate_difference(dof_handler,
                                        solution,
                                        *exact_solution,
                                        norm_per_cell,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::L2_norm,
                                        &p_mask);
      const double error_L2_norm_p =
        VectorTools::compute_global_error(tria,
                                          norm_per_cell,
                                          VectorTools::L2_norm);

      pcout << " " << error_L2_norm_u << " " << error_L2_norm_p;
    }

  pcout << std::endl;
}
