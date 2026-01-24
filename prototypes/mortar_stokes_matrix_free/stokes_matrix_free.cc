// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief StokesOperator: test on rotated hypercube-cylinder geometry.
 */

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/dofs/dof_tools.h>

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

#include "../mortar_functions/generate_grids.h"
#include "../mortar_functions/mortar_operator.h"

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
  const unsigned int n_global_refinements = 2;
  const double       radius               = 0.75;
  const double       outer_radius         = 1.0;
  const double       rotate               = 3.0;
  const double       rotate_pi            = 2 * numbers::PI * rotate / 360.0;
  const bool         rotate_triangulation = true;
  const MPI_Comm     comm                 = MPI_COMM_WORLD;
  const std::string  grid                 = "hyper_cube_with_cylindrical_hole";
  const double       delta_1_scaling      = 0.0001;
  const double       sip_factor           = 10.0;

  try
    {
      ConditionalOStream pcout(std::cout,
                               Utilities::MPI::this_mpi_process(comm) == 0);

      FESystem<dim> fe(FE_Q<dim>(fe_degree), dim + 1);
      MappingQ<dim> mapping_q(mapping_degree);
      QGauss<dim>   quadrature(fe_degree + 1);

      parallel::distributed::Triangulation<dim> tria(comm);
      if (grid == "hyper_cube_with_cylindrical_hole")
        hyper_cube_with_cylindrical_hole(radius,
                                         outer_radius,
                                         rotate_triangulation ? rotate_pi : 0.0,
                                         tria);
      else if (grid == "hyper_cube")
        GridGenerator::hyper_cube(tria, -outer_radius, +outer_radius);
      else if (grid == "hyper_cube_with_cylindrical_hole_with_tolerance")
        hyper_cube_with_cylindrical_hole_with_tolerance(radius,
                                                        outer_radius,
                                                        0.0,
                                                        tria);
      else
        AssertThrow(false, ExcNotImplemented());
      tria.refine_global(n_global_refinements);

      MappingQCache<dim> mapping(mapping_degree);

      if (grid != "hyper_cube_with_cylindrical_hole" || rotate_triangulation)
        mapping.initialize(mapping_q, tria);
      else
        mapping.initialize(
          mapping_q,
          tria,
          [&](const auto &cell, const auto &point) {
            if (cell->center().norm() > radius)
              return point;

            return static_cast<Point<dim>>(
              Physics::Transformations::Rotations::rotation_matrix_2d(
                rotate_pi) *
              point);
          },
          false);

      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      AffineConstraints<double> constraints;
      const IndexSet            locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         locally_relevant_dofs);

      if (grid == "hyper_cube" ||
          grid == "hyper_cube_with_cylindrical_hole_with_tolerance")
        {
          DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
        }
      else if (grid == "hyper_cube_with_cylindrical_hole")
        {
          DoFTools::make_zero_boundary_constraints(dof_handler, 1, constraints);
          DoFTools::make_zero_boundary_constraints(dof_handler, 2, constraints);
          DoFTools::make_zero_boundary_constraints(dof_handler, 3, constraints);
          DoFTools::make_zero_boundary_constraints(dof_handler, 4, constraints);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }
      constraints.close();

      StokesOperator<dim, double> op(
        mapping, dof_handler, constraints, quadrature, delta_1_scaling);

      if (grid == "hyper_cube_with_cylindrical_hole")
        {
          const auto mortar_manager =
            std::make_shared<MortarManagerCircle<dim>>(
              4 * Utilities::pow(2, n_global_refinements + 1),
              radius,
              construct_quadrature(quadrature),
              rotate_pi);

          op.add_coupling(mortar_manager, 0, 5, sip_factor);
        }

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

          FEValues<dim>              fe_values(mapping,
                                  fe,
                                  quadrature,
                                  update_values | update_gradients |
                                    update_JxW_values |
                                    update_quadrature_points);
          FEValuesViews::Vector<dim> velocities(fe_values, 0);
          FEValuesViews::Scalar<dim> pressure(fe_values, dim);

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);

                const double delta_1 =
                  delta_1_scaling * cell->minimum_vertex_distance();

                Vector<double> rhs_local(fe.n_dofs_per_cell());
                std::vector<types::global_dof_index> indices(
                  fe.n_dofs_per_cell());

                cell->get_dof_indices(indices);

                for (const unsigned int q :
                     fe_values.quadrature_point_indices())
                  {
                    const auto JxW   = fe_values.JxW(q);
                    const auto point = fe_values.quadrature_point(q);

                    Tensor<1, dim> source;
                    for (unsigned int d = 0; d < dim; ++d)
                      source[d] = rhs_func->value(point, d);

                    for (const unsigned int i : fe_values.dof_indices())
                      rhs_local(i) +=
                        delta_1 * source * pressure.gradient(i, q) * JxW;
                  }

                constraints.distribute_local_to_global(rhs_local, indices, rhs);
              }
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
          DiagonalMatrix<LinearAlgebra::distributed::Vector<double>>
            preconditioner;
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
              dim + 1,
              DataComponentInterpretation::component_is_part_of_vector);
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
          data_out.write_vtu_in_parallel("out.vtu", MPI_COMM_WORLD);
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
