// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief PoissonOperator: test on rotated hypercube-cylinder geometry.
 */

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
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

  try
    {
      const unsigned int fe_degree            = 3;
      const unsigned int mapping_degree       = 3;
      const unsigned int dim                  = 3;
      const unsigned int n_global_refinements = 1;
      const double       radius               = 1.0;
      const double       rotate               = 3.0;
      const double       rotate_pi = 2 * numbers::PI * rotate / 360.0;
      const bool         rotate_triangulation = false;
      const MPI_Comm     comm                 = MPI_COMM_WORLD;

      ConditionalOStream pcout(std::cout,
                               Utilities::MPI::this_mpi_process(comm) == 0);

      FE_Q<dim>     fe(fe_degree);
      MappingQ<dim> mapping_q(mapping_degree);
      QGauss<dim>   quadrature(fe_degree + 1);

      // generate merged grid
      parallel::distributed::Triangulation<dim> tria(comm);
      hyper_cube_with_cylindrical_hole(radius,
                                       2.0,
                                       rotate_triangulation ? rotate_pi : 0.0,
                                       tria);
      tria.refine_global(n_global_refinements);

      MappingQCache<dim> mapping(mapping_degree);

      if (rotate_triangulation)
        mapping.initialize(mapping_q, tria);
      else
        mapping.initialize(
          mapping_q,
          tria,
          [&](const auto &cell, const auto &point) {
            auto temp = cell->center();
            temp[2]   = 0;

            if (temp.norm() > radius)
              return point;

            Tensor<1, 3, Number> axis;
            axis[2] = 1.;

            return static_cast<Point<dim>>(
              Physics::Transformations::Rotations::rotation_matrix_3d(
                axis, rotate_pi) *
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
      DoFTools::make_zero_boundary_constraints(dof_handler, 1, constraints);
      DoFTools::make_zero_boundary_constraints(dof_handler, 2, constraints);
      DoFTools::make_zero_boundary_constraints(dof_handler, 3, constraints);
      DoFTools::make_zero_boundary_constraints(dof_handler, 4, constraints);
      constraints.close();

      PoissonOperator<dim, 1, double> op(mapping,
                                         dof_handler,
                                         constraints,
                                         quadrature);

      const auto mortar_manager = std::make_shared<MortarManagerCircle<dim>>(
        std::vector<unsigned int>{8 * Utilities::pow(2, n_global_refinements),
                                  1 * Utilities::pow(2, n_global_refinements)},
        std::vector<double>{radius, radius},
        construct_quadrature(quadrature),
        rotate_pi);

      op.add_coupling(mortar_manager, 0, 5);

      LinearAlgebra::distributed::Vector<double> rhs, solution;
      op.initialize_dof_vector(rhs);
      op.initialize_dof_vector(solution);

      VectorTools::create_right_hand_side(mapping,
                                          dof_handler,
                                          quadrature,
                                          Functions::ConstantFunction<dim>(1.0),
                                          rhs,
                                          constraints);

      ReductionControl        reduction_control(10000, 1e-20, 1e-6);
      SolverGMRES<VectorType> solver(reduction_control);

      if (false)
        {
          // 0) no preconditioner
          PreconditionIdentity preconditioner;
          solution = 0.0;
          solver.solve(op, solution, rhs, preconditioner);
          pcout << reduction_control.last_step() << std::endl;
        }
      if (false)
        {
          // 1) with preconditioner: inverse diagonal
          DiagonalMatrix<LinearAlgebra::distributed::Vector<double>>
            preconditioner;
          op.compute_inverse_diagonal(preconditioner.get_vector());
          solution = 0.0;
          solver.solve(op, solution, rhs, preconditioner);
          pcout << reduction_control.last_step() << std::endl;
        }
      if (false)
        {
          // 2) with preconditioner: algebraic multigrid
          TrilinosWrappers::PreconditionILU preconditioner;
          preconditioner.initialize(op.get_system_matrix());
          solution = 0.0;
          solver.solve(op, solution, rhs, preconditioner);
          pcout << reduction_control.last_step() << std::endl;
        }
      if (true)
        {
          // 3) with preconditioner: direct solver
          TrilinosWrappers::SolverDirect preconditioner;
          preconditioner.initialize(op.get_system_matrix());
          solution = 0.0;
          solver.solve(op, solution, rhs, preconditioner);
          pcout << reduction_control.last_step() << std::endl;
        }

      DataOut<dim> data_out;

      DataOutBase::VtkFlags flags;
      flags.write_higher_order_cells = true;
      data_out.set_flags(flags);

      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");

      Vector<double> ranks(tria.n_active_cells());
      ranks = Utilities::MPI::this_mpi_process(comm);
      data_out.add_data_vector(ranks, "ranks");

      data_out.build_patches(
        mapping,
        fe_degree + 1,
        DataOut<dim>::CurvedCellRegion::curved_inner_cells);
      data_out.write_vtu_in_parallel("out.vtu", MPI_COMM_WORLD);
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
