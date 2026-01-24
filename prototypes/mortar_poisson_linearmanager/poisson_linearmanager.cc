// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief PoissonOperator: test on split hypercube geometry.
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

#include <fstream>

#include "../mortar_functions/generate_grids.h"
#include "../mortar_functions/mortar_operator.h"

using namespace dealii;

template <int dim>
class MyMortarManager : public MortarManagerBase<dim>
{
public:
  template <int dim2>
  MyMortarManager(const unsigned int      n_subdivisions,
                  const Quadrature<dim2> &quadrature,
                  const double            shift)
    : MortarManagerBase<dim>(n_subdivisions,
                             1.0 / (2.0 * numbers::PI),
                             quadrature,
                             shift * (2.0 * numbers::PI))
  {}

protected:
  Tensor<1, dim, double>
  get_normal(const Point<dim> &) const override
  {
    return Point<dim>(1.0, 0.0);
  }

  Point<dim>
  from_1D(const double rad) const override
  {
    return Point<dim>(0.5, rad / (2.0 * numbers::PI));
  }

  double
  to_1D(const Point<dim> &face_center) const override
  {
    return (2.0 * numbers::PI) * face_center[1];
  }
};

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
      const unsigned int dim                  = 2;
      const unsigned int n_global_refinements = 2;
      const double       shift                = 0.4;
      const MPI_Comm     comm                 = MPI_COMM_WORLD;

      ConditionalOStream pcout(std::cout,
                               Utilities::MPI::this_mpi_process(comm) == 0);

      FE_Q<dim>     fe(fe_degree);
      MappingQ<dim> mapping_q(mapping_degree);
      QGauss<dim>   quadrature(fe_degree + 1);

      // generate merged grid
      parallel::distributed::Triangulation<dim> tria(comm);
      split_hyper_cube(tria);

      std::vector<
        GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
        periodic_faces;
      GridTools::collect_periodic_faces(tria, 2, 3, 1, periodic_faces);
      GridTools::collect_periodic_faces(tria, 6, 7, 1, periodic_faces);
      tria.add_periodicity(periodic_faces);

      tria.refine_global(n_global_refinements);

      MappingQCache<dim> mapping(mapping_degree);

      mapping.initialize(
        mapping_q,
        tria,
        [&](const auto &cell, const auto &point) {
          if (cell->center()[0] > 0.5)
            return point;

          auto p = point;
          p[1] += shift;

          return p;
        },
        false);

      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      AffineConstraints<double> constraints;
      const IndexSet            locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         locally_relevant_dofs);
      DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
      DoFTools::make_zero_boundary_constraints(dof_handler, 5, constraints);

      DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
      DoFTools::make_periodicity_constraints(dof_handler, 6, 7, 1, constraints);


      constraints.close();

      PoissonOperator<dim, 1, double> op(mapping,
                                         dof_handler,
                                         constraints,
                                         quadrature);

      const std::shared_ptr<MortarManagerBase<dim>> mortar_manager =
        std::make_shared<MyMortarManager<dim>>(
          Utilities::pow(2, n_global_refinements), quadrature, shift);

      op.add_coupling(mortar_manager, 1, 4);

      LinearAlgebra::distributed::Vector<double> rhs, solution;
      op.initialize_dof_vector(rhs);
      op.initialize_dof_vector(solution);

      VectorTools::create_right_hand_side(
        mapping,
        dof_handler,
        quadrature,
        ScalarFunctionFromFunctionObject<dim>([](const auto &p) {
          return p[0] * std::pow(std::sin(numbers::PI * p[1]), 2.0);
        }),
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

      constraints.distribute(solution);

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
