// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test solves the Laplace equation outside a disk of radius R with
 * the immersed boundary constraints of a scalar field, which impose a flux
 * into the disk on its surface with the sharp-interface stencils. The exact
 * solution is u = A + ln(r/R), where r is the distance to the center of the
 * disk. It is imposed on the boundary of the square domain. With a unit
 * diffusivity D and R = 0.25, the flux into the disk is D du/dr = 4 on its
 * surface, where u = A. Three boundary conditions give this flux:
 * - neumann: the flux q = 4 is imposed, and A = 1;
 * - robin of order 1: the flux alpha u + q is imposed with alpha = 1 and
 *   q = 2, so that A = 2;
 * - robin of order 2: the flux alpha u^2 is imposed with alpha = 1, so that
 *   A = 2.
 * The values of A ensure that a wrong order or a wrong sign of q would change
 * the solution. The problems are solved with Newton's method. The test prints,
 * for two element degrees and three meshes, the number of constrained degrees
 * of freedom, the L2 norm of the error on the fluid cells, and the relative
 * error of the total flux into the disk computed from the residual of the
 * equation on the fluid cells, whose exact value is 2 pi D.
 */

// Lethe
#include <core/ib_particle.h>
#include <core/ib_scalar_constraints.h>
#include <core/parameters.h>
#include <core/vector.h>

// deal.II
#include <deal.II/base/function.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cmath>
#include <string>
#include <vector>

/**
 * @brief Exact solution u = A + ln(r/R) of the Laplace equation outside a
 * disk.
 */
class LogarithmicSolution : public Function<2>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] p_center Center of the disk.
   *
   * @param[in] p_radius Radius of the disk.
   *
   * @param[in] p_surface_value Value A of the solution on the surface of the
   * disk.
   */
  LogarithmicSolution(const Point<2> &p_center,
                      const double    p_radius,
                      const double    p_surface_value)
    : center(p_center)
    , radius(p_radius)
    , surface_value(p_surface_value)
  {}

  /**
   * @brief Return the exact solution at a point.
   *
   * @param[in] point Point at which the solution is evaluated.
   *
   * @return Value of the exact solution.
   */
  double
  value(const Point<2> &point, const unsigned int /*component*/) const override
  {
    return surface_value + std::log(point.distance(center) / radius);
  }

private:
  const Point<2> center;
  const double   radius;
  const double   surface_value;
};

/**
 * @brief Boundary condition imposed on the surface of the disk, and value of
 * the exact solution on this surface.
 */
struct TestCase
{
  /// Name printed in the output.
  std::string name;
  /// Content of the subsection "tracer boundary condition".
  std::string boundary_condition_parameters;
  /// Value A of the exact solution on the surface of the disk.
  double surface_value;
};

void
test(const TestCase    &test_case,
     const unsigned int fe_degree,
     const unsigned int n_refinements)
{
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  // The center of the disk is shifted, so that no vertex of the mesh lies on
  // its surface.
  const Point<2>            center(0.0123, -0.0071);
  const double              radius = 0.25;
  const LogarithmicSolution exact_solution(center,
                                           radius,
                                           test_case.surface_value);

  IBParticle<2> particle;
  particle.initialize_all();
  particle.initialize_shape("sphere", std::vector<double>{radius});
  particle.set_position(center);
  const std::vector<IBParticle<2>> particles{particle};

  parallel::distributed::Triangulation<2> triangulation(mpi_communicator);
  GridGenerator::hyper_cube(triangulation, -1., 1.);
  triangulation.refine_global(n_refinements);

  const FE_Q<2>     fe(fe_degree);
  const MappingQ<2> mapping(1);
  DoFHandler<2>     dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  // Constraints of the Newton correction (homogeneous) and of the solution
  // (exact solution on the boundary of the square)
  AffineConstraints<double> zero_constraints(locally_owned_dofs,
                                             locally_relevant_dofs);
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<2>(), zero_constraints);
  zero_constraints.close();

  AffineConstraints<double> nonzero_constraints(locally_owned_dofs,
                                                locally_relevant_dofs);
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, exact_solution, nonzero_constraints);
  nonzero_constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints, true);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  TrilinosWrappers::SparseMatrix system_matrix;
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  // Initial guess that satisfies the boundary condition on the square
  GlobalVectorType locally_owned_solution(locally_owned_dofs, mpi_communicator);
  nonzero_constraints.distribute(locally_owned_solution);
  GlobalVectorType evaluation_point(locally_owned_dofs,
                                    locally_relevant_dofs,
                                    mpi_communicator);
  evaluation_point = locally_owned_solution;
  evaluation_point.update_ghost_values();

  // Immersed boundary constraints imposing the flux on the surface of the
  // disk, with the default stencil of the solver
  Parameters::IBTracerBoundaryCondition<2> boundary_condition;
  ParameterHandler                         prm;
  boundary_condition.declare_parameters(prm);
  prm.parse_input_from_string("subsection tracer boundary condition\n" +
                              test_case.boundary_condition_parameters +
                              "end\n");
  boundary_condition.parse_parameters(prm);

  const unsigned int     stencil_degree = 2;
  const double           length_ratio   = 4.;
  IBScalarConstraints<2> immersed_boundary_constraints(stencil_degree,
                                                       length_ratio,
                                                       true,
                                                       boundary_condition);
  constexpr double       diffusivity       = 1.;
  constexpr double       inverse_time_step = 0.;
  immersed_boundary_constraints.update(dof_handler,
                                       dof_handler,
                                       mapping,
                                       particles,
                                       zero_constraints,
                                       diffusivity,
                                       inverse_time_step);

  const QGauss<2>           quadrature(fe_degree + 1);
  const unsigned int        n_dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int        n_q_points      = quadrature.size();
  FEValues<2>               fe_values(mapping,
                        fe,
                        quadrature,
                        update_gradients | update_JxW_values);
  FullMatrix<double>        local_matrix(n_dofs_per_cell, n_dofs_per_cell);
  Vector<double>            local_rhs(n_dofs_per_cell);
  std::vector<Tensor<1, 2>> solution_gradients(n_q_points);
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);
  GlobalVectorType system_rhs(locally_owned_dofs, mpi_communicator);
  GlobalVectorType newton_update(locally_owned_dofs, mpi_communicator);

  // Flux into the disk, computed from the residual of the equation on the
  // fluid cells at the last evaluation point
  std::vector<double> particle_fluxes;

  // Newton's method, with a direct solver so that the result does not depend
  // on the convergence of an iterative solver
  constexpr unsigned int max_newton_iterations = 20;
  constexpr double       newton_tolerance      = 1e-10;
  for (unsigned int iteration = 0; iteration < max_newton_iterations;
       ++iteration)
    {
      // Laplace equation on the fluid cells, in Newton form
      system_matrix = 0;
      system_rhs    = 0;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (!cell->is_locally_owned() ||
              !immersed_boundary_constraints.is_fluid_cell(cell))
            continue;

          fe_values.reinit(cell);
          fe_values.get_function_gradients(evaluation_point,
                                           solution_gradients);
          local_matrix = 0;
          local_rhs    = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                  local_matrix(i, j) += fe_values.shape_grad(i, q) *
                                        fe_values.shape_grad(j, q) *
                                        fe_values.JxW(q);
                local_rhs(i) -= fe_values.shape_grad(i, q) *
                                solution_gradients[q] * fe_values.JxW(q);
              }
          cell->get_dof_indices(local_dof_indices);
          zero_constraints.distribute_local_to_global(local_matrix,
                                                      local_rhs,
                                                      local_dof_indices,
                                                      system_matrix,
                                                      system_rhs);
        }
      system_matrix.compress(VectorOperation::add);
      system_rhs.compress(VectorOperation::add);

      particle_fluxes =
        immersed_boundary_constraints.calculate_particle_fluxes(system_rhs);

      immersed_boundary_constraints.assemble_matrix(dof_handler,
                                                    mapping,
                                                    quadrature,
                                                    zero_constraints,
                                                    evaluation_point,
                                                    system_matrix);
      immersed_boundary_constraints.assemble_rhs(dof_handler,
                                                 mapping,
                                                 quadrature,
                                                 zero_constraints,
                                                 evaluation_point,
                                                 system_rhs);

      if (system_rhs.l2_norm() < newton_tolerance)
        break;

      TrilinosWrappers::SolverDirect solver;
      solver.solve(system_matrix, newton_update, system_rhs);
      zero_constraints.distribute(newton_update);

      locally_owned_solution += newton_update;
      evaluation_point = locally_owned_solution;
      evaluation_point.update_ghost_values();
    }

  // L2 norm of the error on the fluid cells
  const QGauss<2>     error_quadrature(fe_degree + 2);
  FEValues<2>         fe_values_error(mapping,
                              fe,
                              error_quadrature,
                              update_values | update_quadrature_points |
                                update_JxW_values);
  std::vector<double> solution_values(error_quadrature.size());
  double              squared_error = 0.;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() ||
          !immersed_boundary_constraints.is_fluid_cell(cell))
        continue;

      fe_values_error.reinit(cell);
      fe_values_error.get_function_values(evaluation_point, solution_values);
      for (unsigned int q = 0; q < error_quadrature.size(); ++q)
        {
          const double error =
            solution_values[q] -
            exact_solution.value(fe_values_error.quadrature_point(q), 0);
          squared_error += error * error * fe_values_error.JxW(q);
        }
    }
  const double error_norm =
    std::sqrt(Utilities::MPI::sum(squared_error, mpi_communicator));

  const double exact_flux = 2. * numbers::PI * diffusivity;
  const double flux_relative_error =
    std::abs(particle_fluxes[0] - exact_flux) / exact_flux;

  deallog << test_case.name << ", degree " << fe_degree << ", refinements "
          << n_refinements << ": "
          << immersed_boundary_constraints.get_n_interface_dofs()
          << " constrained degrees of freedom ("
          << immersed_boundary_constraints.get_n_dofs_without_stencil()
          << " without stencil), L2 error " << error_norm
          << ", flux relative error " << flux_relative_error << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();
      deallog.precision(3);

      const std::vector<TestCase> test_cases{
        {"Neumann",
         "  set type = neumann\n"
         "  subsection flux\n"
         "    set Function expression = 4\n"
         "  end\n",
         1.},
        {"Robin of order 1",
         "  set type               = robin\n"
         "  set reaction constant  = 1\n"
         "  set reaction order     = 1\n"
         "  subsection flux\n"
         "    set Function expression = 2\n"
         "  end\n",
         2.},
        {"Robin of order 2",
         "  set type               = robin\n"
         "  set reaction constant  = 1\n"
         "  set reaction order     = 2\n",
         2.}};

      for (const TestCase &test_case : test_cases)
        for (unsigned int fe_degree = 1; fe_degree <= 2; ++fe_degree)
          for (unsigned int n_refinements = 4; n_refinements <= 6;
               ++n_refinements)
            test(test_case, fe_degree, n_refinements);
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
