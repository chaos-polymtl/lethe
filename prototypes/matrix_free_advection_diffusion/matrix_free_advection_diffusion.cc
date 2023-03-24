/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 - 2022 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------*/


#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// We use a parameter file for all the settings
struct Settings
{
  bool
  try_parse(const std::string &prm_filename);

  enum PreconditionerType
  {
    amg,
    gmg,
    ilu
  };

  enum GeometryType
  {
    hyperball,
    hypercube
  };

  enum SourceTermType
  {
    zero,
    mms
  };

  enum ProblemType
  {
    boundary_layer,
    double_glazing
  };

  PreconditionerType preconditioner;
  GeometryType       geometry;
  SourceTermType     source_term;
  ProblemType        problem_type;

  int          dimension;
  unsigned int element_order;
  unsigned int number_of_cycles;
  unsigned int initial_refinement;
  double       peclet_number;
  bool         output;
  std::string  output_name;
  std::string  output_path;
};



bool
Settings::try_parse(const std::string &prm_filename)
{
  ParameterHandler prm;
  prm.declare_entry("dim",
                    "2",
                    Patterns::Integer(),
                    "The problem dimension <2|3>");
  prm.declare_entry("element order",
                    "1",
                    Patterns::Integer(),
                    "Order of FE element <1|2|3>");
  prm.declare_entry("number of cycles",
                    "1",
                    Patterns::Integer(),
                    "Number of cycles <1 up to 9-dim >");
  prm.declare_entry("geometry",
                    "hyperball",
                    Patterns::Selection("hyperball|hypercube"),
                    "Geometry <hyperball|hypercube>");
  prm.declare_entry("initial refinement",
                    "1",
                    Patterns::Integer(),
                    "Global refinement 1st cycle");
  prm.declare_entry("peclet number", "10", Patterns::Double(), "Peclet number");
  prm.declare_entry("output",
                    "true",
                    Patterns::Bool(),
                    "Output vtu files <true|false>");
  prm.declare_entry("output name",
                    "solution",
                    Patterns::FileName(),
                    "Name for vtu files");
  prm.declare_entry("output path",
                    "./",
                    Patterns::FileName(),
                    "Path for vtu output files");
  prm.declare_entry("preconditioner",
                    "AMG",
                    Patterns::Selection("AMG|GMG|ILU"),
                    "GMRES Preconditioner <AMG|GMG|ILU>");
  prm.declare_entry("source term",
                    "zero",
                    Patterns::Selection("zero|mms"),
                    "Source term <zero|mms>");
  prm.declare_entry("problem type",
                    "boundary layer",
                    Patterns::Selection("boundary layer|double glazing"),
                    "Problem type <boundary layer|double glazing>");

  if (prm_filename.size() == 0)
    {
      std::cout
        << "****  Error: No input file provided!\n"
        << "****  Error: Call this program as './matrix_based_non_linear_poisson input.prm\n"
        << '\n'
        << "****  You may want to use one of the input files in this\n"
        << "****  directory, or use the following default values\n"
        << "****  to create an input file:\n";
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        prm.print_parameters(std::cout, ParameterHandler::Text);
      return false;
    }

  try
    {
      prm.parse_input(prm_filename);
    }
  catch (std::exception &e)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        std::cerr << e.what() << std::endl;
      return false;
    }

  if (prm.get("preconditioner") == "AMG")
    this->preconditioner = amg;
  else if (prm.get("preconditioner") == "GMG")
    this->preconditioner = gmg;
  else if (prm.get("preconditioner") == "ILU")
    this->preconditioner = ilu;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("geometry") == "hyperball")
    this->geometry = hyperball;
  else if (prm.get("geometry") == "hypercube")
    this->geometry = hypercube;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("source term") == "zero")
    this->source_term = zero;
  else if (prm.get("source term") == "mms")
    this->source_term = mms;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("problem type") == "boundary layer")
    this->problem_type = boundary_layer;
  else if (prm.get("problem type") == "double glazing")
    this->problem_type = double_glazing;
  else
    AssertThrow(false, ExcNotImplemented());

  this->dimension          = prm.get_integer("dim");
  this->element_order      = prm.get_integer("element order");
  this->number_of_cycles   = prm.get_integer("number of cycles");
  this->initial_refinement = prm.get_integer("initial refinement");
  this->peclet_number      = prm.get_double("peclet number");
  this->output             = prm.get_bool("output");
  this->output_name        = prm.get("output name");
  this->output_path        = prm.get("output path");


  return true;
}

template <int dim>
class SourceTerm : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p,
        const unsigned int /* component */ = 0) const override
  {
    const double coeff  = dim * numbers::PI * numbers::PI;
    double       factor = 1.0;
    for (unsigned int d = 0; d < dim; ++d)
      {
        factor *= std::sin(numbers::PI * p[d]);
      }
    return -std::exp(factor) + coeff * factor;
  }
};

template <int dim>
class AnalyticalSolution : public Function<dim>
{
public:
  AnalyticalSolution(double peclet)
    : Pe(peclet)
  {}

  virtual double
  value(const Point<dim> &p,
        const unsigned int /* component */ = 0) const override
  {
    if (dim == 2)
      {
        // Boundary layer case
        double eps = 1 / Pe;
        return p[0] *
               ((1 - std::exp((p[1] - 1) / eps)) / (1 - std::exp(-2 / eps)));
      }
    return 0;
  }

private:
  double Pe;
};

template <int dim>
class AdvectionField : public TensorFunction<1, dim>
{
public:
  AdvectionField(const Settings::ProblemType &test_case)
    : TensorFunction<1, dim>()
    , test_case(test_case)
  {}

  Tensor<1, dim>
  value(const Point<dim> &p) const override
  {
    (void)p;
    Tensor<1, dim> result;
    if (dim == 2)
      {
        if (test_case == Settings::boundary_layer)
          {
            result[0] = 0.0;
            result[1] = 1.0;
          }
        else if (test_case == Settings::double_glazing)
          {
            result[0] = 2 * p[1] * (1 - p[0] * p[0]);
            result[1] = -2 * p[0] * (1 - p[1] * p[1]);
          }
        else
          {
            result[0] = 1.0;
            result[1] = 0.0;
          }
      }
    else if (dim == 3)
      {
        result[0] = 1.0;
        result[1] = 0.0;
        result[2] = 0.0;
      }
    return result;
  }

private:
  Settings::ProblemType test_case;
};

template <int dim>
class BoundaryFunction : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p,
        const unsigned int /* component */ = 0) const override
  {
    return p[0];
  }
};

template <int dim, int fe_degree, typename number>
class AdvectionDiffusionOperator
  : public MatrixFreeOperators::Base<dim,
                                     LinearAlgebra::distributed::Vector<number>>
{
public:
  using value_type = number;

  using FECellIntegrator =
    FEEvaluation<dim, fe_degree, fe_degree + 1, 1, number>;

  AdvectionDiffusionOperator();

  virtual void
  clear() override;

  void
  reinit_operator_parameters( const Settings &parameters);

  void
  evaluate_newton_step(
    const LinearAlgebra::distributed::Vector<number> &newton_step);

  virtual void
  compute_diagonal() override;

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix(const AffineConstraints<number> &constraints);

private:
  virtual void
  apply_add(
    LinearAlgebra::distributed::Vector<number> &      dst,
    const LinearAlgebra::distributed::Vector<number> &src) const override;

  void
  local_apply(const MatrixFree<dim, number> &                   data,
              LinearAlgebra::distributed::Vector<number> &      dst,
              const LinearAlgebra::distributed::Vector<number> &src,
              const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  local_compute(FECellIntegrator &integrator) const;

  Table<2, VectorizedArray<number>>      nonlinear_values;
  mutable TrilinosWrappers::SparseMatrix system_matrix;
  Settings parameters;
};

template <int dim, int fe_degree, typename number>
AdvectionDiffusionOperator<dim, fe_degree, number>::AdvectionDiffusionOperator()
  : MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>()
{
  nonlinear_values.reinit(0, 0);
  system_matrix.clear();
}


template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::clear()
{
  nonlinear_values.reinit(0, 0);
  MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number>>::
    clear();
}

template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::reinit_operator_parameters(
  const Settings &parameters)
{
  this->parameters = parameters;
}

template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::evaluate_newton_step(
  const LinearAlgebra::distributed::Vector<number> &newton_step)
{
  const unsigned int n_cells = this->data->n_cell_batches();
  FECellIntegrator   phi(*this->data);

  nonlinear_values.reinit(n_cells, phi.n_q_points);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values_plain(newton_step);
      phi.evaluate(EvaluationFlags::values);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          nonlinear_values(cell, q) = 0.0; // std::exp(phi.get_value(q));
        }
    }
}

template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::local_apply(
  const MatrixFree<dim, number> &                   data,
  LinearAlgebra::distributed::Vector<number> &      dst,
  const LinearAlgebra::distributed::Vector<number> &src,
  const std::pair<unsigned int, unsigned int> &     cell_range) const
{
  FECellIntegrator phi(data);

  AdvectionField<dim>   advection_field(parameters.problem_type);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0),
                      phi.get_matrix_free().n_cell_batches());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);

      phi.reinit(cell);

      phi.gather_evaluate(src,
                          EvaluationFlags::values | EvaluationFlags::gradients);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          // Get advection field vector
          Point<dim, VectorizedArray<number>> point_batch =
            phi.quadrature_point(q);
          Tensor<1, dim> advection_vector;
          for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
            {
              Point<dim> single_point;
              for (unsigned int d = 0; d < dim; ++d)
                single_point[d] = point_batch[d][v];
              advection_vector = advection_field.value(single_point);
            }

          phi.submit_value(-nonlinear_values(cell, q) * phi.get_value(q) +
                             advection_vector * phi.get_gradient(q),
                           q);
          phi.submit_gradient(1 / parameters.peclet_number * phi.get_gradient(q), q);
        }

      phi.integrate_scatter(EvaluationFlags::values |
                              EvaluationFlags::gradients,
                            dst);
    }
}

template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::apply_add(
  LinearAlgebra::distributed::Vector<number> &      dst,
  const LinearAlgebra::distributed::Vector<number> &src) const
{
  this->data->cell_loop(&AdvectionDiffusionOperator::local_apply,
                        this,
                        dst,
                        src);
}

template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::local_compute(
  FECellIntegrator &phi) const
{
  AssertDimension(nonlinear_values.size(0),
                  phi.get_matrix_free().n_cell_batches());
  AssertDimension(nonlinear_values.size(1), phi.n_q_points);

  const unsigned int cell = phi.get_current_cell_index();

  phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

  AdvectionField<dim>   advection_field(parameters.problem_type);

  for (unsigned int q = 0; q < phi.n_q_points; ++q)
    {
      // Get advection field vector
      Point<dim, VectorizedArray<number>> point_batch = phi.quadrature_point(q);
      Tensor<1, dim>                      advection_vector;
      for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
        {
          Point<dim> single_point;
          for (unsigned int d = 0; d < dim; ++d)
            single_point[d] = point_batch[d][v];
          advection_vector = advection_field.value(single_point);
        }

      phi.submit_value(-nonlinear_values(cell, q) * phi.get_value(q) +
                         advection_vector * phi.get_gradient(q),
                       q);
      phi.submit_gradient(1 / parameters.peclet_number * phi.get_gradient(q), q);
    }

  phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

template <int dim, int fe_degree, typename number>
void
AdvectionDiffusionOperator<dim, fe_degree, number>::compute_diagonal()
{
  this->inverse_diagonal_entries.reset(
    new DiagonalMatrix<LinearAlgebra::distributed::Vector<number>>());
  LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
    this->inverse_diagonal_entries->get_vector();
  this->data->initialize_dof_vector(inverse_diagonal);

  MatrixFreeTools::compute_diagonal(*this->data,
                                    inverse_diagonal,
                                    &AdvectionDiffusionOperator::local_compute,
                                    this);

  for (auto &diagonal_element : inverse_diagonal)
    {
      diagonal_element =
        (std::abs(diagonal_element) > 1.0e-10) ? (1.0 / diagonal_element) : 1.0;
    }
}

template <int dim, int fe_degree, typename number>
const TrilinosWrappers::SparseMatrix &
AdvectionDiffusionOperator<dim, fe_degree, number>::get_system_matrix(
  const AffineConstraints<number> &constraints)
{
  if (system_matrix.m() == 0 && system_matrix.n() == 0)
    {
      AffineConstraints<number> constraints_copy;
      constraints_copy.copy_from(constraints);

      const auto &dof_handler = this->get_matrix_free()->get_dof_handler();
      const auto &mg_level    = this->get_matrix_free()->get_mg_level();

      TrilinosWrappers::SparsityPattern dsp(
        mg_level != numbers::invalid_unsigned_int ?
          dof_handler.locally_owned_mg_dofs(mg_level) :
          dof_handler.locally_owned_dofs(),
        dof_handler.get_triangulation().get_communicator());

      if (mg_level != numbers::invalid_unsigned_int)
        MGTools::make_sparsity_pattern(dof_handler, dsp, mg_level, constraints);
      else
        DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);


      dsp.compress();
      system_matrix.reinit(dsp);

      MatrixFreeTools::compute_matrix(
        *this->data,
        constraints,
        system_matrix,
        &AdvectionDiffusionOperator::local_compute,
        this);
    }
  return system_matrix;
}

template <int dim, int fe_degree>
class MatrixFreeAdvectionDiffusion
{
public:
  MatrixFreeAdvectionDiffusion(const Settings &parameters);

  void
  run();

private:
  void
  make_grid();

  void
  setup_system();

  void
  setup_gmg();

  void
  evaluate_residual(
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src) const;

  void
  local_evaluate_residual(
    const MatrixFree<dim, double> &                   data,
    LinearAlgebra::distributed::Vector<double> &      dst,
    const LinearAlgebra::distributed::Vector<double> &src,
    const std::pair<unsigned int, unsigned int> &     cell_range) const;

  void
  assemble_rhs();

  double
  compute_residual(const double alpha);

  void
  compute_update();

  void
  solve();

  double
  compute_solution_norm() const;

  double
  compute_l2_error() const;

  void
  output_results(const unsigned int cycle) const;

  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ<dim>                       mapping;

  FE_Q<dim>                 fe;
  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> constraints;
  using SystemMatrixType = AdvectionDiffusionOperator<dim, fe_degree, double>;
  SystemMatrixType  system_matrix;
  MGConstrainedDoFs mg_constrained_dofs;
  using LevelMatrixType = AdvectionDiffusionOperator<dim, fe_degree, double>;
  MGLevelObject<LevelMatrixType>                            mg_matrices;
  MGLevelObject<LinearAlgebra::distributed::Vector<double>> mg_solution;
  MGTransferMatrixFree<dim, double>                         mg_transfer;

  LinearAlgebra::distributed::Vector<double> solution;
  LinearAlgebra::distributed::Vector<double> newton_update;
  LinearAlgebra::distributed::Vector<double> system_rhs;

  unsigned int       linear_iterations;
  ConditionalOStream pcout;
  TimerOutput        computing_timer;

  Settings parameters;
};

template <int dim, int fe_degree>
MatrixFreeAdvectionDiffusion<dim, fe_degree>::MatrixFreeAdvectionDiffusion(
  const Settings &parameters)
  : triangulation(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
  , mapping(fe_degree)
  , fe(fe_degree)
  , dof_handler(triangulation)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , computing_timer(MPI_COMM_WORLD,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , parameters(parameters)
{}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "make grid");

  switch (parameters.geometry)
    {
        case Settings::hyperball: {
          SphericalManifold<dim>                boundary_manifold;
          TransfiniteInterpolationManifold<dim> inner_manifold;

          GridGenerator::hyper_ball(triangulation);

          triangulation.set_all_manifold_ids(1);
          triangulation.set_all_manifold_ids_on_boundary(0);

          triangulation.set_manifold(0, boundary_manifold);

          inner_manifold.initialize(triangulation);
          triangulation.set_manifold(1, inner_manifold);

          break;
        }
        case Settings::hypercube: {
          GridGenerator::hyper_cube(triangulation, -1.0, 1.0, true);
          break;
        }
    }

  triangulation.refine_global(parameters.initial_refinement);
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup system");

  system_matrix.clear();

  system_matrix.reinit_operator_parameters(parameters);

  dof_handler.distribute_dofs(fe);

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);

  // Set homogeneous constraints for the matrix-free operator
  if (parameters.problem_type == Settings::boundary_layer)
    {
      // Create zero BCs for the delta.
      // Left wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      // Right wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      // Top wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      // Bottom wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
    }
  else if (parameters.problem_type == Settings::double_glazing)
    {
      // Left wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      // Right wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      // Top wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
      // Bottom wall
      VectorTools::interpolate_boundary_values(dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(),
                                               constraints);
    }
  constraints.close();

  {
    typename MatrixFree<dim, double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim, double>::AdditionalData::none;
    additional_data.mapping_update_flags =
      (update_values | update_gradients | update_JxW_values |
       update_quadrature_points);
    auto system_mf_storage = std::make_shared<MatrixFree<dim, double>>();
    system_mf_storage->reinit(mapping,
                              dof_handler,
                              constraints,
                              QGauss<1>(fe.degree + 1),
                              additional_data);

    system_matrix.initialize(system_mf_storage);
  }

  system_matrix.initialize_dof_vector(solution);
  system_matrix.initialize_dof_vector(newton_update);
  system_matrix.initialize_dof_vector(system_rhs);

  // Change values of initial newton iteration
  if (parameters.problem_type == Settings::boundary_layer)
    {
      // Set boundary values for the initial newton iteration
      std::map<types::global_dof_index, double> boundary_values_left_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ConstantFunction<dim>(
                                                 -1.0),
                                               boundary_values_left_wall);

      for (auto &boundary_value : boundary_values_left_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;

      std::map<types::global_dof_index, double> boundary_values_right_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ConstantFunction<dim>(
                                                 1.0),
                                               boundary_values_right_wall);

      for (auto &boundary_value : boundary_values_right_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;

      std::map<types::global_dof_index, double> boundary_values_top_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values_top_wall);

      for (auto &boundary_value : boundary_values_top_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;

      std::map<types::global_dof_index, double> boundary_values_bottom_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               2,
                                               BoundaryFunction<dim>(),
                                               boundary_values_bottom_wall);

      for (auto &boundary_value : boundary_values_bottom_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;
    }
  else if (parameters.problem_type == Settings::double_glazing)
    {
      // Set boundary values for the initial newton iteration
      std::map<types::global_dof_index, double> boundary_values_left_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values_left_wall);

      for (auto &boundary_value : boundary_values_left_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;

      std::map<types::global_dof_index, double> boundary_values_right_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               Functions::ConstantFunction<dim>(
                                                 1.0),
                                               boundary_values_right_wall);

      for (auto &boundary_value : boundary_values_right_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;

      std::map<types::global_dof_index, double> boundary_values_top_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values_top_wall);

      for (auto &boundary_value : boundary_values_top_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;

      std::map<types::global_dof_index, double> boundary_values_bottom_wall;
      VectorTools::interpolate_boundary_values(dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(),
                                               boundary_values_bottom_wall);

      for (auto &boundary_value : boundary_values_bottom_wall)
        if (solution.locally_owned_elements().is_element(boundary_value.first))
          solution(boundary_value.first) = boundary_value.second;
    }
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::setup_gmg()
{
  TimerOutput::Scope t(computing_timer, "setup GMG");

  dof_handler.distribute_mg_dofs();

  mg_matrices.clear_elements();

  const unsigned int nlevels = triangulation.n_global_levels();
  mg_matrices.resize(0, nlevels - 1);
  mg_solution.resize(0, nlevels - 1);

  const std::set<types::boundary_id> dirichlet_boundary_ids = {0};
  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary_ids);

  mg_transfer.initialize_constraints(mg_constrained_dofs);
  mg_transfer.build(dof_handler);

  for (unsigned int level = 0; level < nlevels; ++level)
    {
      mg_matrices[level].reinit_operator_parameters(parameters);
      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

      AffineConstraints<double> level_constraints;
      level_constraints.reinit(relevant_dofs);
      level_constraints.add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints.close();

      typename MatrixFree<dim, double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim, double>::AdditionalData::partition_color;
      additional_data.mapping_update_flags =
        (update_values | update_gradients | update_JxW_values |
         update_quadrature_points);
      additional_data.mg_level = level;
      auto mg_mf_storage_level = std::make_shared<MatrixFree<dim, double>>();
      mg_mf_storage_level->reinit(mapping,
                                  dof_handler,
                                  level_constraints,
                                  QGauss<1>(fe.degree + 1),
                                  additional_data);

      mg_matrices[level].initialize(mg_mf_storage_level,
                                    mg_constrained_dofs,
                                    level);
      mg_matrices[level].initialize_dof_vector(mg_solution[level]);
    }
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::evaluate_residual(
  LinearAlgebra::distributed::Vector<double> &      dst,
  const LinearAlgebra::distributed::Vector<double> &src) const
{
  auto matrix_free = system_matrix.get_matrix_free();

  matrix_free->cell_loop(&MatrixFreeAdvectionDiffusion::local_evaluate_residual,
                         this,
                         dst,
                         src,
                         true);
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::local_evaluate_residual(
  const MatrixFree<dim, double> &                   data,
  LinearAlgebra::distributed::Vector<double> &      dst,
  const LinearAlgebra::distributed::Vector<double> &src,
  const std::pair<unsigned int, unsigned int> &     cell_range) const
{
  FEEvaluation<dim, fe_degree, fe_degree + 1, 1, double> phi(data);
  SourceTerm<dim>                                        source_term_function;
  AdvectionField<dim> advection_field(parameters.problem_type);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi.reinit(cell);

      phi.read_dof_values_plain(src);
      phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          VectorizedArray<double> source_value = VectorizedArray<double>(0.0);

          if (parameters.source_term == Settings::mms)
            {
              Point<dim, VectorizedArray<double>> point_batch =
                phi.quadrature_point(q);

              for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
                {
                  Point<dim> single_point;
                  for (unsigned int d = 0; d < dim; ++d)
                    single_point[d] = point_batch[d][v];
                  source_value[v] = source_term_function.value(single_point);
                }
            }

          Point<dim, VectorizedArray<double>> point_batch =
            phi.quadrature_point(q);

          Tensor<1, dim> advection_vector;
          for (unsigned int v = 0; v < VectorizedArray<double>::size(); ++v)
            {
              Point<dim> single_point;
              for (unsigned int d = 0; d < dim; ++d)
                single_point[d] = point_batch[d][v];
              advection_vector = advection_field.value(single_point);
            }

          phi.submit_value(advection_vector * phi.get_gradient(q) -
                             source_value,
                           q);
          phi.submit_gradient(1 / parameters.peclet_number *
                                phi.get_gradient(q),
                              q);
        }

      phi.integrate_scatter(EvaluationFlags::values |
                              EvaluationFlags::gradients,
                            dst);
    }
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::assemble_rhs()
{
  TimerOutput::Scope t(computing_timer, "assemble right hand side");

  evaluate_residual(system_rhs, solution);

  system_rhs *= -1.0;
}

template <int dim, int fe_degree>
double
MatrixFreeAdvectionDiffusion<dim, fe_degree>::compute_residual(
  const double alpha)
{
  TimerOutput::Scope t(computing_timer, "compute residual");

  LinearAlgebra::distributed::Vector<double> residual;
  LinearAlgebra::distributed::Vector<double> evaluation_point;

  system_matrix.initialize_dof_vector(residual);
  system_matrix.initialize_dof_vector(evaluation_point);

  evaluation_point = solution;
  if (alpha > 1e-12)
    {
      evaluation_point.add(alpha, newton_update);
    }

  evaluate_residual(residual, evaluation_point);

  return residual.l2_norm();
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::compute_update()
{
  TimerOutput::Scope t(computing_timer, "compute update");

  solution.update_ghost_values();

  system_matrix.evaluate_newton_step(solution);

  SolverControl solver_control(1000, 1.e-12);
  SolverGMRES<LinearAlgebra::distributed::Vector<double>> gmres(solver_control);

  newton_update = 0.0;

  switch (parameters.preconditioner)
    {
        case Settings::amg: {
          TrilinosWrappers::PreconditionAMG                 preconditioner;
          TrilinosWrappers::PreconditionAMG::AdditionalData data;

          if (fe_degree > 1)
            data.higher_order_elements = true;

          preconditioner.initialize(
            system_matrix.get_system_matrix(constraints), data);
          gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
          break;
        }
        case Settings::gmg: {
          mg_transfer.interpolate_to_mg(dof_handler, mg_solution, solution);
          // Set up smoother
          using SmootherType =
            PreconditionChebyshev<LevelMatrixType,
                                  LinearAlgebra::distributed::Vector<double>>;
          mg::SmootherRelaxation<SmootherType,
                                 LinearAlgebra::distributed::Vector<double>>
                                                               mg_smoother;
          MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
          smoother_data.resize(0, triangulation.n_global_levels() - 1);
          for (unsigned int level = 0; level < triangulation.n_global_levels();
               ++level)
            {
              if (level > 0)
                {
                  smoother_data[level].smoothing_range     = 15.;
                  smoother_data[level].degree              = 4;
                  smoother_data[level].eig_cg_n_iterations = 10;
                }
              else
                {
                  smoother_data[0].smoothing_range = 1e-3;
                  smoother_data[0].degree = numbers::invalid_unsigned_int;
                  smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
                }

              mg_matrices[level].evaluate_newton_step(mg_solution[level]);
              mg_matrices[level].compute_diagonal();
              smoother_data[level].preconditioner =
                mg_matrices[level].get_matrix_diagonal_inverse();
            }
          mg_smoother.initialize(mg_matrices, smoother_data);

          // Set up preconditioned coarse-grid solver
          SolverControl coarse_solver_control(1000, 1e-6, false, false);
          SolverGMRES<LinearAlgebra::distributed::Vector<double>> coarse_solver(
            coarse_solver_control);

          TrilinosWrappers::PreconditionAMG                 precondition_amg;
          TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
          precondition_amg.initialize(
            mg_matrices[0].get_system_matrix(constraints), amg_data);

          MGCoarseGridIterativeSolver<
            LinearAlgebra::distributed::Vector<double>,
            SolverGMRES<LinearAlgebra::distributed::Vector<double>>,
            LevelMatrixType,
            TrilinosWrappers::PreconditionAMG>
            mg_coarse(coarse_solver, mg_matrices[0], precondition_amg);

          // Set up multigrid
          mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_matrix(
            mg_matrices);

          MGLevelObject<
            MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
            mg_interface_matrices;
          mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1);
          for (unsigned int level = 0; level < triangulation.n_global_levels();
               ++level)
            {
              mg_interface_matrices[level].initialize(mg_matrices[level]);
            }
          mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_interface(
            mg_interface_matrices);

          Multigrid<LinearAlgebra::distributed::Vector<double>> mg(
            mg_matrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother);
          mg.set_edge_matrices(mg_interface, mg_interface);


          PreconditionMG<dim,
                         LinearAlgebra::distributed::Vector<double>,
                         MGTransferMatrixFree<dim, double>>
            preconditioner(dof_handler, mg, mg_transfer);

          gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
          break;
        }
        case Settings::ilu: {
          TrilinosWrappers::PreconditionILU                 preconditioner;
          TrilinosWrappers::PreconditionILU::AdditionalData data_ilu;
          preconditioner.initialize(
            system_matrix.get_system_matrix(constraints), data_ilu);

          gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
          break;
        }
      default:
        Assert(
          false,
          ExcMessage(
            "This program supports only AMG, GMG and ILU as preconditioners for the GMRES solver."));
    }

  constraints.distribute(newton_update);

  linear_iterations = solver_control.last_step();

  solution.zero_out_ghost_values();
}


template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::solve()
{
  TimerOutput::Scope t(computing_timer, "solve");

  // tolerances after a few steps.
  const unsigned int itmax = 10;
  const double       TOLf  = 1e-12;
  const double       TOLx  = 1e-10;


  Timer solver_timer;
  solver_timer.start();

  for (unsigned int newton_step = 1; newton_step <= itmax; ++newton_step)
    {
      assemble_rhs();
      compute_update();
      const double ERRx = newton_update.l2_norm();
      const double ERRf = compute_residual(1.0);
      solution.add(1.0, newton_update);

      pcout << "   Nstep " << newton_step << ", errf = " << ERRf
            << ", errx = " << ERRx << ", it = " << linear_iterations
            << std::endl;

      if (ERRf < TOLf || ERRx < TOLx)
        {
          solver_timer.stop();

          pcout << "Convergence step " << newton_step << " value " << ERRf
                << " (used wall time: " << solver_timer.wall_time() << " s)"
                << std::endl;

          break;
        }
      else if (newton_step == itmax)
        {
          solver_timer.stop();
          pcout << "WARNING: No convergence of Newton's method after "
                << newton_step << " steps." << std::endl;

          break;
        }
    }
}

template <int dim, int fe_degree>
double
MatrixFreeAdvectionDiffusion<dim, fe_degree>::compute_solution_norm() const
{
  solution.update_ghost_values();

  Vector<float> norm_per_cell(triangulation.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    Functions::ZeroFunction<dim>(),
                                    norm_per_cell,
                                    QGauss<dim>(fe.degree + 1),
                                    VectorTools::H1_seminorm);

  solution.zero_out_ghost_values();

  return VectorTools::compute_global_error(triangulation,
                                           norm_per_cell,
                                           VectorTools::H1_seminorm);
}

template <int dim, int fe_degree>
double
MatrixFreeAdvectionDiffusion<dim, fe_degree>::compute_l2_error() const
{
  solution.update_ghost_values();

  Vector<float> error_per_cell(triangulation.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    AnalyticalSolution<dim>(
                                      parameters.peclet_number),
                                    error_per_cell,
                                    QGauss<dim>(fe.degree + 1),
                                    VectorTools::L2_norm);

  solution.zero_out_ghost_values();

  return VectorTools::compute_global_error(triangulation,
                                           error_per_cell,
                                           VectorTools::L2_norm);
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::output_results(
  const unsigned int cycle) const
{
  if (triangulation.n_global_active_cells() > 1e6)
    return;

  solution.update_ghost_values();

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");

  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    {
      subdomain(i) = triangulation.locally_owned_subdomain();
    }
  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches(mapping, fe.degree, DataOut<dim>::curved_inner_cells);

  DataOutBase::VtkFlags flags;
  flags.compression_level = DataOutBase::VtkFlags::best_speed;
  data_out.set_flags(flags);

  std::string test_case = "";
  if (parameters.problem_type == Settings::boundary_layer)
    test_case = "_boundary_layer";
  else if (parameters.problem_type == Settings::double_glazing)
    test_case = "_double_glazing";

  data_out.write_vtu_with_pvtu_record(
    parameters.output_path,
    parameters.output_name + std::to_string(dim) + "d_Pe" +
      std::to_string(parameters.peclet_number) + test_case,
    cycle,
    MPI_COMM_WORLD,
    3);

  solution.zero_out_ghost_values();
}

template <int dim, int fe_degree>
void
MatrixFreeAdvectionDiffusion<dim, fe_degree>::run()
{
  {
    const unsigned int n_ranks =
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned int n_vect_doubles = VectorizedArray<double>::size();
    const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

    std::string DAT_header = "START DATE: " + Utilities::System::get_date() +
                             ", TIME: " + Utilities::System::get_time() +
                             ", MATRIX-FREE SOLVER";
    std::string MPI_header = "Running with " + std::to_string(n_ranks) +
                             " MPI process" + (n_ranks > 1 ? "es" : "");
    std::string VEC_header =
      "Vectorization over " + std::to_string(n_vect_doubles) +
      " doubles = " + std::to_string(n_vect_bits) + " bits (" +
      Utilities::System::get_current_vectorization_level() +
      "), VECTORIZATION_LEVEL=" +
      std::to_string(DEAL_II_COMPILER_VECTORIZATION_LEVEL);
    std::string SOL_header     = "Finite element space: " + fe.get_name();
    std::string PRECOND_header = "";
    if (parameters.preconditioner == Settings::amg)
      PRECOND_header = "Preconditioner: AMG";
    else if (parameters.preconditioner == Settings::gmg)
      PRECOND_header = "Preconditioner: GMG";
    else if (parameters.preconditioner == Settings::ilu)
      PRECOND_header = "Preconditioner: ILU";
    std::string GEOMETRY_header = "";
    if (parameters.geometry == Settings::hyperball)
      GEOMETRY_header = "Geometry: hyperball";
    else if (parameters.geometry == Settings::hypercube)
      GEOMETRY_header = "Geometry: hypercube";
    std::string SOURCE_header = "";
    if (parameters.source_term == Settings::zero)
      SOURCE_header = "Source term: zero";
    else if (parameters.source_term == Settings::mms)
      SOURCE_header = "Source term: according to MMS";
    std::string REFINE_header =
      "Initial refinement: " + std::to_string(parameters.initial_refinement);
    std::string CYCLES_header =
      "Total number of cycles: " + std::to_string(parameters.number_of_cycles);
    std::string PECLET_header =
      "Peclet number: " + std::to_string(parameters.peclet_number);
    std::string PROBLEM_TYPE_header = "";
    if (parameters.problem_type == Settings::boundary_layer)
      PROBLEM_TYPE_header = "Test problem: boundary layer";
    else if (parameters.problem_type == Settings::double_glazing)
      PROBLEM_TYPE_header = "Test problem: double glazing";

    pcout << std::string(80, '=') << std::endl;
    pcout << DAT_header << std::endl;
    pcout << std::string(80, '-') << std::endl;

    pcout << MPI_header << std::endl;
    pcout << VEC_header << std::endl;
    pcout << SOL_header << std::endl;
    pcout << PRECOND_header << std::endl;
    pcout << GEOMETRY_header << std::endl;
    pcout << SOURCE_header << std::endl;
    pcout << REFINE_header << std::endl;
    pcout << CYCLES_header << std::endl;
    pcout << PECLET_header << std::endl;
    pcout << PROBLEM_TYPE_header << std::endl;

    pcout << std::string(80, '=') << std::endl;
  }

  for (unsigned int cycle = 0; cycle < parameters.number_of_cycles; ++cycle)
    {
      pcout << std::string(80, '-') << std::endl;
      pcout << "Cycle " << cycle << std::endl;
      pcout << std::string(80, '-') << std::endl;

      if (cycle == 0)
        {
          make_grid();
        }
      else
        {
          triangulation.refine_global(1);
        }

      Timer timer;

      pcout << "Set up system..." << std::endl;
      setup_system();

      if (parameters.preconditioner == Settings::gmg)
        setup_gmg();

      pcout << "   Triangulation: " << triangulation.n_global_active_cells()
            << " cells" << std::endl;
      pcout << "   DoFHandler:    " << dof_handler.n_dofs() << " DoFs"
            << std::endl;

      if (parameters.preconditioner == Settings::gmg)
        for (unsigned int level = 0; level < triangulation.n_global_levels();
             ++level)
          pcout << "   MG Level " << level << ": " << dof_handler.n_dofs(level)
                << " DoFs" << std::endl;

      pcout << std::endl;


      pcout << "Solve using Newton's method..." << std::endl;
      solve();
      pcout << std::endl;


      timer.stop();
      pcout << "Time for setup+solve (CPU/Wall) " << timer.cpu_time() << '/'
            << timer.wall_time() << " s" << std::endl;
      pcout << std::endl;

      const double norm = compute_solution_norm();
      if (parameters.output)
        {
          pcout << "Output results..." << std::endl;
          output_results(cycle);
        }

      pcout << "  H1 seminorm: " << norm << std::endl;
      pcout << std::endl;

      if (parameters.problem_type == Settings::boundary_layer)
        {
          pcout << "  L2 norm: " << compute_l2_error() << std::endl;
        }

      computing_timer.print_summary();
      computing_timer.reset();
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize       mpi_init(argc, argv, 1);
  dealii::Utilities::System::MemoryStats stats;
  dealii::Utilities::System::get_memory_stats(stats);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);


  Settings parameters;
  if (!parameters.try_parse((argc > 1) ? (argv[1]) : ""))
    return 0;

  try
    {
      switch (parameters.dimension)
        {
            case 2: {
              switch (parameters.element_order)
                {
                    case 1: {
                      MatrixFreeAdvectionDiffusion<2, 1>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                    case 2: {
                      MatrixFreeAdvectionDiffusion<2, 2>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                    case 3: {
                      MatrixFreeAdvectionDiffusion<2, 3>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                }
              break;
            }

            case 3: {
              switch (parameters.element_order)
                {
                    case 1: {
                      MatrixFreeAdvectionDiffusion<3, 1>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                    case 2: {
                      MatrixFreeAdvectionDiffusion<3, 2>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                    case 3: {
                      MatrixFreeAdvectionDiffusion<3, 3>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                }
              break;
            }

          default:
            Assert(
              false,
              ExcMessage(
                "This program only works in 2d and 3d and for element orders equal to 1, 2 or 3."));
        }


      pcout << "MEMORY STATS: " << std::endl;
      const auto print = [&pcout](const double value) {
        const auto min_max_avg =
          dealii::Utilities::MPI::min_max_avg(value, MPI_COMM_WORLD);

        pcout << "MIN: " << min_max_avg.min << " MAX: " << min_max_avg.max
              << " AVG: " << min_max_avg.avg << " SUM: " << min_max_avg.sum
              << std::endl;
      };

      pcout << "VmPeak: ";
      print(stats.VmPeak);

      pcout << "VmSize: ";
      print(stats.VmSize);

      pcout << "VmHWM: ";
      print(stats.VmHWM);

      pcout << "VmRSS: ";
      print(stats.VmRSS);

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
