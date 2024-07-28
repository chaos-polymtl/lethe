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
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

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
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/meshworker/mesh_loop.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Define all the parameters that can be specified in the .prm file
struct Settings
{
  bool
  try_parse(const std::string &prm_filename);

  enum PreconditionerType
  {
    amg,
    gmg
  };

  enum GeometryType
  {
    hyperball,
    hypercube,
    hyperrectangle
  };

  enum SourceTermType
  {
    zero,
    mms
  };

  PreconditionerType preconditioner;
  GeometryType       geometry;
  SourceTermType     source_term;

  int          dimension;
  unsigned int element_order;
  unsigned int number_of_cycles;
  unsigned int initial_refinement;
  unsigned int repetitions;
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
                    Patterns::Selection("hyperball|hypercube|hyperrectangle"),
                    "Geometry <hyperball|hypercube|hyperrectangle>");
  prm.declare_entry("initial refinement",
                    "1",
                    Patterns::Integer(),
                    "Global refinement 1st cycle");
  prm.declare_entry("repetitions",
                    "1",
                    Patterns::Integer(),
                    "Repetitions in z direction for the hyper rectangle");
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
                    Patterns::Selection("AMG|GMG"),
                    "Preconditioner <AMG|GMG>");
  prm.declare_entry("source term",
                    "zero",
                    Patterns::Selection("zero|mms"),
                    "Source term <zero|mms>");

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
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("geometry") == "hyperball")
    this->geometry = hyperball;
  else if (prm.get("geometry") == "hypercube")
    this->geometry = hypercube;
  else if (prm.get("geometry") == "hyperrectangle")
    this->geometry = hyperrectangle;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("source term") == "zero")
    this->source_term = zero;
  else if (prm.get("source term") == "mms")
    this->source_term = mms;
  else
    AssertThrow(false, ExcNotImplemented());

  this->dimension          = prm.get_integer("dim");
  this->element_order      = prm.get_integer("element order");
  this->number_of_cycles   = prm.get_integer("number of cycles");
  this->initial_refinement = prm.get_integer("initial refinement");
  this->repetitions        = prm.get_integer("repetitions");
  this->output             = prm.get_bool("output");
  this->output_name        = prm.get("output name");
  this->output_path        = prm.get("output path");

  return true;
}

template <int dim>
class MMSSolution : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p,
        const unsigned int /* component */ = 0) const override
  {
    double val = 1.0;
    for (unsigned int d = 0; d < dim; ++d)
      {
        val *= std::sin(numbers::PI * p[d]);
      }
    return val;
  }
};

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

// Structure created to use the cell worker approach to assemble mg
template <int dim>
struct ScratchData
{
  ScratchData(const Mapping<dim>       &mapping,
              const FiniteElement<dim> &fe,
              const unsigned int        quadrature_degree,
              const UpdateFlags         update_flags)
    : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags)
  {}

  ScratchData(const ScratchData<dim> &scratch_data)
    : fe_values(scratch_data.fe_values.get_mapping(),
                scratch_data.fe_values.get_fe(),
                scratch_data.fe_values.get_quadrature(),
                scratch_data.fe_values.get_update_flags())
  {}

  FEValues<dim>       fe_values;
  std::vector<double> old_solution_values;
};

// Structure created to use the cell worker approach to assemble mg
struct CopyData
{
  unsigned int                         level;
  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  template <class Iterator>
  void
  reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    cell_rhs.reinit(dofs_per_cell);

    local_dof_indices.resize(dofs_per_cell);
    cell->get_active_or_mg_dof_indices(local_dof_indices);
    level = cell->level();
  }
};

// Main class for the nonlinear Poisson problem given by
// -Laplacian(u) = exp(u) using Newton's method, the matrix-based
// approach and zero Dirichlet boundary conditions.
template <int dim, int fe_degree>
class MatrixBasedPoissonProblem
{
public:
  MatrixBasedPoissonProblem(const Settings &parameters);

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
  assemble_rhs();

  void
  assemble_matrix();

  void
  assemble_gmg();

  template <class Iterator>
  void
  cell_worker(const Iterator   &cell,
              ScratchData<dim> &scratch_data,
              CopyData         &copy_data);

  void
  assemble_gmg_meshworker();

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

  using VectorType = TrilinosWrappers::MPI::Vector;
  using MatrixType = TrilinosWrappers::SparseMatrix;

  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ<dim>                       mapping;

  FE_Q<dim>                 fe;
  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> constraints;
  MatrixType                system_matrix;
  SparsityPattern           sparsity_pattern;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  MGLevelObject<MatrixType> mg_matrix;
  MGLevelObject<MatrixType> mg_interface_in;
  MGConstrainedDoFs         mg_constrained_dofs;
  MGLevelObject<VectorType> mg_solution;

  VectorType solution;
  VectorType newton_update;
  VectorType system_rhs;

  unsigned int       linear_iterations;
  ConditionalOStream pcout;
  TimerOutput        computing_timer;
  MPI_Comm           mpi_communicator;

  Settings parameters;
};

template <int dim, int fe_degree>
MatrixBasedPoissonProblem<dim, fe_degree>::MatrixBasedPoissonProblem(
  const Settings &parameters)
  : triangulation(MPI_COMM_WORLD,
                  Triangulation<dim>::limit_level_difference_at_vertices,
                  (parameters.preconditioner == Settings::amg) ?
                    parallel::distributed::Triangulation<dim>::default_setting :
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
  , mapping(fe_degree)
  , fe(fe_degree)
  , dof_handler(triangulation)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , computing_timer(MPI_COMM_WORLD,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , mpi_communicator(MPI_COMM_WORLD)
  , parameters(parameters)
{}


template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "make grid");

  switch (parameters.geometry)
    {
      case Settings::hyperball:
        {
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
      case Settings::hypercube:
        {
          GridGenerator::hyper_cube(triangulation);
          break;
        }
      case Settings::hyperrectangle:
        {
          std::vector<unsigned int> repetitions(dim);
          for (unsigned int i = 0; i < dim - 1; i++)
            {
              repetitions[i] = 1;
            }
          repetitions[dim - 1] = parameters.repetitions;

          GridGenerator::subdivided_hyper_rectangle(
            triangulation,
            repetitions,
            Point<dim>(),
            (dim == 2) ? Point<dim>(1., parameters.repetitions) :
                         Point<dim>(1., 1., parameters.repetitions));
          break;
        }
    }

  triangulation.refine_global(parameters.initial_refinement);
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup system");

  system_matrix.clear();

  dof_handler.distribute_dofs(fe);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  newton_update.reinit(locally_owned_dofs,
                       locally_relevant_dofs,
                       mpi_communicator);

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             dof_handler.locally_owned_dofs(),
                                             mpi_communicator,
                                             locally_relevant_dofs);
  constraints.condense(dsp);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::setup_gmg()
{
  TimerOutput::Scope timing(computing_timer, "setup GMG");

  dof_handler.distribute_mg_dofs();

  mg_constrained_dofs.clear();
  mg_constrained_dofs.initialize(dof_handler);

  const std::set<types::boundary_id> boundary_ids = {types::boundary_id(0)};
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, boundary_ids);

  const unsigned int n_levels = triangulation.n_global_levels();

  mg_matrix.resize(0, n_levels - 1);
  mg_matrix.clear_elements();
  mg_interface_in.resize(0, n_levels - 1);
  mg_interface_in.clear_elements();
  mg_solution.resize(0, n_levels - 1);

  for (unsigned int level = 0; level < n_levels; ++level)
    {
      const IndexSet dof_set =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

      {
        TrilinosWrappers::SparsityPattern dsp(
          dof_handler.locally_owned_mg_dofs(level),
          dof_handler.locally_owned_mg_dofs(level),
          dof_set,
          mpi_communicator);
        MGTools::make_sparsity_pattern(dof_handler, dsp, level);

        dsp.compress();
        mg_matrix[level].reinit(dsp);
      }

      {
        TrilinosWrappers::SparsityPattern dsp(
          dof_handler.locally_owned_mg_dofs(level),
          dof_handler.locally_owned_mg_dofs(level),
          dof_set,
          mpi_communicator);

        MGTools::make_interface_sparsity_pattern(dof_handler,
                                                 mg_constrained_dofs,
                                                 dsp,
                                                 level);
        dsp.compress();
        mg_interface_in[level].reinit(dsp);
      }
      {
        const IndexSet locally_owned_mg_dofs =
          dof_handler.locally_owned_mg_dofs(level);
        const IndexSet locally_relevant_mg_dofs =
          DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

        mg_solution[level].reinit(locally_owned_mg_dofs,
                                  locally_relevant_mg_dofs,
                                  mpi_communicator);
      }
    }
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::assemble_rhs()
{
  TimerOutput::Scope t(computing_timer, "assemble right hand side");

  const QGauss<dim> quadrature_formula(fe.degree + 1);

  system_rhs = 0;
  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values |
                            update_quadrature_points);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
  const unsigned int n_q_points    = fe_values.n_quadrature_points;

  Vector<double>      cell_rhs(dofs_per_cell);
  SourceTerm<dim>     source_term;
  std::vector<double> source_term_values(n_q_points);

  std::vector<double>         newton_step_values(n_q_points);
  std::vector<Tensor<1, dim>> newton_step_gradients(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_rhs = 0.0;

          fe_values.reinit(cell);

          if (parameters.source_term == Settings::mms)
            source_term.value_list(fe_values.get_quadrature_points(),
                                   source_term_values);

          fe_values.get_function_values(solution, newton_step_values);
          fe_values.get_function_gradients(solution, newton_step_gradients);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double nonlinearity = std::exp(newton_step_values[q]);
              const double dx           = fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  if (parameters.source_term == Settings::mms)
                    cell_rhs(i) +=
                      (-grad_phi_i * newton_step_gradients[q] +
                       phi_i * nonlinearity + phi_i * source_term_values[q]) *
                      dx;
                  else
                    cell_rhs(i) += (-grad_phi_i * newton_step_gradients[q] +
                                    phi_i * nonlinearity) *
                                   dx;
                }
            }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_rhs,
                                                 local_dof_indices,
                                                 system_rhs);
        }
    }

  system_rhs.compress(VectorOperation::add);
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::assemble_matrix()
{
  TimerOutput::Scope t(computing_timer, "assemble matrix");

  const QGauss<dim> quadrature_formula(fe.degree + 1);

  system_matrix = 0;

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
  const unsigned int n_q_points    = fe_values.n_quadrature_points;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<double> newton_step_values(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_matrix = 0.0;

          fe_values.reinit(cell);

          fe_values.get_function_values(solution, newton_step_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double nonlinearity = std::exp(newton_step_values[q]);
              const double dx           = fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      const double         phi_j = fe_values.shape_value(j, q);
                      const Tensor<1, dim> grad_phi_j =
                        fe_values.shape_grad(j, q);

                      cell_matrix(i, j) += (grad_phi_i * grad_phi_j -
                                            phi_i * nonlinearity * phi_j) *
                                           dx;
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 system_matrix);
        }
    }

  system_matrix.compress(VectorOperation::add);
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::assemble_gmg()
{
  QGauss<dim> quadrature_formula(fe_degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
  const unsigned int n_q_points    = fe_values.n_quadrature_points;

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<double> newton_step_values(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  std::vector<AffineConstraints<double>> boundary_constraints(
    triangulation.n_global_levels());

  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    {
      const IndexSet dof_set =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
      boundary_constraints[level].reinit(dof_set);
      boundary_constraints[level].add_lines(
        mg_constrained_dofs.get_refinement_edge_indices(level));
      boundary_constraints[level].add_lines(
        mg_constrained_dofs.get_boundary_indices(level));

      boundary_constraints[level].close();
    }

  for (const auto &cell : dof_handler.cell_iterators())
    if (cell->level_subdomain_id() == triangulation.locally_owned_subdomain())
      {
        cell_matrix = 0;

        fe_values.reinit(cell);

        fe_values.get_function_values(solution, newton_step_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            const double nonlinearity = std::exp(newton_step_values[q]);
            const double dx           = fe_values.JxW(q);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                const double         phi_i      = fe_values.shape_value(i, q);
                const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    const double         phi_j = fe_values.shape_value(j, q);
                    const Tensor<1, dim> grad_phi_j =
                      fe_values.shape_grad(j, q);

                    cell_matrix(i, j) +=
                      (grad_phi_i * grad_phi_j - phi_i * nonlinearity * phi_j) *
                      dx;
                  }
              }
          }

        cell->get_mg_dof_indices(local_dof_indices);

        boundary_constraints[cell->level()].distribute_local_to_global(
          cell_matrix, local_dof_indices, mg_matrix[cell->level()]);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            if (mg_constrained_dofs.is_interface_matrix_entry(
                  cell->level(), local_dof_indices[i], local_dof_indices[j]))
              mg_interface_in[cell->level()].add(local_dof_indices[i],
                                                 local_dof_indices[j],
                                                 cell_matrix(i, j));
      }

  for (unsigned int i = 0; i < triangulation.n_global_levels(); ++i)
    {
      mg_matrix[i].compress(VectorOperation::add);
      mg_interface_in[i].compress(VectorOperation::add);
    }
}

template <int dim, int fe_degree>
template <class Iterator>
void
MatrixBasedPoissonProblem<dim, fe_degree>::cell_worker(
  const Iterator   &cell,
  ScratchData<dim> &scratch_data,
  CopyData         &copy_data)
{
  FEValues<dim> &fe_values = scratch_data.fe_values;
  fe_values.reinit(cell);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
  const unsigned int n_q_points    = fe_values.n_quadrature_points;

  copy_data.reinit(cell, dofs_per_cell);

  std::vector<double> newton_step_values(n_q_points);
  VectorType          old_solution = solution;
  fe_values.get_function_values(old_solution, newton_step_values);

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double nonlinearity = std::exp(newton_step_values[q]);
      const double dx           = fe_values.JxW(q);

      for (unsigned int i = 0; i < dofs_per_cell; i++)
        {
          const double         phi_i      = fe_values.shape_value(i, q);
          const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            {
              const double         phi_j      = fe_values.shape_value(j, q);
              const Tensor<1, dim> grad_phi_j = fe_values.shape_grad(j, q);

              copy_data.cell_matrix(i, j) +=
                (grad_phi_i * grad_phi_j - phi_i * nonlinearity * phi_j) * dx;
            }
        }
    }
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::assemble_gmg_meshworker()
{
  std::vector<AffineConstraints<double>> boundary_constraints(
    triangulation.n_global_levels());

  for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level)
    {
      const IndexSet dof_set =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
      boundary_constraints[level].reinit(dof_set);
      boundary_constraints[level].add_lines(
        mg_constrained_dofs.get_refinement_edge_indices(level));
      boundary_constraints[level].add_lines(
        mg_constrained_dofs.get_boundary_indices(level));

      boundary_constraints[level].close();
    }
  auto cell_worker =
    [&](const typename DoFHandler<dim>::level_cell_iterator &cell,
        ScratchData<dim>                                    &scratch_data,
        CopyData                                            &copy_data) {
      this->cell_worker(cell, scratch_data, copy_data);
    };

  auto copier = [&](const CopyData &cd) {
    boundary_constraints[cd.level].distribute_local_to_global(
      cd.cell_matrix, cd.local_dof_indices, mg_matrix[cd.level]);

    const unsigned int dofs_per_cell = cd.local_dof_indices.size();
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        if (mg_constrained_dofs.is_interface_matrix_entry(
              cd.level, cd.local_dof_indices[i], cd.local_dof_indices[j]))
          {
            mg_interface_in[cd.level].add(cd.local_dof_indices[i],
                                          cd.local_dof_indices[j],
                                          cd.cell_matrix(i, j));
          }
  };

  using CellFilter =
    FilteredIterator<typename DoFHandler<dim>::level_cell_iterator>;
  const unsigned int n_gauss_points = fe_degree + 1;

  ScratchData<dim> scratch_data(mapping,
                                fe,
                                n_gauss_points,
                                update_values | update_gradients |
                                  update_JxW_values | update_quadrature_points);

  MeshWorker::mesh_loop(CellFilter(IteratorFilters::LocallyOwnedLevelCell(),
                                   dof_handler.begin_mg()),
                        CellFilter(IteratorFilters::LocallyOwnedLevelCell(),
                                   dof_handler.end_mg()),
                        cell_worker,
                        copier,
                        scratch_data,
                        CopyData(),
                        MeshWorker::assemble_own_cells);
}

template <int dim, int fe_degree>
double
MatrixBasedPoissonProblem<dim, fe_degree>::compute_residual(const double alpha)
{
  TimerOutput::Scope t(computing_timer, "compute residual");

  VectorType residual;
  VectorType evaluation_point(system_rhs);
  VectorType local_newton_update(system_rhs);
  local_newton_update = newton_update;
  VectorType local_evaluation_point;

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  residual.reinit(locally_owned_dofs, mpi_communicator);
  local_evaluation_point.reinit(locally_owned_dofs,
                                locally_relevant_dofs,
                                mpi_communicator);

  evaluation_point = solution;

  if (alpha > 1e-12)
    {
      evaluation_point.add(alpha, local_newton_update);
    }

  local_evaluation_point = solution;

  const QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values |
                            update_quadrature_points);

  const unsigned int dofs_per_cell = fe_values.dofs_per_cell;
  const unsigned int n_q_points    = fe_values.n_quadrature_points;

  Vector<double>      cell_residual(dofs_per_cell);
  SourceTerm<dim>     source_term;
  std::vector<double> source_term_values(n_q_points);

  std::vector<double>         values(n_q_points);
  std::vector<Tensor<1, dim>> gradients(n_q_points);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell_residual = 0.0;

          fe_values.reinit(cell);

          if (parameters.source_term == Settings::mms)
            source_term.value_list(fe_values.get_quadrature_points(),
                                   source_term_values);

          fe_values.get_function_values(local_evaluation_point, values);
          fe_values.get_function_gradients(local_evaluation_point, gradients);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double nonlinearity = std::exp(values[q]);
              const double dx           = fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  if (parameters.source_term == Settings::mms)
                    cell_residual(i) +=
                      (grad_phi_i * gradients[q] - phi_i * nonlinearity -
                       phi_i * source_term_values[q]) *
                      dx;
                  else
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

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::compute_update()
{
  TimerOutput::Scope t(computing_timer, "compute update");

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs, mpi_communicator);

  SolverControl              solver_control(100, 1.e-12);
  TrilinosWrappers::SolverCG cg(solver_control);

  switch (parameters.preconditioner)
    {
      case Settings::amg:
        {
          TrilinosWrappers::PreconditionAMG                 preconditioner;
          TrilinosWrappers::PreconditionAMG::AdditionalData data;

          if (fe_degree > 1)
            data.higher_order_elements = true;

          preconditioner.initialize(system_matrix, data);

          cg.solve(system_matrix,
                   completely_distributed_solution,
                   system_rhs,
                   preconditioner);
          break;
        }
      case Settings::gmg:
        {
          MGTransferPrebuilt<VectorType> mg_transfer(mg_constrained_dofs);
          mg_transfer.build(dof_handler);

          mg_transfer.copy_to_mg(dof_handler, mg_solution, solution);

          assemble_gmg();
          // Second option: uncomment to use the meshworker
          // assemble_gmg_meshworker();

          SolverControl        coarse_solver_control(1000, 1e-12, false, false);
          SolverCG<VectorType> coarse_solver(coarse_solver_control);

          PreconditionIdentity identity;
          MGCoarseGridIterativeSolver<VectorType,
                                      SolverCG<VectorType>,
                                      MatrixType,
                                      PreconditionIdentity>
            coarse_grid_solver(coarse_solver, mg_matrix[0], identity);

          using Smoother = TrilinosWrappers::PreconditionChebyshev;
          MGSmootherPrecondition<MatrixType, Smoother, VectorType> smoother;
          MGLevelObject<typename Smoother::AdditionalData> smoother_data;
          smoother_data.resize(0, triangulation.n_global_levels() - 1);

          for (unsigned int level = 0; level < triangulation.n_global_levels();
               ++level)
            {
              if (level > 0)
                {
                  smoother_data[level].degree = 4;
                }
              else
                {
                  smoother_data[0].degree = 2;
                }
            }
          smoother.initialize(mg_matrix, smoother_data);

          mg::Matrix<VectorType> mg_m(mg_matrix);
          mg::Matrix<VectorType> mg_in(mg_interface_in);
          mg::Matrix<VectorType> mg_out(mg_interface_in);

          Multigrid<VectorType> mg(
            mg_m, coarse_grid_solver, mg_transfer, smoother, smoother);
          mg.set_edge_matrices(mg_out, mg_in);

          PreconditionMG<dim, VectorType, MGTransferPrebuilt<VectorType>>
            preconditioner(dof_handler, mg, mg_transfer);

          SolverCG<VectorType> solver(solver_control);
          solver.solve(system_matrix,
                       completely_distributed_solution,
                       system_rhs,
                       preconditioner);
          break;
        }
      default:
        Assert(false,
               ExcMessage(
                 "This program supports only AMG and GMG as preconditioner."));
    }

  constraints.distribute(completely_distributed_solution);

  linear_iterations = solver_control.last_step();

  newton_update = completely_distributed_solution;
}


template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::solve()
{
  TimerOutput::Scope t(computing_timer, "solve");

  const unsigned int itmax = 10;
  const double       TOLf  = 1e-12;
  const double       TOLx  = 1e-10;

  Timer solver_timer;
  solver_timer.start();

  for (unsigned int newton_step = 1; newton_step <= itmax; ++newton_step)
    {
      assemble_matrix();
      assemble_rhs();

      compute_update();

      // Get vectors without any ghosts
      VectorType local_newton_update(system_rhs);
      local_newton_update = newton_update;
      VectorType local_solution(system_rhs);
      local_solution = solution;

      const double ERRx = local_newton_update.l2_norm();
      const double ERRf = compute_residual(1.0);

      local_solution.add(1.0, local_newton_update);

      solution = local_solution;

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
MatrixBasedPoissonProblem<dim, fe_degree>::compute_solution_norm() const
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

  return VectorTools::compute_global_error(triangulation,
                                           norm_per_cell,
                                           VectorTools::H1_seminorm);
}

template <int dim, int fe_degree>
double
MatrixBasedPoissonProblem<dim, fe_degree>::compute_l2_error() const
{
  solution.update_ghost_values();

  Vector<float> error_per_cell(triangulation.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    MMSSolution<dim>(),
                                    error_per_cell,
                                    QGauss<dim>(fe.degree + 1),
                                    VectorTools::L2_norm);

  return VectorTools::compute_global_error(triangulation,
                                           error_per_cell,
                                           VectorTools::L2_norm);
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::output_results(
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
  data_out.write_vtu_with_pvtu_record(parameters.output_path,
                                      parameters.output_name +
                                        std::to_string(dim) + "d",
                                      cycle,
                                      MPI_COMM_WORLD,
                                      3);
}

template <int dim, int fe_degree>
void
MatrixBasedPoissonProblem<dim, fe_degree>::run()
{
  {
    const unsigned int n_ranks =
      Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned int n_vect_doubles = VectorizedArray<double>::size();
    const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

    std::string DAT_header = "START DATE: " + Utilities::System::get_date() +
                             ", TIME: " + Utilities::System::get_time() +
                             ", MATRIX-BASED SOLVER";
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
    std::string GEOMETRY_header = "";
    if (parameters.geometry == Settings::hyperball)
      GEOMETRY_header = "Geometry: hyperball";
    else if (parameters.geometry == Settings::hypercube)
      GEOMETRY_header = "Geometry: hypercube";
    else if (parameters.geometry == Settings::hyperrectangle)
      GEOMETRY_header = "Geometry: hyper rectangle";
    std::string SOURCE_header = "";
    if (parameters.source_term == Settings::zero)
      SOURCE_header = "Source term: zero";
    else if (parameters.source_term == Settings::mms)
      SOURCE_header = "Source term: according to MMS";
    std::string REFINE_header =
      "Initial refinement: " + std::to_string(parameters.initial_refinement);
    std::string REPETITIONS_header =
      "Repetitions in z: " + std::to_string(parameters.repetitions);
    std::string CYCLES_header =
      "Total number of cycles: " + std::to_string(parameters.number_of_cycles);

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
    pcout << REPETITIONS_header << std::endl;
    pcout << CYCLES_header << std::endl;

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
        {
          for (unsigned int level = 0; level < triangulation.n_global_levels();
               ++level)
            pcout << "   MG Level " << level << ": "
                  << dof_handler.n_dofs(level) << " DoFs" << std::endl;
        }
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

      if (parameters.source_term == Settings::mms)
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
          case 2:
            {
              switch (parameters.element_order)
                {
                  case 1:
                    {
                      MatrixBasedPoissonProblem<2, 1>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 2:
                    {
                      MatrixBasedPoissonProblem<2, 2>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 3:
                    {
                      MatrixBasedPoissonProblem<2, 3>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                }
              break;
            }

          case 3:
            {
              switch (parameters.element_order)
                {
                  case 1:
                    {
                      MatrixBasedPoissonProblem<3, 1>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 2:
                    {
                      MatrixBasedPoissonProblem<3, 2>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 3:
                    {
                      MatrixBasedPoissonProblem<3, 3>
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
