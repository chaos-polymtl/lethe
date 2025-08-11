// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/refinement.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/smoothness_estimator.h>
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
    amg
  };

  enum GeometryType
  {
    hyperball,
    hypercube,
    hyperrectangle,
    hypercube_with_hole
  };

  enum SourceTermType
  {
    zero,
    constant
  };

  PreconditionerType preconditioner;
  GeometryType       geometry;
  SourceTermType     source_term;

  int          dimension;
  unsigned int element_order;
  unsigned int max_element_order;
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
  prm.declare_entry("max element order",
                    "3",
                    Patterns::Integer(),
                    "Maximum refinement order of FE element <1|2|3|4|5|6|7>");
  prm.declare_entry("number of cycles",
                    "1",
                    Patterns::Integer(),
                    "Number of cycles <1 up to 9-dim >");
  prm.declare_entry(
    "geometry",
    "hyperball",
    Patterns::Selection(
      "hyperball|hypercube|hyperrectangle|hypercube_with_hole"),
    "Geometry <hyperball|hypercube|hyperrectangle|hypercube_with_hole>");
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
                    Patterns::Selection("AMG"),
                    "Preconditioner <AMG>");
  prm.declare_entry("source term",
                    "zero",
                    Patterns::Selection("zero|constant"),
                    "Source term <zero|constant>");

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
#if DEAL_II_VERSION_GTE(9, 7, 0)
        prm.print_parameters(std::cout, ParameterHandler::DefaultStyle);
#else
        prm.print_parameters(std::cout, ParameterHandler::Text);
#endif
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
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("geometry") == "hyperball")
    this->geometry = hyperball;
  else if (prm.get("geometry") == "hypercube")
    this->geometry = hypercube;
  else if (prm.get("geometry") == "hyperrectangle")
    this->geometry = hyperrectangle;
  else if (prm.get("geometry") == "hypercube_with_hole")
    this->geometry = hypercube_with_hole;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("source term") == "zero")
    this->source_term = zero;
  else if (prm.get("source term") == "constant")
    this->source_term = constant;
  else
    AssertThrow(false, ExcNotImplemented());

  this->dimension          = prm.get_integer("dim");
  this->element_order      = prm.get_integer("element order");
  this->max_element_order  = prm.get_integer("max element order");
  this->number_of_cycles   = prm.get_integer("number of cycles");
  this->initial_refinement = prm.get_integer("initial refinement");
  this->repetitions        = prm.get_integer("repetitions");
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
    double product = 1.;
    for (unsigned int d = 0; d < dim; ++d)
      {
        product *= (p[d] + 1);
      }
    return product;
  }
};

// Main class for the nonlinear Poisson problem given by
// -Laplacian(u) = 1 using Newton's method, the matrix-based
// approach and zero Dirichlet boundary conditions.
template <int dim, int fe_degree>
class HPMatrixBasedPoissonProblem
{
public:
  HPMatrixBasedPoissonProblem(const Settings &parameters);

  void
  run();

private:
  void
  make_grid();

  void
  setup_system();

  void
  assemble_rhs();

  void
  assemble_matrix();

  double
  compute_residual(const double alpha);

  void
  compute_update();

  void
  solve();

  void
  hp_refine();

  void
  output_results(const unsigned int cycle) const;

  using VectorType = TrilinosWrappers::MPI::Vector;
  using MatrixType = TrilinosWrappers::SparseMatrix;

  parallel::distributed::Triangulation<dim> triangulation;
  const MappingQ<dim>                       mapping;

  hp::MappingCollection<dim>                  mapping_collection;
  hp::FECollection<dim>                       fe_collection;
  hp::QCollection<dim>                        quadrature_collection;
  hp::QCollection<dim - 1>                    face_quadrature_collection;
  DoFHandler<dim>                             dof_handler;
  AffineConstraints<double>                   constraints;
  MatrixType                                  system_matrix;
  SparsityPattern                             sparsity_pattern;
  std::unique_ptr<parallel::CellWeights<dim>> cell_weights;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  VectorType solution;
  VectorType newton_update;
  VectorType system_rhs;
  VectorType residual_to_output;

  unsigned int       linear_iterations;
  ConditionalOStream pcout;
  TimerOutput        computing_timer;
  MPI_Comm           mpi_communicator;

  Settings parameters;
};

template <int dim, int fe_degree>
HPMatrixBasedPoissonProblem<dim, fe_degree>::HPMatrixBasedPoissonProblem(
  const Settings &parameters)
  : triangulation(MPI_COMM_WORLD)
  , mapping(fe_degree)
  , dof_handler(triangulation)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , computing_timer(MPI_COMM_WORLD,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , mpi_communicator(MPI_COMM_WORLD)
  , parameters(parameters)
{
  for (unsigned int degree = parameters.element_order;
       degree <= parameters.max_element_order;
       ++degree)
    {
      fe_collection.push_back(FE_Q<dim>(degree));
      quadrature_collection.push_back(QGauss<dim>(degree + 1));
      face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
    }

  cell_weights = std::make_unique<parallel::CellWeights<dim>>(
    dof_handler, parallel::CellWeights<dim>::ndofs_weighting({1, 1}));

  const unsigned int min_fe_index = parameters.element_order;
  triangulation.signals.post_p4est_refinement.connect(
    [&, min_fe_index]() {
      const parallel::distributed::TemporarilyMatchRefineFlags<dim>
        refine_modifier(triangulation);
      hp::Refinement::limit_p_level_difference(dof_handler,
                                               parameters.element_order,
                                               /*contains=*/min_fe_index);
    },
    boost::signals2::at_front);
}


template <int dim, int fe_degree>
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::make_grid()
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
          GridGenerator::subdivided_hyper_cube(triangulation, 1.0);
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
      case Settings::hypercube_with_hole:
        {
          GridGenerator::subdivided_hyper_cube(triangulation, 4, -1., 1.);

          std::set<typename Triangulation<dim>::active_cell_iterator>
            cells_to_remove;
          for (const auto &cell : triangulation.active_cell_iterators())
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                 ++v)
              if (cell->vertex(v).square() < .1)
                cells_to_remove.insert(cell);

          GridGenerator::create_triangulation_with_removed_cells(
            triangulation, cells_to_remove, triangulation);
        }
    }

  const unsigned int min_fe_index = parameters.element_order;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      cell->set_active_fe_index(min_fe_index);

  triangulation.refine_global(parameters.initial_refinement);
}

template <int dim, int fe_degree>
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup system");

  system_matrix.clear();

  dof_handler.distribute_dofs(fe_collection);

  locally_owned_dofs    = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  newton_update.reinit(locally_owned_dofs,
                       locally_relevant_dofs,
                       mpi_communicator);

  constraints.clear();
  constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  constraints.make_consistent_in_parallel(locally_owned_dofs,
                                          locally_relevant_dofs,
                                          mpi_communicator);

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
HPMatrixBasedPoissonProblem<dim, fe_degree>::assemble_rhs()
{
  TimerOutput::Scope t(computing_timer, "assemble right hand side");

  system_rhs = 0;

  hp::FEValues<dim> hp_fe_values(fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

  Vector<double>  cell_rhs;
  SourceTerm<dim> source_term;

  std::vector<types::global_dof_index> local_dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          cell_rhs.reinit(dofs_per_cell);

          cell_rhs = 0.0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          const unsigned int  n_q_points = fe_values.n_quadrature_points;
          std::vector<double> source_term_values(n_q_points);

          std::vector<double>         newton_step_values(n_q_points);
          std::vector<Tensor<1, dim>> newton_step_gradients(n_q_points);

          if (parameters.source_term == Settings::constant)
            source_term.value_list(fe_values.get_quadrature_points(),
                                   source_term_values);

          fe_values.get_function_values(solution, newton_step_values);
          fe_values.get_function_gradients(solution, newton_step_gradients);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double nonlinearity = 0.;
              const double dx           = fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  if (parameters.source_term == Settings::constant)
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
          local_dof_indices.resize(dofs_per_cell);
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
HPMatrixBasedPoissonProblem<dim, fe_degree>::assemble_matrix()
{
  TimerOutput::Scope t(computing_timer, "assemble matrix");

  system_matrix = 0;

  hp::FEValues<dim> hp_fe_values(fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

  FullMatrix<double> cell_matrix;


  std::vector<types::global_dof_index> local_dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

          cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
          cell_matrix = 0.0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          const unsigned int n_q_points = fe_values.n_quadrature_points;

          std::vector<double> newton_step_values(n_q_points);

          fe_values.get_function_values(solution, newton_step_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double nonlinearity = 0.;
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
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix,
                                                 local_dof_indices,
                                                 system_matrix);
        }
    }

  system_matrix.compress(VectorOperation::add);
}

template <int dim, int fe_degree>
double
HPMatrixBasedPoissonProblem<dim, fe_degree>::compute_residual(
  const double alpha)
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

  hp::FEValues<dim> hp_fe_values(fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                   update_JxW_values |
                                   update_quadrature_points);


  Vector<double>  cell_residual;
  SourceTerm<dim> source_term;

  std::vector<types::global_dof_index> local_dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          cell_residual.reinit(dofs_per_cell);
          cell_residual = 0.0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();


          const unsigned int n_q_points = fe_values.n_quadrature_points;

          std::vector<double> source_term_values(n_q_points);

          std::vector<double>         values(n_q_points);
          std::vector<Tensor<1, dim>> gradients(n_q_points);

          if (parameters.source_term == Settings::constant)
            source_term.value_list(fe_values.get_quadrature_points(),
                                   source_term_values);

          fe_values.get_function_values(local_evaluation_point, values);
          fe_values.get_function_gradients(local_evaluation_point, gradients);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double nonlinearity = 0.;
              const double dx           = fe_values.JxW(q);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  if (parameters.source_term == Settings::constant)
                    cell_residual(i) +=
                      (grad_phi_i * gradients[q] - phi_i * nonlinearity -
                       phi_i * source_term_values[q]) *
                      dx;
                  else
                    cell_residual(i) +=
                      (grad_phi_i * gradients[q] - phi_i * nonlinearity) * dx;
                }
            }
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_residual,
                                                 local_dof_indices,
                                                 residual);
        }
    }

  residual.compress(VectorOperation::add);
  residual.update_ghost_values();

  residual_to_output = residual;

  return residual.l2_norm();
}

template <int dim, int fe_degree>
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::compute_update()
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
      default:
        Assert(false,
               ExcMessage("This program supports only AMG preconditioner."));
    }

  constraints.distribute(completely_distributed_solution);

  linear_iterations = solver_control.last_step();

  newton_update = completely_distributed_solution;
}


template <int dim, int fe_degree>
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::solve()
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
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::hp_refine()
{
  TimerOutput::Scope t(computing_timer, "hp refinement");

  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate(
    dof_handler,
    face_quadrature_collection,
    std::map<types::boundary_id, const Function<dim> *>(),
    solution,
    estimated_error_per_cell);

  Vector<float> smoothness_indicators(triangulation.n_active_cells());
  FESeries::Fourier<dim, dim> fourier =
    SmoothnessEstimator::Fourier::default_fe_series(fe_collection);
  SmoothnessEstimator::Fourier::coefficient_decay(fourier,
                                                  dof_handler,
                                                  solution,
                                                  smoothness_indicators);

  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
    triangulation, estimated_error_per_cell, 0.3, 0.03);

  hp::Refinement::p_adaptivity_from_relative_threshold(dof_handler,
                                                       smoothness_indicators,
                                                       0.2,
                                                       0.2);

  hp::Refinement::choose_p_over_h(dof_handler);

  triangulation.prepare_coarsening_and_refinement();
  hp::Refinement::limit_p_level_difference(dof_handler);

  triangulation.execute_coarsening_and_refinement();
}

template <int dim, int fe_degree>
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::output_results(
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

  // Output the finite element degree
  Vector<float> fe_degrees(triangulation.n_active_cells());
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_degrees(cell->active_cell_index()) =
            fe_collection[cell->active_fe_index()].degree;
        }
    }

  data_out.add_data_vector(fe_degrees, "fe_degree");

  data_out.add_data_vector(residual_to_output, "residual");

  data_out.build_patches();

  DataOutBase::VtkFlags flags;
  flags.compression_level = DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(flags);

  const std::string filename =
    (parameters.output_path + parameters.output_name + "-" +
     std::to_string(cycle) + ".vtu");
  data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD);
}

template <int dim, int fe_degree>
void
HPMatrixBasedPoissonProblem<dim, fe_degree>::run()
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
    std::string PRECOND_header = "";
    if (parameters.preconditioner == Settings::amg)
      PRECOND_header = "Preconditioner: AMG";
    std::string GEOMETRY_header = "";
    if (parameters.geometry == Settings::hyperball)
      GEOMETRY_header = "Geometry: hyperball";
    else if (parameters.geometry == Settings::hypercube)
      GEOMETRY_header = "Geometry: hypercube";
    else if (parameters.geometry == Settings::hyperrectangle)
      GEOMETRY_header = "Geometry: hyper rectangle";
    else if (parameters.geometry == Settings::hypercube_with_hole)
      GEOMETRY_header = "Geometry: hypercube with hole";
    std::string SOURCE_header = "";
    if (parameters.source_term == Settings::zero)
      SOURCE_header = "Source term: zero";
    else if (parameters.source_term == Settings::constant)
      SOURCE_header = "Source term: 1";
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

      Timer timer;

      pcout << "Set up system..." << std::endl;
      setup_system();

      pcout << "   Triangulation: " << triangulation.n_global_active_cells()
            << " cells" << std::endl;
      pcout << "   DoFHandler:    " << dof_handler.n_dofs() << " DoFs"
            << std::endl;

      pcout << std::endl;

      pcout << "Solve using Newton's method..." << std::endl;
      solve();
      pcout << std::endl;

      timer.stop();
      pcout << "Time for setup+solve (CPU/Wall) " << timer.cpu_time() << '/'
            << timer.wall_time() << " s" << std::endl;
      pcout << std::endl;

      if (parameters.output)
        {
          pcout << "Output results..." << std::endl;
          output_results(cycle);
          pcout << std::endl;
        }

      pcout << "Refine mesh using hp refinement..." << std::endl;
      hp_refine();
      pcout << std::endl;

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
                      HPMatrixBasedPoissonProblem<2, 1>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 2:
                    {
                      HPMatrixBasedPoissonProblem<2, 2>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 3:
                    {
                      HPMatrixBasedPoissonProblem<2, 3>
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
                      HPMatrixBasedPoissonProblem<3, 1>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 2:
                    {
                      HPMatrixBasedPoissonProblem<3, 2>
                        non_linear_poisson_problem(parameters);
                      non_linear_poisson_problem.run();

                      break;
                    }
                  case 3:
                    {
                      HPMatrixBasedPoissonProblem<3, 3>
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
