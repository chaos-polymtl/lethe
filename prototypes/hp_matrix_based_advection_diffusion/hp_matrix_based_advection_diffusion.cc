// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/fe/mapping_q.h>

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
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/smoothness_estimator.h>
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
    ilu
  };

  enum GeometryType
  {
    hyperball,
    hypercube,
    hypercube_with_hole,
    hyperrectangle
  };

  enum SourceTermType
  {
    zero,
    mms
  };

  enum ProblemType
  {
    boundary_layer,
    double_glazing,
    boundary_layer_with_hole
  };

  PreconditionerType preconditioner;
  GeometryType       geometry;
  SourceTermType     source_term;
  ProblemType        problem_type;

  int          dimension;
  unsigned int element_order;
  unsigned int number_of_cycles;
  unsigned int initial_refinement;
  unsigned int repetitions;
  double       peclet_number;
  bool         stabilization;
  bool         nonlinearity;
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
  prm.declare_entry(
    "geometry",
    "hyperball",
    Patterns::Selection(
      "hyperball|hypercube|hypercube with hole|hyperrectangle"),
    "Geometry <hyperball|hypercube|hypercube with hole|hyperrectangle>");
  prm.declare_entry("initial refinement",
                    "1",
                    Patterns::Integer(),
                    "Global refinement 1st cycle");
  prm.declare_entry("repetitions",
                    "1",
                    Patterns::Integer(),
                    "Repetitions in one direction for the hyperrectangle");
  prm.declare_entry("peclet number", "10", Patterns::Double(), "Peclet number");
  prm.declare_entry("stabilization",
                    "false",
                    Patterns::Bool(),
                    "Enable stabilization <true|false>");
  prm.declare_entry("nonlinearity",
                    "false",
                    Patterns::Bool(),
                    "Add exponential nonlinearity <true|false>");
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
                    Patterns::Selection("AMG|ILU"),
                    "GMRES Preconditioner <AMG|ILU>");
  prm.declare_entry("source term",
                    "zero",
                    Patterns::Selection("zero|mms"),
                    "Source term <zero|mms>");
  prm.declare_entry(
    "problem type",
    "boundary layer",
    Patterns::Selection(
      "boundary layer|double glazing|boundary layer with hole"),
    "Problem type <boundary layer|double glazing|boundary layer with hole>");

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
  else if (prm.get("preconditioner") == "ILU")
    this->preconditioner = ilu;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("geometry") == "hyperball")
    this->geometry = hyperball;
  else if (prm.get("geometry") == "hypercube")
    this->geometry = hypercube;
  else if (prm.get("geometry") == "hypercube with hole")
    this->geometry = hypercube_with_hole;
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

  if (prm.get("problem type") == "boundary layer")
    this->problem_type = boundary_layer;
  else if (prm.get("problem type") == "double glazing")
    this->problem_type = double_glazing;
  else if (prm.get("problem type") == "boundary layer with hole")
    this->problem_type = boundary_layer_with_hole;
  else
    AssertThrow(false, ExcNotImplemented());

  this->dimension          = prm.get_integer("dim");
  this->element_order      = prm.get_integer("element order");
  this->number_of_cycles   = prm.get_integer("number of cycles");
  this->initial_refinement = prm.get_integer("initial refinement");
  this->repetitions        = prm.get_integer("repetitions");
  this->peclet_number      = prm.get_double("peclet number");
  this->stabilization      = prm.get_bool("stabilization");
  this->nonlinearity       = prm.get_bool("nonlinearity");
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

// Only available for the boundary layer case in 2D
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
    else if (test_case == Settings::boundary_layer_with_hole)
      {
        result[0] = -std::sin(numbers::PI / 6);
        result[1] = std::cos(numbers::PI / 6);
      }
    else
      {
        result[0] = 1.0;
        result[1] = 0.0;
      }

    if (dim == 3)
      {
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

// For hypercube with hole case
template <int dim>
class BoundaryValues : public Function<dim>
{
public:
  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  virtual void
  value_list(const std::vector<Point<dim>> &points,
             std::vector<double>           &values,
             const unsigned int             component = 0) const override;
};

template <int dim>
double
BoundaryValues<dim>::value(const Point<dim>  &p,
                           const unsigned int component) const
{
  Assert(component == 0, ExcIndexRange(component, 0, 1));
  (void)component;

  // Set boundary to 1 if $x=1$, or if $x>0.5$ and $y=-1$.
  if (std::fabs(p[0] - 1) < 1e-8 || (std::fabs(p[1] + 1) < 1e-8 && p[0] >= 0.5))
    {
      return 1.0;
    }
  else
    {
      return 0.0;
    }
}

template <int dim>
void
BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points,
                                std::vector<double>           &values,
                                const unsigned int             component) const
{
  AssertDimension(values.size(), points.size());

  for (unsigned int i = 0; i < points.size(); ++i)
    values[i] = BoundaryValues<dim>::value(points[i], component);
}

// Main class for the advection-diffusion problem given by
// −ϵ∆u + w · ∇u = f using Newton's method, the matrix-based
// approach and different test problems.
template <int dim, int fe_degree>
class MatrixBasedAdvectionDiffusion
{
public:
  MatrixBasedAdvectionDiffusion(const Settings &parameters);

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
  hp_refine();

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

  hp::FECollection<dim>    fe_collection;
  hp::QCollection<dim>     quadrature_collection;
  hp::QCollection<dim - 1> face_quadrature_collection;

  DoFHandler<dim>           dof_handler;
  AffineConstraints<double> zero_constraints;
  AffineConstraints<double> nonzero_constraints;
  MatrixType                system_matrix;
  SparsityPattern           sparsity_pattern;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

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
MatrixBasedAdvectionDiffusion<dim, fe_degree>::MatrixBasedAdvectionDiffusion(
  const Settings &parameters)
  : triangulation(MPI_COMM_WORLD,
                  Triangulation<dim>::limit_level_difference_at_vertices,
                  (parameters.preconditioner == Settings::amg) ?
                    parallel::distributed::Triangulation<dim>::default_setting :
                    parallel::distributed::Triangulation<
                      dim>::construct_multigrid_hierarchy)
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
  for (unsigned int degree = 1; degree <= 4; ++degree)
    {
      fe_collection.push_back(FE_Q<dim>(degree));
      quadrature_collection.push_back(QGauss<dim>(degree + 1));
      face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1));
    }
}


template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::make_grid()
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
          GridGenerator::hyper_cube(triangulation, -1.0, 1.0, true);
          break;
        }
      case Settings::hypercube_with_hole:
        {
          GridGenerator::hyper_cube_with_cylindrical_hole(triangulation,
                                                          0.3,
                                                          1.0);

          const SphericalManifold<dim> manifold_description(Point<dim>(0, 0));
          triangulation.set_manifold(1, manifold_description);
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
            (dim == 2) ? Point<dim>(-1., -1.) : Point<dim>(-1., -1., -1.),
            // Point<dim>(),
            (dim == 2) ? Point<dim>(1., parameters.repetitions) :
                         Point<dim>(1., 1., parameters.repetitions),
            true);
          break;
        }
    }

  triangulation.refine_global(parameters.initial_refinement);
}



template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::setup_system()
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

  zero_constraints.clear();
  zero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

  nonzero_constraints.clear();
  nonzero_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

  if (parameters.problem_type == Settings::boundary_layer)
    {
      // Create zero BCs for the delta.
      // Left wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Right wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Top wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Bottom wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);

      // Nonzero constraints
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Functions::ConstantFunction<dim>(
                                                 -1.0),
                                               nonzero_constraints);

      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               Functions::ConstantFunction<dim>(
                                                 1.0),
                                               nonzero_constraints);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               nonzero_constraints);
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 2, BoundaryFunction<dim>(), nonzero_constraints);
    }
  else if (parameters.problem_type == Settings::double_glazing)
    {
      // Left wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Right wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Top wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Bottom wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);

      // Nonzero constraints
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Functions::ConstantFunction<dim>(
                                                 -1.0),
                                               nonzero_constraints);

      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               Functions::ConstantFunction<dim>(
                                                 1.0),
                                               nonzero_constraints);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(),
                                               nonzero_constraints);
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 2, BoundaryFunction<dim>(), nonzero_constraints);
    }
  else if (parameters.problem_type == Settings::boundary_layer_with_hole)
    {
      // Left wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);
      // Right wall
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(),
                                               zero_constraints);

      // Nonzero constraints
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 0, BoundaryValues<dim>(), nonzero_constraints);

      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 1, BoundaryValues<dim>(), nonzero_constraints);
    }

  DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);
  zero_constraints.close();

  DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
  nonzero_constraints.close();


  zero_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               mpi_communicator);

  nonzero_constraints.make_consistent_in_parallel(locally_owned_dofs,
                                                  locally_relevant_dofs,
                                                  mpi_communicator);

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(dof_handler, dsp, zero_constraints, false);
  SparsityTools::distribute_sparsity_pattern(dsp,
                                             dof_handler.locally_owned_dofs(),
                                             mpi_communicator,
                                             locally_relevant_dofs);
  zero_constraints.condense(dsp);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  // Instantiante local_solution to get solution without ghost elements
  VectorType local_solution(system_rhs);
  local_solution = solution;
  nonzero_constraints.distribute(local_solution);
  solution = local_solution;
}

template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::assemble_rhs()
{
  TimerOutput::Scope t(computing_timer, "assemble right hand side");

  system_rhs = 0;
  hp::FEValues<dim> hp_fe_values(fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                   update_hessians | update_JxW_values |
                                   update_quadrature_points);

  Vector<double>  cell_rhs;
  SourceTerm<dim> source_term;

  std::vector<types::global_dof_index> local_dof_indices;

  AdvectionField<dim> advection_field(parameters.problem_type);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          cell_rhs.reinit(dofs_per_cell);
          cell_rhs = 0.0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          const unsigned int n_q_points = fe_values.n_quadrature_points;

          std::vector<double>         source_term_values(n_q_points);
          std::vector<double>         newton_step_values(n_q_points);
          std::vector<Tensor<1, dim>> newton_step_gradients(n_q_points);
          std::vector<double>         newton_step_laplacians(n_q_points);
          std::vector<Tensor<1, dim>> advection_term_values(n_q_points);

          if (parameters.source_term == Settings::mms)
            source_term.value_list(fe_values.get_quadrature_points(),
                                   source_term_values);

          advection_field.value_list(fe_values.get_quadrature_points(),
                                     advection_term_values);

          fe_values.get_function_values(solution, newton_step_values);
          fe_values.get_function_gradients(solution, newton_step_gradients);

          if (parameters.stabilization)
            fe_values.get_function_laplacians(solution, newton_step_laplacians);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double dx = fe_values.JxW(q);
              const double nonlinearity =
                (parameters.nonlinearity ?
                   1e-01 * std::exp(newton_step_values[q]) :
                   0.0);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  if (parameters.source_term == Settings::mms)
                    cell_rhs(i) +=
                      (-1 / parameters.peclet_number * grad_phi_i *
                         newton_step_gradients[q] -
                       advection_term_values[q] * newton_step_gradients[q] *
                         phi_i +
                       phi_i * source_term_values[q] + phi_i * nonlinearity) *
                      dx;
                  else
                    cell_rhs(i) += (-1 / parameters.peclet_number * grad_phi_i *
                                      newton_step_gradients[q] -
                                    advection_term_values[q] *
                                      newton_step_gradients[q] * phi_i +
                                    phi_i * nonlinearity) *
                                   dx;

                  if (parameters.stabilization)
                    {
                      double h = 0;
                      if (dim == 2)
                        {
                          h = std::sqrt(4. * cell->measure() / M_PI) /
                              parameters.element_order;
                        }
                      else if (dim == 3)
                        {
                          h = std::pow(6 * cell->measure() / M_PI, 1. / 3.) /
                              parameters.element_order;
                        }

                      double tau = std::pow(
                        std::pow(1 / (parameters.peclet_number * h * h), 2) +
                          std::pow(2 * advection_term_values[q].norm() / h, 2),
                        -0.5);

                      cell_rhs(i) -=
                        ((-1 / parameters.peclet_number *
                          newton_step_laplacians[q]) +
                         (advection_term_values[q] * newton_step_gradients[q]) -
                         nonlinearity) *
                        (tau * (advection_term_values[q] * grad_phi_i)) * dx;
                    }
                }
            }
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
          zero_constraints.distribute_local_to_global(cell_rhs,
                                                      local_dof_indices,
                                                      system_rhs);
        }
    }

  system_rhs.compress(VectorOperation::add);
}

template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::assemble_matrix()
{
  TimerOutput::Scope t(computing_timer, "assemble matrix");

  system_matrix = 0;

  hp::FEValues<dim> hp_fe_values(fe_collection,
                                 quadrature_collection,
                                 update_values | update_gradients |
                                   update_hessians | update_JxW_values |
                                   update_quadrature_points);

  FullMatrix<double> cell_matrix;

  std::vector<types::global_dof_index> local_dof_indices;

  AdvectionField<dim> advection_field(parameters.problem_type);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
          cell_matrix = 0.0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          const unsigned int  n_q_points = fe_values.n_quadrature_points;
          std::vector<double> newton_step_values(n_q_points);
          std::vector<Tensor<1, dim>> advection_term_values(n_q_points);

          fe_values.get_function_values(solution, newton_step_values);

          advection_field.value_list(fe_values.get_quadrature_points(),
                                     advection_term_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double dx = fe_values.JxW(q);
              const double nonlinearity =
                (parameters.nonlinearity ?
                   1e-01 * std::exp(newton_step_values[q]) :
                   0.0);

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
                        (1 / parameters.peclet_number * grad_phi_i *
                           grad_phi_j +
                         advection_term_values[q] * grad_phi_j * phi_i -
                         phi_i * nonlinearity * phi_j) *
                        dx;

                      if (parameters.stabilization)
                        {
                          double h = 0;
                          if (dim == 2)
                            {
                              h = std::sqrt(4. * cell->measure() / M_PI) /
                                  parameters.element_order;
                            }
                          else if (dim == 3)
                            {
                              h =
                                std::pow(6 * cell->measure() / M_PI, 1. / 3.) /
                                parameters.element_order;
                            }

                          double tau = std::pow(
                            std::pow(1 / (parameters.peclet_number * h * h),
                                     2) +
                              std::pow(2 * advection_term_values[q].norm() / h,
                                       2),
                            -0.5);

                          auto shape_hessian_j = fe_values.shape_hessian(j, q);
                          auto shape_laplacian_j = trace(shape_hessian_j);
                          cell_matrix(i, j) +=
                            ((-1 / parameters.peclet_number *
                              shape_laplacian_j) +
                             (advection_term_values[q] * grad_phi_j) +
                             (phi_j * nonlinearity)) *
                            (tau * (advection_term_values[q] * grad_phi_i)) *
                            dx;
                        }
                    }
                }
            }
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
          zero_constraints.distribute_local_to_global(cell_matrix,
                                                      local_dof_indices,
                                                      system_matrix);
        }
    }

  system_matrix.compress(VectorOperation::add);
}

template <int dim, int fe_degree>
double
MatrixBasedAdvectionDiffusion<dim, fe_degree>::compute_residual(
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
                                   update_hessians | update_JxW_values |
                                   update_quadrature_points);

  Vector<double>  cell_residual;
  SourceTerm<dim> source_term;

  std::vector<types::global_dof_index> local_dof_indices;

  AdvectionField<dim> advection_field(parameters.problem_type);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
          cell_residual.reinit(dofs_per_cell);
          cell_residual = 0.0;

          hp_fe_values.reinit(cell);

          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

          const unsigned int  n_q_points = fe_values.n_quadrature_points;
          std::vector<double> source_term_values(n_q_points);
          std::vector<double> values(n_q_points);
          std::vector<Tensor<1, dim>> gradients(n_q_points);
          std::vector<double>         laplacians(n_q_points);
          std::vector<Tensor<1, dim>> advection_term_values(n_q_points);

          if (parameters.source_term == Settings::mms)
            source_term.value_list(fe_values.get_quadrature_points(),
                                   source_term_values);

          advection_field.value_list(fe_values.get_quadrature_points(),
                                     advection_term_values);

          fe_values.get_function_values(local_evaluation_point, values);
          fe_values.get_function_gradients(local_evaluation_point, gradients);
          if (parameters.stabilization)
            fe_values.get_function_laplacians(local_evaluation_point,
                                              laplacians);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double dx = fe_values.JxW(q);
              const double nonlinearity =
                (parameters.nonlinearity ? 1e-01 * std::exp(values[q]) : 0.0);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  const double         phi_i      = fe_values.shape_value(i, q);
                  const Tensor<1, dim> grad_phi_i = fe_values.shape_grad(i, q);

                  if (parameters.source_term == Settings::mms)
                    cell_residual(i) +=
                      (1 / parameters.peclet_number * grad_phi_i *
                         gradients[q] +
                       advection_term_values[q] * gradients[q] * phi_i +
                       phi_i * source_term_values[q] - phi_i * nonlinearity) *
                      dx;
                  else
                    cell_residual(i) +=
                      (1 / parameters.peclet_number * grad_phi_i *
                         gradients[q] +
                       advection_term_values[q] * gradients[q] * phi_i -
                       phi_i * nonlinearity) *
                      dx;

                  if (parameters.stabilization)
                    {
                      double h = 0;
                      if (dim == 2)
                        {
                          h = std::sqrt(4. * cell->measure() / M_PI) /
                              parameters.element_order;
                        }
                      else if (dim == 3)
                        {
                          h = std::pow(6 * cell->measure() / M_PI, 1. / 3.) /
                              parameters.element_order;
                        }

                      double tau = std::pow(
                        std::pow(1 / (parameters.peclet_number * h * h), 2) +
                          std::pow(2 * advection_term_values[q].norm() / h, 2),
                        -0.5);

                      cell_residual(i) +=
                        ((-1 / parameters.peclet_number * laplacians[q]) +
                         (advection_term_values[q] * gradients[q]) -
                         nonlinearity) *
                        (tau * (advection_term_values[q] * grad_phi_i)) * dx;
                    }
                }
            }
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);
          zero_constraints.distribute_local_to_global(cell_residual,
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
MatrixBasedAdvectionDiffusion<dim, fe_degree>::compute_update()
{
  TimerOutput::Scope t(computing_timer, "compute update");

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    locally_owned_dofs, mpi_communicator);

  SolverControl solver_control(1000, 1.e-4 * system_rhs.l2_norm(), true, true);
  TrilinosWrappers::SolverGMRES gmres(solver_control);

  switch (parameters.preconditioner)
    {
      case Settings::amg:
        {
          TrilinosWrappers::PreconditionAMG                 preconditioner;
          TrilinosWrappers::PreconditionAMG::AdditionalData data;

          if (fe_degree > 1)
            data.higher_order_elements = true;

          data.elliptic      = false;
          data.smoother_type = "Jacobi";

          preconditioner.initialize(system_matrix, data);

          gmres.solve(system_matrix,
                      completely_distributed_solution,
                      system_rhs,
                      preconditioner);
          break;
        }
      case Settings::ilu:
        {
          TrilinosWrappers::PreconditionILU                 preconditioner;
          TrilinosWrappers::PreconditionILU::AdditionalData data_ilu;
          preconditioner.initialize(system_matrix, data_ilu);

          gmres.solve(system_matrix,
                      completely_distributed_solution,
                      system_rhs,
                      preconditioner);
          break;
        }
      default:
        Assert(false,
               ExcMessage(
                 "This program supports only AMG and ILU as preconditioner."));
    }

  zero_constraints.distribute(completely_distributed_solution);

  linear_iterations = solver_control.last_step();

  newton_update = completely_distributed_solution;
}


template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::hp_refine()
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
MatrixBasedAdvectionDiffusion<dim, fe_degree>::solve()
{
  TimerOutput::Scope t(computing_timer, "solve");

  const unsigned int itmax = 14;
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
MatrixBasedAdvectionDiffusion<dim, fe_degree>::compute_solution_norm() const
{
  solution.update_ghost_values();

  Vector<float> norm_per_cell(triangulation.n_active_cells());

  return VectorTools::compute_global_error(triangulation,
                                           norm_per_cell,
                                           VectorTools::H1_seminorm);
}

template <int dim, int fe_degree>
double
MatrixBasedAdvectionDiffusion<dim, fe_degree>::compute_l2_error() const
{
  solution.update_ghost_values();

  Vector<float> error_per_cell(triangulation.n_active_cells());

  return VectorTools::compute_global_error(triangulation,
                                           error_per_cell,
                                           VectorTools::L2_norm);
}

template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::output_results(
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

  data_out.build_patches();

  DataOutBase::VtkFlags flags;
  flags.compression_level = DataOutBase::CompressionLevel::best_speed;
  data_out.set_flags(flags);

  std::string test_case = "";
  if (parameters.problem_type == Settings::boundary_layer)
    test_case = "_boundary_layer";
  else if (parameters.problem_type == Settings::double_glazing)
    test_case = "_double_glazing";

  if (parameters.stabilization)
    data_out.write_vtu_with_pvtu_record(
      parameters.output_path,
      parameters.output_name + std::to_string(dim) + "d_Pe" +
        std::to_string(parameters.peclet_number) + test_case + +"_stabilized",
      cycle,
      MPI_COMM_WORLD,
      3);
  else
    data_out.write_vtu_with_pvtu_record(
      parameters.output_path,
      parameters.output_name + std::to_string(dim) + "d_Pe" +
        std::to_string(parameters.peclet_number) + test_case,
      cycle,
      MPI_COMM_WORLD,
      3);
}

template <int dim, int fe_degree>
void
MatrixBasedAdvectionDiffusion<dim, fe_degree>::run()
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
    std::string SOL_header     = "Finite element space: "; // + fe.get_name();
    std::string PRECOND_header = "";
    if (parameters.preconditioner == Settings::amg)
      PRECOND_header = "Preconditioner: AMG";
    else if (parameters.preconditioner == Settings::ilu)
      PRECOND_header = "Preconditioner: ILU";
    std::string GEOMETRY_header = "";
    if (parameters.geometry == Settings::hyperball)
      GEOMETRY_header = "Geometry: hyperball";
    else if (parameters.geometry == Settings::hypercube)
      GEOMETRY_header = "Geometry: hypercube";
    else if (parameters.geometry == Settings::hypercube_with_hole)
      GEOMETRY_header = "Geometry: hypercube with hole";
    else if (parameters.geometry == Settings::hyperrectangle)
      GEOMETRY_header = "Geometry: hyperrectangle";
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
    else if (parameters.problem_type == Settings::boundary_layer_with_hole)
      PROBLEM_TYPE_header = "Test problem: boundary layer with hole";
    std::string STABILIZATION_header = "";
    if (parameters.stabilization)
      STABILIZATION_header = "Stabilization: true";
    else
      STABILIZATION_header = "Stabilization: false";
    std::string NONLINEARITY_header = "";
    if (parameters.nonlinearity)
      NONLINEARITY_header = "Nonlinearity: true";
    else
      NONLINEARITY_header = "Nonlinearity: false";
    std::string REPETITIONS_header =
      "Repetitions in one direction: " + std::to_string(parameters.repetitions);

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
    pcout << STABILIZATION_header << std::endl;
    pcout << NONLINEARITY_header << std::endl;
    pcout << REPETITIONS_header << std::endl;

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

      if (parameters.problem_type == Settings::boundary_layer &&
          parameters.geometry == Settings::hypercube)
        {
          pcout << "  L2 norm: " << compute_l2_error() << std::endl;
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
                      MatrixBasedAdvectionDiffusion<2, 1>
                        advection_diffusion_problem(parameters);
                      advection_diffusion_problem.run();

                      break;
                    }
                  case 2:
                    {
                      MatrixBasedAdvectionDiffusion<2, 2>
                        advection_diffusion_problem(parameters);
                      advection_diffusion_problem.run();

                      break;
                    }
                  case 3:
                    {
                      MatrixBasedAdvectionDiffusion<2, 3>
                        advection_diffusion_problem(parameters);
                      advection_diffusion_problem.run();

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
                      MatrixBasedAdvectionDiffusion<3, 1>
                        advection_diffusion_problem(parameters);
                      advection_diffusion_problem.run();

                      break;
                    }
                  case 2:
                    {
                      MatrixBasedAdvectionDiffusion<3, 2>
                        advection_diffusion_problem(parameters);
                      advection_diffusion_problem.run();

                      break;
                    }
                  case 3:
                    {
                      MatrixBasedAdvectionDiffusion<3, 3>
                        advection_diffusion_problem(parameters);
                      advection_diffusion_problem.run();

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
