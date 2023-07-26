/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2023 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
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
#include <deal.II/fe/fe_system.h>
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
#include <deal.II/lac/sparsity_tools.h>
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
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
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
    gcmg,
    lsmg,
    ilu,
    none
  };

  enum GeometryType
  {
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
                    "Order of FE element for u and p <1|2|3>");
  prm.declare_entry("number of cycles",
                    "1",
                    Patterns::Integer(),
                    "Number of cycles <1 up to 9-dim >");
  prm.declare_entry("geometry",
                    "hypercube",
                    Patterns::Selection("hypercube|hyperrectangle"),
                    "Geometry <hypercube|hyperrectangle>");
  prm.declare_entry("initial refinement",
                    "1",
                    Patterns::Integer(),
                    "Global refinement 1st cycle");
  prm.declare_entry("repetitions",
                    "1",
                    Patterns::Integer(),
                    "Repetitions in one direction for the hyperrectangle");
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
                    Patterns::Selection("AMG|LSMG|GCMG|ILU|none"),
                    "GMRES Preconditioner <AMG|LSMG|GCMG|ILU|none>");
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
  else if (prm.get("preconditioner") == "GCMG")
    this->preconditioner = gcmg;
  else if (prm.get("preconditioner") == "LSMG")
    this->preconditioner = lsmg;
  else if (prm.get("preconditioner") == "ILU")
    this->preconditioner = ilu;
  else if (prm.get("preconditioner") == "none")
    this->preconditioner = none;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("geometry") == "hypercube")
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

// Function for analytical solution in the case of MMS verification
template <int dim>
class AnalyticalSolution : public Function<dim>
{
public:
  AnalyticalSolution()
    : Function<dim>(dim + 1)
  {}

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &value) const override;
};

template <int dim>
void
AnalyticalSolution<dim>::vector_value(const Point<dim> &p,
                                      Vector<double> &  values) const
{
  AssertDimension(values.size(), dim + 1);

  // 2D Stokes analytical solution for the Rayleigh-Kothe-ish vortex
  // u=sin(a*x)*sin(a*x)*cos(a*y)*sin(a*y)
  // v=-cos(a*x)*sin(a*x)*sin(a*y)*sin(a*y)
  // p=sin(a*x)*sin(a*y)
  // with a = PI or any multiple of PI
  const double a = numbers::PI;
  const double x = p[0];
  const double y = p[1];

  values(0) =
    std::sin(a * x) * std::sin(a * x) * std::cos(a * y) * std::sin(a * y);
  values(1) =
    -std::cos(a * x) * std::sin(a * x) * std::sin(a * y) * std::sin(a * y);
  values(2) = std::sin(a * x) * std::sin(a * y);
}

// Function for the full source term
template <int dim>
class FullSourceTerm : public Function<dim>
{
public:
  FullSourceTerm()
  {}

  virtual double
  value(const Point<dim> &p, const unsigned int component = 0) const override;

  template <typename number>
  number
  value(const Point<dim, number> &p, const unsigned int component = 0) const;
};

template <int dim>
template <typename number>
number
FullSourceTerm<dim>::value(const Point<dim, number> &p,
                           const unsigned int        component) const
{
  const double a = numbers::PI;
  const double x = p[0];
  const double y = p[1];
  if (component == 0)
    {
      return (2 * a * a *
                (std::sin(a * x) * std::sin(a * x) -
                 std::cos(a * x) * std::cos(a * x)) *
                std::sin(a * y) * std::cos(a * y) +
              4 * a * a * std::sin(a * x) * std::sin(a * x) * std::sin(a * y) *
                std::cos(a * y) +
              a * std::sin(a * y) * std::cos(a * x));
    }
  else if (component == 1)
    {
      return (-2 * a * a *
                (std::sin(a * y) * std::sin(a * y) -
                 std::cos(a * y) * std::cos(a * y)) *
                std::sin(a * x) * std::cos(a * x) -
              4 * a * a * std::sin(a * x) * std::sin(a * y) * std::sin(a * y) *
                std::cos(a * x) +
              a * std::sin(a * x) * std::cos(a * y));
    }
  else
    return 0;
}

template <int dim>
double
FullSourceTerm<dim>::value(const Point<dim> & p,
                           const unsigned int component) const
{
  return value<double>(p, component);
}

// Matrix-free helper function
template <int dim, typename Number>
VectorizedArray<Number>
evaluate_function(const Function<dim> &                      function,
                  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  VectorizedArray<Number> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];
      result[v] = function.value(p);
    }
  return result;
}

// Matrix-free helper function
template <int dim, typename Number, int components>
Tensor<1, components, VectorizedArray<Number>>
evaluate_function(const Function<dim> &                      function,
                  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  Tensor<1, components, VectorizedArray<Number>> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (unsigned int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];
      for (unsigned int d = 0; d < components; ++d)
        result[d][v] = function.value(p, d);
    }
  return result;
}

// Matrix-free differential operator for a vector-valued problem.
// It "replaces" the traditional assemble_matrix() function.
template <int dim, typename number>
class StokesOperator : public Subscriptor
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;

  using VectorType = LinearAlgebra::distributed::Vector<number>;

  using value_type = number;

  using size_type = VectorizedArray<number>;

  StokesOperator();

  StokesOperator(const MappingQ<dim> &            mapping,
                 const DoFHandler<dim> &          dof_handler,
                 const AffineConstraints<number> &constraints,
                 const QGauss<1> &                quadrature,
                 const Settings &                 parameters,
                 const unsigned int               mg_level);

  void
  reinit(const MappingQ<dim> &            mapping,
         const DoFHandler<dim> &          dof_handler,
         const AffineConstraints<number> &constraints,
         const QGauss<1> &                quadrature,
         const Settings &                 parameters,
         const unsigned int               mg_level);

  void
  compute_stabilization_parameter();

  types::global_dof_index
  m() const;

  number
  el(unsigned int, unsigned int) const;

  void
  clear();

  void
  initialize_dof_vector(VectorType &vec) const;

  const std::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner() const;

  void
  vmult(VectorType &dst, const VectorType &src) const;

  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  void
  vmult_interface_down(VectorType &dst, VectorType const &src) const;

  void
  vmult_interface_up(VectorType &dst, VectorType const &src) const;

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const;

  const MatrixFree<dim, number> &
  get_system_matrix_free() const;

  const AlignedVector<VectorizedArray<number>>
  get_stabilization_parameter() const;

  void
  compute_inverse_diagonal(VectorType &diagonal) const;

private:
  void
  do_cell_integral_local(FECellIntegrator &integrator) const;

  void
  do_cell_integral_range(
    const MatrixFree<dim, number> &              matrix_free,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &range) const;

  static IndexSet
  get_refinement_edges(const MatrixFree<dim, number> &matrix_free);

private:
  MatrixFree<dim, number>                matrix_free;
  AffineConstraints<number>              constraints;
  mutable TrilinosWrappers::SparseMatrix system_matrix;
  Settings                               parameters;
  AlignedVector<VectorizedArray<number>> stabilization_parameter;

  // Variables needed for local smoothing
  std::vector<unsigned int>                      constrained_indices;
  mutable std::vector<std::pair<number, number>> constrained_values;
  std::vector<unsigned int>                      edge_constrained_indices;
  bool has_edge_constrained_indices = false;
  mutable std::vector<std::pair<number, number>> edge_constrained_values;
  std::vector<bool>                              edge_constrained_cell;
};

template <int dim, typename number>
StokesOperator<dim, number>::StokesOperator()
{}

template <int dim, typename number>
StokesOperator<dim, number>::StokesOperator(
  const MappingQ<dim> &            mapping,
  const DoFHandler<dim> &          dof_handler,
  const AffineConstraints<number> &constraints,
  const QGauss<1> &                quadrature,
  const Settings &                 parameters,
  const unsigned int               mg_level)
{
  this->reinit(
    mapping, dof_handler, constraints, quadrature, parameters, mg_level);
}

template <int dim, typename number>
void
StokesOperator<dim, number>::reinit(
  const MappingQ<dim> &            mapping,
  const DoFHandler<dim> &          dof_handler,
  const AffineConstraints<number> &constraints,
  const QGauss<1> &                quadrature,
  const Settings &                 parameters,
  const unsigned int               mg_level)
{
  this->system_matrix.clear();
  this->constraints.copy_from(constraints);

  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    (update_values | update_gradients | update_JxW_values |
     update_quadrature_points | update_hessians);
  additional_data.mg_level = mg_level;

  matrix_free.reinit(
    mapping, dof_handler, constraints, quadrature, additional_data);

  this->parameters = parameters;

  this->compute_stabilization_parameter();

  // Next lines required for local smoothing
  constrained_indices.clear();
  for (auto i : this->matrix_free.get_constrained_dofs())
    constrained_indices.push_back(i);
  constrained_values.resize(constrained_indices.size());

  if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
    {
      std::vector<types::global_dof_index> interface_indices;
      IndexSet                             refinement_edge_indices;
      refinement_edge_indices = get_refinement_edges(this->matrix_free);
      refinement_edge_indices.fill_index_vector(interface_indices);

      edge_constrained_indices.clear();
      edge_constrained_indices.reserve(interface_indices.size());
      edge_constrained_values.resize(interface_indices.size());
      const IndexSet &locally_owned =
        this->matrix_free.get_dof_handler().locally_owned_mg_dofs(
          this->matrix_free.get_mg_level());
      for (unsigned int i = 0; i < interface_indices.size(); ++i)
        if (locally_owned.is_element(interface_indices[i]))
          edge_constrained_indices.push_back(
            locally_owned.index_within_set(interface_indices[i]));

      this->has_edge_constrained_indices =
        Utilities::MPI::max(edge_constrained_indices.size(),
                            dof_handler.get_communicator()) > 0;

      if (this->has_edge_constrained_indices)
        {
          edge_constrained_cell.resize(matrix_free.n_cell_batches(), false);

          VectorType temp;
          matrix_free.initialize_dof_vector(temp);

          for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
            temp.local_element(edge_constrained_indices[i]) = 1.0;

          temp.update_ghost_values();

          FECellIntegrator integrator(matrix_free);

          for (unsigned int cell = 0; cell < matrix_free.n_cell_batches();
               ++cell)
            {
              integrator.reinit(cell);
              integrator.read_dof_values(temp);

              for (unsigned int i = 0; i < integrator.dofs_per_cell; ++i)
                if ((integrator.begin_dof_values()[i] ==
                     VectorizedArray<number>()) == false)
                  {
                    edge_constrained_cell[cell] = true;
                    break;
                  }
            }
        }
    }
}

template <int dim, typename number>
void
StokesOperator<dim, number>::compute_stabilization_parameter()
{
  unsigned int n_cells =
    matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
  stabilization_parameter.resize(n_cells);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (auto lane = 0u;
           lane < matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          double h_k = matrix_free.get_cell_iterator(cell, lane)->measure();

          if (dim == 2)
            stabilization_parameter[cell][lane] =
              std::sqrt(4. * h_k / M_PI) * std::sqrt(4. * h_k / M_PI);
          else if (dim == 3)
            stabilization_parameter[cell][lane] =
              std::pow(6 * h_k / M_PI, 1. / 3.) *
              std::pow(6 * h_k / M_PI, 1. / 3.);
        }
    }
}

template <int dim, typename number>
types::global_dof_index
StokesOperator<dim, number>::m() const
{
  if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
    return this->matrix_free.get_dof_handler().n_dofs(
      this->matrix_free.get_mg_level());
  else
    return this->matrix_free.get_dof_handler().n_dofs();
}

template <int dim, typename number>
number
StokesOperator<dim, number>::el(unsigned int, unsigned int) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}

template <int dim, typename number>
void
StokesOperator<dim, number>::clear()
{}

// The matrix free object initialized the vector accordingly
template <int dim, typename number>
void
StokesOperator<dim, number>::initialize_dof_vector(VectorType &vec) const
{
  matrix_free.initialize_dof_vector(vec);
}

template <int dim, typename number>
const std::shared_ptr<const Utilities::MPI::Partitioner> &
StokesOperator<dim, number>::get_vector_partitioner() const
{
  return matrix_free.get_vector_partitioner();
}

// Performs an operator evaluation by looping over all cells
// and evluating the integrals in the matrix-free way
template <int dim, typename number>
void
StokesOperator<dim, number>::vmult(VectorType &dst, const VectorType &src) const
{
  // save values for edge constrained dofs and set them to 0 in src vector
  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    {
      edge_constrained_values[i] = std::pair<number, number>(
        src.local_element(edge_constrained_indices[i]),
        dst.local_element(edge_constrained_indices[i]));

      const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
        .local_element(edge_constrained_indices[i]) = 0.;
    }

  this->matrix_free.cell_loop(
    &StokesOperator::do_cell_integral_range, this, dst, src, true);

  // set constrained dofs as the sum of current dst value and src value
  for (unsigned int i = 0; i < constrained_indices.size(); ++i)
    dst.local_element(constrained_indices[i]) =
      src.local_element(constrained_indices[i]);

  // restoring edge constrained dofs in src and dst
  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    {
      const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
        .local_element(edge_constrained_indices[i]) =
        edge_constrained_values[i].first;
      dst.local_element(edge_constrained_indices[i]) =
        edge_constrained_values[i].first;
    }
}

// Performs the transposed operator evaluation. Since we have
// non-symmetric matrices this is different from the vmult call.
// TODO: implement this correctly (let's start with something symetric)
template <int dim, typename number>
void
StokesOperator<dim, number>::Tvmult(VectorType &      dst,
                                    const VectorType &src) const
{
  this->vmult(dst, src);
}

template <int dim, typename number>
void
StokesOperator<dim, number>::vmult_interface_down(VectorType &      dst,
                                                  VectorType const &src) const
{
  this->matrix_free.cell_loop(
    &StokesOperator::do_cell_integral_range, this, dst, src, true);

  // set constrained dofs as the sum of current dst value and src value
  for (unsigned int i = 0; i < constrained_indices.size(); ++i)
    dst.local_element(constrained_indices[i]) =
      src.local_element(constrained_indices[i]);
}

template <int dim, typename number>
void
StokesOperator<dim, number>::vmult_interface_up(VectorType &      dst,
                                                VectorType const &src) const
{
  if (has_edge_constrained_indices == false)
    {
      dst = number(0.);
      return;
    }

  dst = 0.0;

  // make a copy of src vector and set everything to 0 except edge
  // constrained dofs
  VectorType src_cpy;
  src_cpy.reinit(src, /*omit_zeroing_entries=*/false);

  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    src_cpy.local_element(edge_constrained_indices[i]) =
      src.local_element(edge_constrained_indices[i]);

  // do loop with copy of src
  this->matrix_free.cell_loop(
    &StokesOperator::do_cell_integral_range, this, dst, src_cpy, false);
}


// Computes the matrix efficiently using the optimized compute_matrix call.
// This function is usually only used to build the system matrix for the coarse
// grid level in the multigrid algorithm.
template <int dim, typename number>
const TrilinosWrappers::SparseMatrix &
StokesOperator<dim, number>::get_system_matrix() const
{
  if (system_matrix.m() == 0 && system_matrix.n() == 0)
    {
      const auto &dof_handler = this->matrix_free.get_dof_handler();

      TrilinosWrappers::SparsityPattern dsp(
        this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int ?
          dof_handler.locally_owned_mg_dofs(this->matrix_free.get_mg_level()) :
          dof_handler.locally_owned_dofs(),
        dof_handler.get_triangulation().get_communicator());

      if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
        MGTools::make_sparsity_pattern(dof_handler,
                                       dsp,
                                       this->matrix_free.get_mg_level(),
                                       this->constraints);
      else
        DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints);

      dsp.compress();
      system_matrix.reinit(dsp);

      MatrixFreeTools::compute_matrix(matrix_free,
                                      constraints,
                                      system_matrix,
                                      &StokesOperator::do_cell_integral_local,
                                      this);
    }
  return this->system_matrix;
}

template <int dim, typename number>
const MatrixFree<dim, number> &
StokesOperator<dim, number>::get_system_matrix_free() const
{
  return this->matrix_free;
}

template <int dim, typename number>
const AlignedVector<VectorizedArray<number>>
StokesOperator<dim, number>::get_stabilization_parameter() const
{
  return this->stabilization_parameter;
}

// The diagonal of the matrix is needed, e.g., for smoothers used in
// multigrid. Thereforem it is computed by a sequence of operator
// evaluations to unit basis vectors using the optimized compute_diagonal
// call.
template <int dim, typename number>
void
StokesOperator<dim, number>::compute_inverse_diagonal(
  VectorType &diagonal) const
{
  this->matrix_free.initialize_dof_vector(diagonal);
  MatrixFreeTools::compute_diagonal(matrix_free,
                                    diagonal,
                                    &StokesOperator::do_cell_integral_local,
                                    this);

  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    diagonal.local_element(edge_constrained_indices[i]) = 0.0;

  for (auto &i : diagonal)
    i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0;
}

template <int dim, typename number>
void
StokesOperator<dim, number>::do_cell_integral_local(
  FECellIntegrator &integrator) const
{
  using FECellIntegratorType =
    FEEvaluation<dim, -1, 0, dim + 1, double, VectorizedArray<double>>;

  FullSourceTerm<dim> source_term_function;

  integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                      EvaluationFlags::hessians);

  auto tau = integrator.read_cell_data(stabilization_parameter);

  for (unsigned int q = 0; q < integrator.n_q_points; ++q)
    {
      // Evaluate source term function
      Tensor<1, dim + 1, VectorizedArray<double>> source_value;

      if (parameters.source_term == Settings::mms)
        {
          Point<dim, VectorizedArray<double>> point_batch =
            integrator.quadrature_point(q);
          source_value =
            evaluate_function<dim, double, dim + 1>(source_term_function,
                                                    point_batch);
        }

      // Gather the original value/gradient
      typename FECellIntegratorType::value_type value = integrator.get_value(q);
      typename FECellIntegratorType::gradient_type gradient =
        integrator.get_gradient(q);
      typename FECellIntegratorType::gradient_type hessian_diagonal =
        integrator.get_hessian_diagonal(q);

      // Result value/gradient we will use
      typename FECellIntegratorType::value_type    value_result;
      typename FECellIntegratorType::gradient_type gradient_result;

      // Assemble -nabla u + nabla p = 0 for the first 3 components
      // The corresponding weak form is nabla v * nabla u  - p nabla \cdot v = 0
      // ; Assemble q div(u) = 0 for the last component

      for (unsigned int i = 0; i < dim; ++i)
        {
          gradient_result[i] = gradient[i];
          gradient_result[i][i] += -value[dim];

          value_result[dim] += gradient[i][i];
        }

      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            gradient_result[dim][i] += -tau * hessian_diagonal[i][k];
        }
      gradient_result[dim] += tau * gradient[dim];

      integrator.submit_gradient(gradient_result, q);
      integrator.submit_value(value_result, q);
    }

  integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

template <int dim, typename number>
void
StokesOperator<dim, number>::do_cell_integral_range(
  const MatrixFree<dim, number> &              matrix_free,
  VectorType &                                 dst,
  const VectorType &                           src,
  const std::pair<unsigned int, unsigned int> &range) const
{
  FECellIntegrator integrator(matrix_free, range);

  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      integrator.reinit(cell);

      integrator.read_dof_values(src);

      do_cell_integral_local(integrator);

      integrator.distribute_local_to_global(dst);
    }
}

template <int dim, typename number>
IndexSet
StokesOperator<dim, number>::get_refinement_edges(
  const MatrixFree<dim, number> &matrix_free)
{
  const unsigned int level = matrix_free.get_mg_level();

  std::vector<IndexSet> refinement_edge_indices;
  refinement_edge_indices.clear();
  const unsigned int nlevels =
    matrix_free.get_dof_handler().get_triangulation().n_global_levels();
  refinement_edge_indices.resize(nlevels);
  for (unsigned int l = 0; l < nlevels; l++)
    refinement_edge_indices[l] =
      IndexSet(matrix_free.get_dof_handler().n_dofs(l));

  MGTools::extract_inner_interface_dofs(matrix_free.get_dof_handler(),
                                        refinement_edge_indices);
  return refinement_edge_indices[level];
}

// Multigrid global coarsening algorithm and parameters
struct MultigridParameters
{
  struct
  {
    std::string  type            = "gmres_with_amg";
    unsigned int maxiter         = 2000;
    double       abstol          = 1e-14;
    double       reltol          = 1e-4;
    unsigned int smoother_sweeps = 1;
    unsigned int n_cycles        = 1;
    std::string  smoother_type   = "ILU";
  } coarse_solver;

  struct
  {
    std::string  type         = "relaxation";
    unsigned int n_iterations = 10;
    double       relaxation   = 0.5;
  } smoother;
};

template <typename VectorType,
          int dim,
          typename SystemMatrixType,
          typename LevelMatrixType,
          typename MGTransferType>
static void
gcmg_solve(SolverControl &            solver_control,
           VectorType &               dst,
           const VectorType &         src,
           const MultigridParameters &mg_data,
           const DoFHandler<dim> &    dof,
           const SystemMatrixType &   fine_matrix,
           const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices,
           const MGTransferType &                                 mg_transfer)
{
  AssertThrow(mg_data.coarse_solver.type == "gmres_with_amg",
              ExcNotImplemented());
  AssertThrow(mg_data.smoother.type == "relaxation", ExcNotImplemented());

  const unsigned int min_level = mg_matrices.min_level();
  const unsigned int max_level = mg_matrices.max_level();

  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType =
    PreconditionRelaxation<LevelMatrixType, SmootherPreconditionerType>;
  using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;

  mg::Matrix<VectorType> mg_matrix(mg_matrices);

  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(min_level,
                                                                     max_level);

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>();
      mg_matrices[level]->compute_inverse_diagonal(
        smoother_data[level].preconditioner->get_vector());
      smoother_data[level].n_iterations = mg_data.smoother.n_iterations;
      smoother_data[level].relaxation   = mg_data.smoother.relaxation;
    }

  MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType> mg_smoother;
  mg_smoother.initialize(mg_matrices, smoother_data);

  ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter,
                                              mg_data.coarse_solver.abstol,
                                              mg_data.coarse_solver.reltol,
                                              false,
                                              false);
  SolverGMRES<VectorType> coarse_grid_solver(coarse_grid_solver_control);

  std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  TrilinosWrappers::PreconditionAMG                 precondition_amg;
  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps;
  amg_data.n_cycles        = mg_data.coarse_solver.n_cycles;
  amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str();

  precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(),
                              amg_data);

  mg_coarse =
    std::make_unique<MGCoarseGridIterativeSolver<VectorType,
                                                 SolverGMRES<VectorType>,
                                                 LevelMatrixType,
                                                 decltype(precondition_amg)>>(
      coarse_grid_solver, *mg_matrices[min_level], precondition_amg);

  Multigrid<VectorType> mg(
    mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother);

  PreconditionerType preconditioner(dof, mg, mg_transfer);

  SolverGMRES<VectorType>(solver_control)
    .solve(fine_matrix, dst, src, preconditioner);
}


template <typename VectorType, typename OperatorType, int dim>
void
solve_with_gcmg(SolverControl &            solver_control,
                const OperatorType &       system_matrix,
                VectorType &               dst,
                const VectorType &         src,
                const MultigridParameters &mg_data,
                const MappingQ<dim>        mapping,
                const DoFHandler<dim> &    dof_handler,
                const QGauss<1> &          quadrature,
                const Settings &           parameters)
{
  MGLevelObject<DoFHandler<dim>>                     dof_handlers;
  MGLevelObject<std::unique_ptr<OperatorType>>       operators;
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;

  std::vector<std::shared_ptr<const Triangulation<dim>>>
    coarse_grid_triangulations;

  coarse_grid_triangulations =
    MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      dof_handler.get_triangulation());

  const unsigned int n_h_levels = coarse_grid_triangulations.size();

  unsigned int minlevel = 0;
  unsigned int maxlevel = n_h_levels - 1;

  dof_handlers.resize(minlevel, maxlevel);
  operators.resize(minlevel, maxlevel);
  transfers.resize(minlevel, maxlevel);

  for (unsigned int l = minlevel; l <= maxlevel; ++l)
    {
      dof_handlers[l].reinit(*coarse_grid_triangulations[l]);
      dof_handlers[l].distribute_dofs(dof_handler.get_fe());
    }

  MGLevelObject<AffineConstraints<typename VectorType::value_type>> constraints(
    minlevel, maxlevel);

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      const auto &dof_handler = dof_handlers[level];
      auto &      constraint  = constraints[level];

      const IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraint.reinit(locally_relevant_dofs);

      DoFTools::make_hanging_node_constraints(dof_handler, constraint);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               0,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraint);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               1,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraint);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               2,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraint);
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               3,
                                               Functions::ZeroFunction<dim>(
                                                 dim + 1),
                                               constraint);

      constraint.close();

      operators[level] =
        std::make_unique<OperatorType>(mapping,
                                       dof_handler,
                                       constraint,
                                       quadrature,
                                       parameters,
                                       numbers::invalid_unsigned_int);
    }

  for (unsigned int level = minlevel; level < maxlevel; ++level)
    transfers[level + 1].reinit(dof_handlers[level + 1],
                                dof_handlers[level],
                                constraints[level + 1],
                                constraints[level]);

  MGTransferGlobalCoarsening<dim, VectorType> transfer(
    transfers,
    [&](const auto l, auto &vec) { operators[l]->initialize_dof_vector(vec); });

  ConditionalOStream pcout(std::cout,
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    pcout << "   MG Level " << level << ": " << dof_handlers[level].n_dofs()
          << " DoFs, "
          << coarse_grid_triangulations[level]->n_global_active_cells()
          << " cells" << std::endl;

  gcmg_solve(solver_control,
             dst,
             src,
             mg_data,
             dof_handler,
             system_matrix,
             operators,
             transfer);
}

template <typename VectorType,
          int dim,
          typename SystemMatrixType,
          typename LevelMatrixType,
          typename MGTransferType>
static void
lsmg_solve(
  SolverControl &                                        solver_control,
  VectorType &                                           dst,
  const VectorType &                                     src,
  const MultigridParameters &                            mg_data,
  const DoFHandler<dim> &                                dof,
  const SystemMatrixType &                               fine_matrix,
  const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices,
  const MGTransferType &                                 mg_transfer,
  const MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
    &ls_mg_interface_in,
  const MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
    &ls_mg_interface_out,
  const MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
    &ls_mg_matrices)
{
  AssertThrow(mg_data.coarse_solver.type == "gmres_with_amg",
              ExcNotImplemented());
  AssertThrow(mg_data.smoother.type == "relaxation", ExcNotImplemented());

  const unsigned int min_level = mg_matrices.min_level();
  const unsigned int max_level = mg_matrices.max_level();

  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType =
    PreconditionRelaxation<LevelMatrixType, SmootherPreconditionerType>;
  using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>;

  mg::Matrix<VectorType> mg_matrix(ls_mg_matrices);

  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(min_level,
                                                                     max_level);

  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>();
      mg_matrices[level]->compute_inverse_diagonal(
        smoother_data[level].preconditioner->get_vector());
      smoother_data[level].n_iterations = mg_data.smoother.n_iterations;
      smoother_data[level].relaxation   = mg_data.smoother.relaxation;
    }

  MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType> mg_smoother;
  mg_smoother.initialize(mg_matrices, smoother_data);

  ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter,
                                              mg_data.coarse_solver.abstol,
                                              mg_data.coarse_solver.reltol,
                                              false,
                                              false);
  SolverGMRES<VectorType> coarse_grid_solver(coarse_grid_solver_control);

  std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  TrilinosWrappers::PreconditionAMG                 precondition_amg;
  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps;
  amg_data.n_cycles        = mg_data.coarse_solver.n_cycles;
  amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str();

  precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(),
                              amg_data);

  mg_coarse =
    std::make_unique<MGCoarseGridIterativeSolver<VectorType,
                                                 SolverGMRES<VectorType>,
                                                 LevelMatrixType,
                                                 decltype(precondition_amg)>>(
      coarse_grid_solver, *mg_matrices[min_level], precondition_amg);

  // Additional step for local smoothing: interface matrices
  mg::Matrix<LinearAlgebra::distributed::Vector<double>> mg_interface_matrix_in(
    ls_mg_interface_in);
  mg::Matrix<LinearAlgebra::distributed::Vector<double>>
    mg_interface_matrix_out(ls_mg_interface_out);

  Multigrid<VectorType> mg(
    mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother);

  if (dof.get_triangulation().has_hanging_nodes())
    mg.set_edge_matrices(mg_interface_matrix_in, mg_interface_matrix_out);

  PreconditionerType preconditioner(dof, mg, mg_transfer);

  // PreconditionIdentity identity;
  SolverGMRES<VectorType>(solver_control)
    .solve(fine_matrix, dst, src, preconditioner);
}

template <typename VectorType, typename OperatorType, int dim>
void
solve_with_lsmg(SolverControl &            solver_control,
                const OperatorType &       system_matrix,
                VectorType &               dst,
                const VectorType &         src,
                const MultigridParameters &mg_data,
                const MappingQ<dim>        mapping,
                const DoFHandler<dim> &    dof_handler,
                const QGauss<1> &          quadrature,
                const Settings &           parameters)
{
  MGLevelObject<std::unique_ptr<OperatorType>> operators;
  MGTransferMatrixFree<dim, double>            mg_transfer;
  MGLevelObject<VectorType>                    mg_solution;
  MGLevelObject<AffineConstraints<double>>     level_constraints;
  MGConstrainedDoFs                            mg_constrained_dofs;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_operators;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_in;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_out;

  const unsigned int n_h_levels =
    dof_handler.get_triangulation().n_global_levels();

  unsigned int minlevel = 0;
  unsigned int maxlevel = n_h_levels - 1;

  operators.resize(0, n_h_levels - 1);
  mg_solution.resize(0, n_h_levels - 1);
  level_constraints.resize(0, n_h_levels - 1);
  ls_mg_interface_in.resize(0, n_h_levels - 1);
  ls_mg_interface_out.resize(0, n_h_levels - 1);
  ls_mg_operators.resize(0, n_h_levels - 1);

  const std::set<types::boundary_id> dirichlet_boundary_ids = {0, 1, 2, 3};
  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler,
                                                     dirichlet_boundary_ids);

  std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>> partitioners(
    dof_handler.get_triangulation().n_global_levels());

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    {
      const IndexSet relevant_dofs =
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

      // AffineConstraints<double> level_constraints;
      level_constraints[level].reinit(relevant_dofs);
      level_constraints[level].add_lines(
        mg_constrained_dofs.get_boundary_indices(level));
      level_constraints[level].close();

      operators[level] =
        std::make_unique<OperatorType>(mapping,
                                       dof_handler,
                                       level_constraints[level],
                                       quadrature,
                                       parameters,
                                       level);

      operators[level]->initialize_dof_vector(mg_solution[level]);

      ls_mg_operators[level].initialize(*operators[level]);
      ls_mg_interface_in[level].initialize(*operators[level]);
      ls_mg_interface_out[level].initialize(*operators[level]);

      partitioners[level] = operators[level]->get_vector_partitioner();
    }

  mg_transfer.initialize_constraints(mg_constrained_dofs);
  mg_transfer.build(dof_handler, partitioners);
  mg_transfer.interpolate_to_mg(dof_handler, mg_solution, dst);


  ConditionalOStream pcout(std::cout,
                           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                            0));

  for (unsigned int level = minlevel; level <= maxlevel; ++level)
    pcout << "   MG Level " << level << ": " << dof_handler.n_dofs(level)
          << " DoFs, " << dof_handler.get_triangulation().n_cells(level)
          << " cells" << std::endl;

  lsmg_solve(solver_control,
             dst,
             src,
             mg_data,
             dof_handler,
             system_matrix,
             operators,
             mg_transfer,
             ls_mg_interface_in,
             ls_mg_interface_out,
             ls_mg_operators);
}

// Main class to solve the stokes vector-valued problem given by
// -∇^2 u + ∇p = f_1 and ∇ · u = f_2 using Newton's method
// and the matrix-free approach.
template <int dim>
class MatrixFreeStokes
{
public:
  MatrixFreeStokes(const Settings &parameters);

  void
  run();

private:
  void
  make_grid();

  void
  refine_grid();

  void
  setup_system();

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

  void
  compute_solution_norm() const;

  void
  compute_l2_error() const;

  void
  output_results(const unsigned int cycle) const;


  parallel::distributed::Triangulation<dim> triangulation;

  Settings parameters;

  const unsigned int element_order;

  const MappingQ<dim> mapping;

  FESystem<dim>   fe_system;
  DoFHandler<dim> dof_handler;

  AffineConstraints<double> constraints;
  using SystemMatrixType = StokesOperator<dim, double>;
  SystemMatrixType  system_matrix;
  MGConstrainedDoFs mg_constrained_dofs;
  using LevelMatrixType = StokesOperator<dim, double>;

  TrilinosWrappers::SparseMatrix mb_system_matrix;
  TrilinosWrappers::SparseMatrix mb_system_matrix_coarse;

  LinearAlgebra::distributed::Vector<double> solution;
  LinearAlgebra::distributed::Vector<double> newton_update;
  LinearAlgebra::distributed::Vector<double> system_rhs;

  unsigned int       linear_iterations;
  ConditionalOStream pcout;
  TimerOutput        computing_timer;
};

template <int dim>
MatrixFreeStokes<dim>::MatrixFreeStokes(const Settings &parameters)
  : triangulation(
      MPI_COMM_WORLD,
      Triangulation<dim>::limit_level_difference_at_vertices,
      parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy)
  , parameters(parameters)
  , element_order(parameters.element_order)
  , mapping(parameters.element_order)
  , fe_system(FE_Q<dim>(parameters.element_order), dim + 1)
  , dof_handler(triangulation)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , computing_timer(MPI_COMM_WORLD,
                    pcout,
                    TimerOutput::never,
                    TimerOutput::wall_times)
{}

template <int dim>
void
MatrixFreeStokes<dim>::make_grid()
{
  TimerOutput::Scope t(computing_timer, "make grid");

  switch (parameters.geometry)
    {
        case Settings::hypercube: {
          GridGenerator::hyper_cube(triangulation, -1.0, 1.0, true);
          break;
        }
        case Settings::hyperrectangle: {
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
            (dim == 2) ? Point<dim>(1., parameters.repetitions) :
                         Point<dim>(1., 1., parameters.repetitions),
            true);
          break;
        }
    }

  triangulation.refine_global(parameters.initial_refinement);
}

// To test locally-refined meshes
template <int dim>
void
MatrixFreeStokes<dim>::refine_grid()
{
  unsigned int n_refinements = 2;

  for (unsigned int i = 1; i < n_refinements; i++)
    {
      for (auto cell : triangulation.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool flag = true;
            for (int d = 0; d < dim; d++)
              if (cell->center()[d] > 0.0)
                flag = false;
            if (flag)
              cell->set_refine_flag();
          }
      triangulation.execute_coarsening_and_refinement();
    }
}

template <int dim>
void
MatrixFreeStokes<dim>::setup_system()
{
  TimerOutput::Scope t(computing_timer, "setup system");

  dof_handler.distribute_dofs(fe_system);

  if (parameters.preconditioner == Settings::lsmg)
    dof_handler.distribute_mg_dofs();

  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  // Set homogeneous constraints for the matrix-free operator
  // Create zero BCs for the delta.
  // Left wall
  VectorTools::interpolate_boundary_values(
    dof_handler, 0, Functions::ZeroFunction<dim>(dim + 1), constraints);

  VectorTools::interpolate_boundary_values(
    dof_handler, 1, Functions::ZeroFunction<dim>(dim + 1), constraints);

  VectorTools::interpolate_boundary_values(
    dof_handler, 2, Functions::ZeroFunction<dim>(dim + 1), constraints);

  VectorTools::interpolate_boundary_values(
    dof_handler, 3, Functions::ZeroFunction<dim>(dim + 1), constraints);

  constraints.close();

  system_matrix.reinit(mapping,
                       dof_handler,
                       constraints,
                       QGauss<1>(element_order + 1),
                       parameters,
                       numbers::invalid_unsigned_int);

  system_matrix.initialize_dof_vector(solution);
  system_matrix.initialize_dof_vector(newton_update);
  system_matrix.initialize_dof_vector(system_rhs);
}

template <int dim>
void
MatrixFreeStokes<dim>::evaluate_residual(
  LinearAlgebra::distributed::Vector<double> &      dst,
  const LinearAlgebra::distributed::Vector<double> &src) const
{
  system_matrix.get_system_matrix_free().cell_loop(
    &MatrixFreeStokes::local_evaluate_residual, this, dst, src, true);
}

template <int dim>
void
MatrixFreeStokes<dim>::local_evaluate_residual(
  const MatrixFree<dim, double> &                   data,
  LinearAlgebra::distributed::Vector<double> &      dst,
  const LinearAlgebra::distributed::Vector<double> &src,
  const std::pair<unsigned int, unsigned int> &     cell_range) const
{
  using FECellIntegratorType =
    FEEvaluation<dim, -1, 0, dim + 1, double, VectorizedArray<double>>;

  FECellIntegratorType phi(data);

  FullSourceTerm<dim> source_term_function;

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values_plain(src);
      phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                   EvaluationFlags::hessians);
      auto tau =
        phi.read_cell_data(system_matrix.get_stabilization_parameter());

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          // Evaluate source term function
          Tensor<1, dim + 1, VectorizedArray<double>> source_value;

          if (parameters.source_term == Settings::mms)
            {
              Point<dim, VectorizedArray<double>> point_batch =
                phi.quadrature_point(q);
              source_value =
                evaluate_function<dim, double, dim + 1>(source_term_function,
                                                        point_batch);
            }

          // Gather the original value/gradient
          typename FECellIntegratorType::value_type    value = phi.get_value(q);
          typename FECellIntegratorType::gradient_type gradient =
            phi.get_gradient(q);
          typename FECellIntegratorType::gradient_type hessian_diagonal =
            phi.get_hessian_diagonal(q);

          // Result value/gradient we will use
          typename FECellIntegratorType::value_type    value_result;
          typename FECellIntegratorType::gradient_type gradient_result;

          // Assemble -nabla^2 u + nabla p = 0 for the first 3 components
          // The corresponding weak form is nabla v * nabla u  - p nabla \cdot v
          // = 0 ; Assemble q div(u) = 0 for the last component
          for (unsigned int i = 0; i < dim; ++i)
            {
              gradient_result[i] = gradient[i];
              gradient_result[i][i] += -value[dim];
              value_result[i] = -source_value[i];

              value_result[dim] += gradient[i][i];
            }

          for (unsigned int i = 0; i < dim; ++i)
            {
              gradient_result[dim][i] += -tau * source_value[i];
              for (unsigned int k = 0; k < dim; ++k)
                gradient_result[dim][i] += -tau * hessian_diagonal[i][k];
            }
          gradient_result[dim] += tau * gradient[dim];

          phi.submit_gradient(gradient_result, q);
          phi.submit_value(value_result, q);
        }

      phi.integrate_scatter(EvaluationFlags::values |
                              EvaluationFlags::gradients,
                            dst);
    }
}

template <int dim>
void
MatrixFreeStokes<dim>::assemble_rhs()
{
  TimerOutput::Scope t(computing_timer, "assemble right hand side");

  evaluate_residual(system_rhs, solution);

  system_rhs *= -1.0;
}

template <int dim>
double
MatrixFreeStokes<dim>::compute_residual(const double alpha)
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

template <int dim>
void
MatrixFreeStokes<dim>::compute_update()
{
  TimerOutput::Scope t(computing_timer, "compute update");

  solution.update_ghost_values();

  SolverControl solver_control(10000,
                               std::max(1.e-4 * system_rhs.l2_norm(), 1e-14),
                               true,
                               true);
  SolverGMRES<LinearAlgebra::distributed::Vector<double>>::AdditionalData
    gmres_parameters;
  gmres_parameters.max_n_tmp_vectors = 1000;
  SolverGMRES<LinearAlgebra::distributed::Vector<double>> gmres(
    solver_control, gmres_parameters);

  newton_update = 0.0;

  switch (parameters.preconditioner)
    {
        case Settings::amg: {
          TrilinosWrappers::PreconditionAMG                 preconditioner;
          TrilinosWrappers::PreconditionAMG::AdditionalData data;

          if (parameters.element_order > 1)
            data.higher_order_elements = true;

          data.elliptic      = false;
          data.smoother_type = "Jacobi";

          preconditioner.initialize(system_matrix.get_system_matrix(), data);
          gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
          break;
        }
        case Settings::gcmg: {
          MultigridParameters mg_data;

          solve_with_gcmg(solver_control,
                          system_matrix,
                          newton_update,
                          system_rhs,
                          mg_data,
                          mapping,
                          dof_handler,
                          QGauss<1>(element_order + 1),
                          parameters);
          break;
        }
        case Settings::lsmg: {
          MultigridParameters mg_data;

          solve_with_lsmg(solver_control,
                          system_matrix,
                          newton_update,
                          system_rhs,
                          mg_data,
                          mapping,
                          dof_handler,
                          QGauss<1>(element_order + 1),
                          parameters);
          break;
        }
        case Settings::ilu: {
          TrilinosWrappers::PreconditionILU                 preconditioner;
          TrilinosWrappers::PreconditionILU::AdditionalData data_ilu;
          preconditioner.initialize(system_matrix.get_system_matrix(),
                                    data_ilu);

          gmres.solve(system_matrix, newton_update, system_rhs, preconditioner);
          break;
        }
        case Settings::none: {
          PreconditionIdentity preconditioner;
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


template <int dim>
void
MatrixFreeStokes<dim>::solve()
{
  TimerOutput::Scope t(computing_timer, "solve");

  const unsigned int itmax = 14;
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

template <int dim>
void
MatrixFreeStokes<dim>::compute_solution_norm() const
{
  solution.update_ghost_values();

  const ComponentSelectFunction<dim> p_mask(dim, dim + 1);
  const ComponentSelectFunction<dim> u_mask(std::make_pair(0, dim), dim + 1);

  Vector<float> norm_per_cell(triangulation.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    Functions::ZeroFunction<dim>(dim + 1),
                                    norm_per_cell,
                                    QGauss<dim>(element_order + 1),
                                    VectorTools::H1_seminorm,
                                    &u_mask);

  solution.zero_out_ghost_values();

  const double u_h1_norm =
    VectorTools::compute_global_error(triangulation,
                                      norm_per_cell,
                                      VectorTools::H1_seminorm);

  solution.update_ghost_values();

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    Functions::ZeroFunction<dim>(dim + 1),
                                    norm_per_cell,
                                    QGauss<dim>(element_order + 1),
                                    VectorTools::H1_seminorm,
                                    &p_mask);

  solution.zero_out_ghost_values();

  const double p_h1_norm =
    VectorTools::compute_global_error(triangulation,
                                      norm_per_cell,
                                      VectorTools::H1_seminorm);

  pcout << "  u H1 seminorm: " << u_h1_norm << std::endl;
  pcout << "  p H1 seminorm: " << p_h1_norm << std::endl;
  pcout << std::endl;
}

template <int dim>
void
MatrixFreeStokes<dim>::compute_l2_error() const
{
  solution.update_ghost_values();

  const ComponentSelectFunction<dim> p_mask(dim, dim + 1);
  const ComponentSelectFunction<dim> u_mask(std::make_pair(0, dim), dim + 1);

  Vector<float> error_per_cell(triangulation.n_active_cells());

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    AnalyticalSolution<dim>(),
                                    error_per_cell,
                                    QGauss<dim>(element_order + 1),
                                    VectorTools::L2_norm,
                                    &u_mask);

  solution.zero_out_ghost_values();

  const double u_l2_error =
    VectorTools::compute_global_error(triangulation,
                                      error_per_cell,
                                      VectorTools::L2_norm);

  solution.update_ghost_values();

  VectorTools::integrate_difference(mapping,
                                    dof_handler,
                                    solution,
                                    AnalyticalSolution<dim>(),
                                    error_per_cell,
                                    QGauss<dim>(element_order + 1),
                                    VectorTools::L2_norm,
                                    &p_mask);

  solution.zero_out_ghost_values();

  const double p_l2_error =
    VectorTools::compute_global_error(triangulation,
                                      error_per_cell,
                                      VectorTools::L2_norm);

  pcout << "  u L2 norm error: " << u_l2_error << std::endl;
  pcout << "  p L2 norm error: " << p_l2_error << std::endl;
}

template <int dim>
void
MatrixFreeStokes<dim>::output_results(const unsigned int cycle) const
{
  if (triangulation.n_global_active_cells() > 1e6)
    return;

  std::vector<std::string> solution_names(dim, "u");
  solution_names.push_back("p");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    dealii_data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  dealii_data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution,
                           solution_names,
                           DataOut<dim>::type_dof_data,
                           dealii_data_component_interpretation);

  Vector<float> subdomain(triangulation.n_active_cells());
  subdomain = triangulation.locally_owned_subdomain();

  data_out.add_data_vector(subdomain, "subdomain");

  data_out.build_patches();

  DataOutBase::VtkFlags flags;
  flags.compression_level = DataOutBase::VtkFlags::best_speed;
  data_out.set_flags(flags);

  data_out.write_vtu_with_pvtu_record(parameters.output_path,
                                      parameters.output_name +
                                        std::to_string(dim),
                                      cycle,
                                      MPI_COMM_WORLD,
                                      3);
}

template <int dim>
void
MatrixFreeStokes<dim>::run()
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
    std::string SOL_header = "Finite element space: " + fe_system.get_name();
    std::string PRECOND_header = "";
    if (parameters.preconditioner == Settings::amg)
      PRECOND_header = "Preconditioner: AMG";
    else if (parameters.preconditioner == Settings::gcmg)
      PRECOND_header = "Preconditioner: GCMG";
    else if (parameters.preconditioner == Settings::lsmg)
      PRECOND_header = "Preconditioner: LSMG";
    else if (parameters.preconditioner == Settings::ilu)
      PRECOND_header = "Preconditioner: ILU";
    else if (parameters.preconditioner == Settings::none)
      PRECOND_header = "Preconditioner: none";
    std::string GEOMETRY_header = "";
    if (parameters.geometry == Settings::hypercube)
      GEOMETRY_header = "Geometry: hypercube";
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
      else
        {
          triangulation.refine_global(1);
          // To test locally refined mesh for local smoothing:
          // refine_grid();
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
        }

      // compute_solution_norm();

      compute_l2_error();

      computing_timer.print_summary();
      computing_timer.reset();
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

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
              MatrixFreeStokes<2> vector_valued_problem(parameters);
              vector_valued_problem.run();

              break;
            }

            case 3: {
              MatrixFreeStokes<3> vector_valued_problem(parameters);
              vector_valued_problem.run();
              break;
            }

          default:
            Assert(
              false,
              ExcMessage(
                "This program only works in 2d and 3d and for element orders equal to 1, 2 or 3."));
        }
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
