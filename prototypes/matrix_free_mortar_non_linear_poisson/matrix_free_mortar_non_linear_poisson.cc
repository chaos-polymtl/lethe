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
#include <deal.II/lac/vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_remote_evaluation.h>
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

// Define all the parameters that can be specified in the .prm file
struct Settings
{
  bool
  try_parse(const std::string &prm_filename);

  enum PreconditionerType
  {
    gmg
  };

  enum GeometryType
  {
    hypercube,
    hyperrectangle,
    cocircle,
    corectangle
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
                    "cocircle",
                    Patterns::Selection("cocircle|corectangle"),
                    "Geometry <cocircle|corectangle>");
  prm.declare_entry("initial refinement",
                    "1",
                    Patterns::Integer(),
                    "Global refinement 1st cycle");
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
                    "GMG",
                    Patterns::Selection("GMG"),
                    "Preconditioner <GMG>");
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

  if (prm.get("preconditioner") == "GMG")
    this->preconditioner = gmg;
  else
    AssertThrow(false, ExcNotImplemented());

  if (prm.get("geometry") == "hypercube")
    this->geometry = hypercube;
  else if (prm.get("geometry") == "hyperrectangle")
    this->geometry = hyperrectangle;
  else if (prm.get("geometry") == "cocircle")
    this->geometry = cocircle;
  else if (prm.get("geometry") == "corectangle")
    this->geometry = corectangle;
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


using namespace dealii;

template <int dim,
          typename number,
          typename VectorizedArrayType = VectorizedArray<number>>
class PoissonOperator
{
public:
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, 1, number, VectorizedArrayType>;
  using FEFaceIntegrator =
    FEFaceEvaluation<dim, -1, 0, 1, number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<number>;

  PoissonOperator(
    const MatrixFree<dim, double, VectorizedArrayType> &matrix_free,
    const std::vector<unsigned int>                    &non_matching_faces)
    : matrix_free(matrix_free)
    , panalty_factor(
        compute_pentaly_factor(matrix_free.get_dof_handler().get_fe().degree,
                               1.0))
  {
    // store all boundary faces in one set
    for (const auto &face_pair : non_matching_faces)
      faces.insert(face_pair);


    std::vector<
      std::pair<types::boundary_id, std::function<std::vector<bool>()>>>
      non_matching_faces_marked_vertices;

    for (const auto &nm_face : non_matching_faces)
      {
        auto marked_vertices = [&]() {
          const auto &tria = matrix_free.get_dof_handler().get_triangulation();

          std::vector<bool> mask(tria.n_vertices(), true);

          for (const auto &cell : tria.active_cell_iterators())
            for (auto const &f : cell->face_indices())
              if (cell->face(f)->at_boundary() &&
                  cell->face(f)->boundary_id() == nm_face)
                for (const auto v : cell->vertex_indices())
                  mask[cell->vertex_index(v)] = false;

          return mask;
        };

        non_matching_faces_marked_vertices.emplace_back(
          std::make_pair(nm_face, marked_vertices));
      }

    phi_r_comm =
      Utilities::compute_remote_communicator_faces_point_to_point_interpolation(
        matrix_free, non_matching_faces_marked_vertices);

    phi_r_cache =
      std::make_shared<FERemoteEvaluation<dim, 1, VectorizedArrayType>>(
        phi_r_comm, matrix_free.get_dof_handler());

    phi_r_sigma_cache =
      std::make_shared<FERemoteEvaluation<dim, 1, VectorizedArrayType>>(
        phi_r_comm, matrix_free.get_dof_handler().get_triangulation());

    compute_penalty_parameters();
  }

  void
  initialize_dof_vector(VectorType &vec)
  {
    matrix_free.initialize_dof_vector(vec);
  }

  void
  rhs(VectorType &vec) const
  {
    const int dummy = 0;

    matrix_free.template cell_loop<VectorType, int>(
      [&](const auto &data, auto &dst, const auto &, const auto cells) {
        FECellIntegrator phi(data);
        for (unsigned int cell = cells.first; cell < cells.second; ++cell)
          {
            phi.reinit(cell);
            for (unsigned int q = 0; q < phi.n_q_points; ++q)
              phi.submit_value(1.0, q);

            phi.integrate_scatter(EvaluationFlags::values, dst);
          }
      },
      vec,
      dummy,
      true);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    phi_r_cache->gather_evaluate(src,
                                 EvaluationFlags::values |
                                   EvaluationFlags::gradients);

    matrix_free.loop(
      &PoissonOperator<dim, number, VectorizedArrayType>::do_vmult_cell,
      &PoissonOperator<dim, number, VectorizedArrayType>::do_vmult_face,
      &PoissonOperator<dim, number, VectorizedArrayType>::do_vmult_boundary,
      this,
      dst,
      src,
      true);
  }

  void
  do_vmult_cell(const MatrixFree<dim, number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src, EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          phi.submit_gradient(phi.get_gradient(q), q);
        phi.integrate_scatter(EvaluationFlags::gradients, dst);
      }
  }

  void
  do_vmult_face(const MatrixFree<dim, number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &face_range) const
  {
    (void)data;
    (void)dst;
    (void)src;
    (void)face_range;
  }

  void
  do_vmult_boundary(
    const MatrixFree<dim, number>               &data,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &face_range) const
  {
    FEFaceIntegrator phi_m(data, true);

    auto phi_r       = phi_r_cache->get_data_accessor();
    auto phi_r_sigma = phi_r_sigma_cache->get_data_accessor();

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        if (is_internal_face(face) == false)
          continue; // nothing to do

        phi_m.reinit(face);

        phi_m.gather_evaluate(src,
                              EvaluationFlags::values |
                                EvaluationFlags::gradients);

        phi_r.reinit(face);
        phi_r_sigma.reinit(face);

        const auto sigma_m = phi_m.read_cell_data(array_penalty_parameter);

        for (unsigned int q = 0; q < phi_m.n_q_points; ++q)
          {
            const auto value_m = phi_m.get_value(q);
            const auto value_p = phi_r.get_value(q);

            const auto gradient_m = phi_m.get_gradient(q);
            const auto gradient_p = phi_r.get_gradient(q);

            const auto sigma_p = phi_r_sigma.get_value(q);
            const auto sigma   = std::max(sigma_m, sigma_p) * panalty_factor;

            const auto jump_value = (value_m - value_p) * 0.5;
            const auto avg_gradient =
              phi_m.get_normal_vector(q) * (gradient_m + gradient_p) * 0.5;

            phi_m.submit_normal_derivative(-jump_value, q);
            phi_m.submit_value(jump_value * sigma * 2.0 - avg_gradient, q);
          }
        phi_m.integrate_scatter(EvaluationFlags::values |
                                  EvaluationFlags::gradients,
                                dst);
      }
  }

  void
  do_residual_cell(
    const MatrixFree<dim, number>               &data,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.gather_evaluate(src,
                            EvaluationFlags::values |
                              EvaluationFlags::gradients);
        for (unsigned int q = 0; q < phi.n_q_points; ++q)
          {
            phi.submit_gradient(phi.get_gradient(q), q);
            phi.submit_value(-1.0, q);
          }
        phi.integrate_scatter(EvaluationFlags::values |
                                EvaluationFlags::gradients,
                              dst);
      }
  }

  void
  evaluate_residual(VectorType &dst, const VectorType &src) const
  {
    phi_r_cache->gather_evaluate(src,
                                 EvaluationFlags::values |
                                   EvaluationFlags::gradients);

    matrix_free.loop(
      &PoissonOperator<dim, number, VectorizedArrayType>::do_residual_cell,
      &PoissonOperator<dim, number, VectorizedArrayType>::do_vmult_face,
      &PoissonOperator<dim, number, VectorizedArrayType>::do_vmult_boundary,
      this,
      dst,
      src,
      true);
  }

private:
  bool
  is_internal_face(const unsigned int face) const
  {
    return faces.find(matrix_free.get_boundary_id(face)) != faces.end();
  }

  void
  compute_penalty_parameters()
  {
    // step 1) compute penalty parameter of each cell
    const unsigned int n_cells =
      matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
    array_penalty_parameter.resize(n_cells);

    dealii::Mapping<dim> const &mapping =
      *matrix_free.get_mapping_info().mapping;
    dealii::FiniteElement<dim> const &fe =
      matrix_free.get_dof_handler().get_fe();
    const unsigned int degree = fe.degree;

    dealii::QGauss<dim>   quadrature(degree + 1);
    dealii::FEValues<dim> fe_values(mapping,
                                    fe,
                                    quadrature,
                                    dealii::update_JxW_values);

    dealii::QGauss<dim - 1>   face_quadrature(degree + 1);
    dealii::FEFaceValues<dim> fe_face_values(mapping,
                                             fe,
                                             face_quadrature,
                                             dealii::update_JxW_values);

    Vector<number> array_penalty_parameter_scalar(
      matrix_free.get_dof_handler().get_triangulation().n_active_cells());

    for (unsigned int i = 0; i < n_cells; ++i)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(i);
           ++v)
        {
          typename dealii::DoFHandler<dim>::cell_iterator cell =
            matrix_free.get_cell_iterator(i, v);
          fe_values.reinit(cell);

          number volume = 0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            volume += fe_values.JxW(q);

          number surface_area = 0;
          for (const auto f : cell->face_indices())
            {
              fe_face_values.reinit(cell, f);
              const number factor =
                (cell->at_boundary(f) && !cell->has_periodic_neighbor(f)) ? 1. :
                                                                            0.5;
              for (unsigned int q = 0; q < face_quadrature.size(); ++q)
                surface_area += fe_face_values.JxW(q) * factor;
            }

          array_penalty_parameter[i][v] = surface_area / volume;
          array_penalty_parameter_scalar[cell->active_cell_index()] =
            surface_area / volume;
        }


    dealii::FEEvaluation<dim, -1, 0, 1, number> fe_eval(matrix_free, 0, 0);

    phi_r_sigma_cache->gather_evaluate(array_penalty_parameter_scalar,
                                       EvaluationFlags::values);
  }

  static number
  compute_pentaly_factor(const unsigned int degree, const number factor)
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  const MatrixFree<dim, number, VectorizedArrayType> &matrix_free;

  FERemoteEvaluationCommunicator<dim> phi_r_comm;
  mutable std::shared_ptr<FERemoteEvaluation<dim, 1, VectorizedArrayType>>
    phi_r_cache;
  mutable std::shared_ptr<FERemoteEvaluation<dim, 1, VectorizedArrayType>>
    phi_r_sigma_cache;

  const double                               panalty_factor;
  dealii::AlignedVector<VectorizedArrayType> array_penalty_parameter;

  std::set<unsigned int> faces;
};


// Main class to solve the Mortar Non-Linear Poisson problem
template <int dim>
class MatrixFreeMortarNonLinearPoisson
{
public:
  MatrixFreeMortarNonLinearPoisson(const Settings &parameters)
    : parameters(parameters)
    , fe_q(parameters.element_order)
    , quad(parameters.element_order + 1)
    , tria(MPI_COMM_WORLD)
    , dof_handler(tria){};

  void
  solve();

private:
  void
  make_grid();
  //
  // void
  // refine_grid();
  //
  void
  setup_system();

  // void
  // evaluate_residual(
  //   LinearAlgebra::distributed::Vector<double> &      dst,
  //   const LinearAlgebra::distributed::Vector<double> &src) const;

  // void
  // local_evaluate_residual(
  //   const MatrixFree<dim, double> &                   data,
  //   LinearAlgebra::distributed::Vector<double> &      dst,
  //   const LinearAlgebra::distributed::Vector<double> &src,
  //   const std::pair<unsigned int, unsigned int> &     cell_range) const;
  //
  // void
  // assemble_rhs();
  //
  // double
  // compute_residual(const double alpha);
  //
  // void
  // compute_update();

  // void
  // compute_solution_norm() const;
  //
  // void
  // compute_l2_error() const;
  //
  // void
  // output_results(const unsigned int cycle) const;
private:
  Settings parameters;

  const MappingQ1<dim>                      mapping;
  const FE_Q<dim>                           fe_q;
  const QGauss<dim>                         quad;
  parallel::distributed::Triangulation<dim> tria;
  DoFHandler<dim>                           dof_handler;


  std::vector<std::pair<unsigned int, unsigned int>> face_pairs;
  std::vector<unsigned int>                          nm_face_pairs;
};

template <int dim>
void
MatrixFreeMortarNonLinearPoisson<dim>::make_grid()
{
  if (parameters.geometry == Settings::GeometryType::corectangle)
    {
      Triangulation<dim> tria_0, tria_1;
      GridGenerator::subdivided_hyper_rectangle(
        tria_0, {7, 7}, {0.0, 0.0}, {1.0, 1.0}, true);

      GridGenerator::subdivided_hyper_rectangle(
        tria_1, {6, 3}, {1.0, 0.0}, {3.0, 1.0}, true);

      for (const auto &face : tria_1.active_face_iterators())
        if (face->at_boundary())
          face->set_boundary_id(face->boundary_id() + 2 * dim);

      GridGenerator::merge_triangulations(
        tria_0, tria_1, tria, 0., false, true);
      AssertDimension(tria_0.n_vertices() + tria_1.n_vertices(),
                      tria.n_vertices());

      tria.refine_global(parameters.initial_refinement);

      face_pairs.emplace_back(1, 2 * dim);
      face_pairs.emplace_back(2 * dim, 1);

      nm_face_pairs.emplace_back(1);
      nm_face_pairs.emplace_back(2 * dim);
    }
  else if (parameters.geometry == Settings::GeometryType::cocircle)
    {
      // Generate two grids of two non-matching circles
      const double r1_i = 0.25;
      const double r1_o = 0.5;
      const double r2_i = 0.5;
      const double r2_o = 1;


      Triangulation<dim> circle_one;
      GridGenerator::hyper_shell(circle_one,
                                 (dim == 2) ? Point<dim>(0, 0) :
                                              Point<dim>(0, 0, 0),
                                 r1_i,
                                 r1_o,
                                 0,
                                 true);
      Triangulation<dim> circle_two;
      GridGenerator::hyper_shell(circle_two,
                                 (dim == 2) ? Point<dim>(0, 0) :
                                              Point<dim>(0, 0, 0),
                                 r2_i,
                                 r2_o,
                                 6,
                                 true);

      // shift boundary id of circle two
      for (const auto &face : circle_two.active_face_iterators())
        if (face->at_boundary())
          face->set_boundary_id(face->boundary_id() + 2);

      GridGenerator::merge_triangulations(
        circle_one, circle_two, tria, 0, true, true);
      tria.set_manifold(0,
                        SphericalManifold<dim>(
                          (dim == 2) ? Point<dim>(0, 0) : Point<dim>(0, 0, 0)));

      tria.refine_global(parameters.initial_refinement);

      face_pairs.emplace_back(1, 2);
      face_pairs.emplace_back(2, 1);

      nm_face_pairs.emplace_back(1);
      nm_face_pairs.emplace_back(2);
    }
  else
    {
      throw(std::runtime_error("Not implemented"));
    }
}



template <int dim>
void
MatrixFreeMortarNonLinearPoisson<dim>::solve()
{
  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int fe_degree            = parameters.element_order;
  const unsigned int n_global_refinements = parameters.initial_refinement;

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // Make triangulation
  make_grid();

  // create DoFHandler
  dof_handler.distribute_dofs(fe_q);

  // create MatrixFree
  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
  data.mapping_update_flags =
    update_quadrature_points | update_gradients | update_values;
  data.mapping_update_flags_boundary_faces = data.mapping_update_flags;
  data.mapping_update_flags_inner_faces    = data.mapping_update_flags;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  AffineConstraints<Number>                    constraints;

  for (unsigned int d = 0; d < 4 * dim; ++d)
    if (face_pairs.size() == 0)
      {
        DoFTools::make_zero_boundary_constraints(dof_handler, d, constraints);
      }
    else
      {
        for (const auto &face_pair : face_pairs)
          if (d != face_pair.first && d != face_pair.second)
            DoFTools::make_zero_boundary_constraints(dof_handler,
                                                     d,
                                                     constraints);
      }

  constraints.close();

  matrix_free.reinit(mapping, dof_handler, constraints, quad, data);

  pcout << "Statistics:" << std::endl;
  pcout << " - n cells: " << tria.n_global_active_cells() << std::endl;
  pcout << " - n dof:   " << dof_handler.n_dofs() << std::endl;
  pcout << std::endl;

  std::ofstream mesh_file("mesh.vtu");
  GridOut().write_vtu(tria, mesh_file);

  PoissonOperator<dim, Number, VectorizedArrayType> op(matrix_free,
                                                       nm_face_pairs);

  VectorType rhs, solution, newton_update, residual;

  op.initialize_dof_vector(rhs);
  op.initialize_dof_vector(solution);
  op.initialize_dof_vector(newton_update);
  op.initialize_dof_vector(residual);


  op.rhs(rhs);
  rhs.zero_out_ghost_values();

  constraints.set_zero(rhs);

  const unsigned int itmax = 10;
  const double       TOLf  = 1e-12;
  const double       TOLx  = 1e-10;

  Timer solver_timer;
  solver_timer.start();

  op.evaluate_residual(residual,solution);
  pcout << "Initial norm of the residual is: " << residual.l2_norm() << std::endl;

  for (unsigned int newton_step = 1; newton_step <= itmax; ++newton_step)
    {
      ReductionControl reduction_control(10000, 1e-20, 1e-12);

      // note: we need to use GMRES, since the system is non-symmetrical
      SolverGMRES<VectorType> solver(reduction_control);
      op.evaluate_residual(residual, solution);
      // Multiply by -1 to have J(x) dx = - R(x)
      residual*=-1;

      solver.solve(op, newton_update, residual, PreconditionIdentity());

      unsigned int  linear_iterations = reduction_control.last_step();

      const double ERRx = newton_update.l2_norm();
      solution.add(1.0, newton_update);
      op.evaluate_residual(residual,solution);
      const double ERRf = residual.l2_norm();

       pcout << "   Nstep " << newton_step << ", errf = " << ERRf
            << ", errx = " << ERRx << ", it = " << linear_iterations
            << std::endl;

       if (ERRf < TOLf || ERRx < TOLx)
        {
          pcout << "Convergence step " << newton_step << " value " << ERRf
                << std::endl;
          break;
        }
    }

  // output result
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution_poisson");
  data_out.add_data_vector(residual, "residual");
  Vector<float> ranks(tria.n_active_cells());
  ranks = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  data_out.add_data_vector(ranks, "ranks");
  data_out.build_patches();
  data_out.write_vtu_in_parallel("solution_poisson.vtu", MPI_COMM_WORLD);
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

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
              MatrixFreeMortarNonLinearPoisson<2> problem(parameters);
              problem.solve();

              break;
            }

          case 3:
            {
              MatrixFreeMortarNonLinearPoisson<3> problem(parameters);
              problem.solve();
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
