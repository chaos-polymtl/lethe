// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_tools.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

using namespace dealii;

template <int dim, int spacedim>
std::tuple<std::vector<Point<spacedim>>, std::vector<double>>
compute_quadrature(double             rad_00,
                   double             rad_01,
                   double             rad_10,
                   double             rad_11,
                   const double       radius,
                   const unsigned int n_quadrature_points,
                   const unsigned int mapping_degree = 10)
{
  if (rad_00 > rad_10) // normalize
    {
      std::swap(rad_00, rad_10);
      std::swap(rad_01, rad_11);
    }

  if (rad_01 <= rad_10)
    return {}; // no cut

  const double left_rad  = std::max(rad_00, rad_10);
  const double right_rad = std::min(rad_01, rad_11);

  std::vector<Point<spacedim>> points;
  points.emplace_back(std::cos(left_rad) * radius, std::sin(left_rad) * radius);
  points.emplace_back(std::cos(right_rad) * radius,
                      std::sin(right_rad) * radius);

  std::vector<CellData<dim>> cells;
  CellData<dim>              cell(2);
  cell.vertices[0] = 0;
  cell.vertices[1] = 1;
  cells.emplace_back(cell);

  Triangulation<dim, spacedim> tria;
  tria.create_triangulation(points, cells, {});

  tria.set_manifold(0, SphericalManifold<dim, spacedim>());
  tria.set_all_manifold_ids(0);

  MappingQ<dim, spacedim>   mapping(mapping_degree);
  QGauss<dim>               quadrature(n_quadrature_points);
  FE_Nothing<dim, spacedim> fe;

  FEValues<dim, spacedim> phi(mapping,
                              fe,
                              quadrature,
                              update_quadrature_points | update_JxW_values);
  phi.reinit(tria.begin());

  return {phi.get_quadrature_points(), phi.get_JxW_values()};
}

template <int structdim, int dim, int spacedim>
std::tuple<std::vector<Point<spacedim>>, std::vector<double>>
compute_quadrature(
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &cell_0,
  const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &cell_1,
  const double                                                radius,
  const unsigned int n_quadrature_points,
  const unsigned int mapping_degree = 10)
{
  AssertDimension(structdim, 1);
  AssertDimension(spacedim, 2);

  // helper function to 1) guarantee that all segments have the same orientation
  // and 2) periodicities are handled correctly -> split up into two segment
  const auto create_sections = [](const auto &cell) {
    const auto point_to_rad = [](const auto point) {
      const auto temp = std::atan(std::abs(point[1]) / std::abs(point[0]));

      if (point[1] >= 0.0)
        {
          if (point[0] >= 0.0)
            return temp;
          else
            return numbers::PI - temp;
        }
      else
        {
          if (point[0] >= 0.0)
            return 2 * numbers::PI - temp;
          else
            return numbers::PI + temp;
        }
    };

    const auto cross_product = [](const auto &p0, const auto &p1) {
      Tensor<1, 3, double> t0;
      Tensor<1, 3, double> t1;

      for (unsigned int i = 0; i < 2; ++i)
        {
          t0[i] = p0[i];
          t1[i] = p1[i];
        }

      return cross_product_3d(t0, t1)[2];
    };

    double rad_0 = point_to_rad(cell->vertex(0));
    double rad_1 = point_to_rad(cell->vertex(1));

    if (cross_product(cell->vertex(0), cell->vertex(1)) < 0.0)
      std::swap(rad_0, rad_1);

    std::vector<std::pair<double, double>> sections;

    if (rad_0 < rad_1) // normal
      {
        sections.emplace_back(rad_0, rad_1);
      }
    else // periodic
      {
        sections.emplace_back(rad_0, 2 * numbers::PI);
        sections.emplace_back(0, rad_1);
      }

    return sections;
  };

  const auto subsections_0 = create_sections(cell_0);
  const auto subsections_1 = create_sections(cell_1);

  std::vector<Point<spacedim>> all_points;
  std::vector<double>          all_JxWs;

  for (const auto &subsection_0 : subsections_0)
    for (const auto &subsection_1 : subsections_1)
      {
        const auto [points, JxWs] =
          compute_quadrature<structdim, spacedim>(subsection_0.first,
                                                  subsection_0.second,
                                                  subsection_1.first,
                                                  subsection_1.second,
                                                  radius,
                                                  n_quadrature_points,
                                                  mapping_degree);

        all_points.insert(all_points.end(), points.begin(), points.end());
        all_JxWs.insert(all_JxWs.end(), JxWs.begin(), JxWs.end());
      }

  return {all_points, all_JxWs};
}

template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class PoissonOperator : public Subscriptor
{
public:
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  PoissonOperator(const MappingQ<dim>             &mapping,
                  const DoFHandler<dim>           &dof_handler,
                  const AffineConstraints<Number> &constraints,
                  const Quadrature<dim>           &quadrature)
    : mapping(mapping)
    , dof_handler(dof_handler)
    , constraints(constraints)
    , quadrature(quadrature)
    , panalty_factor(compute_pentaly_factor(dof_handler.get_fe().degree, 1.0))
    , valid_system(false)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;
    data.mapping_update_flags_boundary_faces = data.mapping_update_flags;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    const auto &tria = dof_handler.get_triangulation();

    const double radius = 0.5; // TODO

    const auto compute_penalty_parameter = [&](const auto &cell) {
      const dealii::FiniteElement<dim> &fe =
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

      fe_values.reinit(cell);

      Number volume = 0;
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        volume += fe_values.JxW(q);

      Number surface_area = 0;
      for (const auto f : cell->face_indices())
        {
          fe_face_values.reinit(cell, f);
          const Number factor =
            (cell->at_boundary(f) && !cell->has_periodic_neighbor(f)) ? 1. :
                                                                        0.5;
          for (unsigned int q = 0; q < face_quadrature.size(); ++q)
            surface_area += fe_face_values.JxW(q) * factor;
        }

      return surface_area / volume;
    };

    for (const auto &cell_0 : tria.active_cell_iterators())
      for (const auto &face_0 : cell_0->face_iterators())
        if (face_0->boundary_id() == 1)
          for (const auto &cell_1 : tria.active_cell_iterators())
            for (const auto &face_1 : cell_1->face_iterators())
              if (face_1->boundary_id() == 2)
                {
                  const auto [points, JxWs] =
                    compute_quadrature<dim - 1, dim, dim>(face_0,
                                                          face_1,
                                                          radius,
                                                          quadrature.size());

                  std::vector<Point<dim>> points_0(points.size());
                  mapping.transform_points_real_to_unit_cell(cell_0,
                                                             points,
                                                             points_0);

                  std::vector<Point<dim>> points_1(points.size());
                  mapping.transform_points_real_to_unit_cell(cell_1,
                                                             points,
                                                             points_1);

                  std::vector<Tensor<1, dim, Number>> normal;

                  for (const auto &p : points)
                    normal.emplace_back(p / p.norm());

                  const Number penalty_parameter =
                    std::min(compute_penalty_parameter(cell_0),
                             compute_penalty_parameter(cell_1));


                  all_intersections.emplace_back(JxWs,
                                                 cell_0,
                                                 points_0,
                                                 cell_1,
                                                 points_1,
                                                 normal,
                                                 penalty_parameter);
                }
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
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
    matrix_free.cell_loop(
      &PoissonOperator<dim, Number, VectorizedArrayType>::do_vmult_cell,
      this,
      dst,
      src,
      true);

    FEPointEvaluation<1, dim, dim, double> phi_m(
      mapping, dof_handler.get_fe(), update_values | update_gradients);
    FEPointEvaluation<1, dim, dim, double> phi_p(
      mapping, dof_handler.get_fe(), update_values | update_gradients);

    for (const auto &[JxWs,
                      cell_0,
                      points_0,
                      cell_1,
                      points_1,
                      normals,
                      penalty_parameter] : all_intersections)
      {
        const auto cell_dof_0 = cell_0->as_dof_handler_iterator(dof_handler);
        const auto cell_dof_1 = cell_1->as_dof_handler_iterator(dof_handler);

        phi_m.reinit(cell_0, points_0);
        phi_p.reinit(cell_1, points_1);

        // gather
        Vector<double> buffer_0;
        buffer_0.reinit(cell_dof_0->get_fe().n_dofs_per_cell());
        cell_dof_0->get_dof_values(src, buffer_0);

        Vector<double> buffer_1;
        buffer_1.reinit(cell_dof_1->get_fe().n_dofs_per_cell());
        cell_dof_1->get_dof_values(src, buffer_1);

        // evaluate
        phi_m.evaluate(buffer_0, EvaluationFlags::values);
        phi_p.evaluate(buffer_1, EvaluationFlags::values);

        // quadrature loop
        for (const auto q : phi_m.quadrature_point_indices())
          {
            const auto normal = normals[q];
            const auto JxW    = JxWs[q];

            const auto value_m = phi_m.get_value(q);
            const auto value_p = phi_p.get_value(q);

            const auto gradient_m = phi_m.get_gradient(q);
            const auto gradient_p = phi_p.get_gradient(q);

            const auto jump_value = (value_m - value_p) * 0.5 * JxW;
            const auto avg_gradient =
              normal * (gradient_m + gradient_p) * 0.5 * JxW;

            const double sigma = penalty_parameter * panalty_factor;

            phi_m.submit_gradient(-jump_value * normal, q);
            phi_m.submit_value(+jump_value * sigma * 2.0 - avg_gradient, q);

            phi_p.submit_gradient(-jump_value * normal, q);
            phi_p.submit_value(-jump_value * sigma * 2.0 + avg_gradient, q);
          }

        // integrate
        phi_m.test_and_sum(buffer_0, EvaluationFlags::values);
        phi_p.test_and_sum(buffer_1, EvaluationFlags::values);

        // scatter
        cell_dof_0->distribute_local_to_global(buffer_0, dst);
        cell_dof_1->distribute_local_to_global(buffer_1, dst);
      }
  }

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const
  {
    initialize_system_matrix();

    return system_matrix;
  }

  void
  initialize_system_matrix() const
  {
    const auto &dof_handler = matrix_free.get_dof_handler();
    const auto &constraints = matrix_free.get_affine_constraints();

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(this->matrix_free.get_mg_level() !=
                       numbers::invalid_unsigned_int ?
                     dof_handler.locally_owned_mg_dofs(
                       this->matrix_free.get_mg_level()) :
                     dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
          MGTools::make_sparsity_pattern(dof_handler,
                                         dsp,
                                         this->matrix_free.get_mg_level(),
                                         constraints);
        else
          DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);

        // for coupling term
        for (const auto &[JxWs,
                          cell_0,
                          points_0,
                          cell_1,
                          points_1,
                          normals,
                          penalty_parameter] : all_intersections)
          {
            const auto cell_dof_0 =
              cell_0->as_dof_handler_iterator(dof_handler);
            const auto cell_dof_1 =
              cell_1->as_dof_handler_iterator(dof_handler);

            std::vector<types::global_dof_index> dof_indices_0(
              cell_dof_0->get_fe().n_dofs_per_cell());
            std::vector<types::global_dof_index> dof_indices_1(
              cell_dof_1->get_fe().n_dofs_per_cell());

            cell_dof_0->get_dof_indices(dof_indices_0);
            cell_dof_1->get_dof_indices(dof_indices_1);

            constraints.add_entries_local_to_global(dof_indices_0,
                                                    dof_indices_1,
                                                    dsp);
            constraints.add_entries_local_to_global(dof_indices_1,
                                                    dof_indices_0,
                                                    dsp);
          }

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::compute_matrix(
          matrix_free,
          constraints,
          system_matrix,
          &PoissonOperator<dim, Number, VectorizedArrayType>::
            do_vmult_cell_single,
          this);

        // for coupling term
        FEPointEvaluation<1, dim, dim, double> phi_m(
          mapping, dof_handler.get_fe(), update_values | update_gradients);
        FEPointEvaluation<1, dim, dim, double> phi_p(
          mapping, dof_handler.get_fe(), update_values | update_gradients);

        for (const auto &[JxWs,
                          cell_0,
                          points_0,
                          cell_1,
                          points_1,
                          normals,
                          penalty_parameter] : all_intersections)
          {
            const auto cell_dof_0 =
              cell_0->as_dof_handler_iterator(dof_handler);
            const auto cell_dof_1 =
              cell_1->as_dof_handler_iterator(dof_handler);

            phi_m.reinit(cell_0, points_0);
            phi_p.reinit(cell_1, points_1);

            const unsigned int n_dofs_per_cell_0 =
              cell_dof_0->get_fe().n_dofs_per_cell();
            const unsigned int n_dofs_per_cell_1 =
              cell_dof_1->get_fe().n_dofs_per_cell();

            Vector<double> buffer_0;
            buffer_0.reinit(n_dofs_per_cell_0);

            Vector<double> buffer_1;
            buffer_1.reinit(n_dofs_per_cell_1);

            std::vector<types::global_dof_index> dof_indices_0(
              cell_dof_0->get_fe().n_dofs_per_cell());
            std::vector<types::global_dof_index> dof_indices_1(
              cell_dof_0->get_fe().n_dofs_per_cell());

            cell_dof_0->get_dof_indices(dof_indices_0);
            cell_dof_1->get_dof_indices(dof_indices_1);

            // evaluate
            phi_m.evaluate(buffer_0, EvaluationFlags::values);
            phi_p.evaluate(buffer_1, EvaluationFlags::values);

            for (unsigned int b = 0; b < 2; ++b)
              {
                FullMatrix<double> matrix_0;
                FullMatrix<double> matrix_1;

                if (b == 0)
                  {
                    matrix_0.reinit(n_dofs_per_cell_0, n_dofs_per_cell_0);
                    matrix_1.reinit(n_dofs_per_cell_1, n_dofs_per_cell_0);
                  }
                else
                  {
                    matrix_0.reinit(n_dofs_per_cell_0, n_dofs_per_cell_1);
                    matrix_1.reinit(n_dofs_per_cell_1, n_dofs_per_cell_1);
                  }

                for (unsigned int j = 0;
                     j < (b == 0 ? n_dofs_per_cell_0 : n_dofs_per_cell_1);
                     ++j)
                  {
                    buffer_0 = 0.0;
                    buffer_1 = 0.0;

                    if (b == 0)
                      buffer_0[j] = 1.0;
                    else
                      buffer_1[j] = 1.0;

                    phi_m.evaluate(buffer_0, EvaluationFlags::values);
                    phi_p.evaluate(buffer_1, EvaluationFlags::values);

                    // quadrature loop
                    for (const auto q : phi_m.quadrature_point_indices())
                      {
                        const auto normal = normals[q];
                        const auto JxW    = JxWs[q];

                        const auto value_m = phi_m.get_value(q);
                        const auto value_p = phi_p.get_value(q);

                        const auto gradient_m = phi_m.get_gradient(q);
                        const auto gradient_p = phi_p.get_gradient(q);

                        const auto jump_value = (value_m - value_p) * 0.5 * JxW;
                        const auto avg_gradient =
                          normal * (gradient_m + gradient_p) * 0.5 * JxW;

                        const double sigma = penalty_parameter * panalty_factor;

                        phi_m.submit_gradient(-jump_value * normal, q);
                        phi_m.submit_value(+jump_value * sigma * 2.0 -
                                             avg_gradient,
                                           q);

                        phi_p.submit_gradient(-jump_value * normal, q);
                        phi_p.submit_value(-jump_value * sigma * 2.0 +
                                             avg_gradient,
                                           q);
                      }

                    // integrate
                    phi_m.test_and_sum(buffer_0, EvaluationFlags::values);
                    phi_p.test_and_sum(buffer_1, EvaluationFlags::values);

                    for (unsigned int i = 0; i < n_dofs_per_cell_0; ++i)
                      matrix_0[i][j] = buffer_0[i];
                    for (unsigned int i = 0; i < n_dofs_per_cell_1; ++i)
                      matrix_1[i][j] = buffer_1[i];
                  }


                if (b == 0)
                  {
                    constraints.distribute_local_to_global(matrix_0,
                                                           dof_indices_0,
                                                           system_matrix);
                    constraints.distribute_local_to_global(matrix_1,
                                                           dof_indices_1,
                                                           dof_indices_0,
                                                           system_matrix);
                  }
                else
                  {
                    constraints.distribute_local_to_global(matrix_0,
                                                           dof_indices_0,
                                                           dof_indices_1,
                                                           system_matrix);
                    constraints.distribute_local_to_global(matrix_1,
                                                           dof_indices_1,
                                                           system_matrix);
                  }
              }
          }

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }

private:
  const MappingQ<dim>             &mapping;
  const DoFHandler<dim>           &dof_handler;
  const AffineConstraints<Number> &constraints;
  const Quadrature<dim>           &quadrature;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;

  std::vector<std::tuple<std::vector<double>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         std::vector<Tensor<1, dim, Number>>,
                         Number>>
    all_intersections;

  void
  do_vmult_cell(const MatrixFree<dim, Number>               &data,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FECellIntegrator phi(data);

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        phi.reinit(cell);
        phi.read_dof_values(src);

        do_vmult_cell_single(phi);
        phi.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(FECellIntegrator &phi) const
  {
    phi.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);

    phi.integrate(EvaluationFlags::gradients);
  }

  static Number
  compute_pentaly_factor(const unsigned int degree, const Number factor)
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  mutable TrilinosWrappers::SparseMatrix system_matrix;
  const Number                           panalty_factor;
  mutable bool                           valid_system;
};

template <int dim>
void
create_triangulation(Triangulation<dim> &tria)
{
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
                    SphericalManifold<dim>((dim == 2) ? Point<dim>(0, 0) :
                                                        Point<dim>(0, 0, 0)));
}

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

  FE_Q<dim>     fe(fe_degree);
  MappingQ<dim> mapping(mapping_degree);
  QGauss<dim>   quadrature(fe_degree + 1);

  Triangulation<dim> tria;
  create_triangulation(tria);
  tria.refine_global(n_global_refinements);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  DoFTools::make_zero_boundary_constraints(dof_handler, 0, constraints);
  DoFTools::make_zero_boundary_constraints(dof_handler, 3, constraints);
  constraints.close();

  PoissonOperator<dim, double> op(mapping,
                                  dof_handler,
                                  constraints,
                                  quadrature);

  LinearAlgebra::distributed::Vector<double> rhs, solution;
  op.initialize_dof_vector(rhs);
  op.initialize_dof_vector(solution);

  op.rhs(rhs);

  ReductionControl reduction_control(10000, 1e-20, 1e-4);

  TrilinosWrappers::PreconditionAMG
    preconditioner; // TrilinosWrappers::SolverDirect
  preconditioner.initialize(op.get_system_matrix());

  SolverGMRES<VectorType> solver(reduction_control);
  solver.solve(op, solution, rhs, preconditioner);

  std::cout << reduction_control.last_step() << std::endl;

  DataOut<dim> data_out;

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(mapping,
                         fe_degree + 1,
                         DataOut<dim>::CurvedCellRegion::curved_inner_cells);
  data_out.write_vtu_in_parallel("poisson_dg.vtu", MPI_COMM_WORLD);
}
