// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test CouplingOperatorBase on rotated hypercube-cylinder geometry
 * by implementing SIPG and comparing it to a regular DG simulation.
 */

#include <deal.II/base/conditional_ostream.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
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

#include "./tests.h"

using namespace dealii;


template <int dim, typename Number>
class MyCouplingOperator : public CouplingOperatorBase<dim, Number>
{
public:
  using FEPointIntegrator = FEPointEvaluation<1, dim, dim, Number>;
  using value_type        = typename FEPointIntegrator::value_type;

  MyCouplingOperator(const Mapping<dim>              &mapping,
                     const DoFHandler<dim>           &dof_handler,
                     const AffineConstraints<Number> &constraints,
                     const Quadrature<dim>            quadrature,
                     const unsigned int               n_subdivisions,
                     const double                     radius,
                     const double                     rotate_pi,
                     const unsigned int               bid_0,
                     const unsigned int               bid_1,
                     const double                     sip_factor = 1.0)
    : CouplingOperatorBase<dim, Number>(mapping,
                                        dof_handler,
                                        constraints,
                                        quadrature,
                                        n_subdivisions,
                                        1 /*n components*/,
                                        2 /*n data points*/,
                                        radius,
                                        rotate_pi,
                                        bid_0,
                                        bid_1,
                                        sip_factor,
                                        get_relevant_dof_indices(
                                          dof_handler.get_fe()),
                                        0.0 /*TODO*/)
    , fe_sub(dof_handler.get_fe().base_element(
               dof_handler.get_fe().component_to_base_index(0).first),
             1)
    , phi_m(mapping, fe_sub, update_values | update_gradients)
  {}

  static std::vector<unsigned int>
  get_relevant_dof_indices(const FiniteElement<dim> &fe)
  {
    std::vector<unsigned int> result;

    for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
      result.push_back(i);

    AssertDimension(fe.n_dofs_per_cell(), result.size());

    return result;
  }

  void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const override
  {
    this->phi_m.reinit(cell, points);
  }

  void
  local_evaluate(const Vector<Number> &buffer,
                 const unsigned int    ptr_q,
                 const unsigned int    q_stride,
                 Number               *all_value_m) const override
  {
    this->phi_m.evaluate(buffer,
                         EvaluationFlags::values | EvaluationFlags::gradients);

    for (const auto q : this->phi_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        const auto normal     = this->all_normals[q_index];
        const auto value_m    = this->phi_m.get_value(q);
        const auto gradient_m = contract(this->phi_m.get_gradient(q), normal);

        BufferRW<Number> buffer_m(all_value_m, q * 2 * q_stride);

        buffer_m.write(value_m);
        buffer_m.write(gradient_m);
      }
  }

  void
  local_integrate(Vector<Number>    &buffer,
                  const unsigned int ptr_q,
                  const unsigned int q_stride,
                  Number            *all_value_m,
                  Number            *all_value_p) const override
  {
    for (const auto q : this->phi_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        BufferRW<Number> buffer_m(all_value_m, q * 2 * q_stride);
        BufferRW<Number> buffer_p(all_value_p, q * 2 * q_stride);

        const auto value_m           = buffer_m.template read<value_type>();
        const auto value_p           = buffer_p.template read<value_type>();
        const auto normal_gradient_m = buffer_m.template read<value_type>();
        const auto normal_gradient_p = buffer_p.template read<value_type>();

        const auto JxW = this->all_weights[q_index];

        const auto normal = this->all_normals[q_index];

        const auto value_jump = (value_m - value_p);
        const auto gradient_normal_avg =
          (normal_gradient_m - normal_gradient_p) * 0.5;

        typename FEPointIntegrator::value_type p_value_result    = {};
        typename FEPointIntegrator::value_type p_gradient_result = {};

        p_value_result += value_jump;
        p_value_result -= gradient_normal_avg;
        p_gradient_result -= value_jump;


        phi_m.submit_gradient(outer(p_gradient_result, normal) * 0.5 * JxW, q);
        phi_m.submit_value(p_value_result * JxW, q);
      }

    this->phi_m.test_and_sum(buffer,
                             EvaluationFlags::values |
                               EvaluationFlags::gradients);
  }

  const FESystem<dim>       fe_sub;
  mutable FEPointIntegrator phi_m;
};


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;

  const unsigned int fe_degree            = 1;
  const unsigned int mapping_degree       = 5;
  const unsigned int dim                  = 2;
  const unsigned int n_global_refinements = 0;
  const double       radius               = 0.75;
  const double       outer_radius         = 1.0;
  const MPI_Comm     comm                 = MPI_COMM_WORLD;

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(comm) == 0);

  FE_DGQ<dim>   fe(fe_degree);
  MappingQ<dim> mapping(mapping_degree);
  QGauss<dim>   quadrature(fe_degree + 1);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tria(comm);
      hyper_cube_with_cylindrical_hole_with_tolerance(radius,
                                                      outer_radius,
                                                      0.0,
                                                      tria);

      for (const auto &cell : tria.active_cell_iterators())
        if (!cell->is_artificial())
          if (cell->center().distance({}) > radius)
            cell->set_material_id(1);

      tria.refine_global(n_global_refinements);

      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      AffineConstraints<double> constraints;
      const IndexSet            locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         locally_relevant_dofs);

      typename MatrixFree<dim, Number>::AdditionalData data;
      data.mapping_update_flags =
        update_quadrature_points | update_gradients | update_values;
      data.mapping_update_flags_inner_faces =
        update_values | update_gradients | update_JxW_values;
      data.mapping_update_flags_boundary_faces =
        update_values | update_gradients | update_JxW_values;

      MatrixFree<dim, Number> matrix_free;
      matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);


      using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, 1, Number>;

      TrilinosWrappers::SparsityPattern dsp;

      dsp.reinit(dof_handler.locally_owned_dofs(),
                 dof_handler.get_communicator());

      DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints);

      dsp.compress();

      TrilinosWrappers::SparseMatrix system_matrix;
      system_matrix.reinit(dsp);


      MatrixFreeTools::internal::
        ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
          data_cell;

      MatrixFreeTools::internal::
        ComputeMatrixScratchData<dim, VectorizedArray<Number>, true>
          data_face;

      data_face.dof_numbers               = {0, 0};
      data_face.quad_numbers              = {0, 0};
      data_face.n_components              = {1, 1};
      data_face.first_selected_components = {0, 0};
      data_face.batch_type                = {1, 2};

      data_face
        .op_create = [&](const std::pair<unsigned int, unsigned int> &range) {
        std::vector<
          std::unique_ptr<FEEvaluationData<dim, VectorizedArray<Number>, true>>>
          phi;

        phi.emplace_back(std::make_unique<FEFaceIntegrator>(
          matrix_free, range, true, 0, 0, 0));

        phi.emplace_back(std::make_unique<FEFaceIntegrator>(
          matrix_free, range, false, 0, 0, 0));

        return phi;
      };

      data_face.op_reinit = [](auto &phi, const unsigned batch) {
        static_cast<FEFaceIntegrator &>(*phi[0]).reinit(batch);
        static_cast<FEFaceIntegrator &>(*phi[1]).reinit(batch);
      };

      data_face.op_compute = [&](auto &phi) {
        auto &phi_m = static_cast<FEFaceIntegrator &>(*phi[0]);
        auto &phi_p = static_cast<FEFaceIntegrator &>(*phi[1]);


        VectorizedArrayType mask = 0.0;

        const unsigned int face = phi_m.get_current_cell_index();
        for (unsigned int v = 0;
             v < matrix_free.n_active_entries_per_face_batch(face);
             ++v)
          if (matrix_free.get_face_iterator(face, v, true)
                .first->material_id() !=
              matrix_free.get_face_iterator(face, v, false)
                .first->material_id())
            mask[v] = 1.0;

        phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
        phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

        for (const auto q : phi_m.quadrature_point_indices())
          {
            const auto normal = phi_m.normal_vector(q);

            const auto value_jump = phi_m.get_value(q) - phi_p.get_value(q);
            const auto gradient_normal_avg =
              (phi_m.get_gradient(q) + phi_p.get_gradient(q)) * normal * 0.5;

            typename FEFaceIntegrator::value_type p_value_result    = {};
            typename FEFaceIntegrator::value_type p_gradient_result = {};

            p_value_result += value_jump;
            p_value_result -= gradient_normal_avg;
            p_gradient_result -= value_jump;

            phi_m.submit_gradient(outer(p_gradient_result, normal) * 0.5 * mask,
                                  q);
            phi_p.submit_gradient(outer(p_gradient_result, normal) * 0.5 * mask,
                                  q);
            phi_m.submit_value(p_value_result * mask, q);
            phi_p.submit_value(-p_value_result * mask, q);
          }

        phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
        phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
      };

      MatrixFreeTools::internal::
        ComputeMatrixScratchData<dim, VectorizedArray<Number>, true>
          data_boundary;

      MatrixFreeTools::internal::compute_matrix(matrix_free,
                                                constraints,
                                                data_cell,
                                                data_face,
                                                data_boundary,
                                                system_matrix);

      system_matrix.compress(VectorOperation::add);

      pcout << system_matrix.frobenius_norm() << std::endl;

      if (false)
        system_matrix.print(std::cout);
    }

  if (true)
    {
      parallel::distributed::Triangulation<dim> tria(comm);
      hyper_cube_with_cylindrical_hole(radius, outer_radius, 0.0, tria);

      tria.refine_global(n_global_refinements);

      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      AffineConstraints<double> constraints;
      const IndexSet            locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         locally_relevant_dofs);

      MyCouplingOperator<dim, double> coupling_operator(
        mapping,
        dof_handler,
        constraints,
        quadrature,
        4 * Utilities::pow(2, n_global_refinements + 1),
        radius,
        0,
        0,
        5);


      TrilinosWrappers::SparsityPattern dsp;

      dsp.reinit(dof_handler.locally_owned_dofs(),
                 dof_handler.get_communicator());

      DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints);
      coupling_operator.add_sparsity_pattern_entries(dsp);

      dsp.compress();

      TrilinosWrappers::SparseMatrix system_matrix;
      system_matrix.reinit(dsp);

      coupling_operator.add_system_matrix_entries(system_matrix);

      system_matrix.compress(VectorOperation::add);

      pcout << system_matrix.frobenius_norm() << std::endl;

      if (false)
        system_matrix.print(std::cout);
    }
}
