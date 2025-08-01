// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test CouplingOperator on rotated hypershell-hypershell geometry
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
class MyCouplingOperator : public CouplingEvaluationBase<dim, Number>
{
public:
  using FEPointIntegrator = FEPointEvaluation<1, dim, dim, Number>;
  using value_type        = typename FEPointIntegrator::value_type;

  MyCouplingOperator(const Mapping<dim>    &mapping,
                     const DoFHandler<dim> &dof_handler)
    : fe_sub(dof_handler.get_fe().base_element(
               dof_handler.get_fe().component_to_base_index(0).first),
             1)
    , phi_m(mapping, fe_sub, update_values)
  {
    for (unsigned int i = 0; i < dof_handler.get_fe().n_dofs_per_cell(); ++i)
      relevant_dof_indices.push_back(i);
  }

  unsigned int
  data_size() const override
  {
    return 1;
  }

  const std::vector<unsigned int> &
  get_relevant_dof_indices() const override
  {
    return relevant_dof_indices;
  }

  void
  local_reinit(const typename Triangulation<dim>::cell_iterator &cell,
               const ArrayView<const Point<dim, Number>> &points) const override
  {
    this->phi_m.reinit(cell, points);
  }

  void
  local_evaluate(const CouplingEvaluationData<dim, Number> &data,
                 const Vector<Number>                      &buffer,
                 const unsigned int                         ptr_q,
                 const unsigned int                         q_stride,
                 Number *all_value_m) const override
  {
    (void)data;
    (void)ptr_q;

    this->phi_m.evaluate(buffer, EvaluationFlags::values);

    for (const auto q : this->phi_m.quadrature_point_indices())
      {
        const auto value_m = this->phi_m.get_value(q);

        BufferRW<Number> buffer_m(all_value_m, q * 1 * q_stride);

        buffer_m.write(value_m);
      }
  }

  void
  local_integrate(const CouplingEvaluationData<dim, Number> &data,
                  Vector<Number>                            &buffer,
                  const unsigned int                         ptr_q,
                  const unsigned int                         q_stride,
                  Number                                    *all_value_m,
                  Number *all_value_p) const override
  {
    for (const auto q : this->phi_m.quadrature_point_indices())
      {
        const unsigned int q_index = ptr_q + q;

        BufferRW<Number> buffer_m(all_value_m, q * 1 * q_stride);
        BufferRW<Number> buffer_p(all_value_p, q * 1 * q_stride);

        const auto value_m = buffer_m.template read<value_type>();
        const auto value_p = buffer_p.template read<value_type>();
        const auto JxW     = data.all_weights[q_index];

        phi_m.submit_value((value_m - value_p) * JxW, q);
      }

    this->phi_m.test_and_sum(buffer, EvaluationFlags::values);
  }

  const FESystem<dim>       fe_sub;
  mutable FEPointIntegrator phi_m;

  std::vector<unsigned int> relevant_dof_indices;
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
  const double       radius               = 1.0;
  const MPI_Comm     comm                 = MPI_COMM_WORLD;

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(comm) == 0);

  FE_DGQ<dim>   fe(fe_degree);
  MappingQ<dim> mapping(mapping_degree);
  QGauss<dim>   quadrature(fe_degree + 1);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tria(comm);
      hyper_shell_with_hyper_shell(radius, tria, 1e-7);

      for (const auto &cell : tria.active_cell_iterators())
        if (!cell->is_artificial())
          if (cell->center().distance({}) > 0.5 * radius)
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
                 dof_handler.get_mpi_communicator());

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

        phi_m.evaluate(EvaluationFlags::values);
        phi_p.evaluate(EvaluationFlags::values);

        for (const auto q : phi_m.quadrature_point_indices())
          {
            const auto jump = phi_m.get_value(q) - phi_p.get_value(q);

            phi_m.submit_value(jump * mask, q);
            phi_p.submit_value(-jump * mask, q);
          }


        phi_m.integrate(EvaluationFlags::values);
        phi_p.integrate(EvaluationFlags::values);
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
      hyper_shell_with_hyper_shell(radius, tria, 0.0);

      tria.refine_global(n_global_refinements);

      DoFHandler<dim> dof_handler(tria);
      dof_handler.distribute_dofs(fe);

      AffineConstraints<double> constraints;
      const IndexSet            locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(dof_handler);
      constraints.reinit(dof_handler.locally_owned_dofs(),
                         locally_relevant_dofs);

      const auto mortar_manager = std::make_shared<MortarManagerCircle<dim>>(
        6 * Utilities::pow(2, n_global_refinements),
        quadrature,
        0.5 * radius,
        0.0);

      const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
        coupling_evaluator =
          std::make_shared<MyCouplingOperator<dim, Number>>(mapping,
                                                            dof_handler);

      CouplingOperator<dim, double> coupling_operator(mapping,
                                                      dof_handler,
                                                      constraints,
                                                      coupling_evaluator,
                                                      mortar_manager,
                                                      1,
                                                      2,
                                                      1.0);


      TrilinosWrappers::SparsityPattern dsp;

      dsp.reinit(dof_handler.locally_owned_dofs(),
                 dof_handler.get_mpi_communicator());

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
