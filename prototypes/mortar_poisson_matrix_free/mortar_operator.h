// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#include <deal.II/base/mpi_noncontiguous_partitioner.h>
#include <deal.II/base/mpi_noncontiguous_partitioner.templates.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/tools.h>

#include "./auxiliary_functions.h"
#include "./mortar_evaluation.h"
#include "./mortar_manager.h"

using namespace dealii;


/*-------------- Operator base -------------------------------*/
/**
 * @brief Base class for the Coupling Operator Base
 */
template <int dim, typename Number>
class CouplingOperator
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] bid_m Boundary ID of the face whose outwards-pointing
   *   normal shows in the same direction as the normal provided by
   *   @p mortar_manager.
   * @param[in] bid_p Boundary ID of the face whose outwards-pointing
   *   normal shows in the opposite direction as the normal provided by
   *   @p mortar_manager.
   */
  CouplingOperator(
    const Mapping<dim>                                        &mapping,
    const DoFHandler<dim>                                     &dof_handler,
    const AffineConstraints<Number>                           &constraints,
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator,
    const std::shared_ptr<MortarManagerBase<dim>>              mortar_manager,
    const unsigned int                                         bid_m,
    const unsigned int                                         bid_p,
    const double                                               sip_factor)

    : mapping(mapping)
    , dof_handler(dof_handler)
    , bid_m(bid_m)
    , bid_p(bid_p)
    , evaluator(evaluator)
    , mortar_manager(mortar_manager)
  {
    this->q_data_size          = evaluator->data_size();
    this->relevant_dof_indices = evaluator->get_relevant_dof_indices();
    this->n_dofs_per_cell      = this->relevant_dof_indices.size();

    data.penalty_factor =
      compute_penalty_factor(dof_handler.get_fe().degree, sip_factor);

    // Number of cells
    const unsigned int n_sub_cells = mortar_manager->get_n_total_mortars();

    std::vector<types::global_dof_index> is_local_cell;
    std::vector<types::global_dof_index> is_ghost_cell;

#ifdef DEBUG
    std::vector<double> vec_local_cells(n_sub_cells * 2, 0.0);
    std::vector<double> vec_ghost_cells(n_sub_cells * 2, 0.0);
#endif

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto face_no : cell->face_indices())
          if ((cell->face(face_no)->boundary_id() == bid_m) ||
              (cell->face(face_no)->boundary_id() == bid_p))
            {
              const auto face = cell->face(face_no);

              const auto center = get_face_center(cell, face);

              // Indices of mortars on face of cell.
              const auto indices = mortar_manager->get_mortar_indices(
                center, cell->face(face_no)->boundary_id() == bid_m);

              const auto local_dofs = this->get_dof_indices(cell);

              // Loop over the mortar indices at the rotor-stator
              // interface. The logic of local (rotor) and ghost (stator) is the
              // same as in the previous loop.
              for (unsigned int ii = 0; ii < indices.size(); ++ii)
                {
                  const unsigned int i        = indices[ii];
                  unsigned int       id_local = 0, id_ghost = 0;

                  if (face->boundary_id() == bid_m)
                    {
                      id_local = i;
                      id_ghost = i + n_sub_cells;
                    }
                  else if (face->boundary_id() == bid_p)
                    {
                      id_local = i + n_sub_cells;
                      id_ghost = i;
                    }

#ifdef DEBUG
                  vec_local_cells[id_local] += 1.0;
                  vec_ghost_cells[id_ghost] += 1.0;
#endif
                  is_local_cell.emplace_back(id_local);
                  is_ghost_cell.emplace_back(id_ghost);

                  for (const auto l_dof : local_dofs)
                    dof_indices.emplace_back(l_dof);
                }

              // Weights of quadrature points
              const auto weights = mortar_manager->get_weights(
                get_face_center(cell, face),
                cell->face(face_no)->boundary_id() == bid_m);
              data.all_weights.insert(data.all_weights.end(),
                                      weights.begin(),
                                      weights.end());

              // Normals of quadrature points
              if constexpr (dim == 3)
                {
                  const auto points = mortar_manager->get_points(
                    get_face_center(cell, face),
                    cell->face(face_no)->boundary_id() == bid_m);
                  std::vector<Point<dim, Number>> points_ref(points.size());
                  mapping.transform_points_real_to_unit_cell(cell,
                                                             points,
                                                             points_ref);
                  all_points_ref.insert(all_points_ref.end(),
                                        points_ref.begin(),
                                        points_ref.end());

                  std::vector<Point<dim - 1>> quad;

                  for (const auto p : points_ref)
                    {
                      Point<dim - 1> temp;
                      for (int i = 0, j = 0; i < dim; ++i)
                        if ((face_no / 2) != static_cast<unsigned int>(i))
                          temp[j++] = p[i];

                      if ((dim == 3) && ((face_no / 2) == 1))
                        std::swap(temp[0], temp[1]);

                      quad.emplace_back(temp);
                    }

                  FEFaceValues<dim> fe_face_values(mapping,
                                                   cell->get_fe(),
                                                   quad,
                                                   update_normal_vectors);

                  fe_face_values.reinit(cell, face_no);

                  data.all_normals.insert(
                    data.all_normals.end(),
                    fe_face_values.get_normal_vectors().begin(),
                    fe_face_values.get_normal_vectors().end());
                }
              else
                {
                  auto normals = mortar_manager->get_normals(
                    get_face_center(cell, face),
                    cell->face(face_no)->boundary_id() == bid_m);
                  if (face->boundary_id() == bid_p)
                    for (auto &normal : normals)
                      normal *= -1.0;
                  data.all_normals.insert(data.all_normals.end(),
                                          normals.begin(),
                                          normals.end());

                  if constexpr (dim == 1)
                    {
                      if (face_no == 0)
                        all_points_ref.emplace_back(0.0);
                      else if (face_no == 1)
                        all_points_ref.emplace_back(1.0);
                      else
                        AssertThrow(false, ExcNotImplemented());
                    }
                  else if constexpr (dim == 2)
                    {
                      auto points = mortar_manager->get_points_ref(
                        get_face_center(cell, face),
                        cell->face(face_no)->boundary_id() == bid_m);

                      const bool flip =
                        (face->vertex(0)[0] * face->vertex(1)[1] -
                         face->vertex(0)[1] * face->vertex(1)[0]) < 0.0;

                      if (flip)
                        for (auto &p : points)
                          p[0] = 1.0 - p[0];

                      if (face_no / 2 == 0)
                        {
                          for (auto &p : points)
                            all_points_ref.emplace_back(face_no % 2, p[0]);
                        }
                      else if (face_no / 2 == 1)
                        {
                          for (auto &p : points)
                            all_points_ref.emplace_back(p[0], face_no % 2);
                        }
                      else
                        {
                          AssertThrow(false, ExcNotImplemented());
                        }
                    }
                  else if constexpr (dim == 3)
                    {
                      AssertThrow(false, ExcNotImplemented()); // TODO
                    }
                  else
                    AssertThrow(false, ExcNotImplemented());
                }

              // Penalty parameter
              const Number penalty_parameter = compute_penalty_parameter(cell);

              // Store penalty parameter for all quadrature points
              for (unsigned int i = 0; i < mortar_manager->get_n_points(); ++i)
                data.all_penalty_parameter.emplace_back(penalty_parameter);
            }

#ifdef DEBUG
    Utilities::MPI::sum(vec_local_cells,
                        dof_handler.get_mpi_communicator(),
                        vec_local_cells);
    Utilities::MPI::sum(vec_ghost_cells,
                        dof_handler.get_mpi_communicator(),
                        vec_ghost_cells);

    std::set<unsigned int> vec_local_cells_0;
    std::set<unsigned int> vec_local_cells_2;

    if (Utilities::MPI::this_mpi_process(dof_handler.get_mpi_communicator()) ==
        0)
      {
        for (unsigned int i = 0; i < vec_local_cells.size(); ++i)
          {
            if (vec_local_cells[i] == 0)
              vec_local_cells_0.insert(i);
            if (vec_local_cells[i] > 1)
              vec_local_cells_2.insert(i);
          }

        std::set<unsigned int> vec_ghost_cells_0;
        std::set<unsigned int> vec_ghost_cells_2;

        for (unsigned int i = 0; i < vec_ghost_cells.size(); ++i)
          {
            if (vec_ghost_cells[i] == 0)
              vec_ghost_cells_0.insert(i);
            if (vec_ghost_cells[i] > 1)
              vec_ghost_cells_2.insert(i);
          }

        if (!(vec_local_cells_0.empty() && vec_local_cells_2.empty() &&
              vec_ghost_cells_0.empty() && vec_ghost_cells_2.empty()))
          {
            std::cout << "CouplingOperator mortar matching failed:"
                      << std::endl;

            if (!vec_local_cells_0.empty())
              {
                std::cout << " - some local cells are not owned: ";
                for (const auto &i : vec_local_cells_0)
                  std::cout << i << " ";
                std::cout << std::endl;
              }
            if (!vec_local_cells_2.empty())
              {
                std::cout << " - some local cells are owned multiple times: ";
                for (const auto &i : vec_local_cells_2)
                  std::cout << i << " ";
                std::cout << std::endl;
              }
            if (!vec_ghost_cells_0.empty())
              {
                std::cout << " - some ghost cells are not owned: ";
                for (const auto &i : vec_ghost_cells_0)
                  std::cout << i << " ";
                std::cout << std::endl;
              }
            if (!vec_ghost_cells_2.empty())
              {
                std::cout << " - some ghost cells are owned multiple times: ";
                for (const auto &i : vec_ghost_cells_2)
                  std::cout << i << " ";
                std::cout << std::endl;
              }

            MPI_Barrier(dof_handler.get_mpi_communicator());

            AssertThrow(false, ExcInternalError());
          }
      }
#endif

    // Setup communication
    partitioner.reinit(is_local_cell,
                       is_ghost_cell,
                       dof_handler.get_mpi_communicator());

    // Finalized penalty parameters
    const unsigned n_q_points =
      mortar_manager->get_n_points() / mortar_manager->get_n_mortars();
    std::vector<Number> all_penalty_parameter_ghost(
      data.all_penalty_parameter.size());
    partitioner.export_to_ghosted_array<Number, 0>(data.all_penalty_parameter,
                                                   all_penalty_parameter_ghost,
                                                   n_q_points);
    for (unsigned int i = 0; i < data.all_penalty_parameter.size(); ++i)
      data.all_penalty_parameter[i] =
        std::max(data.all_penalty_parameter[i], all_penalty_parameter_ghost[i]);

    // Finalize DoF indices and update constraints
    dof_indices_ghost.resize(dof_indices.size());
    partitioner.export_to_ghosted_array<types::global_dof_index, 0>(
      dof_indices, dof_indices_ghost, n_dofs_per_cell);

    {
      auto locally_owned_dofs = constraints.get_locally_owned_indices();
      auto constraints_to_make_consistent = constraints.get_local_lines();


      for (unsigned int i = 0; i < dof_indices.size(); ++i)
        {
          constraints_to_make_consistent.add_index(dof_indices[i]);
          constraints_to_make_consistent.add_index(dof_indices_ghost[i]);
        }

      constraints_extended.reinit(locally_owned_dofs,
                                  constraints_to_make_consistent);
      constraints_extended.merge(
        constraints,
        AffineConstraints<Number>::MergeConflictBehavior::no_conflicts_allowed,
        true);

      constraints_extended.make_consistent_in_parallel(
        locally_owned_dofs,
        constraints_to_make_consistent,
        dof_handler.get_mpi_communicator());

      partitioner_extended =
        std::make_shared<const Utilities::MPI::Partitioner>(
          locally_owned_dofs,
          constraints_extended.get_local_lines(),
          dof_handler.get_mpi_communicator());
    }
  }

  /**
   * @brief Return object containing problem constraints
   *
   * @return AffineConstraints
   */
  const AffineConstraints<Number> &
  get_affine_constraints() const
  {
    return constraints_extended;
  }

  /**
   * @brief Add matrix-vector multiplication
   *
   * @param[in, out] dst Destination vector holding the result
   * @param[in] src Input source vector
   */
  template <typename VectorType>
  void
  vmult_add(VectorType &dst, const VectorType &src) const
  {
    VectorType dst_internal;

    if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
      {
        dst_internal.reinit(this->partitioner_extended->locally_owned_range(),
                            this->partitioner_extended->get_mpi_communicator());
      }
    else
      dst_internal.reinit(this->partitioner_extended);

    VectorType src_internal;
    if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
      src_internal.reinit(this->partitioner_extended->locally_owned_range(),
                          this->partitioner_extended->ghost_indices(),
                          this->partitioner_extended->get_mpi_communicator());
    else
      src_internal.reinit(this->partitioner_extended);
    src_internal = src;
    src_internal.update_ghost_values();

    // 1) Evaluate
    unsigned int ptr_q = 0;

    Vector<Number> buffer;

    std::vector<Number> all_values_local(data.all_normals.size() * q_data_size);
    std::vector<Number> all_values_ghost(data.all_normals.size() * q_data_size);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto &face : cell->face_iterators())
          if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
            {
              // Quadrature points at the current face. Note: we process
              // all mortars of the face together here.
              const unsigned int n_q_points = mortar_manager->get_n_points();

              evaluator->local_reinit(cell,
                                      ArrayView<const Point<dim, Number>>(
                                        all_points_ref.data() + ptr_q,
                                        n_q_points));

              buffer.reinit(n_dofs_per_cell);

              const auto local_dofs = this->get_dof_indices(cell);

              for (unsigned int i = 0; i < local_dofs.size(); ++i)
                buffer[i] = src_internal[local_dofs[i]];

              evaluator->local_evaluate(data,
                                        buffer,
                                        ptr_q,
                                        1,
                                        all_values_local.data() +
                                          ptr_q * q_data_size);

              ptr_q += n_q_points;
            }

    // 2) Communicate
    const unsigned n_q_points =
      mortar_manager->get_n_points() / mortar_manager->get_n_mortars();
    partitioner.export_to_ghosted_array<Number, 0>(
      ArrayView<const Number>(reinterpret_cast<Number *>(
                                all_values_local.data()),
                              all_values_local.size()),
      ArrayView<Number>(reinterpret_cast<Number *>(all_values_ghost.data()),
                        all_values_ghost.size()),
      n_q_points * q_data_size);

    // 3) Integrate
    ptr_q = 0;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto &face : cell->face_iterators())
          if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
            {
              // Quadrature points at the current face. Note: we process
              // all mortars of the face together here.
              const unsigned int total_n_q_points =
                mortar_manager->get_n_points();

              evaluator->local_reinit(cell,
                                      ArrayView<const Point<dim, Number>>(
                                        all_points_ref.data() + ptr_q,
                                        total_n_q_points));

              buffer.reinit(n_dofs_per_cell);
              evaluator->local_integrate(data,
                                         buffer,
                                         ptr_q,
                                         1,
                                         all_values_local.data() +
                                           ptr_q * q_data_size,
                                         all_values_ghost.data() +
                                           ptr_q * q_data_size);

              const auto local_dofs = this->get_dof_indices(cell);
              constraints_extended.distribute_local_to_global(buffer,
                                                              local_dofs,
                                                              dst_internal);

              ptr_q += total_n_q_points;
            }

    dst_internal.compress(VectorOperation::add);
    dst.add(1.0, dst_internal);
  }

  /**
   * @brief Add mortar coupling terms in diagonal entries
   *
   * @param[in, out] diagonal Matrix diagonal
   */
  template <typename VectorType>
  void
  add_diagonal_entries(VectorType &diagonal) const
  {
    VectorType diagonal_internal;

    if constexpr (std::is_same_v<VectorType, TrilinosWrappers::MPI::Vector>)
      {
        diagonal_internal.reinit(
          this->partitioner_extended->locally_owned_range(),
          this->partitioner_extended->get_mpi_communicator());
      }
    else
      diagonal_internal.reinit(this->partitioner_extended);

    unsigned int ptr_q = 0;

    Vector<Number>      buffer, diagonal_local;
    std::vector<Number> all_values_local, all_values_ghost;

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto &face : cell->face_iterators())
          if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
            {
              // Quadrature points at the current face. Note: we process
              // all mortars of the face together here.
              const unsigned int n_q_points = mortar_manager->get_n_points();

              evaluator->local_reinit(cell,
                                      ArrayView<const Point<dim, Number>>(
                                        all_points_ref.data() + ptr_q,
                                        n_q_points));

              buffer.reinit(n_dofs_per_cell);
              diagonal_local.reinit(n_dofs_per_cell);
              all_values_local.resize(n_q_points * q_data_size);
              all_values_ghost.resize(n_q_points * q_data_size);

              for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                {
                  // Create i-th basis vector
                  for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                    buffer[j] = static_cast<Number>(i == j);

                  // Interpolate i-th basis vector to the quadrature points
                  evaluator->local_evaluate(
                    data, buffer, ptr_q, 1, all_values_local.data());

                  buffer.reinit(n_dofs_per_cell);

                  // integrate the coupling terms of the mortar using the
                  // interpolated information from the cell
                  evaluator->local_integrate(data,
                                             buffer,
                                             ptr_q,
                                             1,
                                             all_values_local.data(),
                                             all_values_ghost.data());

                  diagonal_local[i] = buffer[i];
                }

              const auto local_dofs = this->get_dof_indices(cell);
              constraints_extended.distribute_local_to_global(
                diagonal_local, local_dofs, diagonal_internal);

              ptr_q += n_q_points;
            }

    diagonal_internal.compress(VectorOperation::add);
    diagonal.add(1.0, diagonal_internal);
  }


  /**
   * @brief Add mortar coupling terms in the sparsity pattern
   *
   * @param[in, out] dsp Dynamic Sparsity Pattern object
   */
  void
  add_sparsity_pattern_entries(SparsityPatternBase &dsp) const
  {
    for (unsigned int i = 0; i < dof_indices.size(); i += n_dofs_per_cell)
      {
        std::vector<types::global_dof_index> a(dof_indices.begin() + i,
                                               dof_indices.begin() + i +
                                                 n_dofs_per_cell);
        std::vector<types::global_dof_index> b(dof_indices_ghost.begin() + i,
                                               dof_indices_ghost.begin() + i +
                                                 n_dofs_per_cell);

        constraints_extended.add_entries_local_to_global(a, b, dsp);
        constraints_extended.add_entries_local_to_global(b, a, dsp);
      }
  }

  /**
   * @brief Add mortar coupling terms in the system matrix
   *
   * @param[in, out] system_matrix System matrix
   */
  void
  add_system_matrix_entries(TrilinosWrappers::SparseMatrix &system_matrix) const
  {
    std::vector<Number> all_values_local(data.all_normals.size() *
                                         n_dofs_per_cell * q_data_size);
    std::vector<Number> all_values_ghost(data.all_normals.size() *
                                         n_dofs_per_cell * q_data_size);

    unsigned int ptr_q = 0;

    Vector<Number> buffer;

    // 1) Evaluate
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto &face : cell->face_iterators())
          if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
            {
              // Quadrature points at the current face. Note: we process
              // all mortars of the face together here.
              const unsigned int n_q_points = mortar_manager->get_n_points();

              // Initialize coupling evaluator with the current cell and
              // the relevant quadrature points
              evaluator->local_reinit(cell,
                                      ArrayView<const Point<dim, Number>>(
                                        all_points_ref.data() + ptr_q,
                                        n_q_points));

              // Initialize buffer to store information at dof level
              buffer.reinit(n_dofs_per_cell);

              for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                {
                  // Create i-th basis vector
                  for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                    buffer[j] = static_cast<Number>(i == j);

                  // Interpolate i-th basis vector to the quadrature points
                  evaluator->local_evaluate(data,
                                            buffer,
                                            ptr_q,
                                            n_dofs_per_cell,
                                            all_values_local.data() +
                                              (ptr_q * n_dofs_per_cell + i) *
                                                q_data_size);
                }

              ptr_q += n_q_points;
            }

    const unsigned n_q_points =
      mortar_manager->get_n_points() / mortar_manager->get_n_mortars();

    // 2) Communicate: export data from local to ghost side
    partitioner.export_to_ghosted_array<Number, 0>(
      ArrayView<const Number>(reinterpret_cast<Number *>(
                                all_values_local.data()),
                              all_values_local.size()),
      ArrayView<Number>(reinterpret_cast<Number *>(all_values_ghost.data()),
                        all_values_ghost.size()),
      n_dofs_per_cell * n_q_points * q_data_size);


    ptr_q                 = 0;
    unsigned int ptr_dofs = 0;

    // 3) Integrate
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        for (const auto &face : cell->face_iterators())
          if ((face->boundary_id() == bid_m) || (face->boundary_id() == bid_p))
            {
              // Number of mortars attached to the current cell (i.e. 1 for
              // aligned rotor-stator meshes and 2 for non-aligned case)
              const unsigned int n_mortars = mortar_manager->get_n_mortars();

              // Loop over mortars
              for (unsigned int m = 0; m < n_mortars; ++m)
                {
                  evaluator->local_reinit(cell,
                                          ArrayView<const Point<dim, Number>>(
                                            all_points_ref.data() + ptr_q,
                                            n_q_points));

                  // Loop over local and ghost cells attached the mortar
                  for (unsigned int b = 0; b < 2; ++b)
                    {
                      FullMatrix<Number> cell_matrix(n_dofs_per_cell,
                                                     n_dofs_per_cell);

                      // Loop over cell dofs and integrate the coupling terms
                      // of the mortar using the interpolated information from
                      // the cell
                      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                        {
                          buffer.reinit(n_dofs_per_cell);
                          if (b == 0) // local cell
                            evaluator->local_integrate(
                              data,
                              buffer,
                              ptr_q,
                              n_dofs_per_cell,
                              all_values_local.data() +
                                (ptr_q * n_dofs_per_cell + i) * q_data_size,
                              nullptr);
                          else // ghost cell
                            evaluator->local_integrate(
                              data,
                              buffer,
                              ptr_q,
                              n_dofs_per_cell,
                              nullptr,
                              all_values_ghost.data() +
                                (ptr_q * n_dofs_per_cell + i) * q_data_size);

                          // Copy data from buffer to cell matrix
                          for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                            cell_matrix[j][i] = buffer[j];
                        }

                      // Vector of local dof indices from the cell in the
                      // negative
                      // ('mortar') side
                      std::vector<types::global_dof_index> local_dof_indices(
                        dof_indices.begin() + ptr_dofs,
                        dof_indices.begin() + ptr_dofs + n_dofs_per_cell);

                      if (b == 0) // local cell -> local-local block
                        {
                          constraints_extended.distribute_local_to_global(
                            cell_matrix, local_dof_indices, system_matrix);
                        }
                      else // ghost cell -> local-ghost block
                        {
                          std::vector<types::global_dof_index>
                            local_dof_indices_ghost(dof_indices_ghost.begin() +
                                                      ptr_dofs,
                                                    dof_indices_ghost.begin() +
                                                      ptr_dofs +
                                                      n_dofs_per_cell);

                          constraints_extended.distribute_local_to_global(
                            cell_matrix,
                            local_dof_indices,
                            local_dof_indices_ghost,
                            system_matrix);
                        }
                    }

                  ptr_dofs += n_dofs_per_cell;

                  ptr_q += n_q_points;
                }
            }

    AssertDimension(ptr_q, data.all_normals.size());
    AssertDimension(ptr_dofs, dof_indices.size());
  }


private:
  /**
   * @brief Compute penalty factor used in weak imposition of coupling at the rotor-stator interface
   *
   * @param[in] degree Polynomial degree of the FE approximation
   * @param[in] factor Penalty factor (akin to penalty factor in SIPG)
   *
   * @return penalty factor value
   * penalty_factor = (degree + 1)^2
   */
  static Number
  compute_penalty_factor(const unsigned int degree, const Number factor)
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  /**
   * @brief Compute penalty parameter in a cell
   *
   * @param[in] cell Cell iterator
   * @return Penalty parameter
   *
   * @return penalty parameter value from SIPG method
   * penalty_parameter = (A(∂Ω_e \ Γ_h)/2 + A(∂Ω_e ∩ Γ_h))/V(Ω_e)
   */
  Number
  compute_penalty_parameter(
    const typename Triangulation<dim>::cell_iterator &cell) const
  {
    const unsigned int degree = dof_handler.get_fe().degree;

    FE_Nothing<dim> fe_nothing;

    dealii::QGauss<dim>   quadrature(degree + 1);
    dealii::FEValues<dim> fe_values(mapping,
                                    fe_nothing,
                                    quadrature,
                                    dealii::update_JxW_values);

    dealii::QGauss<dim - 1>   face_quadrature(degree + 1);
    dealii::FEFaceValues<dim> fe_face_values(mapping,
                                             fe_nothing,
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
          (cell->at_boundary(f) && !cell->has_periodic_neighbor(f) &&
           (cell->face(f)->boundary_id() != bid_m &&
            cell->face(f)->boundary_id() != bid_p)) ?
            1. :
            0.5;

        for (unsigned int q = 0; q < face_quadrature.size(); ++q)
          surface_area += fe_face_values.JxW(q) * factor;
      }

    return surface_area / volume;
  }


  /**
   * @brief Returns angle of a point (cell center)
   *
   * @param[in] cell Cell iterator
   * @param[in] face Face iterator
   *
   * @return Angle in radians
   */
  Point<dim>
  get_face_center(const typename Triangulation<dim>::cell_iterator &cell,
                  const typename Triangulation<dim>::face_iterator &face) const
  {
    return mapping.transform_unit_to_real_cell(
      cell, MappingQ1<dim>().transform_real_to_unit_cell(cell, face->center()));
  }

  /**
   * @brief Returns dof indices
   *
   * @param[in] cell Cell iterator
   */
  std::vector<types::global_dof_index>
  get_dof_indices(
    const typename DoFHandler<dim>::active_cell_iterator &cell) const
  {
    std::vector<types::global_dof_index> local_dofs_all(
      dof_handler.get_fe().n_dofs_per_cell());
    cell->get_dof_indices(local_dofs_all);

    std::vector<types::global_dof_index> local_dofs(n_dofs_per_cell);

    for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
      local_dofs[i] = local_dofs_all[relevant_dof_indices[i]];

    return local_dofs;
  }


  /// Mapping of the domain
  const Mapping<dim> &mapping;
  /// DoFHandler associated to the triangulation
  const DoFHandler<dim> &dof_handler;

  std::vector<std::tuple<std::vector<double>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         typename Triangulation<dim>::active_cell_iterator,
                         std::vector<Point<dim>>,
                         std::vector<Tensor<1, dim, Number>>,
                         Number>>
    all_intersections;

protected:
  /// Number of data points per quadrature point
  unsigned int q_data_size;

  /// Boundary ID of the inner domain (rotor)
  const unsigned int bid_m;
  /// Boundary ID of the outer domain (stator)
  const unsigned int bid_p;

  /// List of relevant DoF indices per cell
  std::vector<unsigned int> relevant_dof_indices;
  /// Number of DoFs per cell
  unsigned int n_dofs_per_cell;

  Utilities::MPI::NoncontiguousPartitioner partitioner;

  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> dof_indices_ghost;

  /// Vectors storing information at quadrature points for all cells at the
  /// rotor-stator interface
  std::vector<Point<dim, Number>>     all_points_ref;
  CouplingEvaluationData<dim, Number> data;

  /// Constraints extended according to mortar entries
  AffineConstraints<Number>                          constraints_extended;
  std::shared_ptr<const Utilities::MPI::Partitioner> partitioner_extended;

  std::shared_ptr<CouplingEvaluationBase<dim, Number>> evaluator;
  std::shared_ptr<MortarManagerBase<dim>>              mortar_manager;
};


/**
 * @brief Base class for the operator base
 */
template <int dim,
          int n_components,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class OperatorBase : public Subscriptor
{
public:
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  OperatorBase() = default;

  OperatorBase(const Mapping<dim>              &mapping,
               const DoFHandler<dim>           &dof_handler,
               const AffineConstraints<Number> &constraints,
               const Quadrature<dim>           &quadrature)
  {
    reinit(mapping, dof_handler, constraints, quadrature);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;
  }

  /**
   * @brief Create coupling operator
   */
  void
  add_coupling(const std::shared_ptr<MortarManagerBase<dim>> mortar_manager,
               const unsigned int                            bid_0,
               const unsigned int                            bid_1,
               const double                                  sip_factor = 1.0)
  {
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationSIPG<dim, n_components, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler());

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.cell_loop(
      &OperatorBase<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      {
        // apply constraints
        // TODO: only apply relevant constraints
        const auto &constraints = coupling_operator->get_affine_constraints();
        constraints.distribute(const_cast<VectorType &>(src));
        src.update_ghost_values();

        coupling_operator->vmult_add(dst, src);

        constraints.set_zero(const_cast<VectorType &>(src));
        constraints.set_zero(dst);
      }

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    matrix_free.initialize_dof_vector(diagonal);

    MatrixFreeTools::compute_diagonal(
      matrix_free,
      diagonal,
      &OperatorBase<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell_single,
      this);

    // add coupling terms
    if (coupling_operator)
      coupling_operator->add_diagonal_entries(diagonal);

    for (auto &i : diagonal)
      i = (i != 0.0) ? (1.0 / i) : 1.0;
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

    const auto constraints = (coupling_operator != nullptr) ?
                               (&coupling_operator->get_affine_constraints()) :
                               (&matrix_free.get_affine_constraints());

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::compute_matrix(
          matrix_free,
          *constraints,
          system_matrix,
          &OperatorBase<dim, n_components, Number, VectorizedArrayType>::
            do_vmult_cell_single,
          this);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


protected:
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

  virtual void
  do_vmult_cell_single(FECellIntegrator &phi) const = 0;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;

private:
};


template <int dim,
          int n_components,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class PoissonOperator
  : public OperatorBase<dim, n_components, Number, VectorizedArrayType>
{
public:
  PoissonOperator() = default;

  PoissonOperator(const MappingQ<dim>             &mapping,
                  const DoFHandler<dim>           &dof_handler,
                  const AffineConstraints<Number> &constraints,
                  const Quadrature<dim>           &quadrature)
    : OperatorBase<dim, n_components, Number, VectorizedArrayType>(mapping,
                                                                   dof_handler,
                                                                   constraints,
                                                                   quadrature)
  {}

  void
  do_vmult_cell_single(
    typename OperatorBase<dim, n_components, Number, VectorizedArrayType>::
      FECellIntegrator &phi) const override
  {
    phi.evaluate(EvaluationFlags::gradients);

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      phi.submit_gradient(phi.get_gradient(q), q);

    phi.integrate(EvaluationFlags::gradients);
  }
};


template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class StokesOperator
  : public OperatorBase<dim, dim + 1, Number, VectorizedArrayType>
{
public:
  using BaseClass = OperatorBase<dim, dim + 1, Number, VectorizedArrayType>;
  using FECellIntegrator = typename BaseClass::FECellIntegrator;

  StokesOperator(const MappingQ<dim>             &mapping,
                 const DoFHandler<dim>           &dof_handler,
                 const AffineConstraints<Number> &constraints,
                 const Quadrature<dim>           &quadrature,
                 const double                     delta_1_scaling)
    : BaseClass(mapping, dof_handler, constraints, quadrature)
    , delta_1_scaling(delta_1_scaling)
  {}

  void
  do_vmult_cell_single(FECellIntegrator &phi) const override
  {
    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        delta_1_scaling *
        this->matrix_free.get_cell_iterator(cell, v)->minimum_vertex_distance();

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        typename FECellIntegrator::value_type    value_result    = {};
        typename FECellIntegrator::gradient_type gradient_result = {};

        const auto value    = phi.get_value(q);
        const auto gradient = phi.get_gradient(q);

        const VectorizedArray<Number>                 p_value = value[dim];
        const Tensor<1, dim, VectorizedArray<Number>> p_gradient =
          gradient[dim];

        Tensor<2, dim, VectorizedArray<Number>> u_gradient;

        for (int d = 0; d < dim; ++d)
          u_gradient[d] = gradient[d];

        // a)     (ε(v), 2νε(u))
        if (true)
          symm_scalar_product_add(gradient_result,
                                  u_gradient,
                                  VectorizedArrayType(2.0));
        else
          for (int d = 0; d < dim; ++d)
            gradient_result[d] += u_gradient[d];

        // b)   - (div(v), p)
        for (int d = 0; d < dim; ++d)
          gradient_result[d][d] -= p_value;

        // c)     (q, div(u))
        for (int d = 0; d < dim; ++d)
          value_result[dim] -= u_gradient[d][d];

        // d) δ_1 (∇q, ∇p)
        gradient_result[dim] = -delta_1 * p_gradient;

        phi.submit_value(value_result, q);
        phi.submit_gradient(gradient_result, q);
      }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

private:
  const double delta_1_scaling;
};


/**
 * @brief Base class for the operator base
 */
template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class GeneralStokesOperator : public Subscriptor
{
public:
  using FECellIntegratorU =
    FEEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType>;
  using FECellIntegratorP =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;


  GeneralStokesOperator(const MappingQ<dim>             &mapping,
                        const DoFHandler<dim>           &dof_handler,
                        const AffineConstraints<Number> &constraints,
                        const Quadrature<dim>           &quadrature,
                        const double                     delta_1_scaling,
                        const bool weak_velocity_divergence_term = false)
    : delta_1_scaling(delta_1_scaling)
    , weak_velocity_divergence_term(weak_velocity_divergence_term)
  {
    reinit(mapping, dof_handler, constraints, quadrature);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;
  }


  /**
   * @brief Create coupling operator
   *
   * @param[in] n_subdivisions Number of cells at the interface between inner
   * and outer domains
   * @param[in] radius Radius at the interface between inner and outer domains
   * @param[in] rotate_pi Rotation angle for the inner domain
   * @param[in] bid_0 Boundary ID of inner domain (rotor)
   * @param[in] bid_1 Boundary ID of outer domain (stator)
   * @param[in] sip_factor Penalty factor (akin to symmetric interior penalty
   * factor in SIPG)
   */
  void
  add_coupling(const std::shared_ptr<MortarManagerBase<dim>> mortar_manager,
               const unsigned int                            bid_0,
               const unsigned int                            bid_1,
               const double                                  sip_factor = 1.0)
  {
    const bool is_p_disc = matrix_free.get_dof_handler()
                             .get_fe()
                             .base_element(matrix_free.get_dof_handler()
                                             .get_fe()
                                             .component_to_base_index(dim)
                                             .first)
                             .n_dofs_per_vertex() == 0;

    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationStokes<dim, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler(),
          !is_p_disc,
          weak_velocity_divergence_term);

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.cell_loop(
      &GeneralStokesOperator<dim, Number, VectorizedArrayType>::do_vmult_cell,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      coupling_operator->vmult_add(dst, src);

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    matrix_free.initialize_dof_vector(diagonal);

    MatrixFreeTools::internal::
      ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
        data_cell;

    data_cell.dof_numbers               = {0, 0};
    data_cell.quad_numbers              = {0, 0};
    data_cell.n_components              = {dim, 1};
    data_cell.first_selected_components = {0, dim};
    data_cell.batch_type                = {0, 0};

    data_cell
      .op_create = [&](const std::pair<unsigned int, unsigned int> &range) {
      std::vector<
        std::unique_ptr<FEEvaluationData<dim, VectorizedArray<Number>, false>>>
        phi;

      phi.emplace_back(
        std::make_unique<FECellIntegratorU>(matrix_free, range, 0, 0, 0));

      phi.emplace_back(
        std::make_unique<FECellIntegratorP>(matrix_free, range, 0, 0, dim));

      return phi;
    };

    data_cell.op_reinit = [](auto &phi, const unsigned batch) {
      static_cast<FECellIntegratorU &>(*phi[0]).reinit(batch);
      static_cast<FECellIntegratorP &>(*phi[1]).reinit(batch);
    };

    data_cell.op_compute = [&](auto &phi) {
      auto &phi_0 = static_cast<FECellIntegratorU &>(*phi[0]);
      auto &phi_1 = static_cast<FECellIntegratorP &>(*phi[1]);

      do_vmult_cell_single(phi_0, phi_1);
    };

    std::vector<VectorType *> diagonal_global_components(1);
    diagonal_global_components[0] = &diagonal;

    MatrixFreeTools::internal::compute_diagonal(
      matrix_free, data_cell, {}, {}, diagonal, diagonal_global_components);

    // add coupling terms
    if (coupling_operator)
      coupling_operator->add_diagonal_entries(diagonal);

    for (auto &i : diagonal)
      i = (i != 0.0) ? (1.0 / i) : 1.0;
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

    auto constraints = &matrix_free.get_affine_constraints();

    AffineConstraints<Number> affine_constraints_tmp;

    if (coupling_operator)
      {
        affine_constraints_tmp.copy_from(
          coupling_operator->get_affine_constraints());

        affine_constraints_tmp.close();
      }

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
            data_cell;

        data_cell.dof_numbers               = {0, 0};
        data_cell.quad_numbers              = {0, 0};
        data_cell.n_components              = {dim, 1};
        data_cell.first_selected_components = {0, dim};
        data_cell.batch_type                = {0, 0};

        data_cell.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, false>>>
              phi;

            phi.emplace_back(
              std::make_unique<FECellIntegratorU>(matrix_free, range, 0, 0, 0));

            phi.emplace_back(std::make_unique<FECellIntegratorP>(
              matrix_free, range, 0, 0, dim));

            return phi;
          };

        data_cell.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FECellIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FECellIntegratorP &>(*phi[1]).reinit(batch);
        };

        data_cell.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FECellIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FECellIntegratorP &>(*phi[1]);

          do_vmult_cell_single(phi_0, phi_1);
        };

        MatrixFreeTools::internal::compute_matrix(
          matrix_free, *constraints, data_cell, {}, {}, system_matrix);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


private:
  void
  do_vmult_cell(const MatrixFree<dim, Number>               &matrix_free,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegratorU integrator_u(matrix_free, range, 0, 0, 0);
    FECellIntegratorP integrator_p(matrix_free, range, 0, 0, dim);
    for (unsigned cell = range.first; cell < range.second; ++cell)
      {
        integrator_u.reinit(cell);
        integrator_p.reinit(cell);

        integrator_u.read_dof_values(src);
        integrator_p.read_dof_values(src);
        do_vmult_cell_single(integrator_u, integrator_p);
        integrator_u.distribute_local_to_global(dst);
        integrator_p.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(FECellIntegratorU &phi_u, FECellIntegratorP &phi_p) const
  {
    phi_u.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi_u.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        delta_1_scaling *
        this->matrix_free.get_cell_iterator(cell, v)->minimum_vertex_distance();

    for (unsigned int q = 0; q < phi_u.n_q_points; ++q)
      {
        typename FECellIntegratorP::value_type    p_value_result    = {};
        typename FECellIntegratorP::gradient_type p_gradient_result = {};
        typename FECellIntegratorU::value_type    u_value_result    = {};
        typename FECellIntegratorU::gradient_type u_gradient_result = {};

        const auto p_value    = phi_p.get_value(q);
        const auto p_gradient = phi_p.get_gradient(q);

        const auto u_value    = phi_u.get_value(q);
        const auto u_gradient = phi_u.get_gradient(q);

        // (ε(v), 2νε(u))
        if (false)
          {
            symm_scalar_product_add(u_gradient_result,
                                    u_gradient,
                                    VectorizedArrayType(2.0));
          }
        else
          {
            u_gradient_result = u_gradient;
          }

        // - (div(v), p)
        if constexpr (dim == 1)
          u_gradient_result[0] -= p_value;
        else
          for (int d = 0; d < dim; ++d)
            u_gradient_result[d][d] -= p_value;

        if (weak_velocity_divergence_term)
          {
            // - (∇q, u)
            if constexpr (dim == 1)
              p_gradient_result[0] -= u_value;
            else
              p_gradient_result -= u_value;
          }
        else
          {
            // + (q, div(u))
            if constexpr (dim == 1)
              p_value_result += u_gradient[0];
            else
              for (int d = 0; d < dim; ++d)
                p_value_result += u_gradient[d][d];
          }

        // δ_1 (∇q, ∇p)
        if (delta_1_scaling != 0.0)
          p_gradient_result += delta_1 * p_gradient;

        phi_p.submit_value(p_value_result, q);
        phi_p.submit_gradient(p_gradient_result, q);

        phi_u.submit_value(u_value_result, q);
        phi_u.submit_gradient(u_gradient_result, q);
      }

    phi_u.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;

  const double delta_1_scaling;
  const bool   weak_velocity_divergence_term;
};



/**
 * @brief Base class for the operator base
 */
template <int dim,
          int n_components,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class PoissonOperatorDG : public Subscriptor
{
public:
  using FECellIntegrator =
    FEEvaluation<dim, -1, 0, n_components, Number, VectorizedArrayType>;
  using FEFaceIntegrator =
    FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  PoissonOperatorDG() = default;

  PoissonOperatorDG(const Mapping<dim>              &mapping,
                    const DoFHandler<dim>           &dof_handler,
                    const AffineConstraints<Number> &constraints,
                    const Quadrature<dim>           &quadrature)
  {
    reinit(mapping, dof_handler, constraints, quadrature);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_JxW_values;
    data.mapping_update_flags_inner_faces =
      update_values | update_gradients | update_JxW_values;
    data.mapping_update_flags_boundary_faces =
      update_values | update_gradients | update_JxW_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;

    compute_penalty_parameters();

    panalty_factor = compute_pentaly_factor(dof_handler.get_fe().degree, 1.0);
  }

  /**
   * @brief Create coupling operator
   */
  void
  add_coupling(const unsigned int n_subdivisions,
               const double       radius,
               const double       rotate_pi,
               const unsigned int bid_0,
               const unsigned int bid_1,
               const double       sip_factor = 1.0)
  {
    const std::shared_ptr<MortarManagerBase<dim>> mortar_manager =
      std::make_shared<MortarManagerCircle<dim>>(n_subdivisions,
                                                 radius,
                                                 matrix_free.get_quadrature(),
                                                 rotate_pi);

    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationSIPG<dim, n_components, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler());

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.loop(
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_face,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_boundary,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      coupling_operator->vmult_add(dst, src);

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    matrix_free.initialize_dof_vector(diagonal);

    MatrixFreeTools::compute_diagonal(
      matrix_free,
      diagonal,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        do_vmult_cell_single,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_face_cell,
      &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
        local_apply_boundary_cell,
      this);

    // add coupling terms
    if (coupling_operator)
      coupling_operator->add_diagonal_entries(diagonal);

    for (auto &i : diagonal)
      i = (i != 0.0) ? (1.0 / i) : 1.0;
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

    const auto constraints = (coupling_operator != nullptr) ?
                               (&coupling_operator->get_affine_constraints()) :
                               (&matrix_free.get_affine_constraints());

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::compute_matrix(
          matrix_free,
          *constraints,
          system_matrix,
          &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
            do_vmult_cell_single,
          &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
            local_apply_face_cell,
          &PoissonOperatorDG<dim, n_components, Number, VectorizedArrayType>::
            local_apply_boundary_cell,
          this);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


protected:
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
  local_apply_face(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegrator phi_m(data, true);
    FEFaceIntegrator phi_p(data, false);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_m.reinit(face);
        phi_p.reinit(face);
        phi_m.read_dof_values(src);
        phi_p.read_dof_values(src);
        local_apply_face_cell(phi_m, phi_p);
        phi_m.distribute_local_to_global(dst);
        phi_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_boundary(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegrator phi(data, true);
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi.reinit(face);
        phi.read_dof_values(src);
        local_apply_boundary_cell(phi);
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

  void
  local_apply_face_cell(FEFaceIntegrator &phi_m, FEFaceIntegrator &phi_p) const
  {
    phi_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto sigma = std::max(phi_m.read_cell_data(penalty_parameters),
                                phi_p.read_cell_data(penalty_parameters)) *
                       panalty_factor;

    for (const auto q : phi_m.quadrature_point_indices())
      {
        const auto average_value =
          (phi_m.get_value(q) - phi_p.get_value(q)) * 0.5;
        const auto average_valgrad =
          average_value * 2. * sigma -
          (phi_m.get_normal_derivative(q) + phi_p.get_normal_derivative(q)) *
            0.5;

        phi_m.submit_normal_derivative(-average_value, q);
        phi_p.submit_normal_derivative(-average_value, q);
        phi_m.submit_value(average_valgrad, q);
        phi_p.submit_value(-average_valgrad, q);
      }

    phi_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  void
  local_apply_boundary_cell(FEFaceIntegrator &phi) const
  {
    if ((phi.boundary_id() == 0 || phi.boundary_id() == 5))
      {
        const VectorizedArrayType zero = 0.0;

        for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
          phi.begin_dof_values()[i] = zero;
        return;
      }

    phi.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto sigma = phi.read_cell_data(penalty_parameters) * panalty_factor;

    for (const auto q : phi.quadrature_point_indices())
      {
        const auto average_value = phi.get_value(q);
        const auto average_valgrad =
          average_value * sigma * 2.0 - phi.get_normal_derivative(q);

        phi.submit_normal_derivative(-average_value, q);
        phi.submit_value(average_valgrad, q);
      }

    phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  void
  compute_penalty_parameters()
  {
    const unsigned int n_cells =
      matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
    penalty_parameters.resize(n_cells);

    const auto &mapping = *matrix_free.get_mapping_info().mapping;
    const auto &fe      = matrix_free.get_dof_handler().get_fe();

    FEValues<dim> fe_values(mapping,
                            fe,
                            matrix_free.get_quadrature(),
                            update_JxW_values);

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     matrix_free.get_face_quadrature(),
                                     update_JxW_values);

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(cell);
           ++v)
        {
          const auto dealii_cell = matrix_free.get_cell_iterator(cell, v);
          fe_values.reinit(dealii_cell);

          // compute cell volume
          Number volume = 0.0;
          for (const auto q : fe_values.quadrature_point_indices())
            volume += fe_values.JxW(q);

          // compute surface area
          Number surface_area = 0.0;
          for (const auto f : dealii_cell->face_indices())
            {
              fe_face_values.reinit(dealii_cell, f);

              const Number factor = (dealii_cell->at_boundary(f) &&
                                     !dealii_cell->has_periodic_neighbor(f)) ?
                                      1. :
                                      0.5;

              for (const auto q : fe_face_values.quadrature_point_indices())
                surface_area += fe_face_values.JxW(q) * factor;
            }

          penalty_parameters[cell][v] = surface_area / volume;
        }
  }

  Number
  compute_pentaly_factor(const unsigned int degree, const Number factor) const
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  AlignedVector<VectorizedArrayType> penalty_parameters;
  VectorizedArrayType                panalty_factor;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;
};


/**
 * @brief Base class for the operator base
 */
template <int dim,
          typename Number,
          typename VectorizedArrayType = VectorizedArray<Number>>
class GeneralStokesOperatorDG : public Subscriptor
{
public:
  using FECellIntegratorU =
    FEEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType>;
  using FEFaceIntegratorU =
    FEFaceEvaluation<dim, -1, 0, dim, Number, VectorizedArrayType>;
  using FECellIntegratorP =
    FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;
  using FEFaceIntegratorP =
    FEFaceEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType>;

  using VectorType = LinearAlgebra::distributed::Vector<Number>;


  GeneralStokesOperatorDG(const MappingQ<dim>             &mapping,
                          const DoFHandler<dim>           &dof_handler,
                          const AffineConstraints<Number> &constraints,
                          const Quadrature<dim>           &quadrature,
                          const double                     sip_factor = 1.0,
                          const bool   weak_pressure_gradient_term    = true,
                          const bool   weak_velocity_divergence_term  = true,
                          const double delta_1_scaling                = 0.0)
    : weak_pressure_gradient_term(weak_pressure_gradient_term)
    , weak_velocity_divergence_term(weak_velocity_divergence_term)
    , delta_1_scaling(delta_1_scaling)
  {
    reinit(mapping, dof_handler, constraints, quadrature, sip_factor);
  }

  void
  reinit(const Mapping<dim>              &mapping,
         const DoFHandler<dim>           &dof_handler,
         const AffineConstraints<Number> &constraints,
         const Quadrature<dim>           &quadrature,
         const double                     sip_factor = 1.0)
  {
    typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points | update_gradients | update_values;
    data.mapping_update_flags_inner_faces =
      update_values | update_gradients | update_JxW_values;
    data.mapping_update_flags_boundary_faces =
      update_values | update_gradients | update_JxW_values;

    matrix_free.reinit(mapping, dof_handler, constraints, quadrature, data);

    valid_system = false;

    this->sip_factor = sip_factor;

    compute_penalty_parameters();

    panalty_factor =
      compute_pentaly_factor(dof_handler.get_fe().degree, sip_factor);
  }

  /**
   * @brief Create coupling operator
   */
  void
  add_coupling(const std::shared_ptr<MortarManagerBase<dim>> mortar_manager,
               const unsigned int                            bid_0,
               const unsigned int                            bid_1)
  {
    const std::shared_ptr<CouplingEvaluationBase<dim, Number>>
      coupling_evaluator =
        std::make_shared<CouplingEvaluationStokes<dim, Number>>(
          *matrix_free.get_mapping_info().mapping,
          matrix_free.get_dof_handler(),
          weak_pressure_gradient_term,
          weak_velocity_divergence_term);

    coupling_operator = std::make_shared<CouplingOperator<dim, Number>>(
      *matrix_free.get_mapping_info().mapping,
      matrix_free.get_dof_handler(),
      matrix_free.get_affine_constraints(),
      coupling_evaluator,
      mortar_manager,
      bid_0,
      bid_1,
      sip_factor);

    coupling_bids.insert(bid_0);
    coupling_bids.insert(bid_1);

    compute_penalty_parameters();
  }

  virtual types::global_dof_index
  m() const
  {
    if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
      return this->matrix_free.get_dof_handler().n_dofs(
        this->matrix_free.get_mg_level());
    else
      return this->matrix_free.get_dof_handler().n_dofs();
  }

  Number
  el(unsigned int, unsigned int) const
  {
    DEAL_II_NOT_IMPLEMENTED();
    return 0;
  }

  void
  initialize_dof_vector(VectorType &dst) const
  {
    matrix_free.initialize_dof_vector(dst);
  }

  void
  vmult(VectorType &dst, const VectorType &src) const
  {
    src.update_ghost_values();

    matrix_free.loop(
      &GeneralStokesOperatorDG<dim, Number, VectorizedArrayType>::do_vmult_cell,
      &GeneralStokesOperatorDG<dim, Number, VectorizedArrayType>::
        local_apply_face,
      &GeneralStokesOperatorDG<dim, Number, VectorizedArrayType>::
        local_apply_boundary,
      this,
      dst,
      src,
      true);

    // apply coupling terms
    if (coupling_operator)
      {
        // apply constraints
        // TODO: only apply relevant constraints
        const auto &constraints = coupling_operator->get_affine_constraints();
        constraints.distribute(const_cast<VectorType &>(src));
        src.update_ghost_values();

        coupling_operator->vmult_add(dst, src);

        constraints.set_zero(const_cast<VectorType &>(src));
        constraints.set_zero(dst);
      }

    src.zero_out_ghost_values();
  }

  void
  Tvmult(VectorType &dst, const VectorType &src) const
  {
    vmult(dst, src);
  }

  void
  compute_inverse_diagonal(VectorType &diagonal) const
  {
    (void)diagonal;

    AssertThrow(false, ExcNotImplemented());
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

    auto constraints = &matrix_free.get_affine_constraints();

    AffineConstraints<Number> affine_constraints_tmp;

    if (coupling_operator)
      {
        affine_constraints_tmp.copy_from(
          coupling_operator->get_affine_constraints());
        affine_constraints_tmp.close();
      }

    if (system_matrix.m() == 0 || system_matrix.n() == 0)
      {
        system_matrix.clear();

        TrilinosWrappers::SparsityPattern dsp;

        dsp.reinit(dof_handler.locally_owned_dofs(),
                   dof_handler.get_mpi_communicator());

        DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, *constraints);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_sparsity_pattern_entries(dsp);

        dsp.compress();

        system_matrix.reinit(dsp);
      }

    if (this->valid_system == false)
      {
        system_matrix = 0.0;

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, false>
            data_cell;

        data_cell.dof_numbers               = {0, 0};
        data_cell.quad_numbers              = {0, 0};
        data_cell.n_components              = {dim, 1};
        data_cell.first_selected_components = {0, dim};
        data_cell.batch_type                = {0, 0};

        data_cell.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, false>>>
              phi;

            phi.emplace_back(
              std::make_unique<FECellIntegratorU>(matrix_free, range, 0, 0, 0));

            phi.emplace_back(std::make_unique<FECellIntegratorP>(
              matrix_free, range, 0, 0, dim));

            return phi;
          };

        data_cell.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FECellIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FECellIntegratorP &>(*phi[1]).reinit(batch);
        };

        data_cell.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FECellIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FECellIntegratorP &>(*phi[1]);

          do_vmult_cell_single(phi_0, phi_1);
        };

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, true>
            data_face;

        data_face.dof_numbers               = {0, 0, 0, 0};
        data_face.quad_numbers              = {0, 0, 0, 0};
        data_face.n_components              = {dim, dim, 1, 1};
        data_face.first_selected_components = {0, 0, dim, dim};
        data_face.batch_type                = {1, 2, 1, 2};

        data_face.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, true>>>
              phi;

            phi.emplace_back(std::make_unique<FEFaceIntegratorU>(
              matrix_free, range, true, 0, 0, 0));

            phi.emplace_back(std::make_unique<FEFaceIntegratorU>(
              matrix_free, range, false, 0, 0, 0));

            phi.emplace_back(std::make_unique<FEFaceIntegratorP>(
              matrix_free, range, true, 0, 0, dim));

            phi.emplace_back(std::make_unique<FEFaceIntegratorP>(
              matrix_free, range, false, 0, 0, dim));

            return phi;
          };

        data_face.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FEFaceIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FEFaceIntegratorU &>(*phi[1]).reinit(batch);
          static_cast<FEFaceIntegratorP &>(*phi[2]).reinit(batch);
          static_cast<FEFaceIntegratorP &>(*phi[3]).reinit(batch);
        };

        data_face.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FEFaceIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FEFaceIntegratorU &>(*phi[1]);
          auto &phi_2 = static_cast<FEFaceIntegratorP &>(*phi[2]);
          auto &phi_3 = static_cast<FEFaceIntegratorP &>(*phi[3]);

          local_apply_face_cell(phi_0, phi_1, phi_2, phi_3);
        };

        MatrixFreeTools::internal::
          ComputeMatrixScratchData<dim, VectorizedArray<Number>, true>
            data_boundary;

        data_boundary.dof_numbers               = {0, 0};
        data_boundary.quad_numbers              = {0, 0};
        data_boundary.n_components              = {dim, 1};
        data_boundary.first_selected_components = {0, dim};
        data_boundary.batch_type                = {1, 1};

        data_boundary.op_create =
          [&](const std::pair<unsigned int, unsigned int> &range) {
            std::vector<std::unique_ptr<
              FEEvaluationData<dim, VectorizedArray<Number>, true>>>
              phi;

            phi.emplace_back(std::make_unique<FEFaceIntegratorU>(
              matrix_free, range, true, 0, 0, 0));
            phi.emplace_back(std::make_unique<FEFaceIntegratorP>(
              matrix_free, range, true, 0, 0, dim));

            return phi;
          };

        data_boundary.op_reinit = [](auto &phi, const unsigned batch) {
          static_cast<FEFaceIntegratorU &>(*phi[0]).reinit(batch);
          static_cast<FEFaceIntegratorP &>(*phi[1]).reinit(batch);
        };

        data_boundary.op_compute = [&](auto &phi) {
          auto &phi_0 = static_cast<FEFaceIntegratorU &>(*phi[0]);
          auto &phi_1 = static_cast<FEFaceIntegratorP &>(*phi[1]);

          local_apply_boundary_cell(phi_0, phi_1);
        };

        MatrixFreeTools::internal::compute_matrix(matrix_free,
                                                  *constraints,
                                                  data_cell,
                                                  data_face,
                                                  data_boundary,
                                                  system_matrix);

        // apply coupling terms
        if (coupling_operator)
          coupling_operator->add_system_matrix_entries(system_matrix);

        system_matrix.compress(VectorOperation::add);

        this->valid_system = true;
      }
  }


private:
  void
  do_vmult_cell(const MatrixFree<dim, Number>               &matrix_free,
                VectorType                                  &dst,
                const VectorType                            &src,
                const std::pair<unsigned int, unsigned int> &range) const
  {
    FECellIntegratorU integrator_u(matrix_free, range, 0, 0, 0);
    FECellIntegratorP integrator_p(matrix_free, range, 0, 0, dim);
    for (unsigned cell = range.first; cell < range.second; ++cell)
      {
        integrator_u.reinit(cell);
        integrator_p.reinit(cell);

        integrator_u.read_dof_values(src);
        integrator_p.read_dof_values(src);
        do_vmult_cell_single(integrator_u, integrator_p);
        integrator_u.distribute_local_to_global(dst);
        integrator_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_face(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegratorU phi_u_m(data, face_range, true, 0, 0, 0);
    FEFaceIntegratorU phi_u_p(data, face_range, false, 0, 0, 0);
    FEFaceIntegratorP phi_p_m(data, face_range, true, 0, 0, dim);
    FEFaceIntegratorP phi_p_p(data, face_range, false, 0, 0, dim);

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_u_m.reinit(face);
        phi_u_p.reinit(face);
        phi_p_m.reinit(face);
        phi_p_p.reinit(face);
        phi_u_m.read_dof_values(src);
        phi_u_p.read_dof_values(src);
        phi_p_m.read_dof_values(src);
        phi_p_p.read_dof_values(src);
        local_apply_face_cell(phi_u_m, phi_u_p, phi_p_m, phi_p_p);
        phi_u_m.distribute_local_to_global(dst);
        phi_u_p.distribute_local_to_global(dst);
        phi_p_m.distribute_local_to_global(dst);
        phi_p_p.distribute_local_to_global(dst);
      }
  }

  void
  local_apply_boundary(
    const MatrixFree<dim, Number, VectorizedArrayType> &data,
    VectorType                                         &dst,
    const VectorType                                   &src,
    const std::pair<unsigned int, unsigned int>        &face_range) const
  {
    FEFaceIntegratorU phi_u(data, face_range, true, 0, 0, 0);
    FEFaceIntegratorP phi_p(data, face_range, true, 0, 0, dim);
    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        phi_u.reinit(face);
        phi_p.reinit(face);
        phi_u.read_dof_values(src);
        phi_p.read_dof_values(src);
        local_apply_boundary_cell(phi_u, phi_p);
        phi_u.distribute_local_to_global(dst);
        phi_p.distribute_local_to_global(dst);
      }
  }

  void
  do_vmult_cell_single(FECellIntegratorU &phi_u, FECellIntegratorP &phi_p) const
  {
    phi_u.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

    const auto cell = phi_u.get_current_cell_index();

    VectorizedArrayType delta_1;
    for (unsigned int v = 0;
         v < this->matrix_free.n_active_entries_per_cell_batch(cell);
         ++v)
      delta_1[v] =
        delta_1_scaling *
        this->matrix_free.get_cell_iterator(cell, v)->minimum_vertex_distance();

    for (unsigned int q = 0; q < phi_u.n_q_points; ++q)
      {
        typename FECellIntegratorP::value_type    p_value_result    = {};
        typename FECellIntegratorP::gradient_type p_gradient_result = {};
        typename FECellIntegratorU::value_type    u_value_result    = {};
        typename FECellIntegratorU::gradient_type u_gradient_result = {};

        const auto p_value    = phi_p.get_value(q);
        const auto p_gradient = phi_p.get_gradient(q);
        const auto u_value    = phi_u.get_value(q);
        const auto u_gradient = phi_u.get_gradient(q);

        if (true /*Laplace term*/)
          {
            // (∇v, ∇u)
            u_gradient_result += u_gradient;
          }

        if (weak_pressure_gradient_term)
          {
            // - (div(v), p)
            for (int d = 0; d < dim; ++d)
              u_gradient_result[d][d] -= p_value;
          }
        else
          {
            // + (v, ∇p)
            u_value_result += p_gradient;
          }

        if (weak_velocity_divergence_term)
          {
            // - (∇q, u)
            p_gradient_result -= u_value * vel_div_sign;
          }
        else
          {
            // + (q, div(u))
            for (int d = 0; d < dim; ++d)
              p_value_result += u_gradient[d][d] * vel_div_sign;
          }

        // δ_1 (∇q, ∇p)
        if (delta_1_scaling != 0.0)
          p_gradient_result += delta_1 * p_gradient;

        phi_p.submit_value(p_value_result, q);
        phi_p.submit_gradient(p_gradient_result, q);

        phi_u.submit_value(u_value_result, q);
        phi_u.submit_gradient(u_gradient_result, q);
      }

    phi_u.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
  }

  void
  local_apply_face_cell(FEFaceIntegratorU &phi_u_m,
                        FEFaceIntegratorU &phi_u_p,
                        FEFaceIntegratorP &phi_p_m,
                        FEFaceIntegratorP &phi_p_p) const
  {
    phi_u_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_u_p.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.evaluate(EvaluationFlags::values);
    phi_p_p.evaluate(EvaluationFlags::values);

    const auto sigma = std::max(phi_u_m.read_cell_data(penalty_parameters),
                                phi_u_p.read_cell_data(penalty_parameters)) *
                       panalty_factor;

    VectorizedArrayType mask = 1.0;

    const unsigned int face = phi_u_m.get_current_cell_index();
    for (unsigned int v = 0;
         v < matrix_free.n_active_entries_per_face_batch(face);
         ++v)
      if (matrix_free.get_face_iterator(face, v, true).first->material_id() !=
          matrix_free.get_face_iterator(face, v, false).first->material_id())
        mask[v] = 0.0;

    for (const auto q : phi_u_m.quadrature_point_indices())
      {
        const auto u_value_avg =
          (phi_u_m.get_value(q) + phi_u_p.get_value(q)) * 0.5;
        const auto u_value_jump = phi_u_m.get_value(q) - phi_u_p.get_value(q);
        const auto u_gradient_avg =
          (phi_u_m.get_gradient(q) + phi_u_p.get_gradient(q)) * 0.5;
        const auto p_value_avg =
          (phi_p_m.get_value(q) + phi_p_p.get_value(q)) * 0.5;
        const auto normal = phi_u_m.normal_vector(q);

        typename FECellIntegratorU::value_type u_normal_gradient_avg_result =
          {};
        typename FECellIntegratorU::value_type u_value_jump_result = {};
        typename FECellIntegratorP::value_type p_value_jump_result = {};

        if (true /*Laplace term*/)
          {
            // - (n avg(∇v), jump(u))
            u_normal_gradient_avg_result -= u_value_jump;

            // - (jump(v), avg(∇u) n)
            u_value_jump_result -= u_gradient_avg * normal;

            // + (jump(v), σ jump(u))
            u_value_jump_result += sigma * u_value_jump;
          }

        if (weak_pressure_gradient_term)
          {
            // + (jump(v), avg(p) n)
            u_value_jump_result += p_value_avg * normal;
          }
        else
          {
            // nothing to do
          }

        if (weak_velocity_divergence_term)
          {
            // + (jump(q), avg(u) n)
            p_value_jump_result += u_value_avg * normal * vel_div_sign;
          }
        else
          {
            // - (avg(q), jump(u) n)
            p_value_jump_result -= 0.5 * u_value_jump * normal;
          }

        phi_u_m.submit_normal_derivative(u_normal_gradient_avg_result * 0.5, q);
        phi_u_p.submit_normal_derivative(u_normal_gradient_avg_result * 0.5, q);
        phi_u_m.submit_value(u_value_jump_result, q);
        phi_u_p.submit_value(-u_value_jump_result, q);
        phi_p_m.submit_value(p_value_jump_result, q);
        phi_p_p.submit_value(-p_value_jump_result, q);
      }

    phi_u_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_u_p.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.integrate(EvaluationFlags::values);
    phi_p_p.integrate(EvaluationFlags::values);
  }

  void
  local_apply_boundary_cell(FEFaceIntegratorU &phi_u_m,
                            FEFaceIntegratorP &phi_p_m) const
  {
    if (coupling_bids.find(phi_u_m.boundary_id()) != coupling_bids.end())
      {
        for (unsigned int i = 0; i < phi_u_m.dofs_per_cell; ++i)
          phi_u_m.begin_dof_values()[i] = 0.0;
        for (unsigned int i = 0; i < phi_p_m.dofs_per_cell; ++i)
          phi_p_m.begin_dof_values()[i] = 0.0;
        return;
      }

    phi_u_m.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.evaluate(EvaluationFlags::values);

    const auto sigma =
      phi_u_m.read_cell_data(penalty_parameters) * panalty_factor;

    for (const auto q : phi_u_m.quadrature_point_indices())
      {
        const auto u_value_avg =
          (phi_u_m.get_value(q) - phi_u_m.get_value(q)) * 0.5;
        const auto u_value_jump = phi_u_m.get_value(q) + phi_u_m.get_value(q);
        const auto u_gradient_avg =
          (phi_u_m.get_gradient(q) + phi_u_m.get_gradient(q)) * 0.5;
        const auto p_value_avg =
          (phi_p_m.get_value(q) + phi_p_m.get_value(q)) * 0.5;
        const auto normal = phi_u_m.normal_vector(q);

        typename FECellIntegratorU::value_type u_normal_gradient_avg_result =
          {};
        typename FECellIntegratorU::value_type u_value_jump_result = {};
        typename FECellIntegratorP::value_type p_value_jump_result = {};

        if (true /*Laplace term*/)
          {
            // - (n avg(∇v), jump(u))
            u_normal_gradient_avg_result -= u_value_jump;

            // - (jump(v), avg(∇u) n)
            u_value_jump_result -= u_gradient_avg * normal;

            // + (jump(v), σ jump(u))
            u_value_jump_result += sigma * u_value_jump;
          }

        if (weak_pressure_gradient_term)
          {
            // + (jump(v), avg(p) n)
            u_value_jump_result += p_value_avg * normal;
          }
        else
          {
            // nothing to do
          }

        if (weak_velocity_divergence_term)
          {
            // + (jump(q), avg(u) n)
            p_value_jump_result += u_value_avg * normal * vel_div_sign;
          }
        else
          {
            // nothing to do
          }

        phi_u_m.submit_normal_derivative(u_normal_gradient_avg_result * 0.5, q);
        phi_u_m.submit_value(u_value_jump_result, q);
        phi_p_m.submit_value(p_value_jump_result, q);
      }

    phi_u_m.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
    phi_p_m.integrate(EvaluationFlags::values);
  }

  void
  compute_penalty_parameters()
  {
    const unsigned int n_cells =
      matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
    penalty_parameters.resize(n_cells);

    const auto &mapping = *matrix_free.get_mapping_info().mapping;
    const auto &fe      = matrix_free.get_dof_handler().get_fe();

    FEValues<dim> fe_values(mapping,
                            fe,
                            matrix_free.get_quadrature(),
                            update_JxW_values);

    FEFaceValues<dim> fe_face_values(mapping,
                                     fe,
                                     matrix_free.get_face_quadrature(),
                                     update_JxW_values);

    for (unsigned int cell = 0; cell < n_cells; ++cell)
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(cell);
           ++v)
        {
          const auto dealii_cell = matrix_free.get_cell_iterator(cell, v);
          fe_values.reinit(dealii_cell);

          // compute cell volume
          Number volume = 0.0;
          for (const auto q : fe_values.quadrature_point_indices())
            volume += fe_values.JxW(q);

          // compute surface area
          Number surface_area = 0.0;
          for (const auto f : dealii_cell->face_indices())
            {
              fe_face_values.reinit(dealii_cell, f);

              const Number factor =
                (dealii_cell->at_boundary(f) &&
                 !dealii_cell->has_periodic_neighbor(f) &&
                 (coupling_bids.find(dealii_cell->face(f)->boundary_id()) ==
                  coupling_bids.end())) ?
                  1. :
                  0.5;

              for (const auto q : fe_face_values.quadrature_point_indices())
                surface_area += fe_face_values.JxW(q) * factor;
            }

          penalty_parameters[cell][v] = surface_area / volume;
        }
  }

  Number
  compute_pentaly_factor(const unsigned int degree, const Number factor) const
  {
    return factor * (degree + 1.0) * (degree + 1.0);
  }

  const bool   weak_pressure_gradient_term;
  const bool   weak_velocity_divergence_term;
  const double delta_1_scaling;

  const double vel_div_sign = +1.0;

  mutable double sip_factor;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  mutable TrilinosWrappers::SparseMatrix       system_matrix;
  mutable bool                                 valid_system;

  AlignedVector<VectorizedArrayType> penalty_parameters;
  VectorizedArrayType                panalty_factor;

  std::shared_ptr<CouplingOperator<dim, Number>> coupling_operator;

  std::set<unsigned int> coupling_bids;
};
