// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/multiphysics.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/fluid_dynamics_matrix_free.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_matrix_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#if DEAL_II_VERSION_GTE(9, 8, 0)
#  include <deal.II/multigrid/mg_transfer_matrix_free.templates.h>
#else
#  include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>
#endif
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

/**
 * @brief A base class of preconditioners used by the smoother.
 */
template <typename VectorType>
class PreconditionBase : public EnableObserverPointer
{
public:
  /**
   * @brief Constructor.
   */
  PreconditionBase()
    : comm(MPI_COMM_WORLD)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(comm) == 0)
    , timer(pcout, TimerOutput::never, TimerOutput::wall_times)
  {}

  /**
   * @brief Apply preconditioner.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  virtual void
  vmult(VectorType &dst, const VectorType &src) const = 0;

  /**
   * @brief Print timers.
   */
  void
  timer_print() const
  {
    timer.print_wall_time_statistics(comm);
  }

  /**
   * @brief Reset timers.
   */
  void
  timer_reset() const
  {
    timer.reset();
  }

protected:
  /// Communicator used for timer output.
  const MPI_Comm comm;

  /// Stream used for timer output.
  ConditionalOStream pcout;

  /// Timer.
  mutable TimerOutput timer;
};

/**
 * @brief A wrapper class around a preconditioner to be
 * able to handle vector types that are not supported
 * by the preconditioner.
 */
template <typename VectorType, typename VectorTypePrecondition>
class PreconditionAdapter : public PreconditionBase<VectorType>
{
public:
  /**
   * Constructor.
   *
   * @param[in] preconditioner Preconditioner.
   */
  template <typename PreconditionType>
  PreconditionAdapter(const std::shared_ptr<PreconditionType> &preconditioner)
    : fu_preconditioner([preconditioner](VectorTypePrecondition       &dst,
                                         const VectorTypePrecondition &src) {
      preconditioner->vmult(dst, src);
    })
  {}

  /**
   * @brief Apply preconditioner.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const override
  {
    if constexpr (std::is_same_v<VectorType, VectorTypePrecondition>)
      {
        // vector is supported: nothing to do
        fu_preconditioner(dst, src);
      }
    else
      {
        // vector is not supported: copy vectors
        src_ = src;
        dst_ = dst;
        fu_preconditioner(dst_, src_);
        dst = dst_;
      }
  }

  /**
   * @brief Apply preconditioner.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  template <
    typename U = VectorTypePrecondition,
    std::enable_if_t<!std::is_same_v<VectorType, U>, VectorType> * = nullptr>
  void
  vmult(U &dst, const U &src) const
  {
    fu_preconditioner(dst, src);
  }

private:
  /// Preconditioner.
  const std::function<void(VectorTypePrecondition &,
                           const VectorTypePrecondition &)>
    fu_preconditioner;

  /// Source vector supported by the preconditioner.
  mutable VectorTypePrecondition src_;

  /// Destination vector supported by the preconditioner.
  mutable VectorTypePrecondition dst_;
};

/**
 * @brief A wrapper class around DiagonalMatrix class from deal.II.
 */
template <typename VectorType>
class MyDiagonalMatrix : public PreconditionBase<VectorType>
{
public:
  /**
   * @brief Constructor.

   * @param[in] diagonal Vector containing the diagonal of the matrix.
   */
  MyDiagonalMatrix(const VectorType &diagonal)
    : diagonal_matrix(diagonal)
  {}

  /**
   * @brief Apply preconditioner.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const override
  {
    diagonal_matrix.vmult(dst, src);
  }

private:
  /// DiagonalMatrix class from deal.II.
  DiagonalMatrix<VectorType> diagonal_matrix;
};

/**
 * @brief An Additive Schwarz preconditioner.
 */
template <typename VectorType>
class PreconditionASM : public PreconditionBase<VectorType>
{
private:
  /// Weighting type.
  enum class WeightingType
  {
    none,
    left,
    right,
    symm
  };

public:
  /// Value type.
  using Number = typename VectorType::value_type;

  /**
   * @brief Constructor.
   */
  PreconditionASM()
    : weighting_type(WeightingType::left)
  {}

  /**
   * @brief Initialize inverses of blocks.
   */
  template <int dim,
            typename GlobalSparseMatrixType,
            typename GlobalSparsityPattern,
            typename Number>
  void
  initialize(
    const DoFHandler<dim>           &dof_handler,
    const GlobalSparseMatrixType    &global_sparse_matrix,
    const GlobalSparsityPattern     &global_sparsity_pattern,
    const AffineConstraints<Number> &constraints,
    const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
    const unsigned int                                        mg_level)
  {
    this->timer.enter_subsection("asm::indices");

    std::vector<std::vector<types::global_dof_index>> patches;

    const auto add_indices = [&](const auto &local_dof_indices) {
      std::vector<types::global_dof_index> local_dof_indices_temp;

      for (const auto i : local_dof_indices)
        if (!constraints.is_constrained(i))
          local_dof_indices_temp.emplace_back(i);

      if (!local_dof_indices_temp.empty())
        patches.push_back(local_dof_indices_temp);
    };

    if (mg_level == numbers::invalid_unsigned_int)
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->is_locally_owned() == false)
            continue;

          std::vector<types::global_dof_index> local_dof_indices(
            cell->get_fe().n_dofs_per_cell());
          cell->get_dof_indices(local_dof_indices);

          add_indices(local_dof_indices);
        }
    else
      for (const auto &cell : dof_handler.mg_cell_iterators_on_level(mg_level))
        {
          if (cell->is_locally_owned_on_level() == false)
            continue;

          std::vector<types::global_dof_index> local_dof_indices(
            cell->get_fe().n_dofs_per_cell());
          cell->get_mg_dof_indices(local_dof_indices);

          add_indices(local_dof_indices);
        }

    const IndexSet locally_owned_dofs =
      (mg_level == numbers::invalid_unsigned_int) ?
        dof_handler.locally_owned_dofs() :
        dof_handler.locally_owned_mg_dofs(mg_level);

    IndexSet ghost_dofs(locally_owned_dofs.size());
    for (const auto &indices : patches)
      ghost_dofs.add_indices(indices.begin(), indices.end());

    std::shared_ptr<const Utilities::MPI::Partitioner> partition =
      std::make_shared<Utilities::MPI::Partitioner>(locally_owned_dofs,
                                                    ghost_dofs,
                                                    MPI_COMM_WORLD);

    if (dealii::internal::is_partitioner_contained(partition, partitioner))
      {
        partition = partitioner;
      }
    else
      {
        src_internal.reinit(partition);
        dst_internal.reinit(partition);
      }

    // convert indices to local ones
    this->patches.resize(patches.size());
    for (unsigned int c = 0; c < patches.size(); ++c)
      {
        this->patches[c].resize(0);
        this->patches[c].reserve(patches[c].size());

        for (const auto &i : patches[c])
          this->patches[c].emplace_back(partition->global_to_local(i));
      }

    this->timer.leave_subsection("asm::indices");
    this->timer.enter_subsection("asm::restrict");

    std::vector<FullMatrix<Number>> blocks;

    SparseMatrixTools::restrict_to_full_matrices(global_sparse_matrix,
                                                 global_sparsity_pattern,
                                                 patches,
                                                 blocks);

    this->timer.leave_subsection("asm::restrict");
    this->timer.enter_subsection("asm::invert");

    this->blocks.resize(blocks.size());

    for (unsigned int b = 0; b < blocks.size(); ++b)
      {
        this->blocks[b] =
          LAPACKFullMatrix<Number>(blocks[b].m(), blocks[b].n());
        this->blocks[b] = blocks[b];
        this->blocks[b].compute_lu_factorization();
      }

    this->timer.leave_subsection("asm::invert");

    if (weighting_type != WeightingType::none)
      {
        this->timer.enter_subsection("asm::weight");
        Vector<Number> vector_weights;
        weights.reinit(partition);

        for (const auto &patch : patches)
          {
            const unsigned int dofs_per_cell = patch.size();
            vector_weights.reinit(dofs_per_cell);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              vector_weights[i] = 1.0;

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              weights[patch[i]] += vector_weights[i];
          }

        weights.compress(VectorOperation::add);
        for (auto &i : weights)
          i = (weighting_type == WeightingType::symm) ? std::sqrt(1.0 / i) :
                                                        (1.0 / i);
        weights.update_ghost_values();

        this->timer.leave_subsection("asm::weight");
      }
  }

  /**
   * @brief Apply preconditioner.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const override
  {
    this->timer.enter_subsection("asm::vmult");

    const auto &src_ptr = (src_internal.size() != 0) ? src_internal : src;
    auto       &dst_ptr = (dst_internal.size() != 0) ? dst_internal : dst;

    dst_ptr = 0.0;
    if (src_internal.size() != 0)
      src_internal.copy_locally_owned_data_from(src);
    src_ptr.update_ghost_values();

    Vector<Number> vector_src, vector_dst, vector_weights;

    for (unsigned int c = 0; c < patches.size(); ++c)
      {
        const unsigned int dofs_per_cell = patches[c].size();

        vector_src.reinit(dofs_per_cell);
        vector_dst.reinit(dofs_per_cell);
        if (weighting_type != WeightingType::none)
          vector_weights.reinit(dofs_per_cell);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          vector_src[i] = src_ptr.local_element(patches[c][i]);

        if (weighting_type != WeightingType::none)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            vector_weights[i] = weights.local_element(patches[c][i]);

        if (weighting_type == WeightingType::symm ||
            weighting_type == WeightingType::right)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            vector_src[i] *= vector_weights[i];

        blocks[c].solve(vector_src);
        vector_dst = vector_src;

        if (weighting_type == WeightingType::symm ||
            weighting_type == WeightingType::left)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            vector_dst[i] *= vector_weights[i];

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          dst_ptr.local_element(patches[c][i]) += vector_dst[i];
      }

    src_ptr.zero_out_ghost_values();
    dst_ptr.compress(VectorOperation::add);

    if (dst_internal.size() != 0)
      dst.copy_locally_owned_data_from(dst_internal);

    this->timer.leave_subsection("asm::vmult");
  }

private:
  /// DoF indices of patches.
  std::vector<std::vector<unsigned int>> patches;

  /// Inverse of patch matrices.
  std::vector<LAPACKFullMatrix<Number>> blocks;

  /// Weights.
  mutable VectorType weights;

  /// Weighting type.
  const WeightingType weighting_type;

  /// Internal destination vector with correct ghosting.
  mutable VectorType dst_internal;

  /// Internal source vector with correct ghosting.
  mutable VectorType src_internal;
};



namespace dealii
{
  /**
   * Coarse grid solver using a preconditioner only. This is a little wrapper,
   * transforming a preconditioner into a coarse grid solver.
   */
  template <class VectorType, class PreconditionerType>
  class MGCoarseGridApplyPreconditioner : public MGCoarseGridBase<VectorType>
  {
  public:
    /**
     * Default constructor.
     */
    MGCoarseGridApplyPreconditioner();

    /**
     * Constructor. Store a pointer to the preconditioner for later use.
     */
    MGCoarseGridApplyPreconditioner(const PreconditionerType &precondition);

    /**
     * Clear the pointer.
     */
    void
    clear();

    /**
     * Initialize new data.
     */
    void
    initialize(const PreconditionerType &precondition);

    /**
     * Implementation of the abstract function.
     */
    virtual void
    operator()(const unsigned int level,
               VectorType        &dst,
               const VectorType  &src) const override;

  private:
    /**
     * Reference to the preconditioner.
     */
    ObserverPointer<
      const PreconditionerType,
      MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>>
      preconditioner;
  };



  template <class VectorType, class PreconditionerType>
  MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>::
    MGCoarseGridApplyPreconditioner()
    : preconditioner(0, typeid(*this).name())
  {}



  template <class VectorType, class PreconditionerType>
  MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>::
    MGCoarseGridApplyPreconditioner(const PreconditionerType &preconditioner)
    : preconditioner(&preconditioner, typeid(*this).name())
  {}



  template <class VectorType, class PreconditionerType>
  void
  MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>::initialize(
    const PreconditionerType &preconditioner_)
  {
    preconditioner = &preconditioner_;
  }



  template <class VectorType, class PreconditionerType>
  void
  MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>::clear()
  {
    preconditioner = 0;
  }



  template <class VectorType, class PreconditionerType>
  void
  MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>::operator()(
    const unsigned int /*level*/,
    VectorType       &dst,
    const VectorType &src) const
  {
    preconditioner->vmult(dst, src);
  }
} // namespace dealii

/**
 * A class that wraps MGTransferMatrixFree and MGTransferMF to allow
 * hp-multigrid in the local-smoothing case.
 */
template <int dim, typename Number>
class MGTransferMatrixFreeWrapper
  : public MGTransferBase<LinearAlgebra::distributed::Vector<Number>>
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<Number>;

  MGTransferMatrixFreeWrapper(const unsigned int min_h_level)
    : min_h_level(min_h_level)
    , n_gc_levels(0)
  {}

  void
  build(const MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> &transfers,
        const DoFHandler<dim>                                    &dof_handler,
        const MGConstrainedDoFs &mg_constrained_dofs,
        const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          &external_partitioners_in)
  {
    if (transfers.min_level() != transfers.max_level())
      {
        this->n_gc_levels = transfers.max_level() - transfers.min_level();

        std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          external_partitioners;

        for (unsigned int l = transfers.min_level(); l <= transfers.max_level();
             ++l)
          external_partitioners.emplace_back(external_partitioners_in[l]);

        gc.initialize_two_level_transfers(transfers);
        gc.build(external_partitioners);
      }

    ls.initialize_constraints(mg_constrained_dofs);

    std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
      external_partitioners(min_h_level + external_partitioners_in.size() -
                            n_gc_levels);

    for (unsigned int i = n_gc_levels; i < external_partitioners_in.size(); ++i)
      external_partitioners[min_h_level + i - n_gc_levels] =
        external_partitioners_in[i];

    ls.build(dof_handler, external_partitioners);
  }

  template <class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &dof_handler,
             MGLevelObject<VectorType>       &dst,
             const InVector                  &src) const
  {
    if ((n_gc_levels > 0) && (dst.max_level() <= n_gc_levels))
      {
        gc.copy_to_mg(dof_handler, dst, src);
        return;
      }

    MGLevelObject<VectorType> dst_gc(0, n_gc_levels);
    MGLevelObject<VectorType> dst_ls(min_h_level,
                                     dst.max_level() + min_h_level -
                                       n_gc_levels);

    for (unsigned int l = dst.min_level(); l <= dst.max_level(); ++l)
      {
        if (l <= n_gc_levels)
          dst_gc[l] = dst[l];

        if (l >= n_gc_levels)
          dst_ls[min_h_level + (l - n_gc_levels)] = dst[l];
      }

    ls.copy_to_mg(dof_handler, dst_ls, src);

    if (n_gc_levels > 0)
      gc.copy_to_mg(dof_handler, dst_gc, dst_ls[dst_ls.min_level()]);

    for (unsigned int l = dst.min_level(); l <= dst.max_level(); ++l)
      if (l < n_gc_levels)
        dst[l] = dst_gc[l];
      else
        dst[l] = dst_ls[min_h_level + (l - n_gc_levels)];
  }

  template <class OutVector, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &dof_handler,
               OutVector                       &dst,
               const MGLevelObject<VectorType> &src) const
  {
    if ((n_gc_levels > 0) && (src.max_level() <= n_gc_levels))
      {
        gc.copy_from_mg(dof_handler, dst, src);
        return;
      }

    MGLevelObject<VectorType> src_ls(min_h_level,
                                     src.max_level() + min_h_level -
                                       n_gc_levels);
    for (unsigned int l = n_gc_levels; l <= src.max_level(); ++l)
      src_ls[min_h_level + (l - n_gc_levels)] = src[l];

    ls.copy_from_mg(dof_handler, dst, src_ls);
  }

  template <typename InVectorType>
  void
  interpolate_to_mg(const DoFHandler<dim>     &dof_handler,
                    MGLevelObject<VectorType> &dst,
                    const InVectorType        &src) const
  {
    if ((n_gc_levels > 0) && (dst.max_level() <= n_gc_levels))
      {
        gc.interpolate_to_mg(dof_handler, dst, src);
        return;
      }

    MGLevelObject<VectorType> dst_gc(0, n_gc_levels);
    MGLevelObject<VectorType> dst_ls(min_h_level,
                                     dst.max_level() + min_h_level -
                                       n_gc_levels);

    for (unsigned int l = dst.min_level(); l <= dst.max_level(); ++l)
      {
        if (l <= n_gc_levels)
          dst_gc[l] = dst[l];

        if (l >= n_gc_levels)
          dst_ls[min_h_level + (l - n_gc_levels)] = dst[l];
      }

    ls.interpolate_to_mg(dof_handler, dst_ls, src);

    if (n_gc_levels > 0)
      gc.interpolate_to_mg(dof_handler, dst_gc, dst_ls[dst_ls.min_level()]);

    for (unsigned int l = dst.min_level(); l <= dst.max_level(); ++l)
      if (l < n_gc_levels)
        dst[l] = dst_gc[l];
      else
        dst[l] = dst_ls[min_h_level + (l - n_gc_levels)];
  }

  void
  prolongate(const unsigned int to_level,
             VectorType        &dst,
             const VectorType  &src) const override
  {
    if (to_level <= n_gc_levels)
      gc.prolongate(to_level, dst, src);
    else
      ls.prolongate(to_level + min_h_level - n_gc_levels, dst, src);
  }

  void
  restrict_and_add(const unsigned int from_level,
                   VectorType        &dst,
                   const VectorType  &src) const override
  {
    if (from_level <= n_gc_levels)
      gc.restrict_and_add(from_level, dst, src);
    else
      ls.restrict_and_add(from_level + min_h_level - n_gc_levels, dst, src);
  }

private:
  const unsigned int min_h_level;
  unsigned int       n_gc_levels;

  MGTransferMatrixFree<dim, Number>           ls;
  MGTransferGlobalCoarsening<dim, VectorType> gc;
};

template <int dim>
MFNavierStokesPreconditionGMGBase<dim>::MFNavierStokesPreconditionGMGBase(
  const SimulationParameters<dim> &simulation_parameters,
  const DoFHandler<dim>           &dof_handler,
  const DoFHandler<dim>           &dof_handler_fe_q_iso_q1)
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , simulation_parameters(simulation_parameters)
  , dof_handler(dof_handler)
  , mg_setup_timer(this->pcout, TimerOutput::never, TimerOutput::wall_times)
  , mg_vmult_timer(this->pcout, TimerOutput::never, TimerOutput::wall_times)
{
  (void)dof_handler_fe_q_iso_q1; // TODO: remove
}

template <int dim>
void
MFNavierStokesPreconditionGMGBase<dim>::reinit(
  const std::shared_ptr<Mapping<dim>>              &mapping,
  const std::shared_ptr<Quadrature<dim>>           &cell_quadrature,
  const std::shared_ptr<Function<dim>>              forcing_function,
  const std::shared_ptr<SimulationControl>         &simulation_control,
  const std::shared_ptr<PhysicalPropertiesManager> &physical_properties_manager,
  const std::shared_ptr<FESystem<dim>>              fe)
{
  AssertThrow(*cell_quadrature ==
                QGauss<dim>(this->dof_handler.get_fe().degree + 1),
              ExcNotImplemented());

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      const auto &triangulation = dof_handler.get_triangulation();

      // Define maximum and minimum level according to triangulation
      const unsigned int n_h_levels  = triangulation.n_global_levels();
      unsigned int       min_h_level = 0;
      const unsigned int max_h_level = n_h_levels - 1;

      // If multigrid number of levels or minimum number of cells in level are
      // specified, change the min level and print levels information

      int mg_min_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_min_level;

      int mg_int_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_int_level;

      int mg_level_min_cells =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_level_min_cells;


      std::vector<unsigned int> n_cells_on_levels(
        triangulation.n_global_levels(), 0);

      for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
        for (const auto &cell : triangulation.cell_iterators_on_level(l))
          if (cell->is_locally_owned_on_level())
            n_cells_on_levels[l]++;

      Utilities::MPI::sum(n_cells_on_levels,
                          this->dof_handler.get_mpi_communicator(),
                          n_cells_on_levels);

      if (mg_min_level != -1)
        {
          AssertThrow(
            mg_min_level <= static_cast<int>(MGTools::max_level_for_coarse_mesh(
                              triangulation)),
            ExcMessage(std::string(
              "The maximum level allowed for the coarse mesh (mg min level) is: " +
              std::to_string(
                MGTools::max_level_for_coarse_mesh(triangulation)) +
              ".")));

          min_h_level = mg_min_level;
        }
      else if (mg_level_min_cells != -1)
        {
          AssertThrow(
            mg_level_min_cells <=
              static_cast<int>(
                n_cells_on_levels[MGTools::max_level_for_coarse_mesh(
                  triangulation)]),
            ExcMessage(
              "The mg level min cells specified are larger than the cells of the finest mg level."));

          for (unsigned int level = min_h_level; level <= max_h_level; ++level)
            if (static_cast<int>(n_cells_on_levels[level]) >=
                mg_level_min_cells)
              {
                min_h_level = level;
                break;
              }
        }

      // p-multigrid
      const auto mg_coarsening_type =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_coarsening_type;

      const auto polynomial_coarsening_sequence =
        MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
          this->dof_handler.get_fe().degree,
          this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_p_coarsening_type);

      std::vector<std::pair<unsigned int, unsigned int>> levels;

      if (mg_coarsening_type ==
          Parameters::LinearSolver::MultigridCoarseningSequenceType::hp)
        {
          // p
          for (const auto i : polynomial_coarsening_sequence)
            levels.emplace_back(min_h_level, i);

          // h
          for (unsigned int i = min_h_level + 1; i <= max_h_level; ++i)
            levels.emplace_back(i, polynomial_coarsening_sequence.back());
        }
      else if (mg_coarsening_type ==
               Parameters::LinearSolver::MultigridCoarseningSequenceType::ph)
        {
          AssertThrow(false, ExcNotImplemented());
        }
      else if (mg_coarsening_type ==
               Parameters::LinearSolver::MultigridCoarseningSequenceType::p)
        {
          AssertThrow(false, ExcNotImplemented());
        }
      else if (mg_coarsening_type ==
               Parameters::LinearSolver::MultigridCoarseningSequenceType::h)
        {
          // h
          for (unsigned int i = min_h_level; i <= max_h_level; ++i)
            levels.emplace_back(i, polynomial_coarsening_sequence.back());
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }

      this->minlevel = 0;
      this->maxlevel = levels.size() - 1;

      this->intlevel = (mg_int_level == -1) ? this->minlevel : mg_int_level;

      std::map<unsigned int, unsigned int> p_map;
      for (const auto &[_, p] : levels)
        p_map.insert({p, p_map.size()});

      dof_handlers.resize(0, p_map.size() - 1);

      for (const auto [p, l] : p_map)
        {
          this->dof_handlers[l].reinit(triangulation);

          // To use elements with linear interpolation for coarse-grid we need
          // to create the min level dof handler with the appropriate element
          // type
          if (this->simulation_parameters.linear_solver
                .at(PhysicsID::fluid_dynamics)
                .mg_use_fe_q_iso_q1 &&
              l == 0)
            {
              AssertThrow(
                mg_coarsening_type ==
                  Parameters::LinearSolver::MultigridCoarseningSequenceType::h,
                ExcNotImplemented());

              const auto points = QGaussLobatto<1>(p + 1).get_points();

              this->dof_handler_fe_q_iso_q1.reinit(triangulation);
              this->dof_handler_fe_q_iso_q1.distribute_dofs(
                FESystem<dim>(FE_Q_iso_Q1<dim>(points), dim + 1));
              this->dof_handler_fe_q_iso_q1.distribute_mg_dofs();
              DoFRenumbering::Cuthill_McKee(this->dof_handler_fe_q_iso_q1);
            }

          this->dof_handlers[l].distribute_dofs(
            FESystem<dim>(FE_Q<dim>(p), dim + 1));
          this->dof_handlers[l].distribute_mg_dofs();
          DoFRenumbering::Cuthill_McKee(this->dof_handlers[l]);
        }

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << std::endl;
          this->pcout << "  -Levels of MG preconditioner:" << std::endl;
          for (unsigned int level = this->minlevel; level <= this->maxlevel;
               ++level)
            this->pcout << "    Level " << level << ": "
                        << this->dof_handlers[p_map[levels[level].second]]
                             .n_dofs(levels[level].first)
                        << " DoFs, " << n_cells_on_levels[levels[level].first]
                        << " cells" << std::endl;
          this->pcout << std::endl;
        }

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity == Parameters::Verbosity::extra_verbose)
        {
          this->pcout << "  -MG vertical communication efficiency: "
                      << MGTools::vertical_communication_efficiency(
                           triangulation)
                      << std::endl;

          this->pcout << "  -MG workload imbalance: "
                      << MGTools::workload_imbalance(triangulation)
                      << std::endl;
        }

      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        partitioners(this->maxlevel - this->minlevel + 1);

      // Local object for constraints of the different levels
      MGLevelObject<AffineConstraints<MGNumber>> level_constraints;

      // Resize all multilevel objects according to level
      this->mg_operators.resize(this->minlevel, this->maxlevel);
      level_constraints.resize(this->minlevel, this->maxlevel);
      this->ls_mg_interface_in.resize(this->minlevel, this->maxlevel);
      this->ls_mg_operators.resize(this->minlevel, this->maxlevel);

      // Fill the level constraints
      this->mg_setup_timer.enter_subsection("Set boundary conditions");

      this->mg_constrained_dofs.clear();
      this->mg_constrained_dofs.resize(dof_handlers.min_level(),
                                       dof_handlers.max_level());

      for (unsigned int l = mg_constrained_dofs.min_level();
           l <= mg_constrained_dofs.max_level();
           ++l)
        this->mg_constrained_dofs[l].initialize(this->dof_handlers[l]);

      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar pressure(dim);

      for (auto const &[id, type] :
           this->simulation_parameters.boundary_conditions.type)
        {
          if (type == BoundaryConditions::BoundaryType::slip)
            {
              std::set<types::boundary_id> no_normal_flux_boundaries;
              no_normal_flux_boundaries.insert(id);
              for (unsigned int level = this->minlevel; level <= this->maxlevel;
                   ++level)
                {
                  const unsigned int p_level = p_map[levels[level].second];

                  AffineConstraints<double> temp_constraints;

                  const IndexSet locally_owned_level_dofs =
                    this->dof_handlers[p_level].locally_owned_mg_dofs(
                      levels[level].first);
                  const IndexSet locally_relevant_level_dofs =
                    DoFTools::extract_locally_relevant_level_dofs(
                      this->dof_handlers[p_level], levels[level].first);
                  temp_constraints.reinit(locally_owned_level_dofs,
                                          locally_relevant_level_dofs);

                  VectorTools::compute_no_normal_flux_constraints_on_level(
                    this->dof_handlers[p_level],
                    0,
                    no_normal_flux_boundaries,
                    temp_constraints,
                    *mapping,
                    this->mg_constrained_dofs[p_level]
                      .get_refinement_edge_indices(levels[level].first),
                    levels[level].first);
                  temp_constraints.close();
                  this->mg_constrained_dofs[p_level].add_user_constraints(
                    levels[level].first, temp_constraints);
                }
            }
          else if (type == BoundaryConditions::BoundaryType::periodic)
            {
              /*already taken into account when mg_constrained_dofs is
               * initialized*/
            }
          else if (type == BoundaryConditions::BoundaryType::periodic_neighbor)
            {
              /*already taken into account when mg_constrained_dofs is
               * initialized*/
            }
          else if (type == BoundaryConditions::BoundaryType::pressure)
            {
              Assert(
                false,
                ExcMessage(
                  "Pressure boundary conditions are not supported by the matrix free application."));
            }
          else if (type == BoundaryConditions::BoundaryType::function_weak)
            {
              /*The function weak boundary condition is implemented in the
               * operators*/
            }
          else if (type == BoundaryConditions::BoundaryType::partial_slip)
            {
              Assert(
                false,
                ExcMessage(
                  "Partial slip boundary conditions are not supported by the matrix free application."));
            }
          else if (type == BoundaryConditions::BoundaryType::outlet)
            {
              /*The directional do-nothing boundary condition is implemented in
               * the operators*/
            }
          else if (type == BoundaryConditions::BoundaryType::noslip ||
                   type == BoundaryConditions::BoundaryType::function)
            {
              std::set<types::boundary_id> dirichlet_boundary_id = {id};

              for (const auto [_, l] : p_map)
                this->mg_constrained_dofs[l].make_zero_boundary_constraints(
                  this->dof_handlers[l],
                  dirichlet_boundary_id,
                  fe->component_mask(velocities));
            }
          else if (type == BoundaryConditions::BoundaryType::none)
            {
              /*Default boundary condition*/
            }
          else
            {
              AssertThrow(
                false,
                ExcMessage(
                  "The boundary condition specified is not supported."));
            }
        }

      this->mg_setup_timer.leave_subsection("Set boundary conditions");

      // Create mg operators for each level and additional operators needed only
      // for local smoothing
      for (unsigned int level = this->minlevel; level <= this->maxlevel;
           ++level)
        {
          const unsigned int p_level = p_map[levels[level].second];

          level_constraints[level].clear();

          const auto &level_dof_handler = this->dof_handlers[p_level];

          const IndexSet owned_dofs =
            level_dof_handler.locally_owned_mg_dofs(levels[level].first);
          const IndexSet relevant_dofs =
            DoFTools::extract_locally_relevant_level_dofs(level_dof_handler,
                                                          levels[level].first);

          level_constraints[level].reinit(owned_dofs, relevant_dofs);

          if (this->simulation_parameters.boundary_conditions
                .fix_pressure_constant &&
              level == this->minlevel)
            {
              types::global_dof_index min_index = numbers::invalid_unsigned_int;

              std::vector<types::global_dof_index> dof_indices;

              // Loop over the cells to identify the min index
              for (const auto &cell :
                   level_dof_handler.mg_cell_iterators_on_level(this->minlevel))
                {
                  if (cell->is_locally_owned_on_level())
                    {
                      const auto &fe = cell->get_fe();

                      dof_indices.resize(fe.n_dofs_per_cell());
                      cell->get_mg_dof_indices(dof_indices);

                      for (unsigned int i = 0; i < dof_indices.size(); ++i)
                        if (fe.system_to_component_index(i).first == dim)
                          min_index = std::min(min_index, dof_indices[i]);
                    }
                }

              // Necessary to find the min across all cores.
              min_index =
                Utilities::MPI::min(min_index,
                                    level_dof_handler.get_mpi_communicator());

              AffineConstraints<double> temp_constraints;
              temp_constraints.reinit(owned_dofs, relevant_dofs);
              if (relevant_dofs.is_element(min_index))
                temp_constraints.add_line(min_index);

              this->mg_constrained_dofs[p_level].add_user_constraints(
                levels[level].first, temp_constraints);
            }

          this->mg_constrained_dofs[p_level].merge_constraints(
            level_constraints[level],
            levels[level].first,
            true,
            false,
            true,
            true);

          this->mg_setup_timer.enter_subsection("Set up operators");

          // Provide appropriate quadrature depending on the type of elements of
          // the level
          Quadrature<dim> quadrature_mg =
            QGauss<dim>(level_dof_handler.get_fe().degree + 1);
          if (this->simulation_parameters.linear_solver
                .at(PhysicsID::fluid_dynamics)
                .mg_use_fe_q_iso_q1 &&
              level == this->minlevel)
            {
              AssertThrow(
                mg_coarsening_type ==
                  Parameters::LinearSolver::MultigridCoarseningSequenceType::h,
                ExcNotImplemented());

              const auto points =
                QGaussLobatto<1>(level_dof_handler.get_fe().degree + 1)
                  .get_points();

              quadrature_mg = QIterated<dim>(QGauss<1>(2), points);
            }

          this->create_level_operator(level);

          this->mg_operators[level]->reinit(
            *mapping,
            (this->simulation_parameters.linear_solver
               .at(PhysicsID::fluid_dynamics)
               .mg_use_fe_q_iso_q1 &&
             level == this->minlevel) ?
              this->dof_handler_fe_q_iso_q1 :
              level_dof_handler,
            level_constraints[level],
            quadrature_mg,
            forcing_function,
            physical_properties_manager,
            this->simulation_parameters.stabilization.stabilization,
            levels[level].first,
            simulation_control,
            this->simulation_parameters.boundary_conditions,
            this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .mg_enable_hessians_jacobian,
            true,
            this->simulation_parameters.mortar_parameters.enable);

          this->ls_mg_operators[level].initialize(*(this->mg_operators)[level]);
          this->ls_mg_interface_in[level].initialize(
            *(this->mg_operators)[level]);

          partitioners[level] =
            this->mg_operators[level]->get_vector_partitioner();

          this->mg_setup_timer.leave_subsection("Set up operators");
        }

      // Create transfer operators
      this->mg_setup_timer.enter_subsection("Create transfer operator");

      this->mg_transfer_ls = std::make_shared<LSTransferType>(min_h_level);

      if (p_map.size() > 1)
        {
          this->transfers.resize(0, p_map.size() - 1);

          for (unsigned int level = transfers.min_level();
               level < transfers.max_level();
               ++level)
            this->transfers[level + 1].reinit(this->dof_handlers[level + 1],
                                              this->dof_handlers[level],
                                              level_constraints[level + 1],
                                              level_constraints[level],
                                              min_h_level,
                                              min_h_level);
        }

      this->mg_transfer_ls->build(
        this->transfers,
        this->dof_handlers[this->dof_handlers.max_level()],
        this->mg_constrained_dofs[this->mg_constrained_dofs.max_level()],
        partitioners);

      this->mg_setup_timer.leave_subsection("Create transfer operator");
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      // Create triangulations
      this->mg_setup_timer.enter_subsection("Create level triangulations");
      this->coarse_grid_triangulations =
        MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
          this->dof_handler.get_triangulation());
      this->mg_setup_timer.leave_subsection("Create level triangulations");

      // Modify the triangulations if multigrid number of levels or minimum
      // number of cells in level are specified
      std::vector<std::shared_ptr<const Triangulation<dim>>> temp;

      int mg_min_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_min_level;

      int mg_int_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_int_level;

      AssertThrow(
        (mg_min_level + 1) <=
          static_cast<int>(this->coarse_grid_triangulations.size()),
        ExcMessage(
          "The mg min level specified is higher than the finest mg level."));

      int mg_level_min_cells =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_level_min_cells;

      AssertThrow(
        mg_level_min_cells <=
          static_cast<int>(this
                             ->coarse_grid_triangulations
                               [this->coarse_grid_triangulations.size() - 1]
                             ->n_global_active_cells()),
        ExcMessage(
          "The mg level min cells specified are larger than the cells of the finest mg level."));

      // find first relevant coarse-grid triangulation
      auto ptr = std::find_if(
        this->coarse_grid_triangulations.begin(),
        this->coarse_grid_triangulations.end() - 1,
        [&mg_min_level, &mg_level_min_cells](const auto &tria) {
          if (mg_min_level != -1) // minimum number of levels
            {
              if ((mg_min_level + 1) <=
                  static_cast<int>(tria->n_global_levels()))
                return true;
            }
          else if (mg_level_min_cells != -1) // minimum number of cells
            {
              if (static_cast<int>(tria->n_global_active_cells()) >=
                  mg_level_min_cells)
                return true;
            }
          else
            {
              return true;
            }
          return false;
        });

      // consider all triangulations from that one
      while (ptr != this->coarse_grid_triangulations.end())
        temp.push_back(*(ptr++));

      this->coarse_grid_triangulations = temp;

      // p-multigrid
      const auto mg_coarsening_type =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_coarsening_type;

      const auto polynomial_coarsening_sequence =
        MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
          this->dof_handler.get_fe().degree,
          this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_p_coarsening_type);

      std::vector<std::pair<unsigned int, unsigned int>> levels;

      if (mg_coarsening_type ==
          Parameters::LinearSolver::MultigridCoarseningSequenceType::hp)
        {
          // p
          for (const auto i : polynomial_coarsening_sequence)
            levels.emplace_back(0, i);

          // h
          for (unsigned int i = 1; i < this->coarse_grid_triangulations.size();
               ++i)
            levels.emplace_back(i, polynomial_coarsening_sequence.back());
        }
      else if (mg_coarsening_type ==
               Parameters::LinearSolver::MultigridCoarseningSequenceType::ph)
        {
          // h
          for (unsigned int i = 0;
               i < this->coarse_grid_triangulations.size() - 1;
               ++i)
            levels.emplace_back(i, polynomial_coarsening_sequence.front());

          // p
          for (const auto i : polynomial_coarsening_sequence)
            levels.emplace_back(this->coarse_grid_triangulations.size() - 1, i);
        }
      else if (mg_coarsening_type ==
               Parameters::LinearSolver::MultigridCoarseningSequenceType::p)
        {
          // p
          for (const auto i : polynomial_coarsening_sequence)
            levels.emplace_back(this->coarse_grid_triangulations.size() - 1, i);
        }
      else if (mg_coarsening_type ==
               Parameters::LinearSolver::MultigridCoarseningSequenceType::h)
        {
          // h
          for (unsigned int i = 0; i < this->coarse_grid_triangulations.size();
               ++i)
            levels.emplace_back(i, polynomial_coarsening_sequence.back());
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
        }

      // Define maximum and minimum level according to triangulations
      this->minlevel = 0;
      this->maxlevel = levels.size() - 1;

      this->intlevel = (mg_int_level == -1) ? this->minlevel : mg_int_level;

      // Local object for constraints of the different levels
      MGLevelObject<AffineConstraints<MGNumber>> constraints;

      // Resize all multilevel objects according to level
      this->mg_operators.resize(this->minlevel, this->maxlevel);
      constraints.resize(this->minlevel, this->maxlevel);
      this->transfers.resize(this->minlevel, this->maxlevel);

      // Distribute DoFs for each level
      this->mg_setup_timer.enter_subsection(
        "Create DoFHandlers and distribute DoFs");
      this->dof_handlers.resize(this->minlevel, this->maxlevel);

      for (unsigned int l = this->minlevel; l <= this->maxlevel; ++l)
        {
          this->dof_handlers[l].reinit(
            *this->coarse_grid_triangulations[levels[l].first]);

          // To use elements with linear interpolation for coarse-grid we need
          // to create the min level dof handler with the appropriate element
          // type
          if (this->simulation_parameters.linear_solver
                .at(PhysicsID::fluid_dynamics)
                .mg_use_fe_q_iso_q1 &&
              l == this->minlevel)
            {
              AssertThrow(
                mg_coarsening_type ==
                  Parameters::LinearSolver::MultigridCoarseningSequenceType::h,
                ExcNotImplemented());

              const auto points =
                QGaussLobatto<1>(levels[l].second + 1).get_points();

              this->dof_handlers[l].distribute_dofs(
                FESystem<dim>(FE_Q_iso_Q1<dim>(points), dim + 1));
              DoFRenumbering::Cuthill_McKee(this->dof_handlers[l]);
            }
          else
            {
              this->dof_handlers[l].distribute_dofs(
                FESystem<dim>(FE_Q<dim>(levels[l].second), dim + 1));
              DoFRenumbering::Cuthill_McKee(this->dof_handlers[l]);
            }
        }

      this->mg_setup_timer.leave_subsection(
        "Create DoFHandlers and distribute DoFs");

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << std::endl;
          this->pcout << "  -Levels of MG preconditioner:" << std::endl;
          for (unsigned int level = this->minlevel; level <= this->maxlevel;
               ++level)
            this->pcout << "    Level " << level << ": "
                        << this->dof_handlers[level].n_dofs() << " DoFs, "
                        << this->coarse_grid_triangulations[levels[level].first]
                             ->n_global_active_cells()
                        << " cells" << std::endl;
          this->pcout << std::endl;
        }

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity == Parameters::Verbosity::extra_verbose)
        {
          this->pcout << "  -MG vertical communication efficiency: "
                      << MGTools::vertical_communication_efficiency(
                           this->coarse_grid_triangulations)
                      << std::endl;

          this->pcout << "  -MG workload imbalance: "
                      << MGTools::workload_imbalance(
                           this->coarse_grid_triangulations)
                      << std::endl;
        }

      // Apply constraints and create mg operators for each level
      for (unsigned int level = this->minlevel; level <= this->maxlevel;
           ++level)
        {
          const auto &level_dof_handler = this->dof_handlers[level];
          auto       &level_constraint  = constraints[level];

          this->mg_setup_timer.enter_subsection("Set boundary conditions");

          level_constraint.clear();

          const IndexSet locally_relevant_dofs =
            DoFTools::extract_locally_relevant_dofs(level_dof_handler);

          level_constraint.reinit(level_dof_handler.locally_owned_dofs(),
                                  locally_relevant_dofs);

          DoFTools::make_hanging_node_constraints(level_dof_handler,
                                                  level_constraint);

          FEValuesExtractors::Vector velocities(0);
          FEValuesExtractors::Scalar pressure(dim);

          for (auto const &[id, type] :
               this->simulation_parameters.boundary_conditions.type)
            {
              if (type == BoundaryConditions::BoundaryType::slip)
                {
                  std::set<types::boundary_id> no_normal_flux_boundaries;
                  no_normal_flux_boundaries.insert(id);

                  VectorTools::compute_no_normal_flux_constraints(
                    level_dof_handler,
                    0,
                    no_normal_flux_boundaries,
                    level_constraint,
                    *mapping);
                }
              else if (type == BoundaryConditions::BoundaryType::periodic)
                {
                  DoFTools::make_periodicity_constraints(
                    level_dof_handler,
                    id,
                    this->simulation_parameters.boundary_conditions
                      .periodic_neighbor_id.at(id),
                    this->simulation_parameters.boundary_conditions
                      .periodic_direction.at(id),
                    level_constraint);
                }
              else if (type ==
                       BoundaryConditions::BoundaryType::periodic_neighbor)
                {
                  /*already taken into account in the previous else if*/
                }
              else if (type == BoundaryConditions::BoundaryType::pressure)
                {
                  Assert(
                    false,
                    ExcMessage(
                      "Pressure boundary conditions are not supported by the matrix free application."));
                }
              else if (type == BoundaryConditions::BoundaryType::function_weak)
                {
                  /*The function weak boundary condition is implemented in
                   * the operators*/
                }
              else if (type == BoundaryConditions::BoundaryType::partial_slip)
                {
                  Assert(
                    false,
                    ExcMessage(
                      "Partial slip boundary conditions are not supported by the matrix free application."));
                }
              else if (type == BoundaryConditions::BoundaryType::outlet)
                {
                  /*The directional do-nothing boundary condition is implemented
                   * in the operators*/
                }
              else if (type == BoundaryConditions::BoundaryType::noslip ||
                       type == BoundaryConditions::BoundaryType::function)
                {
                  VectorTools::interpolate_boundary_values(
                    *mapping,
                    level_dof_handler,
                    id,
                    dealii::Functions::ZeroFunction<dim, MGNumber>(dim + 1),
                    level_constraint,
                    fe->component_mask(velocities));
                }
              else if (type == BoundaryConditions::BoundaryType::none)
                {
                  /*Default boundary condition*/
                }
              else
                {
                  AssertThrow(
                    false,
                    ExcMessage(
                      "The boundary condition specified is not supported."));
                }
            }

          if (this->simulation_parameters.boundary_conditions
                .fix_pressure_constant &&
              level == this->minlevel)
            {
              types::global_dof_index min_index = numbers::invalid_unsigned_int;

              std::vector<types::global_dof_index> dof_indices;

              // Loop over the cells to identify the min index
              for (const auto &cell : level_dof_handler.active_cell_iterators())
                {
                  if (cell->is_locally_owned())
                    {
                      const auto &fe = cell->get_fe();

                      dof_indices.resize(fe.n_dofs_per_cell());
                      cell->get_dof_indices(dof_indices);

                      for (unsigned int i = 0; i < dof_indices.size(); ++i)
                        if (fe.system_to_component_index(i).first == dim)
                          min_index = std::min(min_index, dof_indices[i]);
                    }
                }

              // Necessary to find the min across all cores.
              min_index =
                Utilities::MPI::min(min_index,
                                    this->dof_handler.get_mpi_communicator());

              if (locally_relevant_dofs.is_element(min_index))
                level_constraint.add_line(min_index);
            }

          level_constraint.close();

          this->mg_setup_timer.leave_subsection("Set boundary conditions");

          this->mg_setup_timer.enter_subsection("Set up operators");

          // Provide appropriate quadrature depending on the type of elements of
          // the level
          Quadrature<dim> quadrature_mg = QGauss<dim>(levels[level].second + 1);
          if (this->simulation_parameters.linear_solver
                .at(PhysicsID::fluid_dynamics)
                .mg_use_fe_q_iso_q1 &&
              level == this->minlevel)
            {
              AssertThrow(
                mg_coarsening_type ==
                  Parameters::LinearSolver::MultigridCoarseningSequenceType::h,
                ExcNotImplemented());

              const auto points =
                QGaussLobatto<1>(level_dof_handler.get_fe().degree + 1)
                  .get_points();

              quadrature_mg = QIterated<dim>(QGauss<1>(2), points);
            }

          this->create_level_operator(level);

          this->mg_operators[level]->reinit(
            *mapping,
            level_dof_handler,
            level_constraint,
            quadrature_mg,
            forcing_function,
            physical_properties_manager,
            this->simulation_parameters.stabilization.stabilization,
            numbers::invalid_unsigned_int,
            simulation_control,
            this->simulation_parameters.boundary_conditions,
            this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .mg_enable_hessians_jacobian,
            true,
            this->simulation_parameters.mortar_parameters.enable);

          this->mg_setup_timer.leave_subsection("Set up operators");

          // If enabled, create mortar operators
          if (this->simulation_parameters.mortar_parameters.enable)
            {
              this->mg_setup_timer.enter_subsection("Set up mortar operators");

              std::cout << "\n constraints in level \n";
              level_constraint.print(std::cout);

              // Create mortar manager
              this->mg_operators[level]->mortar_manager_mf =
                std::make_shared<MortarManagerCircle<dim>>(
                  quadrature_mg,
                  *mapping,
                  level_dof_handler,
                  this->simulation_parameters.mortar_parameters);

              // Create mortar coupling evaluator
              this->mg_operators[level]->mortar_coupling_evaluator_mf =
                std::make_shared<NavierStokesCouplingEvaluation<dim, double>>(
                  *mapping,
                  level_dof_handler,
                  physical_properties_manager->get_rheology()
                    ->get_kinematic_viscosity());

              this->mg_operators[level]->mortar_coupling_operator_mf =
                std::make_shared<CouplingOperator<dim, double>>(
                  *mapping,
                  level_dof_handler,
                  level_constraint,
                  this->mg_operators[level]->mortar_coupling_evaluator_mf,
                  this->mg_operators[level]->mortar_manager_mf,
                  this->simulation_parameters.mortar_parameters
                    .rotor_boundary_id,
                  this->simulation_parameters.mortar_parameters
                    .stator_boundary_id,
                  this->simulation_parameters.mortar_parameters.sip_factor);

              this->mg_setup_timer.leave_subsection("Set up mortar operators");
            }
        }

      // Create transfer operators
      this->mg_setup_timer.enter_subsection("Create transfer operator");

      for (unsigned int level = this->minlevel; level < this->maxlevel; ++level)
        this->transfers[level + 1].reinit(this->dof_handlers[level + 1],
                                          this->dof_handlers[level],
                                          constraints[level + 1],
                                          constraints[level]);

      this->mg_transfer_gc = std::make_shared<GCTransferType>(
        this->transfers, [&](const auto l, auto &vec) {
          this->mg_operators[l]->initialize_dof_vector(vec);
        });

      this->mg_setup_timer.leave_subsection("Create transfer operator");
    }
  else
    AssertThrow(false, ExcNotImplemented());

  mg_smoother_preconditioners.resize(this->minlevel, this->maxlevel);
}

template <int dim>
void
MFNavierStokesPreconditionGMGBase<dim>::initialize()
{
  // Create smoother, fill parameters for each level and intialize it
  this->mg_setup_timer.enter_subsection("Set up and initialize smoother");

  this->mg_smoother = std::make_shared<
    MGSmootherPrecondition<OperatorType, SmootherType, MGVectorType>>();

  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
    this->minlevel, this->maxlevel);

  for (unsigned int level = this->minlevel; level <= this->maxlevel; ++level)
    {
      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_preconditioner_type ==
          Parameters::LinearSolver::MultigridSmootherPreconditionerType::
            InverseDiagonal)
        {
          MGVectorType diagonal_vector;
          this->mg_operators[level]->compute_inverse_diagonal(diagonal_vector);
          mg_smoother_preconditioners[level] =
            std::make_shared<MyDiagonalMatrix<MGVectorType>>(diagonal_vector);
        }
      else if (this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .mg_smoother_preconditioner_type ==
               Parameters::LinearSolver::MultigridSmootherPreconditionerType::
                 AdditiveSchwarzMethod)
        {
          if (mg_smoother_preconditioners[level] == nullptr)
            mg_smoother_preconditioners[level] =
              std::make_shared<PreconditionASM<MGVectorType>>();

          dynamic_cast<PreconditionASM<MGVectorType> *>(
            mg_smoother_preconditioners[level].get())
            ->initialize(this->mg_operators[level]
                           ->get_system_matrix_free()
                           .get_dof_handler(),
                         this->mg_operators[level]->get_system_matrix(),
                         this->mg_operators[level]->get_sparsity_pattern(),
                         this->mg_operators[level]
                           ->get_system_matrix_free()
                           .get_affine_constraints(),
                         this->mg_operators[level]->get_vector_partitioner(),
                         this->mg_operators[level]
                           ->get_system_matrix_free()
                           .get_mg_level());
        }

      smoother_data[level].preconditioner = mg_smoother_preconditioners[level];

      smoother_data[level].n_iterations =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_smoother_iterations;

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_eig_estimation)
        {
          // Set relaxation to zero so that eigenvalues are estimated
          // internally
          smoother_data[level].relaxation = 0.0;
          smoother_data[level].smoothing_range =
            this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .eig_estimation_smoothing_range;
          smoother_data[level].eig_cg_n_iterations =
            this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .eig_estimation_cg_n_iterations;
          smoother_data[level].eigenvalue_algorithm =
            SmootherType::AdditionalData::EigenvalueAlgorithm::power_iteration;
          smoother_data[level].constraints.copy_from(
            this->mg_operators[level]
              ->get_system_matrix_free()
              .get_affine_constraints());
        }
      else
        smoother_data[level].relaxation =
          this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_relaxation;
    }

  mg_smoother->initialize(this->mg_operators, smoother_data);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_smoother_eig_estimation &&
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .eig_estimation_verbose != Parameters::Verbosity::quiet)
    {
      // Print eigenvalue estimation for all levels
      for (unsigned int level = this->minlevel; level <= this->maxlevel;
           ++level)
        {
          MGVectorType vec;
          this->mg_operators[level]->initialize_dof_vector(vec);
          const auto evs =
            mg_smoother->smoothers[level].estimate_eigenvalues(vec);

          this->pcout << std::endl;
          this->pcout << "  -Eigenvalue estimation level "
                      << level - this->minlevel << ":" << std::endl;
          this->pcout << "    Relaxation parameter: "
                      << mg_smoother->smoothers[level].get_relaxation()
                      << std::endl;
          this->pcout << "    Minimum eigenvalue: "
                      << evs.min_eigenvalue_estimate << std::endl;
          this->pcout << "    Maximum eigenvalue: "
                      << evs.max_eigenvalue_estimate << std::endl;
          this->pcout << std::endl;
        }
    }

  this->mg_setup_timer.leave_subsection("Set up and initialize smoother");

  // Create coarse-grid GMRES solver and AMG preconditioner
  this->mg_setup_timer.enter_subsection("Create coarse-grid solver");

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_coarse_grid_solver ==
      Parameters::LinearSolver::CoarseGridSolverType::gmres)
    {
      const int max_iterations =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_gmres_max_iterations;
      const double tolerance =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_gmres_tolerance;
      const double reduce =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_gmres_reduce;
      this->coarse_grid_solver_control = std::make_shared<ReductionControl>(
        max_iterations, tolerance, reduce, false, false);
      SolverGMRES<TrilinosVectorType>::AdditionalData solver_parameters;
      solver_parameters.max_n_tmp_vectors =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_gmres_max_krylov_vectors;

      this->coarse_grid_solver =
        std::make_shared<SolverGMRES<TrilinosVectorType>>(
          *this->coarse_grid_solver_control, solver_parameters);

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_gmres_preconditioner ==
          Parameters::LinearSolver::PreconditionerType::amg)
        {
          setup_AMG();

          this->mg_coarse = std::make_shared<MGCoarseGridIterativeSolver<
            MGVectorType,
            SolverGMRES<TrilinosVectorType>,
            TrilinosWrappers::SparseMatrix,
            PreconditionAdapter<MGVectorType, TrilinosVectorType>>>(
            *this->coarse_grid_solver,
            this->mg_operators[this->minlevel]->get_system_matrix(),
            *this->coarse_grid_precondition);
        }
      else if (this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .mg_gmres_preconditioner ==
               Parameters::LinearSolver::PreconditionerType::ilu)
        {
          setup_ILU();

          this->mg_coarse = std::make_shared<MGCoarseGridIterativeSolver<
            MGVectorType,
            SolverGMRES<TrilinosVectorType>,
            TrilinosWrappers::SparseMatrix,
            PreconditionAdapter<MGVectorType, TrilinosVectorType>>>(
            *this->coarse_grid_solver,
            this->mg_operators[this->minlevel]->get_system_matrix(),
            *this->coarse_grid_precondition);
        }
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_solver ==
           Parameters::LinearSolver::CoarseGridSolverType::amg)
    {
      setup_AMG();

      this->mg_coarse = std::make_shared<CoarseGridSolverApply>(
        *this->coarse_grid_precondition);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_solver ==
           Parameters::LinearSolver::CoarseGridSolverType::ilu)
    {
      setup_ILU();

      this->mg_coarse = std::make_shared<CoarseGridSolverApply>(
        *this->coarse_grid_precondition);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_solver ==
           Parameters::LinearSolver::CoarseGridSolverType::direct)
    {
      TrilinosWrappers::SolverDirect::AdditionalData data;
      this->direct_solver_control =
        std::make_shared<SolverControl>(100, 1.e-10);

      auto precondition_direct =
        std::make_shared<TrilinosWrappers::SolverDirect>(
          *this->direct_solver_control, data);

      precondition_direct->initialize(
        this->mg_operators[this->minlevel]->get_system_matrix());

      coarse_grid_precondition =
        std::make_shared<PreconditionAdapter<MGVectorType, TrilinosVectorType>>(
          precondition_direct);

      this->mg_coarse = std::make_shared<CoarseGridSolverApply>(
        *this->coarse_grid_precondition);
    }
  else
    AssertThrow(
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .mg_coarse_grid_solver ==
          Parameters::LinearSolver::CoarseGridSolverType::gmres ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .mg_coarse_grid_solver ==
          Parameters::LinearSolver::CoarseGridSolverType::amg ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .mg_coarse_grid_solver ==
          Parameters::LinearSolver::CoarseGridSolverType::ilu ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .mg_coarse_grid_solver ==
          Parameters::LinearSolver::CoarseGridSolverType::direct,
      ExcMessage(
        "This coarse-grid solver is not supported. Supported options are <gmres|amg|ilu|direct>."));

  this->mg_setup_timer.leave_subsection("Create coarse-grid solver");

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      // Create interface matrices needed for local smoothing in case of
      // local refinement
      this->mg_interface_matrix_in =
        std::make_shared<mg::Matrix<MGVectorType>>(this->ls_mg_interface_in);


      if (this->minlevel != this->intlevel)
        {
          // Create main MG object
          this->mg_intermediate =
            std::make_shared<Multigrid<MGVectorType>>(*this->mg_matrix,
                                                      *this->mg_coarse,
                                                      *this->mg_transfer_ls,
                                                      *this->mg_smoother,
                                                      *this->mg_smoother,
                                                      this->minlevel,
                                                      this->intlevel);

          if (this->dof_handler.get_triangulation().has_hanging_nodes())
            this->mg_intermediate->set_edge_in_matrix(
              *this->mg_interface_matrix_in);

          // Create MG preconditioner
          this->ls_multigrid_preconditioner_intermediate =
            std::make_shared<PreconditionMG<dim, MGVectorType, LSTransferType>>(
              this->dof_handler, *this->mg_intermediate, *this->mg_transfer_ls);

          const int max_iterations = this->simulation_parameters.linear_solver
                                       .at(PhysicsID::fluid_dynamics)
                                       .mg_gmres_max_iterations;
          const double tolerance = this->simulation_parameters.linear_solver
                                     .at(PhysicsID::fluid_dynamics)
                                     .mg_gmres_tolerance;
          const double reduce = this->simulation_parameters.linear_solver
                                  .at(PhysicsID::fluid_dynamics)
                                  .mg_gmres_reduce;

          this->coarse_grid_solver_control_intermediate =
            std::make_shared<ReductionControl>(
              max_iterations, tolerance, reduce, false, false);

          this->coarse_grid_solver_intermediate =
            std::make_shared<SolverGMRES<MGVectorType>>(
              *this->coarse_grid_solver_control_intermediate);

          this->mg_coarse_intermediate =
            std::make_shared<MGCoarseGridIterativeSolver<
              MGVectorType,
              SolverGMRES<MGVectorType>,
              OperatorType,
              PreconditionMG<dim, MGVectorType, LSTransferType>>>(
              *this->coarse_grid_solver_intermediate,
              *this->mg_operators[this->intlevel],
              *this->ls_multigrid_preconditioner_intermediate);
        }

      // Create main MG object
      this->mg = std::make_shared<Multigrid<MGVectorType>>(
        *this->mg_matrix,
        (this->minlevel != this->intlevel) ? (*this->mg_coarse_intermediate) :
                                             (*this->mg_coarse),
        *this->mg_transfer_ls,
        *this->mg_smoother,
        *this->mg_smoother,
        this->intlevel,
        this->maxlevel);

      if (this->dof_handler.get_triangulation().has_hanging_nodes())
        this->mg->set_edge_in_matrix(*this->mg_interface_matrix_in);

      // Create MG preconditioner
      this->ls_multigrid_preconditioner =
        std::make_shared<PreconditionMG<dim, MGVectorType, LSTransferType>>(
          this->dof_handler, *this->mg, *this->mg_transfer_ls);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      if (this->minlevel != this->intlevel)
        {
          // Create main MG object
          this->mg_intermediate =
            std::make_shared<Multigrid<MGVectorType>>(*this->mg_matrix,
                                                      *this->mg_coarse,
                                                      *this->mg_transfer_gc,
                                                      *this->mg_smoother,
                                                      *this->mg_smoother,
                                                      this->minlevel,
                                                      this->intlevel);

          // Create MG preconditioner
          this->gc_multigrid_preconditioner_intermediate =
            std::make_shared<PreconditionMG<dim, MGVectorType, GCTransferType>>(
              this->dof_handler, *this->mg_intermediate, *this->mg_transfer_gc);

          const int max_iterations = this->simulation_parameters.linear_solver
                                       .at(PhysicsID::fluid_dynamics)
                                       .mg_gmres_max_iterations;
          const double tolerance = this->simulation_parameters.linear_solver
                                     .at(PhysicsID::fluid_dynamics)
                                     .mg_gmres_tolerance;
          const double reduce = this->simulation_parameters.linear_solver
                                  .at(PhysicsID::fluid_dynamics)
                                  .mg_gmres_reduce;

          this->coarse_grid_solver_control_intermediate =
            std::make_shared<ReductionControl>(
              max_iterations, tolerance, reduce, false, false);

          this->coarse_grid_solver_intermediate =
            std::make_shared<SolverGMRES<MGVectorType>>(
              *this->coarse_grid_solver_control_intermediate);

          this->mg_coarse_intermediate =
            std::make_shared<MGCoarseGridIterativeSolver<
              MGVectorType,
              SolverGMRES<MGVectorType>,
              OperatorType,
              PreconditionMG<dim, MGVectorType, GCTransferType>>>(
              *this->coarse_grid_solver_intermediate,
              *this->mg_operators[this->intlevel],
              *this->gc_multigrid_preconditioner_intermediate);
        }

      // Create main MG object
      this->mg = std::make_shared<Multigrid<MGVectorType>>(
        *this->mg_matrix,
        (this->minlevel != this->intlevel) ? (*this->mg_coarse_intermediate) :
                                             (*this->mg_coarse),
        *this->mg_transfer_gc,
        *this->mg_smoother,
        *this->mg_smoother,
        this->intlevel,
        this->maxlevel,
        Multigrid<MGVectorType>::Cycle::v_cycle);

      // Create MG preconditioner
      this->gc_multigrid_preconditioner =
        std::make_shared<PreconditionMG<dim, MGVectorType, GCTransferType>>(
          this->dof_handler, *this->mg, *this->mg_transfer_gc);
    }

  // Print detailed timings of multigrid vmult
  const auto create_mg_timer_function = [&](const std::string &label) {
    return [label, this](const bool flag, const unsigned int level) {
      const std::string label_full =
        (label == "") ?
          ("gmg::vmult::level_" + std::to_string(level)) :
          ("gmg::vmult::level_" + std::to_string(level) + "::" + label);

      if (flag)
        this->mg_vmult_timer.enter_subsection(label_full);
      else
        this->mg_vmult_timer.leave_subsection(label_full);
    };
  };

  this->mg->connect_pre_smoother_step(
    create_mg_timer_function("0_pre_smoother_step"));
  this->mg->connect_residual_step(create_mg_timer_function("1_residual_step"));
  this->mg->connect_restriction(create_mg_timer_function("2_restriction"));
  this->mg->connect_coarse_solve(create_mg_timer_function(""));
  this->mg->connect_prolongation(create_mg_timer_function("3_prolongation"));
  this->mg->connect_edge_prolongation(
    create_mg_timer_function("4_edge_prolongation"));
  this->mg->connect_post_smoother_step(
    create_mg_timer_function("5_post_smoother_step"));

  if (mg_intermediate)
    {
      this->mg_intermediate->connect_pre_smoother_step(
        create_mg_timer_function("0_pre_smoother_step"));
      this->mg_intermediate->connect_residual_step(
        create_mg_timer_function("1_residual_step"));
      this->mg_intermediate->connect_restriction(
        create_mg_timer_function("2_restriction"));
      this->mg_intermediate->connect_coarse_solve(create_mg_timer_function(""));
      this->mg_intermediate->connect_prolongation(
        create_mg_timer_function("3_prolongation"));
      this->mg_intermediate->connect_edge_prolongation(
        create_mg_timer_function("4_edge_prolongation"));
      this->mg_intermediate->connect_post_smoother_step(
        create_mg_timer_function("5_post_smoother_step"));
    }

  const auto create_mg_precon_timer_function = [&](const std::string &label) {
    return [label, this](const bool flag) {
      const std::string label_full = "gmg::vmult::" + label;

      if (flag)
        this->mg_vmult_timer.enter_subsection(label_full);
      else
        this->mg_vmult_timer.leave_subsection(label_full);
    };
  };

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      this->ls_multigrid_preconditioner->connect_transfer_to_mg(
        create_mg_precon_timer_function("transfer_to_mg"));
      this->ls_multigrid_preconditioner->connect_transfer_to_global(
        create_mg_precon_timer_function("transfer_to_global"));
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      this->gc_multigrid_preconditioner->connect_transfer_to_mg(
        create_mg_precon_timer_function("transfer_to_mg"));
      this->gc_multigrid_preconditioner->connect_transfer_to_global(
        create_mg_precon_timer_function("transfer_to_global"));
    }
}

template <int dim>
void
MFNavierStokesPreconditionGMGBase<dim>::vmult(VectorType       &dst,
                                              const VectorType &src) const
{
  if (this->ls_multigrid_preconditioner)
    this->ls_multigrid_preconditioner->vmult(dst, src);
  else if (this->gc_multigrid_preconditioner)
    this->gc_multigrid_preconditioner->vmult(dst, src);
  else
    AssertThrow(false, ExcNotImplemented());

  // Save number of coarse grid iterations needed in one vmult
  if (this->coarse_grid_solver_control)
    this->coarse_grid_iterations.emplace_back(
      this->coarse_grid_solver_control->last_step());
  if (this->coarse_grid_solver_control_intermediate)
    this->coarse_grid_iterations.emplace_back(
      this->coarse_grid_solver_control_intermediate->last_step());
}

template <int dim>
void
MFNavierStokesPreconditionGMGBase<dim>::print_relevant_info() const
{
  if (this->coarse_grid_solver_control ||
      this->coarse_grid_solver_control_intermediate)
    {
      if (this->coarse_grid_iterations.empty())
        this->pcout << "  -Coarse grid solver took: 0 iterations" << std::endl;
      else
        {
          unsigned int total = this->coarse_grid_iterations[0];
          this->pcout << "  -Coarse grid solver took: "
                      << this->coarse_grid_iterations[0];
          for (unsigned int i = 1; i < this->coarse_grid_iterations.size(); i++)
            {
              this->pcout << " + " << this->coarse_grid_iterations[i];
              total += this->coarse_grid_iterations[i];
            }
          this->pcout << " = " << total << " iterations" << std::endl;

          this->coarse_grid_iterations.clear();
        }
    }
}

template <int dim>
const MGLevelObject<std::shared_ptr<NavierStokesOperatorBase<
  dim,
  typename MFNavierStokesPreconditionGMGBase<dim>::MGNumber>>> &
MFNavierStokesPreconditionGMGBase<dim>::get_mg_operators() const
{
  return this->mg_operators;
}

template <int dim>
const MGLevelObject<std::shared_ptr<PreconditionBase<
  typename MFNavierStokesPreconditionGMGBase<dim>::MGVectorType>>> &
MFNavierStokesPreconditionGMGBase<dim>::get_mg_smoother_preconditioners() const
{
  return this->mg_smoother_preconditioners;
}

template <int dim>
void
MFNavierStokesPreconditionGMGBase<dim>::setup_AMG()
{
  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;

  // Extract matrix of the minlevel to avoid building it twice
  const TrilinosWrappers::SparseMatrix &min_level_matrix =
    this->mg_operators[this->minlevel]->get_system_matrix();

  if (!this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
         .mg_amg_use_default_parameters)
    {
      amg_data.elliptic = false;
      if (this->dof_handler.get_fe().degree > 1)
        amg_data.higher_order_elements = true;
      amg_data.n_cycles =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_n_cycles;
      amg_data.w_cycle =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_w_cycles;
      amg_data.aggregation_threshold =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_aggregation_threshold;
      amg_data.smoother_sweeps =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_smoother_sweeps;
      amg_data.smoother_overlap =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_smoother_overlap;
      amg_data.output_details = false;
      amg_data.smoother_type  = "ILU";
      amg_data.coarse_type    = "ILU";

      std::vector<std::vector<bool>> constant_modes;
      ComponentMask                  components(dim + 1, true);
      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::lsmg)

        {
          // Constant modes for velocity and pressure
          constant_modes = DoFTools::extract_level_constant_modes(
            this->mg_operators[this->minlevel]
              ->get_system_matrix_free()
              .get_mg_level(),
            this->mg_operators[this->minlevel]
              ->get_system_matrix_free()
              .get_dof_handler(),
            components);
        }
      else if (this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .preconditioner ==
               Parameters::LinearSolver::PreconditionerType::gcmg)
        {
          // Constant modes for velocity and pressure
          constant_modes =
            DoFTools::extract_constant_modes(this->dof_handlers[this->minlevel],
                                             components);
        }

      amg_data.constant_modes = constant_modes;

      Teuchos::ParameterList              parameter_ml;
      std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
      amg_data.set_parameters(parameter_ml,
                              distributed_constant_modes,
                              min_level_matrix);

      const double ilu_fill =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_fill;
      const double ilu_atol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_atol;
      const double ilu_rtol =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .amg_precond_ilu_rtol;
      parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
      parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
      parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

      parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
      parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
      parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);

      auto precondition_amg =
        std::make_shared<TrilinosWrappers::PreconditionAMG>();

      precondition_amg->initialize(min_level_matrix, parameter_ml);

      coarse_grid_precondition =
        std::make_shared<PreconditionAdapter<MGVectorType, TrilinosVectorType>>(
          precondition_amg);
    }
  else
    {
      auto precondition_amg =
        std::make_shared<TrilinosWrappers::PreconditionAMG>();

      coarse_grid_precondition =
        std::make_shared<PreconditionAdapter<MGVectorType, TrilinosVectorType>>(
          precondition_amg);
    }
}

template <int dim>
void
MFNavierStokesPreconditionGMGBase<dim>::setup_ILU()
{
  int current_preconditioner_fill_level =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

  auto precondition_ilu = std::make_shared<TrilinosWrappers::PreconditionILU>();

  precondition_ilu->initialize(
    this->mg_operators[this->minlevel]->get_system_matrix(),
    preconditionerOptions);

  coarse_grid_precondition =
    std::make_shared<PreconditionAdapter<MGVectorType, TrilinosVectorType>>(
      precondition_ilu);
}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::create_level_operator(
  const unsigned int level)
{
  if (this->simulation_parameters.physical_properties_manager
        .is_non_newtonian())
    this->mg_operators[level] = std::make_shared<
      NavierStokesNonNewtonianStabilizedOperator<dim, MGNumber>>();
  else
    this->mg_operators[level] =
      std::make_shared<NavierStokesStabilizedOperator<dim, MGNumber>>();
}

template <int dim>
MFNavierStokesPreconditionGMG<dim>::MFNavierStokesPreconditionGMG(
  const SimulationParameters<dim> &simulation_parameters,
  const DoFHandler<dim>           &dof_handler,
  const DoFHandler<dim>           &dof_handler_fe_q_iso_q1)
  : MFNavierStokesPreconditionGMGBase<dim>(simulation_parameters,
                                           dof_handler,
                                           dof_handler_fe_q_iso_q1)
{}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::initialize(
  const std::shared_ptr<SimulationControl> &simulation_control,
  FlowControl<dim>                         &flow_control,
  const VectorType                         &present_solution,
  const VectorType                         &time_derivative_previous_solutions)
{
  // Local objects for the different levels
  MGLevelObject<MGVectorType> mg_solution(this->minlevel, this->maxlevel);
  MGLevelObject<MGVectorType> mg_time_derivative_previous_solutions(
    this->minlevel, this->maxlevel);

  for (unsigned int level = this->minlevel; level <= this->maxlevel; ++level)
    {
      this->mg_operators[level]->initialize_dof_vector(mg_solution[level]);
      if (is_bdf(simulation_control->get_assembly_method()))
        this->mg_operators[level]->initialize_dof_vector(
          mg_time_derivative_previous_solutions[level]);
    }

  this->mg_setup_timer.enter_subsection("Execute relevant transfers");

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      // Create transfer operator and transfer solution to mg levels

      this->mg_transfer_ls->interpolate_to_mg(this->dof_handler,
                                              mg_solution,
                                              present_solution);

      if (is_bdf(simulation_control->get_assembly_method()))
        this->mg_transfer_ls->interpolate_to_mg(
          this->dof_handler,
          mg_time_derivative_previous_solutions,
          time_derivative_previous_solutions);

      this->mg_matrix =
        std::make_shared<mg::Matrix<MGVectorType>>(this->ls_mg_operators);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      this->mg_transfer_gc->interpolate_to_mg(this->dof_handler,
                                              mg_solution,
                                              present_solution);

      if (is_bdf(simulation_control->get_assembly_method()))
        this->mg_transfer_gc->interpolate_to_mg(
          this->dof_handler,
          mg_time_derivative_previous_solutions,
          time_derivative_previous_solutions);

      this->mg_matrix =
        std::make_shared<mg::Matrix<MGVectorType>>(this->mg_operators);
    }

  this->mg_setup_timer.leave_subsection("Execute relevant transfers");

  // Evaluate non linear terms for all mg operators
  for (unsigned int level = this->minlevel; level <= this->maxlevel; ++level)
    {
      mg_solution[level].update_ghost_values();
      this->mg_operators[level]->evaluate_non_linear_term_and_calculate_tau(
        mg_solution[level]);

      if (is_bdf(simulation_control->get_assembly_method()))
        {
          mg_time_derivative_previous_solutions[level].update_ghost_values();
          this->mg_operators[level]
            ->evaluate_time_derivative_previous_solutions(
              mg_time_derivative_previous_solutions[level]);

          if (this->simulation_parameters.flow_control.enable_flow_control)
            this->mg_operators[level]->update_beta_force(
              flow_control.get_beta());
        }
    }

  MFNavierStokesPreconditionGMGBase<dim>::initialize();
}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::initialize_auxiliary_physics(
  const DoFHandler<dim> &temperature_dof_handler,
  const VectorType      &temperature_present_solution)
{
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      AssertThrow(false, ExcNotImplemented());
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      const unsigned int min_level = this->minlevel;
      const unsigned int max_level = this->maxlevel;

      this->temperature_dof_handlers.resize(min_level, max_level);

      for (unsigned int l = min_level; l <= max_level; l++)
        {
          this->temperature_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->temperature_dof_handlers[l].distribute_dofs(
            temperature_dof_handler.get_fe());
          DoFRenumbering::Cuthill_McKee(this->temperature_dof_handlers[l]);
        }

      this->transfers_temperature.resize(min_level, max_level);

      for (unsigned int l = min_level; l < max_level; l++)
        {
          this->transfers_temperature[l + 1].reinit(
            this->temperature_dof_handlers[l + 1],
            this->temperature_dof_handlers[l],
            {},
            {});
        }

      this->mg_transfer_gc_temperature =
        std::make_shared<GCTransferType>(this->transfers_temperature);

#if DEAL_II_VERSION_GTE(9, 7, 0)
      this->mg_transfer_gc_temperature->build(
        temperature_dof_handler, [&](const auto l, auto &vec) {
          vec.reinit(this->temperature_dof_handlers[l].locally_owned_dofs(),
                     DoFTools::extract_locally_active_dofs(
                       this->temperature_dof_handlers[l]),
                     this->temperature_dof_handlers[l].get_mpi_communicator());
        });
#endif

      MGLevelObject<MGVectorType> mg_temperature_solution(this->minlevel,
                                                          this->maxlevel);

      this->mg_transfer_gc_temperature->interpolate_to_mg(
        temperature_dof_handler,
        mg_temperature_solution,
        temperature_present_solution);

      for (unsigned int l = min_level; l <= max_level; l++)
        {
          mg_temperature_solution[l].update_ghost_values();

          this->mg_operators[l]->compute_buoyancy_term(
            mg_temperature_solution[l], this->temperature_dof_handlers[l]);
        }
    }
}

template <int dim>
FluidDynamicsMatrixFree<dim>::FluidDynamicsMatrixFree(
  SimulationParameters<dim> &nsparam)
  : NavierStokesBase<dim, VectorType, IndexSet>(nsparam)
{
  AssertThrow(
    nsparam.fem_parameters.velocity_order ==
      nsparam.fem_parameters.pressure_order,
    dealii::ExcMessage(
      "Matrix free Navier-Stokes does not support different orders for the velocity and the pressure!"));

  this->fe = std::make_shared<FESystem<dim>>(
    FE_Q<dim>(nsparam.fem_parameters.velocity_order), dim + 1);

  physical_properties_manager = std::make_shared<PhysicalPropertiesManager>(
    this->simulation_parameters.physical_properties_manager);

  if (this->physical_properties_manager->is_non_newtonian())
    {
      system_operator = std::make_shared<
        NavierStokesNonNewtonianStabilizedOperator<dim, double>>();
      AssertThrow(
        this->simulation_parameters.stabilization.stabilization ==
          Parameters::Stabilization::NavierStokesStabilization::pspg_supg,
        dealii::ExcMessage(
          "Matrix free Non-Newtonian Navier-Stokes only supports SUPG/PSPG stabilization."));
    }
  else
    system_operator =
      std::make_shared<NavierStokesStabilizedOperator<dim, double>>();

  if (!this->simulation_parameters.source_term.enable)
    {
      this->forcing_function.reset();
    }
}

template <int dim>
FluidDynamicsMatrixFree<dim>::~FluidDynamicsMatrixFree()
{
  this->dof_handler.clear();
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::solve()
{
  this->computing_timer.enter_subsection("Read mesh and manifolds");

  if (this->simulation_parameters.mortar_parameters.enable)
    {
      read_mesh_and_manifolds_for_stator_and_rotor(
        *this->triangulation,
        this->simulation_parameters.mesh,
        this->simulation_parameters.manifolds_parameters,
        this->simulation_parameters.restart_parameters.restart,
        this->simulation_parameters.boundary_conditions,
        this->simulation_parameters.mortar_parameters);
      // Create and initialize mapping cache
      this->mapping_cache =
        std::make_shared<MappingQCache<dim>>(this->velocity_fem_degree);
      this->rotate_rotor_mapping(true);
    }
  else
    read_mesh_and_manifolds(
      *this->triangulation,
      this->simulation_parameters.mesh,
      this->simulation_parameters.manifolds_parameters,
      this->simulation_parameters.restart_parameters.restart,
      this->simulation_parameters.boundary_conditions);

  this->computing_timer.leave_subsection("Read mesh and manifolds");

  this->setup_dofs();
  this->box_refine_mesh(this->simulation_parameters.restart_parameters.restart);
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

  // Only needed if other physics apart from fluid dynamics are enabled.
  if (this->multiphysics->get_active_physics().size() > 1)
    this->update_multiphysics_time_average_solution();

  while (this->simulation_control->integrate())
    {
      if (this->forcing_function)
        this->forcing_function->set_time(
          this->simulation_control->get_current_time());

      this->update_boundary_conditions();
      this->multiphysics->update_boundary_conditions();

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();
      this->update_mortar_configuration();

      if (!this->simulation_control->is_at_start())
        {
          NavierStokesBase<dim, VectorType, IndexSet>::refine_mesh();
        }

      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->computing_timer.enter_subsection(
            "Calculate time derivative previous solutions");

          calculate_time_derivative_previous_solutions();
          this->time_derivative_previous_solutions.update_ghost_values();
          this->system_operator->evaluate_time_derivative_previous_solutions(
            this->time_derivative_previous_solutions);

          this->computing_timer.leave_subsection(
            "Calculate time derivative previous solutions");

          if (this->simulation_parameters.flow_control.enable_flow_control)
            this->system_operator->update_beta_force(
              this->flow_control.get_beta());
        }

      if (this->multiphysics->get_active_physics().size() > 1)
        {
          update_solutions_for_fluid_dynamics();

          if (this->simulation_parameters.multiphysics.buoyancy_force)
            this->system_operator->compute_buoyancy_term(
              temperature_present_solution,
              *this->multiphysics->get_dof_handler(PhysicsID::heat_transfer));
        }

      this->iterate();
      this->postprocess(false);
      this->finish_time_step();

      if (this->simulation_parameters.timer.type ==
          Parameters::Timer::Type::iteration)
        print_mg_setup_times();
    }

  if (this->simulation_parameters.timer.type == Parameters::Timer::Type::end)
    print_mg_setup_times();

  this->finish_simulation();
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::setup_dofs_fd()
{
  TimerOutput::Scope t(this->computing_timer, "Setup DoFs");

  // Clear the preconditioners
  ilu_preconditioner.reset();
  gmg_preconditioner.reset();

  // Clear matrix free operator
  this->system_operator->clear();

  // Fill the dof handler and initialize vectors
  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      this->dof_handler.distribute_mg_dofs();

      // To use elements with linear interpolation for coarse-grid we need to
      // have another dof handler with the appropriate element type
      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_use_fe_q_iso_q1)
        {
          this->dof_handler_fe_q_iso_q1.reinit(*this->triangulation);

          const auto points =
            QGaussLobatto<1>(this->dof_handler.get_fe().degree + 1)
              .get_points();

          this->dof_handler_fe_q_iso_q1.distribute_dofs(
            FESystem<dim>(FE_Q_iso_Q1<dim>(points), dim + 1));
          this->dof_handler_fe_q_iso_q1.distribute_mg_dofs();
          DoFRenumbering::Cuthill_McKee(this->dof_handler_fe_q_iso_q1);
        }
    }

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  this->locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(this->dof_handler);

  // If enabled, rotate rotor mapping
  this->rotate_rotor_mapping(false);

  // Non-zero constraints
  this->define_non_zero_constraints();

  // Check whether the boundary conditions specified in the parameter file are
  // available for this solver

  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::pressure ||
          type == BoundaryConditions::BoundaryType::partial_slip)
        {
          Assert(
            false,
            ExcMessage(
              "The following boundary conditions are not supported by the lethe-fluid-matrix-free application: pressure and partial slip."));
        }
    }

  // Zero constraints
  this->define_zero_constraints();

  // If enabled, create mortar operators
  this->reinit_mortar_operators_mf();

  // Initialize matrix-free object
  unsigned int mg_level = numbers::invalid_unsigned_int;
  this->system_operator->reinit(
    *this->get_mapping(),
    this->dof_handler,
    this->zero_constraints,
    *this->cell_quadrature,
    this->forcing_function,
    this->physical_properties_manager,
    this->simulation_parameters.stabilization.stabilization,
    mg_level,
    this->simulation_control,
    this->simulation_parameters.boundary_conditions,
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .enable_hessians_jacobian,
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .enable_hessians_residual,
    this->simulation_parameters.mortar_parameters.enable);


  // Initialize vectors using operator
  this->system_operator->initialize_dof_vector(this->present_solution);
  this->system_operator->initialize_dof_vector(this->evaluation_point);
  this->system_operator->initialize_dof_vector(this->newton_update);
  this->system_operator->initialize_dof_vector(this->system_rhs);
  this->system_operator->initialize_dof_vector(this->local_evaluation_point);
  this->system_operator->initialize_dof_vector(
    this->time_derivative_previous_solutions);

  // Initialize vectors of previous solutions
  for (auto &solution : this->previous_solutions)
    {
      this->system_operator->initialize_dof_vector(solution);
    }

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->average_velocities->initialize_vectors(
        this->locally_owned_dofs,
        this->locally_relevant_dofs,
        this->fe->n_dofs_per_vertex(),
        this->mpi_communicator);

      // Intialize Trilinos vectors used to pass average velocities to
      // multiphysics interface
      this->multiphysics_average_velocities.reinit(this->locally_owned_dofs,
                                                   this->locally_relevant_dofs,
                                                   this->mpi_communicator);
    }

  double global_volume =
    GridTools::volume(*this->triangulation, *this->get_mapping());

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;

  // Initialize Trilinos vectors used to pass solution to multiphysics
  // interface
  this->multiphysics_present_solution.reinit(this->locally_owned_dofs,
                                             this->locally_relevant_dofs,
                                             this->mpi_communicator);

  // Pre-allocate memory for the previous solutions using the information
  // of the BDF schemes
  this->multiphysics_previous_solutions.resize(
    this->simulation_control->get_number_of_previous_solution_in_assembly());

  for (auto &solution : this->multiphysics_previous_solutions)
    solution.reinit(this->locally_owned_dofs,
                    this->locally_relevant_dofs,
                    this->mpi_communicator);

  // Provide relevant solution to multiphysics interface only if other physics
  // apart from fluid dynamics are enabled.
  if (this->multiphysics->get_active_physics().size() > 1)
    update_solutions_for_multiphysics();
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::update_mortar_configuration()
{
  if (!this->simulation_parameters.mortar_parameters.enable)
    return;

  TimerOutput::Scope t(this->computing_timer, "Update mortar configuration");

  bool refinement_step;
  if (this->simulation_parameters.mesh_adaptation.refinement_at_frequency)
    refinement_step = this->simulation_control->get_step_number() %
                        this->simulation_parameters.mesh_adaptation.frequency ==
                      0;
  else
    refinement_step = this->simulation_control->get_step_number() == 0;

  // We need to update the mortar operator/evaluator, as well as the sparsity
  // pattern, at every iteration. Since this is already done within
  // setup_dofs(), which is called in refine_mesh(), here we make sure that,
  // when there is no mesh refinement, the mortar information is still updated
  if (this->simulation_control->is_at_start() || !refinement_step ||
      this->simulation_parameters.mesh_adaptation.type ==
        Parameters::MeshAdaptation::Type::none)
    {
      // Clear the preconditioner before the matrix they are associated with is
      // cleared
      gmg_preconditioner.reset();
      ilu_preconditioner.reset();

      // Rotate mapping
      this->rotate_rotor_mapping(false);

      // Non Zero constraints
      this->define_non_zero_constraints();

      // Zero constraints
      this->define_zero_constraints();

      // Create mortar manager, operator, and evaluator
      this->reinit_mortar_operators_mf();
    }
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::reinit_mortar_operators_mf()
{
  if (!this->simulation_parameters.mortar_parameters.enable)
    return;

  TimerOutput::Scope t(this->computing_timer, "Reinit mortar operators");

  // Create mortar manager
  this->system_operator->mortar_manager_mf =
    std::make_shared<MortarManagerCircle<dim>>(
      *this->cell_quadrature,
      *this->get_mapping(),
      this->dof_handler,
      this->simulation_parameters.mortar_parameters);

  // Create mortar coupling evaluator
  this->system_operator->mortar_coupling_evaluator_mf =
    std::make_shared<NavierStokesCouplingEvaluation<dim, double>>(
      *this->get_mapping(),
      this->dof_handler,
      this->physical_properties_manager->get_rheology()
        ->get_kinematic_viscosity());

  this->system_operator->mortar_coupling_operator_mf =
    std::make_shared<CouplingOperator<dim, double>>(
      *this->get_mapping(),
      this->dof_handler,
      this->zero_constraints,
      this->system_operator->mortar_coupling_evaluator_mf,
      this->system_operator->mortar_manager_mf,
      this->simulation_parameters.mortar_parameters.rotor_boundary_id,
      this->simulation_parameters.mortar_parameters.stator_boundary_id,
      this->simulation_parameters.mortar_parameters.sip_factor);
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::set_initial_condition_fd(
  Parameters::FluidDynamicsInitialConditionType initial_condition_type,
  bool                                          restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type ==
           Parameters::FluidDynamicsInitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->present_solution.update_ghost_values();
      this->finish_time_step();
    }
  else if (initial_condition_type ==
           Parameters::FluidDynamicsInitialConditionType::viscous)
    {
      // Set the nodal values to have an initial condition that is adequate
      this->set_nodal_values();
      this->present_solution.update_ghost_values();

      // Get viscosity model
      std::shared_ptr<RheologicalModel> original_viscosity_model =
        this->physical_properties_manager->get_rheology();

      // Temporarily set the rheology to be newtonian with predefined
      // viscosity
      std::shared_ptr<Newtonian> temporary_rheology =
        std::make_shared<Newtonian>(
          this->simulation_parameters.initial_condition->kinematic_viscosity);

      this->physical_properties_manager->set_rheology(temporary_rheology);

      // Solve the problem with the temporary viscosity
      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::steady);
      PhysicsSolver<LinearAlgebra::distributed::Vector<double>>::
        solve_non_linear_system(false);
      this->finish_time_step();

      // Reset original rheology for the system operator
      this->physical_properties_manager->set_rheology(original_viscosity_model);
    }
  else if (initial_condition_type ==
           Parameters::FluidDynamicsInitialConditionType::ramp)
    {
      this->pcout << "*********************************" << std::endl;
      this->pcout << " Initial condition using ramp " << std::endl;
      this->pcout << "*********************************" << std::endl;

      Timer timer(this->mpi_communicator);

      // Set the nodal values to have an initial condition that is adequate
      this->set_nodal_values();
      this->present_solution.update_ghost_values();

      // Create a pointer to the current viscosity model
      std::shared_ptr<RheologicalModel> viscosity_model =
        this->physical_properties_manager->get_rheology();

      // Gather the kinematic viscosity for the simulation
      const double original_viscosity =
        viscosity_model->get_kinematic_viscosity();

      // Gather kinematic viscosity ramp parameters
      const int n_iter_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .n_iter;
      double kinematic_viscosity =
        n_iter_viscosity > 0 ?
          this->simulation_parameters.initial_condition->ramp.ramp_viscosity
            .kinematic_viscosity_init :
          original_viscosity;
      const double alpha_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .alpha;

      // Ramp on kinematic viscosity
      for (int i = 0; i < n_iter_viscosity; ++i)
        {
          this->pcout << std::setprecision(4)
                      << "********* Solution for kinematic viscosity = " +
                           std::to_string(kinematic_viscosity) + " *********"
                      << std::endl;

          // Set the temporary viscosity to the one given by the ramp
          viscosity_model->set_kinematic_viscosity(kinematic_viscosity);

          this->simulation_control->set_assembly_method(
            Parameters::SimulationControl::TimeSteppingMethod::steady);

          // Solve the problem with the temporary viscosity
          PhysicsSolver<LinearAlgebra::distributed::Vector<double>>::
            solve_non_linear_system(false);
          this->finish_time_step();

          // Update the viscosity to the next one in the ramp parameters
          kinematic_viscosity +=
            alpha_viscosity * (original_viscosity - kinematic_viscosity);
        }


      // Reset kinematic viscosity to original value for the system operator
      viscosity_model->set_kinematic_viscosity(original_viscosity);

      timer.stop();

      if (this->simulation_parameters.timer.type !=
          Parameters::Timer::Type::none)
        {
          this->pcout << "*********************************" << std::endl;
          this->pcout << " Time spent in ramp: " << timer.wall_time() << "s"
                      << std::endl;
          this->pcout << "*********************************" << std::endl;
        }
    }
  else
    {
      throw std::runtime_error(
        "Type of initial condition is not supported by MF Navier-Stokes");
    }
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::assemble_system_matrix()
{
  // Required for compilation but not used for matrix free solvers.
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  // Update the precomputed values needed for the evaluation of the residual.
  // This is needed, otherwise the line-search mechanism used in the Newton
  // method might fail even though the Newton step should have been accepted
  // due to a wrong evaluation of the residual and, consequently, a wrong
  // evaluation of the step length.
  this->evaluation_point.update_ghost_values();
  this->system_operator->evaluate_non_linear_term_and_calculate_tau(
    this->evaluation_point);

  this->system_operator->evaluate_residual(this->system_rhs,
                                           this->evaluation_point);

  // If mortar is enabled, add RHS entries
  if (this->simulation_parameters.mortar_parameters.enable)
    {
      this->system_operator->mortar_coupling_operator_mf->vmult_add(
        this->system_rhs, this->evaluation_point);
      this->system_rhs.compress(VectorOperation::add);
    }

  this->system_rhs *= -1.0;

  // Provide residual to simulation control for stopping criterion when using
  // steady bdf
  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::update_multiphysics_time_average_solution()
{
  TimerOutput::Scope t(this->computing_timer,
                       "Update multiphysics average solution");

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      TrilinosWrappers::MPI::Vector temp_average_velocities(
        this->locally_owned_dofs, this->mpi_communicator);
      convert_vector_dealii_to_trilinos(
        temp_average_velocities,
        this->average_velocities->get_average_velocities());
      this->multiphysics_average_velocities = temp_average_velocities;

#ifndef LETHE_USE_LDV
      this->multiphysics->set_time_average_solution(
        PhysicsID::fluid_dynamics, &this->multiphysics_average_velocities);
#else
      this->multiphysics->set_time_average_solution(
        PhysicsID::fluid_dynamics,
        &this->average_velocities->get_average_velocities());
#endif
    }
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::calculate_time_derivative_previous_solutions()
{
  this->time_derivative_previous_solutions = 0;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();

  for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
    {
      this->time_derivative_previous_solutions.add(bdf_coefs[p + 1],
                                                   this->previous_solutions[p]);
    }
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::create_GMG()
{
  gmg_preconditioner = std::make_shared<MFNavierStokesPreconditionGMG<dim>>(
    this->simulation_parameters,
    this->dof_handler,
    this->dof_handler_fe_q_iso_q1);

  gmg_preconditioner->reinit(this->get_mapping(),
                             this->cell_quadrature,
                             this->forcing_function,
                             this->simulation_control,
                             this->physical_properties_manager,
                             this->fe);
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::initialize_GMG()
{
  // Initialize everything related to heat transfer within the MG algorithm
  if (this->simulation_parameters.multiphysics.buoyancy_force)
    dynamic_cast<MFNavierStokesPreconditionGMG<dim> *>(gmg_preconditioner.get())
      ->initialize_auxiliary_physics(*this->multiphysics->get_dof_handler(
                                       PhysicsID::heat_transfer),
                                     this->temperature_present_solution);

  dynamic_cast<MFNavierStokesPreconditionGMG<dim> *>(gmg_preconditioner.get())
    ->initialize(this->simulation_control,
                 this->flow_control,
                 this->present_solution,
                 this->time_derivative_previous_solutions);
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::setup_GMG()
{
  TimerOutput::Scope t(this->computing_timer, "Setup GMG");

  if (!gmg_preconditioner)
    this->create_GMG();

  this->initialize_GMG();
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "Setup ILU");

  int current_preconditioner_fill_level =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_operator->get_system_matrix(),
                                 preconditionerOptions);
}


template <int dim>
void
FluidDynamicsMatrixFree<dim>::print_mg_setup_times()
{
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_verbosity == Parameters::Verbosity::extra_verbose)
    {
      announce_string(this->pcout, "Multigrid setup times");
      this->gmg_preconditioner->mg_setup_timer.print_wall_time_statistics(
        MPI_COMM_WORLD);

      announce_string(this->pcout, "Multigrid vmult times");
      this->gmg_preconditioner->mg_vmult_timer.print_wall_time_statistics(
        MPI_COMM_WORLD);

      announce_string(this->pcout, "System operator times");
      this->system_operator->timer.print_wall_time_statistics(MPI_COMM_WORLD);

      auto mg_operators = this->gmg_preconditioner->get_mg_operators();
      for (unsigned int level = mg_operators.min_level();
           level <= mg_operators.max_level();
           level++)
        {
          announce_string(this->pcout,
                          "Operator level " + std::to_string(level) + " times");
          mg_operators[level]->timer.print_wall_time_statistics(MPI_COMM_WORLD);

          // Reset timer if output is set to every iteration
          mg_operators[level]->timer.reset();
        }

      auto mg_smoother_preconditioners =
        this->gmg_preconditioner->get_mg_smoother_preconditioners();
      for (unsigned int level = mg_operators.min_level();
           level <= mg_operators.max_level();
           level++)
        {
          announce_string(this->pcout,
                          "Preconditioner level " + std::to_string(level) +
                            " times");
          mg_smoother_preconditioners[level]->timer_print();

          // Reset timer if output is set to every iteration
          mg_smoother_preconditioners[level]->timer_reset();
        }

      // Reset timers if output is set to every iteration
      this->gmg_preconditioner->mg_setup_timer.reset();
      this->gmg_preconditioner->mg_vmult_timer.reset();
    }
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::update_solutions_for_multiphysics()
{
  TimerOutput::Scope t(this->computing_timer,
                       "Update solutions for multiphysics");

  // Provide the fluid dynamics dof_handler to the multiphysics interface
  this->multiphysics->set_dof_handler(PhysicsID::fluid_dynamics,
                                      &this->dof_handler);

  // Convert the present solution to multiphysics vector type and provide it
  // to the multiphysics interface
  TrilinosWrappers::MPI::Vector temp_solution(this->locally_owned_dofs,
                                              this->mpi_communicator);

  this->present_solution.update_ghost_values();
  convert_vector_dealii_to_trilinos(temp_solution, this->present_solution);
  multiphysics_present_solution = temp_solution;

#ifndef LETHE_USE_LDV
  this->multiphysics->set_solution(PhysicsID::fluid_dynamics,
                                   &this->multiphysics_present_solution);
#else
  this->multiphysics->set_solution(PhysicsID::fluid_dynamics,
                                   &this->present_solution);
#endif

  // Convert the previous solutions to multiphysics vector type and provide
  // them to the multiphysics interface
  const unsigned int number_of_previous_solutions =
    this->simulation_control->get_number_of_previous_solution_in_assembly();

  std::vector<TrilinosWrappers::MPI::Vector> temp_previous_solutions;

  temp_previous_solutions.resize(number_of_previous_solutions);
  for (auto &solution : temp_previous_solutions)
    solution.reinit(this->locally_owned_dofs, this->mpi_communicator);

  for (unsigned int i = 0; i < number_of_previous_solutions; i++)
    {
      this->previous_solutions[i].update_ghost_values();
      convert_vector_dealii_to_trilinos(temp_previous_solutions[i],
                                        this->previous_solutions[i]);

      this->multiphysics_previous_solutions[i] = temp_previous_solutions[i];
    }

#ifndef LETHE_USE_LDV
  this->multiphysics->set_previous_solutions(
    PhysicsID::fluid_dynamics, &this->multiphysics_previous_solutions);
#else
  this->multiphysics->set_previous_solutions(PhysicsID::fluid_dynamics,
                                             &this->previous_solutions);
#endif
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::update_solutions_for_fluid_dynamics()
{
  TimerOutput::Scope t(this->computing_timer,
                       "Update solutions for fluid dynamics");

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      // Get present solution and dof handler of the heat transfer
      const auto &heat_solution =
        *this->multiphysics->get_solution(PhysicsID::heat_transfer);

#ifndef LETHE_USE_LDV
      const auto &heat_dof_handler =
        *this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      // Copy solution to temporary vector
      TrilinosWrappers::MPI::Vector temp_heat_solution;
      temp_heat_solution.reinit(heat_dof_handler.locally_owned_dofs(),
                                this->mpi_communicator);
      temp_heat_solution = heat_solution;

      // Initialize deal.II vector to store solution
      this->temperature_present_solution.reinit(
        heat_dof_handler.locally_owned_dofs(),
        DoFTools::extract_locally_active_dofs(heat_dof_handler),
        this->mpi_communicator);

      // Perform copy between two vector types
      convert_vector_trilinos_to_dealii(this->temperature_present_solution,
                                        temp_heat_solution);

#else
      this->temperature_present_solution = heat_solution;
#endif
      // Update ghost values for deal.II vector
      this->temperature_present_solution.update_ghost_values();
    }
}
template <int dim>
void
FluidDynamicsMatrixFree<dim>::setup_preconditioner()
{
  this->present_solution.update_ghost_values();

  this->computing_timer.enter_subsection("Evaluate non linear term and tau");

  this->system_operator->evaluate_non_linear_term_and_calculate_tau(
    this->present_solution);

  this->computing_timer.leave_subsection("Evaluate non linear term and tau");

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::ilu)
    setup_ILU();
  else if ((this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .preconditioner ==
            Parameters::LinearSolver::PreconditionerType::lsmg) ||
           (this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .preconditioner ==
            Parameters::LinearSolver::PreconditionerType::gcmg))
    setup_GMG();
  else
    AssertThrow(
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::ilu ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::amg,
      ExcMessage(
        "This linear solver does not support this preconditioner. Only <ilu|lsmg|gcmg> preconditioners are supported."));
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::solve_linear_system(
  const bool initial_step,
  const bool /* renewed_matrix */)
{
  const double absolute_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .minimum_residual;
  const double relative_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .relative_residual;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .solver == Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step, absolute_residual, relative_residual);
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .solver == Parameters::LinearSolver::SolverType::direct)
    solve_system_direct(initial_step, absolute_residual, relative_residual);
  else
    AssertThrow(false, ExcMessage("This solver is not allowed"));
  this->rescale_pressure_dofs_in_newton_update();
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::assemble_L2_projection()
{
  // TODO
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::solve_system_GMRES(const bool   initial_step,
                                                 const double absolute_residual,
                                                 const double relative_residual)
{
  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  SolverGMRES<VectorType>::AdditionalData solver_parameters;

  solver_parameters.max_n_tmp_vectors =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .max_krylov_vectors;
  solver_parameters.right_preconditioning = true;

  SolverGMRES<VectorType> solver(solver_control, solver_parameters);

  this->newton_update = 0.0;

  this->computing_timer.enter_subsection("Solve linear system");

  if ((this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
         .preconditioner ==
       Parameters::LinearSolver::PreconditionerType::lsmg) ||
      (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
         .preconditioner == Parameters::LinearSolver::PreconditionerType::gcmg))
    {
      solver.solve(*(this->system_operator),
                   this->newton_update,
                   this->system_rhs,
                   *(this->gmg_preconditioner));

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity != Parameters::Verbosity::quiet)
        this->gmg_preconditioner->print_relevant_info();
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::ilu)
    solver.solve(*(this->system_operator),
                 this->newton_update,
                 this->system_rhs,
                 *(this->ilu_preconditioner));
  else
    AssertThrow(
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::ilu ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::lsmg ||
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
            .preconditioner ==
          Parameters::LinearSolver::PreconditionerType::gcmg,
      ExcMessage(
        "This linear solver does not support this preconditioner. Only <ilu|lsmg|gcmg> preconditioners are supported."));

  this->computing_timer.leave_subsection("Solve linear system");

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps to reach a residual norm of "
                  << solver_control.last_value() << std::endl;
    }

  this->computing_timer.enter_subsection(
    "Distribute constraints after linear solve");

  constraints_used.distribute(this->newton_update);

  this->computing_timer.leave_subsection(
    "Distribute constraints after linear solve");
}

template <int dim>
void
FluidDynamicsMatrixFree<dim>::solve_system_direct(
  const bool   initial_step,
  const double absolute_residual,
  const double relative_residual)
{
  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               linear_solver_tolerance,
                               true,
                               true);

  TrilinosWrappers::SolverDirect solver(solver_control);

  this->newton_update = 0.0;

  this->computing_timer.enter_subsection("Solve linear system");

  solver.initialize(this->system_operator->get_system_matrix());
  solver.solve(this->newton_update, this->system_rhs);

  this->computing_timer.leave_subsection("Solve linear system");

  this->computing_timer.enter_subsection(
    "Distribute constraints after linear solve");

  constraints_used.distribute(this->newton_update);

  this->computing_timer.leave_subsection(
    "Distribute constraints after linear solve");
}

template class FluidDynamicsMatrixFree<2>;
template class FluidDynamicsMatrixFree<3>;

template class MFNavierStokesPreconditionGMG<2>;
template class MFNavierStokesPreconditionGMG<3>;

template class MFNavierStokesPreconditionGMGBase<2>;
template class MFNavierStokesPreconditionGMGBase<3>;
