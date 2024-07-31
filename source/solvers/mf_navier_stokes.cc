/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------*/

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/multiphysics.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/mf_navier_stokes.h>

#include <deal.II/base/multithread_info.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/vector_tools.h>

/**
 * @brief Helper function that allows to convert deal.II vectors to Trilinos vectors.
 *
 * @tparam number Abstract type for number across the class (i.e., double).
 * @param out Destination TrilinosWrappers::MPI::Vector vector.
 * @param in Source LinearAlgebra::distributed::Vector<number> vector.
 */
template <typename number>
void
convert_vector_dealii_to_trilinos(
  TrilinosWrappers::MPI::Vector                    &out,
  const LinearAlgebra::distributed::Vector<number> &in)
{
  LinearAlgebra::ReadWriteVector<double> rwv(out.locally_owned_elements());
  rwv.import_elements(in, VectorOperation::insert);
  out.import_elements(rwv, VectorOperation::insert);
}

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
    SmartPointer<
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



  namespace internal
  {
    namespace MGCoarseGridApplyPreconditioner
    {
      template <class VectorType,
                class PreconditionerType,
                typename std::enable_if<
                  std::is_same<typename VectorType::value_type, double>::value,
                  VectorType>::type * = nullptr>
      void
      solve(const PreconditionerType preconditioner,
            VectorType              &dst,
            const VectorType        &src)
      {
        // to allow the case that the preconditioner was only set up on a
        // subset of processes
        if (preconditioner != nullptr)
          preconditioner->vmult(dst, src);
      }

      template <class VectorType,
                class PreconditionerType,
                typename std::enable_if<
                  !std::is_same<typename VectorType::value_type, double>::value,
                  VectorType>::type * = nullptr>
      void
      solve(const PreconditionerType preconditioner,
            VectorType              &dst,
            const VectorType        &src)
      {
        LinearAlgebra::distributed::Vector<double> src_;
        LinearAlgebra::distributed::Vector<double> dst_;

        src_ = src;
        dst_ = dst;

        // to allow the case that the preconditioner was only set up on a
        // subset of processes
        if (preconditioner != nullptr)
          preconditioner->vmult(dst_, src_);

        dst = dst_;
      }
    } // namespace MGCoarseGridApplyPreconditioner
  }   // namespace internal


  template <class VectorType, class PreconditionerType>
  void
  MGCoarseGridApplyPreconditioner<VectorType, PreconditionerType>::operator()(
    const unsigned int /*level*/,
    VectorType       &dst,
    const VectorType &src) const
  {
    internal::MGCoarseGridApplyPreconditioner::solve(preconditioner, dst, src);
  }
} // namespace dealii

template <int dim>
MFNavierStokesPreconditionGMG<dim>::MFNavierStokesPreconditionGMG(
  const SimulationParameters<dim>         &simulation_parameters,
  const DoFHandler<dim>                   &dof_handler,
  const DoFHandler<dim>                   &dof_handler_fe_q_iso_q1,
  const std::shared_ptr<Mapping<dim>>     &mapping,
  const std::shared_ptr<Quadrature<dim>>  &cell_quadrature,
  const std::shared_ptr<Function<dim>>     forcing_function,
  const std::shared_ptr<SimulationControl> simulation_control,
  const std::shared_ptr<FESystem<dim>>     fe)
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , simulation_parameters(simulation_parameters)
  , dof_handler(dof_handler)
  , dof_handler_fe_q_iso_q1(dof_handler_fe_q_iso_q1)
  , mg_setup_timer(this->pcout, TimerOutput::never, TimerOutput::wall_times)
  , mg_vmult_timer(this->pcout, TimerOutput::never, TimerOutput::wall_times)
{
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      // Define maximum and minimum level according to triangulation
      const unsigned int n_h_levels =
        this->dof_handler.get_triangulation().n_global_levels();
      this->minlevel = 0;
      this->maxlevel = n_h_levels - 1;

      // If multigrid number of levels or minimum number of cells in level are
      // specified, change the min level fnd print levels information

      int mg_min_level =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_min_level;

      AssertThrow(
        mg_min_level <= static_cast<int>(MGTools::max_level_for_coarse_mesh(
                          this->dof_handler.get_triangulation())),
        ExcMessage(std::string(
          "The maximum level allowed for the coarse mesh (mg min level) is: " +
          std::to_string(MGTools::max_level_for_coarse_mesh(
            this->dof_handler.get_triangulation())) +
          ".")));

      int mg_level_min_cells =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_level_min_cells;


      std::vector<unsigned int> n_cells_on_levels(
        this->dof_handler.get_triangulation().n_global_levels(), 0);

      for (unsigned int l = 0;
           l < this->dof_handler.get_triangulation().n_levels();
           ++l)
        for (const auto &cell :
             this->dof_handler.get_triangulation().cell_iterators_on_level(l))
          if (cell->is_locally_owned_on_level())
            n_cells_on_levels[l]++;

      Utilities::MPI::sum(n_cells_on_levels,
                          this->dof_handler.get_communicator(),
                          n_cells_on_levels);
      AssertThrow(
        mg_level_min_cells <= static_cast<int>(n_cells_on_levels[maxlevel]),
        ExcMessage(
          "The mg level min cells specified are larger than the cells of the finest mg level."));


      if (mg_min_level != -1)
        this->minlevel = mg_min_level;

      if (mg_level_min_cells != -1)
        {
          for (unsigned int level = this->minlevel; level <= this->maxlevel;
               ++level)
            if (static_cast<int>(n_cells_on_levels[level]) >=
                mg_level_min_cells)
              {
                this->minlevel = level;
                break;
              }
        }

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity != Parameters::Verbosity::quiet)
        {
          this->pcout << std::endl;
          this->pcout << "  -Levels of MG preconditioner:" << std::endl;
          for (unsigned int level = this->minlevel; level <= this->maxlevel;
               ++level)
            this->pcout << "    Level " << level - this->minlevel << ": "
                        << this->dof_handler.n_dofs(level) << " DoFs, "
                        << n_cells_on_levels[level] << " cells" << std::endl;
          this->pcout << std::endl;
        }

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_verbosity == Parameters::Verbosity::extra_verbose)
        {
          this->pcout << "  -MG vertical communication efficiency: "
                      << MGTools::vertical_communication_efficiency(
                           this->dof_handler.get_triangulation())
                      << std::endl;

          this->pcout << "  -MG workload imbalance: "
                      << MGTools::workload_imbalance(
                           this->dof_handler.get_triangulation())
                      << std::endl;
        }

      std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
        partitioners(this->dof_handler.get_triangulation().n_global_levels());

      // Local object for constraints of the different levels
      MGLevelObject<AffineConstraints<double>> level_constraints;

      // Resize all multilevel objects according to level
      this->mg_operators.resize(this->minlevel, this->maxlevel);
      level_constraints.resize(this->minlevel, this->maxlevel);
      this->ls_mg_interface_in.resize(this->minlevel, this->maxlevel);
      this->ls_mg_operators.resize(this->minlevel, this->maxlevel);

      // Fill the level constraints
      this->mg_setup_timer.enter_subsection("Set boundary conditions");

      this->mg_constrained_dofs.clear();
      this->mg_constrained_dofs.initialize(this->dof_handler);

      FEValuesExtractors::Vector velocities(0);
      FEValuesExtractors::Scalar pressure(dim);

      for (unsigned int i_bc = 0;
           i_bc < this->simulation_parameters.boundary_conditions.size;
           ++i_bc)
        {
          if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
              BoundaryConditions::BoundaryType::slip)
            {
              std::set<types::boundary_id> no_normal_flux_boundaries;
              no_normal_flux_boundaries.insert(
                this->simulation_parameters.boundary_conditions.id[i_bc]);
              for (unsigned int level = this->minlevel; level <= this->maxlevel;
                   ++level)
                {
                  AffineConstraints<double> temp_constraints;
                  temp_constraints.clear();
                  const IndexSet locally_relevant_level_dofs =
                    DoFTools::extract_locally_relevant_level_dofs(
                      this->dof_handler, level);
                  temp_constraints.reinit(locally_relevant_level_dofs);
                  VectorTools::compute_no_normal_flux_constraints_on_level(
                    this->dof_handler,
                    0,
                    no_normal_flux_boundaries,
                    temp_constraints,
                    *mapping,
                    this->mg_constrained_dofs.get_refinement_edge_indices(
                      level),
                    level);
                  temp_constraints.close();
                  this->mg_constrained_dofs.add_user_constraints(
                    level, temp_constraints);
                }
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::periodic)
            {
              /*already taken into account when mg_constrained_dofs is
               * initialized*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::pressure)
            {
              /*do nothing*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::function_weak)
            {
              /*do nothing*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::partial_slip)
            {
              /*do nothing*/
            }
          else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                   BoundaryConditions::BoundaryType::outlet)
            {
              /*do nothing*/
            }
          else
            {
              std::set<types::boundary_id> dirichlet_boundary_id = {
                this->simulation_parameters.boundary_conditions.id[i_bc]};
              this->mg_constrained_dofs.make_zero_boundary_constraints(
                this->dof_handler,
                dirichlet_boundary_id,
                fe->component_mask(velocities));
            }
        }

      this->mg_setup_timer.leave_subsection("Set boundary conditions");

      // Create mg operators for each level and additional operators needed only
      // for local smoothing
      for (unsigned int level = this->minlevel; level <= this->maxlevel;
           ++level)
        {
          level_constraints[level].clear();

          const IndexSet relevant_dofs =
            DoFTools::extract_locally_relevant_level_dofs(this->dof_handler,
                                                          level);

          level_constraints[level].reinit(relevant_dofs);

#if DEAL_II_VERSION_GTE(9, 6, 0)
          this->mg_constrained_dofs.merge_constraints(
            level_constraints[level], level, true, false, true, true);
#else
          AssertThrow(
            false,
            ExcMessage(
              "The constraints for the lsmg preconditioner require a most recent version of deal.II."));
#endif

          if (this->simulation_parameters.boundary_conditions
                .fix_pressure_constant &&
              level == this->minlevel)
            {
              unsigned int min_index = numbers::invalid_unsigned_int;

              std::vector<types::global_dof_index> dof_indices;

              // Loop over the cells to identify the min index
              for (const auto &cell : this->dof_handler.active_cell_iterators())
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
                                    this->dof_handler.get_communicator());

              if (relevant_dofs.is_element(min_index))
                level_constraints[level].add_line(min_index);
            }

          level_constraints[level].close();

          this->mg_setup_timer.enter_subsection("Set up operators");

          // Provide appropriate quadrature depending on the type of elements of
          // the level
          auto quadrature_mg = *cell_quadrature;
          if (this->simulation_parameters.linear_solver
                .at(PhysicsID::fluid_dynamics)
                .mg_use_fe_q_iso_q1 &&
              level == this->minlevel)
            {
              const auto points =
                QGaussLobatto<1>(this->dof_handler.get_fe().degree + 1)
                  .get_points();

              quadrature_mg = QIterated<dim>(QGauss<1>(2), points);
            }

          this->mg_operators[level] =
            std::make_shared<NavierStokesStabilizedOperator<dim, double>>();

          this->mg_operators[level]->reinit(
            *mapping,
            (this->simulation_parameters.linear_solver
               .at(PhysicsID::fluid_dynamics)
               .mg_use_fe_q_iso_q1 &&
             level == this->minlevel) ?
              this->dof_handler_fe_q_iso_q1 :
              this->dof_handler,
            level_constraints[level],
            quadrature_mg,
            forcing_function,
            this->simulation_parameters.physical_properties_manager
              .get_kinematic_viscosity_scale(),
            this->simulation_parameters.stabilization.stabilization,
            level,
            simulation_control,
            this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .mg_enable_hessians_jacobian,
            true);

          this->ls_mg_operators[level].initialize(*(this->mg_operators)[level]);
          this->ls_mg_interface_in[level].initialize(
            *(this->mg_operators)[level]);

          partitioners[level] =
            this->mg_operators[level]->get_vector_partitioner();

          this->mg_setup_timer.leave_subsection("Set up operators");
        }

      // Create transfer operators
      this->mg_setup_timer.enter_subsection("Create transfer operator");

      this->mg_transfer_ls = std::make_shared<LSTransferType>();

      this->mg_transfer_ls->initialize_constraints(this->mg_constrained_dofs);
      this->mg_transfer_ls->build(this->dof_handler, partitioners);

      this->mg_setup_timer.leave_subsection("Create transfer operator");
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      AssertThrow(*cell_quadrature ==
                    QGauss<dim>(this->dof_handler.get_fe().degree + 1),
                  ExcNotImplemented());

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

      // Local object for constraints of the different levels
      MGLevelObject<AffineConstraints<typename VectorType::value_type>>
        constraints;

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
                QGaussLobatto<1>(this->dof_handler.get_fe().degree + 1)
                  .get_points();

              this->dof_handlers[l].distribute_dofs(
                FESystem<dim>(FE_Q_iso_Q1<dim>(points), dim + 1));
            }
          else
            this->dof_handlers[l].distribute_dofs(
              FESystem<dim>(FE_Q<dim>(levels[l].second), dim + 1));
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
          level_constraint.reinit(locally_relevant_dofs);

          DoFTools::make_hanging_node_constraints(level_dof_handler,
                                                  level_constraint);

          FEValuesExtractors::Vector velocities(0);
          FEValuesExtractors::Scalar pressure(dim);

          for (unsigned int i_bc = 0;
               i_bc < this->simulation_parameters.boundary_conditions.size;
               ++i_bc)
            {
              if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                  BoundaryConditions::BoundaryType::slip)
                {
                  std::set<types::boundary_id> no_normal_flux_boundaries;
                  no_normal_flux_boundaries.insert(
                    this->simulation_parameters.boundary_conditions.id[i_bc]);
                  VectorTools::compute_no_normal_flux_constraints(
                    level_dof_handler,
                    0,
                    no_normal_flux_boundaries,
                    level_constraint,
                    *mapping);
                }
              else if (this->simulation_parameters.boundary_conditions
                         .type[i_bc] ==
                       BoundaryConditions::BoundaryType::periodic)
                {
                  DoFTools::make_periodicity_constraints(
                    level_dof_handler,
                    this->simulation_parameters.boundary_conditions.id[i_bc],
                    this->simulation_parameters.boundary_conditions
                      .periodic_id[i_bc],
                    this->simulation_parameters.boundary_conditions
                      .periodic_direction[i_bc],
                    level_constraint);
                }
              else if (this->simulation_parameters.boundary_conditions
                         .type[i_bc] ==
                       BoundaryConditions::BoundaryType::pressure)
                {
                  /*do nothing*/
                }
              else if (this->simulation_parameters.boundary_conditions
                         .type[i_bc] ==
                       BoundaryConditions::BoundaryType::function_weak)
                {
                  /*do nothing*/
                }
              else if (this->simulation_parameters.boundary_conditions
                         .type[i_bc] ==
                       BoundaryConditions::BoundaryType::partial_slip)
                {
                  /*do nothing*/
                }
              else if (this->simulation_parameters.boundary_conditions
                         .type[i_bc] ==
                       BoundaryConditions::BoundaryType::outlet)
                {
                  /*do nothing*/
                }
              else
                {
                  VectorTools::interpolate_boundary_values(
                    *mapping,
                    level_dof_handler,
                    this->simulation_parameters.boundary_conditions.id[i_bc],
                    dealii::Functions::ZeroFunction<dim>(dim + 1),
                    level_constraint,
                    fe->component_mask(velocities));
                }
            }

          if (this->simulation_parameters.boundary_conditions
                .fix_pressure_constant &&
              level == this->minlevel)
            {
              unsigned int min_index = numbers::invalid_unsigned_int;

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
                                    this->dof_handler.get_communicator());

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
                QGaussLobatto<1>(this->dof_handler.get_fe().degree + 1)
                  .get_points();

              quadrature_mg = QIterated<dim>(QGauss<1>(2), points);
            }

          this->mg_operators[level] =
            std::make_shared<NavierStokesStabilizedOperator<dim, double>>();

          this->mg_operators[level]->reinit(
            *mapping,
            level_dof_handler,
            level_constraint,
            quadrature_mg,
            forcing_function,
            this->simulation_parameters.physical_properties_manager
              .get_kinematic_viscosity_scale(),
            this->simulation_parameters.stabilization.stabilization,
            numbers::invalid_unsigned_int,
            simulation_control,
            this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .mg_enable_hessians_jacobian,
            true);

          this->mg_setup_timer.leave_subsection("Set up operators");
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
}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::initialize(
  const std::shared_ptr<SimulationControl> simulation_control,
  FlowControl<dim>                        &flow_control,
  const VectorType                        &present_solution,
  const VectorType                        &time_derivative_previous_solutions)
{
  // Local objects for the different levels
  MGLevelObject<VectorType> mg_solution(this->minlevel, this->maxlevel);
  MGLevelObject<VectorType> mg_time_derivative_previous_solutions(
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
        std::make_shared<mg::Matrix<VectorType>>(this->ls_mg_operators);
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
        std::make_shared<mg::Matrix<VectorType>>(this->mg_operators);
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

  // Create smoother, fill parameters for each level and intialize it
  this->mg_setup_timer.enter_subsection("Set up and initialize smoother");

  this->mg_smoother = std::make_shared<
    MGSmootherPrecondition<OperatorType, SmootherType, VectorType>>();

  MGLevelObject<typename SmootherType::AdditionalData> smoother_data(
    this->minlevel, this->maxlevel);

  for (unsigned int level = this->minlevel; level <= this->maxlevel; ++level)
    {
      VectorType diagonal_vector;
      this->mg_operators[level]->compute_inverse_diagonal(diagonal_vector);
      smoother_data[level].preconditioner =
        std::make_shared<SmootherPreconditionerType>(diagonal_vector);
      smoother_data[level].n_iterations =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_smoother_iterations;

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_eig_estimation)
        {
#if DEAL_II_VERSION_GTE(9, 6, 0)
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
#else
          AssertThrow(
            false,
            ExcMessage(
              "The estimation of eigenvalues within LSMG requires a version of deal.II >= 9.6.0"));
#endif
        }
      else
        smoother_data[level].relaxation =
          this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_smoother_relaxation;
    }

  mg_smoother->initialize(this->mg_operators, smoother_data);

#if DEAL_II_VERSION_GTE(9, 6, 0)
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_smoother_eig_estimation &&
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .eig_estimation_verbose != Parameters::Verbosity::quiet)
    {
      // Print eigenvalue estimation for all levels
      for (unsigned int level = this->minlevel; level <= this->maxlevel;
           ++level)
        {
          VectorType vec;
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
#else
  AssertThrow(
    false,
    ExcMessage(
      "The estimation of eigenvalues within LSMG requires a version of deal.II >= 9.6.0"));
#endif

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
      SolverGMRES<VectorType>::AdditionalData solver_parameters;
      solver_parameters.max_n_tmp_vectors =
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .mg_gmres_max_krylov_vectors;

      this->coarse_grid_solver = std::make_shared<SolverGMRES<VectorType>>(
        *this->coarse_grid_solver_control, solver_parameters);

      if (this->simulation_parameters.linear_solver
            .at(PhysicsID::fluid_dynamics)
            .mg_gmres_preconditioner ==
          Parameters::LinearSolver::PreconditionerType::amg)
        {
          setup_AMG();

          this->mg_coarse = std::make_shared<
            MGCoarseGridIterativeSolver<VectorType,
                                        SolverGMRES<VectorType>,
                                        OperatorType,
                                        TrilinosWrappers::PreconditionAMG>>(
            *this->coarse_grid_solver,
            *this->mg_operators[this->minlevel],
            *this->precondition_amg);
        }
      else if (this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .mg_gmres_preconditioner ==
               Parameters::LinearSolver::PreconditionerType::ilu)
        {
          setup_ILU();

          this->mg_coarse = std::make_shared<
            MGCoarseGridIterativeSolver<VectorType,
                                        SolverGMRES<VectorType>,
                                        OperatorType,
                                        TrilinosWrappers::PreconditionILU>>(
            *this->coarse_grid_solver,
            *this->mg_operators[this->minlevel],
            *this->precondition_ilu);
        }
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_solver ==
           Parameters::LinearSolver::CoarseGridSolverType::amg)
    {
      setup_AMG();

      this->mg_coarse = std::make_shared<
        MGCoarseGridApplyPreconditioner<VectorType,
                                        TrilinosWrappers::PreconditionAMG>>(
        *this->precondition_amg);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_solver ==
           Parameters::LinearSolver::CoarseGridSolverType::ilu)
    {
      setup_ILU();

      this->mg_coarse = std::make_shared<
        MGCoarseGridApplyPreconditioner<VectorType,
                                        TrilinosWrappers::PreconditionILU>>(
        *this->precondition_ilu);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .mg_coarse_grid_solver ==
           Parameters::LinearSolver::CoarseGridSolverType::direct)
    {
#if DEAL_II_VERSION_GTE(9, 6, 0)
      TrilinosWrappers::SolverDirect::AdditionalData data;
      this->direct_solver_control =
        std::make_shared<SolverControl>(100, 1.e-10);

      this->precondition_direct =
        std::make_shared<TrilinosWrappers::SolverDirect>(
          *this->direct_solver_control, data);

      this->precondition_direct->initialize(
        this->mg_operators[this->minlevel]->get_system_matrix());

      this->mg_coarse = std::make_shared<
        MGCoarseGridApplyPreconditioner<VectorType,
                                        TrilinosWrappers::SolverDirect>>(
        *this->precondition_direct);
#else
      AssertThrow(
        false,
        ExcMessage(
          "The usage of a direct solver as coarse grid solver requires a version of deal.II >= 9.6.0"));
#endif
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
        std::make_shared<mg::Matrix<VectorType>>(this->ls_mg_interface_in);

      // Create main MG object
      this->mg = std::make_shared<Multigrid<VectorType>>(*this->mg_matrix,
                                                         *this->mg_coarse,
                                                         *this->mg_transfer_ls,
                                                         *this->mg_smoother,
                                                         *this->mg_smoother,
                                                         this->minlevel);

      if (this->dof_handler.get_triangulation().has_hanging_nodes())
        this->mg->set_edge_in_matrix(*this->mg_interface_matrix_in);

      // Create MG preconditioner
      this->ls_multigrid_preconditioner =
        std::make_shared<PreconditionMG<dim, VectorType, LSTransferType>>(
          this->dof_handler, *this->mg, *this->mg_transfer_ls);
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      // Create main MG object
      this->mg = std::make_shared<Multigrid<VectorType>>(*this->mg_matrix,
                                                         *this->mg_coarse,
                                                         *this->mg_transfer_gc,
                                                         *this->mg_smoother,
                                                         *this->mg_smoother);

      // Create MG preconditioner
      this->gc_multigrid_preconditioner =
        std::make_shared<PreconditionMG<dim, VectorType, GCTransferType>>(
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
MFNavierStokesPreconditionGMG<dim>::vmult(VectorType       &dst,
                                          const VectorType &src) const
{
  if (this->ls_multigrid_preconditioner)
    this->ls_multigrid_preconditioner->vmult(dst, src);
  else if (this->gc_multigrid_preconditioner)
    this->gc_multigrid_preconditioner->vmult(dst, src);
  else
    AssertThrow(false, ExcNotImplemented());

  // Save number of coarse grid iterations needed in one vmult
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_coarse_grid_solver ==
      Parameters::LinearSolver::CoarseGridSolverType::gmres)
    this->coarse_grid_iterations.emplace_back(
      this->coarse_grid_solver_control->last_step());
}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::print_relevant_info() const
{
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .mg_coarse_grid_solver ==
      Parameters::LinearSolver::CoarseGridSolverType::gmres)
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
const MGLevelObject<std::shared_ptr<NavierStokesOperatorBase<dim, double>>> &
MFNavierStokesPreconditionGMG<dim>::get_mg_operators() const
{
  return this->mg_operators;
}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::setup_AMG()
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
#if DEAL_II_VERSION_GTE(9, 6, 0)
          // Constant modes for velocity and pressure
          DoFTools::extract_level_constant_modes(this->minlevel,
                                                 this->dof_handler,
                                                 components,
                                                 constant_modes);
#else
          AssertThrow(
            false,
            ExcMessage(
              "The extraction of constant modes for the AMG coarse-grid solver requires a version of deal.II >= 9.6.0"));
#endif
        }
      else if (this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .preconditioner ==
               Parameters::LinearSolver::PreconditionerType::gcmg)
        {
          // Constant modes for velocity and pressure
          DoFTools::extract_constant_modes(this->dof_handlers[this->minlevel],
                                           components,
                                           constant_modes);
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

      this->precondition_amg =
        std::make_shared<TrilinosWrappers::PreconditionAMG>();

      this->precondition_amg->initialize(min_level_matrix, parameter_ml);
    }
  else
    {
      this->precondition_amg =
        std::make_shared<TrilinosWrappers::PreconditionAMG>();

      this->precondition_amg->initialize(min_level_matrix, amg_data);
    }
}

template <int dim>
void
MFNavierStokesPreconditionGMG<dim>::setup_ILU()
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

  this->precondition_ilu =
    std::make_shared<TrilinosWrappers::PreconditionILU>();

  this->precondition_ilu->initialize(
    this->mg_operators[this->minlevel]->get_system_matrix(),
    preconditionerOptions);
}

template <int dim>
MFNavierStokesSolver<dim>::MFNavierStokesSolver(
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

  system_operator =
    std::make_shared<NavierStokesStabilizedOperator<dim, double>>();

  if (!this->simulation_parameters.source_term.enable)
    {
      this->forcing_function.reset();
    }
}

template <int dim>
MFNavierStokesSolver<dim>::~MFNavierStokesSolver()
{
  this->dof_handler.clear();
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve()
{
  this->computing_timer.enter_subsection("Read mesh and manifolds");

  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->computing_timer.leave_subsection("Read mesh and manifolds");

  this->setup_dofs();
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

      bool refinement_step;
      if (this->simulation_parameters.mesh_adaptation.refinement_at_frequency)
        refinement_step =
          this->simulation_control->get_step_number() %
            this->simulation_parameters.mesh_adaptation.frequency !=
          0;
      else
        refinement_step = this->simulation_control->get_step_number() == 0;
      if ((refinement_step ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          update_boundary_conditions();
          this->multiphysics->update_boundary_conditions();
        }

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

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

      // Provide the fluid dynamics dof_handler, present solution and previous
      // solution to the multiphysics interface only if other physics
      // apart from fluid dynamics are enabled
      if (this->multiphysics->get_active_physics().size() > 1)
        update_solutions_for_multiphysics();

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
MFNavierStokesSolver<dim>::setup_dofs_fd()
{
  TimerOutput::Scope t(this->computing_timer, "Setup DoFs");

  // Clear the preconditioners
  ilu_preconditioner.reset();
  gmg_preconditioner.reset();

  // Clear matrix free operator
  this->system_operator->clear();

  // Fill the dof handler and initialize vectors
  this->dof_handler.distribute_dofs(*this->fe);

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
        }
    }

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  this->locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(this->dof_handler);

  // Non Zero constraints
  define_non_zero_constraints();

  // Zero constraints
  define_zero_constraints();

  // Initialize matrix-free object
  unsigned int mg_level = numbers::invalid_unsigned_int;
  this->system_operator->reinit(
    *this->mapping,
    this->dof_handler,
    this->zero_constraints,
    *this->cell_quadrature,
    this->forcing_function,
    this->simulation_parameters.physical_properties_manager
      .get_kinematic_viscosity_scale(),
    this->simulation_parameters.stabilization.stabilization,
    mg_level,
    this->simulation_control,
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .enable_hessians_jacobian,
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .enable_hessians_residual);


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

      if (this->simulation_parameters.restart_parameters.checkpoint)
        {
          this->average_velocities->initialize_checkpoint_vectors(
            this->locally_owned_dofs,
            this->locally_relevant_dofs,
            this->mpi_communicator);
        }
    }

  double global_volume =
    GridTools::volume(*this->triangulation, *this->mapping);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;

  // Initialize Trilinos vectors used to pass solution to multiphysics interface
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
MFNavierStokesSolver<dim>::update_boundary_conditions()
{
  if (this->simulation_parameters.boundary_conditions.time_dependent)
    {
      double time = this->simulation_control->get_current_time();
      for (unsigned int i_bc = 0;
           i_bc < this->simulation_parameters.boundary_conditions.size;
           ++i_bc)
        {
          this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
            .u.set_time(time);
          this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
            .v.set_time(time);
          this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
            .w.set_time(time);
          this->simulation_parameters.boundary_conditions
            .bcPressureFunction[i_bc]
            .p.set_time(time);
        }
      define_non_zero_constraints();
      // Distribute constraints
      auto &nonzero_constraints = this->nonzero_constraints;
      nonzero_constraints.distribute(this->local_evaluation_point);
      this->present_solution = this->local_evaluation_point;
    }
}

template <int dim>
void
MFNavierStokesSolver<dim>::set_initial_condition_fd(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->present_solution.update_ghost_values();
      this->finish_time_step();
    }
  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      // Set the nodal values to have an initial condition that is adequate
      this->set_nodal_values();
      this->present_solution.update_ghost_values();

      // Create a pointer to the current viscosity model
      std::shared_ptr<RheologicalModel> viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Keep in memory the initial viscosity
      const double viscosity_end = viscosity_model->get_kinematic_viscosity();

      // Set it to the initial condition viscosity
      viscosity_model->set_kinematic_viscosity(
        this->simulation_parameters.initial_condition->kinematic_viscosity);

      // Set the kinematic viscosity in the system operator to be the
      // temporary viscosity
      this->system_operator->set_kinematic_viscosity(
        this->simulation_parameters.physical_properties_manager
          .get_kinematic_viscosity_scale());

      // Set temporary viscosity in the multigrid operators
      if ((this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::lsmg) ||
          (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg))
        {
          // Create the mg operators if they do not exist to be able
          // to change the viscosity for all of them
          if (!gmg_preconditioner)
            gmg_preconditioner =
              std::make_shared<MFNavierStokesPreconditionGMG<dim>>(
                this->simulation_parameters,
                this->dof_handler,
                this->dof_handler_fe_q_iso_q1,
                this->mapping,
                this->cell_quadrature,
                this->forcing_function,
                this->simulation_control,
                this->fe);

          auto mg_operators = this->gmg_preconditioner->get_mg_operators();
          for (unsigned int level = mg_operators.min_level();
               level <= mg_operators.max_level();
               level++)
            {
              mg_operators[level]->set_kinematic_viscosity(
                this->simulation_parameters.physical_properties_manager
                  .get_kinematic_viscosity_scale());
            }
        }

      // Solve the problem with the temporary viscosity
      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::steady);
      PhysicsSolver<LinearAlgebra::distributed::Vector<double>>::
        solve_non_linear_system(false);
      this->finish_time_step();

      // Set the kinematic viscosity in the system operator to be the original
      // viscosity
      viscosity_model->set_kinematic_viscosity(viscosity_end);

      // Reset kinematic viscosity to simulation parameters
      viscosity_model->set_kinematic_viscosity(viscosity_end);
      this->system_operator->set_kinematic_viscosity(viscosity_end);

      if ((this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::lsmg) ||
          (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg))
        {
          auto mg_operators = this->gmg_preconditioner->get_mg_operators();
          for (unsigned int level = mg_operators.min_level();
               level <= mg_operators.max_level();
               level++)
            {
              mg_operators[level]->set_kinematic_viscosity(viscosity_end);
            }
        }
    }
  else if (initial_condition_type == Parameters::InitialConditionType::ramp)
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
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Gather viscosity final parameters
      const double viscosity_end = viscosity_model->get_kinematic_viscosity();

      // Kinematic viscosity ramp parameters
      const int n_iter_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .n_iter;
      double kinematic_viscosity =
        n_iter_viscosity > 0 ?
          this->simulation_parameters.initial_condition->ramp.ramp_viscosity
            .kinematic_viscosity_init :
          viscosity_end;
      const double alpha_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .alpha;

      viscosity_model->set_kinematic_viscosity(kinematic_viscosity);

      // Ramp on kinematic viscosity
      for (int i = 0; i < n_iter_viscosity; ++i)
        {
          this->pcout << std::setprecision(4)
                      << "********* Solution for kinematic viscosity = " +
                           std::to_string(kinematic_viscosity) + " *********"
                      << std::endl;

          viscosity_model->set_kinematic_viscosity(kinematic_viscosity);

          // Set the kinematic viscosity in the system operator to be the
          // temporary viscosity
          this->system_operator->set_kinematic_viscosity(
            this->simulation_parameters.physical_properties_manager
              .get_kinematic_viscosity_scale());

          // Set temporary viscosity in the multigrid operators
          if ((this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .preconditioner ==
               Parameters::LinearSolver::PreconditionerType::lsmg) ||
              (this->simulation_parameters.linear_solver
                 .at(PhysicsID::fluid_dynamics)
                 .preconditioner ==
               Parameters::LinearSolver::PreconditionerType::gcmg))
            {
              // Create the mg operators if they do not exist to be able
              // to change the viscosity for all of them
              if (!gmg_preconditioner)
                gmg_preconditioner =
                  std::make_shared<MFNavierStokesPreconditionGMG<dim>>(
                    this->simulation_parameters,
                    this->dof_handler,
                    this->dof_handler_fe_q_iso_q1,
                    this->mapping,
                    this->cell_quadrature,
                    this->forcing_function,
                    this->simulation_control,
                    this->fe);

              auto mg_operators = this->gmg_preconditioner->get_mg_operators();
              for (unsigned int level = mg_operators.min_level();
                   level <= mg_operators.max_level();
                   level++)
                {
                  mg_operators[level]->set_kinematic_viscosity(
                    this->simulation_parameters.physical_properties_manager
                      .get_kinematic_viscosity_scale());
                }
            }

          this->simulation_control->set_assembly_method(
            Parameters::SimulationControl::TimeSteppingMethod::steady);
          // Solve the problem with the temporary viscosity
          PhysicsSolver<LinearAlgebra::distributed::Vector<double>>::
            solve_non_linear_system(false);
          this->finish_time_step();

          kinematic_viscosity +=
            alpha_viscosity * (viscosity_end - kinematic_viscosity);
        }
      // Reset kinematic viscosity to simulation parameters
      viscosity_model->set_kinematic_viscosity(viscosity_end);
      this->system_operator->set_kinematic_viscosity(viscosity_end);

      if ((this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::lsmg) ||
          (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg))
        {
          auto mg_operators = this->gmg_preconditioner->get_mg_operators();
          for (unsigned int level = mg_operators.min_level();
               level <= mg_operators.max_level();
               level++)
            {
              mg_operators[level]->set_kinematic_viscosity(viscosity_end);
            }
        }

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
MFNavierStokesSolver<dim>::assemble_system_matrix()
{
  // Required for compilation but not used for matrix free solvers.
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");
}

template <int dim>
void
MFNavierStokesSolver<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_operator->evaluate_residual(this->system_rhs,
                                           this->evaluation_point);

  this->system_rhs *= -1.0;

  // Provide residual to simulation control for stopping criterion when using
  // steady bdf
  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}

template <int dim>
void
MFNavierStokesSolver<dim>::update_multiphysics_time_average_solution()
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

      this->multiphysics->set_time_average_solution(
        PhysicsID::fluid_dynamics, &this->multiphysics_average_velocities);
    }
}

template <int dim>
void
MFNavierStokesSolver<dim>::calculate_time_derivative_previous_solutions()
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
MFNavierStokesSolver<dim>::setup_GMG()
{
  TimerOutput::Scope t(this->computing_timer, "Setup GMG");

  if (!gmg_preconditioner)
    gmg_preconditioner = std::make_shared<MFNavierStokesPreconditionGMG<dim>>(
      this->simulation_parameters,
      this->dof_handler,
      this->dof_handler_fe_q_iso_q1,
      this->mapping,
      this->cell_quadrature,
      this->forcing_function,
      this->simulation_control,
      this->fe);

  gmg_preconditioner->initialize(this->simulation_control,
                                 this->flow_control,
                                 this->present_solution,
                                 this->time_derivative_previous_solutions);
}

template <int dim>
void
MFNavierStokesSolver<dim>::setup_ILU()
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
MFNavierStokesSolver<dim>::print_mg_setup_times()
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

      // Reset timers if output is set to every iteration
      this->gmg_preconditioner->mg_setup_timer.reset();
      this->gmg_preconditioner->mg_vmult_timer.reset();
    }
}

template <int dim>
void
MFNavierStokesSolver<dim>::update_solutions_for_multiphysics()
{
  TimerOutput::Scope t(this->computing_timer,
                       "Update solutions for multiphysics");

  // Provide the fluid dynamics dof_handler to the multiphysics interface
  this->multiphysics->set_dof_handler(PhysicsID::fluid_dynamics,
                                      &this->dof_handler);



  // Convert the present solution to multiphysics vector type and provide it to
  // the multiphysics interface
  TrilinosWrappers::MPI::Vector temp_solution(this->locally_owned_dofs,
                                              this->mpi_communicator);

  this->present_solution.update_ghost_values();
  convert_vector_dealii_to_trilinos(temp_solution, this->present_solution);
  multiphysics_present_solution = temp_solution;

  this->multiphysics->set_solution(PhysicsID::fluid_dynamics,
                                   &this->multiphysics_present_solution);

  // Convert the previous solutions to multiphysics vector type and provide them
  // to the multiphysics interface
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

  this->multiphysics->set_previous_solutions(
    PhysicsID::fluid_dynamics, &this->multiphysics_previous_solutions);
}

template <int dim>
void
MFNavierStokesSolver<dim>::define_non_zero_constraints()
{
  double time = this->simulation_control->get_current_time();
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  // Non-zero constraints
  auto &nonzero_constraints = this->get_nonzero_constraints();
  {
    nonzero_constraints.clear();
    nonzero_constraints.reinit(this->locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);
    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->simulation_parameters.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              nonzero_constraints,
              *this->mapping);
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .u.set_time(time);
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .v.set_time(time);
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .w.set_time(time);
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .u,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .v,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .w),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
              this->simulation_parameters.boundary_conditions
                .periodic_direction[i_bc],
              nonzero_constraints);
          }
      }
  }



  this->establish_solid_domain(true);

  nonzero_constraints.close();
}

template <int dim>
void
MFNavierStokesSolver<dim>::define_zero_constraints()
{
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  this->zero_constraints.clear();
  this->locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(this->dof_handler);
  this->zero_constraints.reinit(this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> no_normal_flux_boundaries;
          no_normal_flux_boundaries.insert(
            this->simulation_parameters.boundary_conditions.id[i_bc]);
          VectorTools::compute_no_normal_flux_constraints(
            this->dof_handler,
            0,
            no_normal_flux_boundaries,
            this->zero_constraints,
            *this->mapping);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            this->simulation_parameters.boundary_conditions.id[i_bc],
            this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
            this->simulation_parameters.boundary_conditions
              .periodic_direction[i_bc],
            this->zero_constraints);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::pressure)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::function_weak)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::partial_slip)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::outlet)
        {
          /*do nothing*/
        }
      else
        {
          VectorTools::interpolate_boundary_values(
            *this->mapping,
            this->dof_handler,
            this->simulation_parameters.boundary_conditions.id[i_bc],
            dealii::Functions::ZeroFunction<dim>(dim + 1),
            this->zero_constraints,
            this->fe->component_mask(velocities));
        }
    }

  this->establish_solid_domain(false);

  this->zero_constraints.close();
}

template <int dim>
void
MFNavierStokesSolver<dim>::setup_preconditioner()
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
MFNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
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
  else
    AssertThrow(false, ExcMessage("This solver is not allowed"));
  this->rescale_pressure_dofs_in_newton_update();
}

template <int dim>
void
MFNavierStokesSolver<dim>::assemble_L2_projection()
{
  // TODO
}

template <int dim>
void
MFNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
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

template class MFNavierStokesSolver<2>;
template class MFNavierStokesSolver<3>;
