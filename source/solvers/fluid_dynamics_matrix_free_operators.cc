// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/fluid_dynamics_matrix_free_operators.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/mg_tools.h>

/**
 * @brief Creates and fills a table that works as bool dof mask object
 * needed for the sparsity pattern and computation of the system matrix
 * of the coarse level in case FE_Q_iso_Q1 elements are used.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam spacedim Space dimension.
 * @param[in] fe Type of finite element used.
 * @param[in] quadrature Quadrature being used.
 * @return  Table<2, bool> Boolean mask with pairs of dof indices.
 */
template <int dim, int spacedim>
const Table<2, bool>
create_bool_dof_mask(const FiniteElement<dim, spacedim> &fe,
                     const Quadrature<dim>              &quadrature)
{
  const auto compute_scalar_bool_dof_mask = [&quadrature](const auto &fe) {
    Table<2, bool>           bool_dof_mask(fe.dofs_per_cell, fe.dofs_per_cell);
    MappingQ1<dim, spacedim> mapping;
    FEValues<dim>            fe_values(mapping, fe, quadrature, update_values);

    Triangulation<dim, spacedim> tria;
    GridGenerator::hyper_cube(tria);

    fe_values.reinit(tria.begin());
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
      for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
        {
          double sum = 0;
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            sum += fe_values.shape_value(i, q) * fe_values.shape_value(j, q);
          if (sum != 0)
            bool_dof_mask(i, j) = true;
        }

    return bool_dof_mask;
  };

  Table<2, bool> bool_dof_mask(fe.dofs_per_cell, fe.dofs_per_cell);

  if (fe.n_components() == 1)
    {
      bool_dof_mask = compute_scalar_bool_dof_mask(fe);
    }
  else
    {
      const auto scalar_bool_dof_mask =
        compute_scalar_bool_dof_mask(fe.base_element(0));

      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
          if (scalar_bool_dof_mask[fe.system_to_component_index(i).second]
                                  [fe.system_to_component_index(j).second])
            bool_dof_mask[i][j] = true;
    }


  return bool_dof_mask;
}



template <int dim, typename number>
NavierStokesOperatorBase<dim, number>::NavierStokesOperatorBase()
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , timer(this->pcout, TimerOutput::never, TimerOutput::wall_times)
{}

template <int dim, typename number>
NavierStokesOperatorBase<dim, number>::NavierStokesOperatorBase(
  const Mapping<dim>                               &mapping,
  const DoFHandler<dim>                            &dof_handler,
  const AffineConstraints<number>                  &constraints,
  const Quadrature<dim>                            &quadrature,
  const std::shared_ptr<Function<dim>>              forcing_function,
  const std::shared_ptr<PhysicalPropertiesManager> &physical_properties_manager,
  const StabilizationType                           stabilization,
  const unsigned int                                mg_level,
  const std::shared_ptr<SimulationControl>         &simulation_control,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const bool                                          &enable_hessians_jacobian,
  const bool                                          &enable_hessians_residual,
  const bool                                          &enable_mortar)
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , timer(this->pcout, TimerOutput::never, TimerOutput::wall_times)
{
  this->reinit(mapping,
               dof_handler,
               constraints,
               quadrature,
               forcing_function,
               physical_properties_manager,
               stabilization,
               mg_level,
               simulation_control,
               boundary_conditions,
               enable_hessians_jacobian,
               enable_hessians_residual,
               enable_mortar);
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::reinit(
  const Mapping<dim>                               &mapping,
  const DoFHandler<dim>                            &dof_handler,
  const AffineConstraints<number>                  &constraints,
  const Quadrature<dim>                            &quadrature,
  const std::shared_ptr<Function<dim>>              forcing_function,
  const std::shared_ptr<PhysicalPropertiesManager> &physical_properties_manager,
  const StabilizationType                           stabilization,
  const unsigned int                                mg_level,
  const std::shared_ptr<SimulationControl>         &simulation_control,
  const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
  const bool                                          &enable_hessians_jacobian,
  const bool                                          &enable_hessians_residual,
  const bool                                          &enable_mortar)
{
  this->enable_face_terms = false;

  this->boundary_conditions = boundary_conditions;

  this->enable_mortar = enable_mortar;

  this->system_matrix.clear();
  this->constraints.copy_from(constraints);

  typename MatrixFree<dim, number>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    (update_values | update_gradients | update_JxW_values |
     update_quadrature_points | update_hessians);
  additional_data.mg_level = mg_level;

  // Loop over all boundary conditions to establish if one of them has face
  // terms, if yes, enable them and add appropriate flags in the matrix-free
  // additional data
  for (auto const &[id, type] : boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::function_weak)
        {
          this->enable_face_terms = true;
          additional_data.mapping_update_flags_boundary_faces =
            update_values | update_gradients | update_quadrature_points |
            update_JxW_values | update_normal_vectors;
        }
      else if (type == BoundaryConditions::BoundaryType::outlet)
        {
          this->enable_face_terms = true;
          additional_data.mapping_update_flags_boundary_faces =
            update_values | update_quadrature_points | update_JxW_values |
            update_normal_vectors;
        }
    }

  matrix_free.reinit(
    mapping, dof_handler, this->constraints, quadrature, additional_data);

  this->fe_degree = dof_handler.get_fe().degree;

  this->forcing_function = forcing_function;

  this->properties_manager = physical_properties_manager;

  if (stabilization ==
        Parameters::Stabilization::NavierStokesStabilization::pspg_supg ||
      stabilization ==
        Parameters::Stabilization::NavierStokesStabilization::gls)
    {
      this->stabilization = stabilization;
    }
  else
    throw std::runtime_error(
      "Only SUPG/PSPG and GLS stabilization is supported at the moment.");

  this->simulation_control = simulation_control;

  this->enable_hessians_jacobian = enable_hessians_jacobian;

  this->enable_hessians_residual = enable_hessians_residual;

  this->compute_element_size();

  this->compute_forcing_term();

  bool_dof_mask = create_bool_dof_mask(dof_handler.get_fe(), quadrature);

  constrained_indices.clear();
  for (auto i : this->matrix_free.get_constrained_dofs())
    constrained_indices.push_back(i);
  constrained_values.resize(constrained_indices.size());

  if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
    {
      std::vector<types::global_dof_index> interface_indices;
      IndexSet                             refinement_edge_indices;
      refinement_edge_indices = get_refinement_edges(this->matrix_free);

      interface_indices = refinement_edge_indices.get_index_vector();

      edge_constrained_indices.clear();
      edge_constrained_indices.reserve(interface_indices.size());
      edge_constrained_values.resize(interface_indices.size());
      const IndexSet &locally_owned =
        this->matrix_free.get_dof_handler().locally_owned_mg_dofs(
          this->matrix_free.get_mg_level());

      for (const auto i : interface_indices)
        if (locally_owned.is_element(i))
          edge_constrained_indices.push_back(locally_owned.index_within_set(i));

      this->has_edge_constrained_indices =
        Utilities::MPI::max(edge_constrained_indices.size(),
                            dof_handler.get_mpi_communicator()) > 0;

      if (this->has_edge_constrained_indices)
        {
          edge_constrained_cell.resize(matrix_free.n_cell_batches(), false);

          VectorType temp;
          matrix_free.initialize_dof_vector(temp);

          for (const auto i : edge_constrained_indices)
            temp.local_element(i) = 1.0;

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
NavierStokesOperatorBase<dim, number>::compute_element_size()
{
  const unsigned int n_cells =
    matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();
  element_size.resize(n_cells);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (auto lane = 0u;
           lane < matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          const double h_k =
            matrix_free.get_cell_iterator(cell, lane)->measure();

          element_size[cell][lane] = compute_cell_diameter<dim>(h_k, fe_degree);
        }
    }
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::compute_forcing_term()
{
  if (this->forcing_function)
    {
      const unsigned int n_cells =
        matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();

      FECellIntegrator integrator(matrix_free);

      forcing_terms.reinit(n_cells, integrator.n_q_points);

      for (unsigned int cell = 0; cell < n_cells; ++cell)
        {
          integrator.reinit(cell);
          for (const auto q : integrator.quadrature_point_indices())
            {
              const Point<dim, VectorizedArray<number>> point_batch =
                integrator.quadrature_point(q);
              forcing_terms(cell, q) =
                evaluate_function<dim, number, dim>(*(this->forcing_function),
                                                    point_batch);
            }
        }
    }
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::compute_buoyancy_term(
  const VectorType      &temperature_solution,
  const DoFHandler<dim> &temperature_dof_handler)
{
  const auto thermal_expansion_model =
    this->properties_manager->get_thermal_expansion();

  const double reference_temperature =
    this->properties_manager->get_reference_temperature();

  const unsigned int n_cells =
    matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();

  FEValues<dim> fe_values(*(this->matrix_free.get_mapping_info().mapping),
                          temperature_dof_handler.get_fe(),
                          this->matrix_free.get_quadrature(),
                          update_values);

  buoyancy_term.reinit(n_cells, fe_values.n_quadrature_points);

  std::vector<number> cell_temperature_solution(fe_values.n_quadrature_points);
  std::vector<double> thermal_expansion(fe_values.n_quadrature_points);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      for (auto lane = 0u;
           lane < matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          fe_values.reinit(
            matrix_free.get_cell_iterator(cell, lane)
              ->as_dof_handler_iterator(temperature_dof_handler));

          fe_values.get_function_values(temperature_solution,
                                        cell_temperature_solution);

          thermal_expansion_model->vector_value({}, thermal_expansion);

          for (const auto q : fe_values.quadrature_point_indices())
            {
              // The following term is precomputed: (1 - β (T -Tref))
              const auto scaling =
                1 - thermal_expansion[q] *
                      (cell_temperature_solution[q] - reference_temperature);

              // Store in table the total buoyancy term given as:
              // (1 - β (T -Tref))g where the g has already been precomputed in
              // compute_forcing_term() function.
              for (unsigned int c = 0; c < dim; ++c)
                buoyancy_term[cell][q][c][lane] =
                  forcing_terms[cell][q][c][lane] * scaling;
            }
        }
    }
}


template <int dim, typename number>
types::global_dof_index
NavierStokesOperatorBase<dim, number>::m() const
{
  if (this->matrix_free.get_mg_level() != numbers::invalid_unsigned_int)
    return this->matrix_free.get_dof_handler().n_dofs(
      this->matrix_free.get_mg_level());
  else
    return this->matrix_free.get_dof_handler().n_dofs();
}

template <int dim, typename number>
number
NavierStokesOperatorBase<dim, number>::el(unsigned int, unsigned int) const
{
  Assert(false, ExcNotImplemented());
  return 0;
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::clear()
{
  matrix_free.clear();
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::initialize_dof_vector(
  VectorType &vec) const
{
  matrix_free.initialize_dof_vector(vec);
}

template <int dim, typename number>
const std::shared_ptr<const Utilities::MPI::Partitioner> &
NavierStokesOperatorBase<dim, number>::get_vector_partitioner() const
{
  return matrix_free.get_vector_partitioner();
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::vmult(VectorType       &dst,
                                             const VectorType &src) const
{
  this->timer.enter_subsection("operator::vmult");

  // save values for edge constrained dofs and set them to 0 in src vector
  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    {
      edge_constrained_values[i] =
        src.local_element(edge_constrained_indices[i]);

      const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
        .local_element(edge_constrained_indices[i]) = 0.;
    }

  // If mortar is enabled, we need to make the ghost values available so that
  // they can be accessed by the mortar operator
  if (this->enable_mortar)
    src.update_ghost_values();

  if (this->enable_face_terms)
    this->matrix_free.loop(
      &NavierStokesOperatorBase::do_cell_integral_range,
      &NavierStokesOperatorBase::do_internal_face_integral_range,
      &NavierStokesOperatorBase::do_boundary_face_integral_range<false>,
      this,
      dst,
      src,
      true);
  else
    this->matrix_free.cell_loop(
      &NavierStokesOperatorBase::do_cell_integral_range, this, dst, src, true);

  // If enabled, add mortar entries
  if (this->enable_mortar)
    {
      this->mortar_coupling_operator_mf->vmult_add(dst, src);
      dst.compress(VectorOperation::add);
      src.zero_out_ghost_values();
    }

  // copy constrained dofs from src to dst (corresponding to diagonal
  // entries with value 1.0)
  for (const auto &constrained_index : constrained_indices)
    dst.local_element(constrained_index) = src.local_element(constrained_index);

  // restoring edge constrained dofs in src and copy the values to dst
  // (corresponding to diagonal entries with value 1.0)
  for (unsigned int i = 0; i < edge_constrained_indices.size(); ++i)
    {
      const_cast<LinearAlgebra::distributed::Vector<number> &>(src)
        .local_element(edge_constrained_indices[i]) =
        edge_constrained_values[i];
      dst.local_element(edge_constrained_indices[i]) =
        edge_constrained_values[i];
    }

  this->timer.leave_subsection("operator::vmult");
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::Tvmult(VectorType       &dst,
                                              const VectorType &src) const
{
  this->vmult(dst, src);
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::vmult_interface_down(
  VectorType       &dst,
  VectorType const &src) const
{
  if (this->enable_face_terms)
    this->matrix_free.loop(
      &NavierStokesOperatorBase::do_cell_integral_range,
      &NavierStokesOperatorBase::do_internal_face_integral_range,
      &NavierStokesOperatorBase::do_boundary_face_integral_range<false>,
      this,
      dst,
      src,
      true);
  else
    this->matrix_free.cell_loop(
      &NavierStokesOperatorBase::do_cell_integral_range, this, dst, src, true);

  // set constrained dofs as the sum of current dst value and src value
  for (const auto i : constrained_indices)
    dst.local_element(i) = src.local_element(i);
}



template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::vmult_interface_up(
  VectorType       &dst,
  VectorType const &src) const
{
  this->timer.enter_subsection("operator::vmult_interface_up");

  if (has_edge_constrained_indices == false)
    {
      dst = number(0.);

      this->timer.leave_subsection("operator::vmult_interface_up");
      return;
    }

  dst = 0.0;

  // make a copy of src vector and set everything to 0 except edge
  // constrained dofs
  VectorType src_cpy;
  src_cpy.reinit(src, /*omit_zeroing_entries=*/false);

  for (const auto i : edge_constrained_indices)
    src_cpy.local_element(i) = src.local_element(i);

  // do loop with copy of src
  if (this->enable_face_terms)
    this->matrix_free.loop(
      &NavierStokesOperatorBase::do_cell_integral_range,
      &NavierStokesOperatorBase::do_internal_face_integral_range,
      &NavierStokesOperatorBase::do_boundary_face_integral_range<false>,
      this,
      dst,
      src_cpy,
      true);
  else
    this->matrix_free.cell_loop(
      &NavierStokesOperatorBase::do_cell_integral_range,
      this,
      dst,
      src_cpy,
      false);

  this->timer.leave_subsection("operator::vmult_interface_up");
}



template <int dim, typename number>
const DynamicSparsityPattern &
NavierStokesOperatorBase<dim, number>::get_sparsity_pattern() const
{
  return dsp;
}



template <int dim, typename number>
const TrilinosWrappers::SparseMatrix &
NavierStokesOperatorBase<dim, number>::get_system_matrix() const
{
  this->timer.enter_subsection("operator::get_system_matrix");

  if (system_matrix.m() == 0 && system_matrix.n() == 0)
    {
      const auto &dof_handler = this->matrix_free.get_dof_handler();

      const unsigned int mg_level = this->matrix_free.get_mg_level();

      IndexSet locally_relevant_dofs;
      IndexSet locally_owned_dofs;

      if (mg_level != numbers::invalid_unsigned_int)
        {
          locally_relevant_dofs =
            DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                          mg_level);
          locally_owned_dofs = dof_handler.locally_owned_mg_dofs(mg_level);
        }
      else
        {
          locally_owned_dofs = dof_handler.locally_owned_dofs();
          locally_relevant_dofs =
            DoFTools::extract_locally_relevant_dofs(dof_handler);
        }

      this->dsp.reinit(locally_relevant_dofs.size(),
                       locally_relevant_dofs.size(),
                       locally_relevant_dofs);

      if (mg_level != numbers::invalid_unsigned_int)
        {
          // the following code does the same as
          // MGTools::make_sparsity_pattern() but also
          // considers bool_dof_mask for FE_Q_iso_Q1
          std::vector<types::global_dof_index> dofs_on_this_cell;

          for (const auto &cell : dof_handler.cell_iterators_on_level(mg_level))
            if (cell->is_locally_owned_on_level())
              {
                const unsigned int dofs_per_cell =
                  dof_handler.get_fe().n_dofs_per_cell();
                dofs_on_this_cell.resize(dofs_per_cell);
                cell->get_mg_dof_indices(dofs_on_this_cell);

                constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                        dsp,
                                                        false,
                                                        bool_dof_mask);
              }
        }
      else
        {
          // the following code does the same as
          // DoFTools::make_sparsity_pattern() but also
          // considers bool_dof_mask for FE_Q_iso_Q1
          std::vector<types::global_dof_index> dofs_on_this_cell;

          for (const auto &cell : dof_handler.active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int dofs_per_cell =
                  cell->get_fe().n_dofs_per_cell();
                dofs_on_this_cell.resize(dofs_per_cell);
                cell->get_dof_indices(dofs_on_this_cell);

                constraints.add_entries_local_to_global(dofs_on_this_cell,
                                                        dsp,
                                                        false,
                                                        bool_dof_mask);
              }
        }

      // If mortar is enabled, add sparsity pattern entries
      if (this->enable_mortar)
        this->mortar_coupling_operator_mf->add_sparsity_pattern_entries(dsp);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs,
        dof_handler.get_triangulation().get_mpi_communicator(),
        locally_relevant_dofs);

      system_matrix.reinit(
        locally_owned_dofs,
        locally_owned_dofs,
        dsp,
        dof_handler.get_triangulation().get_mpi_communicator());
    }

  system_matrix = 0.0;

  // the following code does the same as
  // MatrixFreeTools::compute_matrix() but with
  // minor corrections for FE_Q_iso_Q1

  const auto &lexicographic_numbering =
    matrix_free.get_shape_info().lexicographic_numbering;

  unsigned int cell   = numbers::invalid_unsigned_int;
  unsigned int face   = numbers::invalid_unsigned_int;
  unsigned int column = numbers::invalid_unsigned_int;

  std::function<void(FEFaceIntegrator &)> boundary_function;

  if (enable_face_terms)
    {
      boundary_function = [&](auto &integrator) {
        if (face != integrator.get_current_cell_index())
          {
            face   = integrator.get_current_cell_index();
            column = 0;
          }

        this->do_boundary_face_integral_local<false>(integrator);

        // remove spurious entries for FE_Q_iso_Q1
        for (unsigned int i = 0; i < lexicographic_numbering.size(); ++i)
          {
            if (!bool_dof_mask[lexicographic_numbering[i]]
                              [lexicographic_numbering[column]])
              integrator.begin_dof_values()[i] = 0.0;
          }


        column++;
      };
    }

  MatrixFreeTools::
    compute_matrix<dim, -1, 0, dim + 1, number, VectorizedArray<number>>(
      matrix_free,
      constraints,
      system_matrix,
      [&](auto &integrator) {
        if (cell != integrator.get_current_cell_index())
          {
            cell   = integrator.get_current_cell_index();
            column = 0;
          }

        do_cell_integral_local(integrator);

        // remove spurious entries for FE_Q_iso_Q1
        for (unsigned int i = 0; i < lexicographic_numbering.size(); ++i)
          if (!bool_dof_mask[lexicographic_numbering[i]]
                            [lexicographic_numbering[column]])
            integrator.begin_dof_values()[i] = 0.0;

        column++;
      },
      {},
      boundary_function);

  // If mortar is enabled, add system matrix entries
  if (this->enable_mortar)
    {
      this->mortar_coupling_operator_mf->add_system_matrix_entries(
        system_matrix);
      system_matrix.compress(VectorOperation::add);
    }

  // make sure that diagonal entries related to constrained dofs
  // have a value of 1.0 (note this is consistent to vmult() and
  // compute_inverse_diagonal())
  for (const auto &local_row : constrained_indices)
    {
      const auto global_row =
        get_vector_partitioner()->local_to_global(local_row);
      system_matrix.set(global_row, global_row, 1.0);
    }

  for (const auto &local_row : edge_constrained_indices)
    {
      const auto global_row =
        get_vector_partitioner()->local_to_global(local_row);
      system_matrix.set(global_row, global_row, 1.0);
    }

  system_matrix.compress(VectorOperation::insert);

  this->timer.leave_subsection("operator::get_system_matrix");

  return this->system_matrix;
}

template <int dim, typename number>
const MatrixFree<dim, number> &
NavierStokesOperatorBase<dim, number>::get_system_matrix_free() const
{
  return this->matrix_free;
}

template <int dim, typename number>
const AlignedVector<VectorizedArray<number>>
NavierStokesOperatorBase<dim, number>::get_element_size() const
{
  return this->element_size;
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::
  evaluate_non_linear_term_and_calculate_tau(const VectorType &newton_step)
{
  this->timer.enter_subsection("operator::evaluate_non_linear_term");

  const unsigned int n_cells          = matrix_free.n_cell_batches();
  const unsigned int n_inner_faces    = matrix_free.n_inner_face_batches();
  const unsigned int n_boundary_faces = matrix_free.n_boundary_face_batches();

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  FECellIntegrator integrator(matrix_free);
  FEFaceIntegrator face_integrator(matrix_free, true, 0);


  // 1. Precompute values on cells:

  // Set appropriate size for tables
  nonlinear_previous_values.reinit(n_cells, integrator.n_q_points);
  nonlinear_previous_gradient.reinit(n_cells, integrator.n_q_points);
  nonlinear_previous_hessian_diagonal.reinit(n_cells, integrator.n_q_points);
  stabilization_parameter.reinit(n_cells, integrator.n_q_points);
  stabilization_parameter_lsic.reinit(n_cells, integrator.n_q_points);
  kinematic_viscosity_vector.reinit(n_cells, integrator.n_q_points);
  grad_kinematic_viscosity_shear_rate.reinit(n_cells, integrator.n_q_points);
  kinematic_viscosity_gradient.reinit(n_cells, integrator.n_q_points);
  previous_shear_rate.reinit(n_cells, integrator.n_q_points);
  previous_shear_rate_magnitude.reinit(n_cells, integrator.n_q_points);

  // Define 1/dt if the simulation is transient
  double sdt = 0.0;

  bool transient =
    (is_bdf(this->simulation_control->get_assembly_method())) ? true : false;

  if (transient)
    {
      const auto time_steps_vector =
        this->simulation_control->get_time_steps_vector();
      const double dt = time_steps_vector[0];
      sdt             = 1. / dt;
    }

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      integrator.reinit(cell);
      integrator.read_dof_values_plain(newton_step);

      if (this->enable_hessians_jacobian)
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients |
                            EvaluationFlags::hessians);
      else
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients);

      // Get previously calculated element size needed for tau
      const auto h = integrator.read_cell_data(this->get_element_size());

      for (const auto q : integrator.quadrature_point_indices())
        {
          nonlinear_previous_values(cell, q)   = integrator.get_value(q);
          nonlinear_previous_gradient(cell, q) = integrator.get_gradient(q);

          if (this->enable_hessians_jacobian)
            nonlinear_previous_hessian_diagonal(cell, q) =
              integrator.get_hessian_diagonal(q);

          // Calculate tau
          VectorizedArray<number> u_mag_squared = 1e-12;
          for (unsigned int k = 0; k < dim; ++k)
            u_mag_squared +=
              Utilities::fixed_power<2>(integrator.get_value(q)[k]);

          stabilization_parameter(cell, q) =
            1. / std::sqrt(Utilities::fixed_power<2>(sdt) +
                           4. * u_mag_squared / h / h +
                           9. * Utilities::fixed_power<2>(
                                  4. * kinematic_viscosity / (h * h)));

          stabilization_parameter_lsic(cell, q) =
            std::sqrt(u_mag_squared) * h * 0.5;
        }

      // Compute kinematic viscosity-related entries for non-Newtonian fluids
      // according to the rheological model
      if (this->properties_manager->is_non_newtonian())
        {
          typename FECellIntegrator::gradient_type shear_rate;
          typename FECellIntegrator::value_type    grad_shear_rate = {};
          VectorizedArray<number>                  shear_rate_magnitude;

          for (const auto q : integrator.quadrature_point_indices())
            {
              typename FECellIntegrator::gradient_type gradient =
                integrator.get_gradient(q);

              typename FECellIntegrator::hessian_type hessian =
                integrator.get_hessian(q);

              // Calculate shear rate
              for (unsigned int i = 0; i < dim; ++i)
                {
                  for (unsigned int j = 0; j < dim; ++j)
                    {
                      shear_rate[i][j] = gradient[i][j] + gradient[j][i];
                    }
                }

              // Store shear rate in the previous_shear_rate
              this->previous_shear_rate[cell][q] = shear_rate;

              // Calculate shear rate magnitude
              shear_rate_magnitude = shear_rate * shear_rate;

              shear_rate_magnitude = std::max(sqrt(0.5 * shear_rate_magnitude),
                                              VectorizedArray<number>(1e-12));

              // Store shear rate magnitude in the previous_shear_rate
              this->previous_shear_rate_magnitude[cell][q] =
                shear_rate_magnitude;

              // Compute gradient of shear rate
              // ∂d gamma_dot = 1/(2*gamma_dot)*(∂iuj + ∂jui) * ∂d(∂iuj + ∂jui)
              for (int d = 0; d < dim; ++d)
                {
                  grad_shear_rate[d] = 0.;
                  for (unsigned int i = 0; i < dim; ++i)
                    {
                      for (unsigned int k = 0; k < dim; ++k)
                        {
                          grad_shear_rate[d] +=
                            VectorizedArray<number>(0.5) *
                            (gradient[i][k] + gradient[k][i]) *
                            (hessian[i][d][k] + hessian[k][d][i]) /
                            shear_rate_magnitude;
                        }
                    }
                }

              // Store the shear rate magnitude in the appropriate data
              // structure needed for the set field vector function

              VectorizedArray<number> viscosity;
              VectorizedArray<number> grad_viscosity_shear_rate;

              for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
                {
                  // Get information from physical properties class
                  std::map<field, double> fields;
                  fields[field::shear_rate] = shear_rate_magnitude[v];

                  viscosity[v] =
                    this->properties_manager->get_rheology()->value(fields);

                  grad_viscosity_shear_rate[v] =
                    this->properties_manager->get_rheology()->jacobian(
                      fields, field::shear_rate);
                }

              grad_kinematic_viscosity_shear_rate(cell, q) =
                grad_viscosity_shear_rate;

              kinematic_viscosity_gradient(cell, q) =
                grad_shear_rate * grad_viscosity_shear_rate;

              kinematic_viscosity_vector(cell, q) = viscosity;

              // Recalculate stabilization parameter using kinematic
              // viscosity vector
              VectorizedArray<number> u_mag_squared = 1e-12;
              for (unsigned int k = 0; k < dim; ++k)
                u_mag_squared +=
                  Utilities::fixed_power<2>(integrator.get_value(q)[k]);

              stabilization_parameter(cell, q) =
                1. /
                std::sqrt(
                  Utilities::fixed_power<2>(sdt) + 4. * u_mag_squared / h / h +
                  9. * Utilities::fixed_power<2>(4. * viscosity / (h * h)));
            }
        }
    }

  // 2. Precompute values on faces if required

  if (enable_face_terms)
    {
      effective_beta_face.reinit(n_boundary_faces);
      face_target_velocity.reinit(n_boundary_faces, face_integrator.n_q_points);
      face_nonlinear_previous_values.reinit(n_boundary_faces,
                                            face_integrator.n_q_points);

      // Boundary faces are padded after the internal faces
      // Consequently, we need to offset the indices to avoid storing too
      // much information regarding the boundary faces (e.g. too large
      // vectors)
      for (unsigned int face = n_inner_faces;
           face < n_inner_faces + n_boundary_faces;
           ++face)
        {
          face_integrator.reinit(face);

          // Check if the boundary condition is a weak Dirichlet boundary
          // condition or an outlet boundary condition, otherwise, no terms
          // need to be precomputed for faces
          if (this->boundary_conditions.type.at(
                face_integrator.boundary_id()) !=
                BoundaryConditions::BoundaryType::function_weak &&
              this->boundary_conditions.type.at(
                face_integrator.boundary_id()) !=
                BoundaryConditions::BoundaryType::outlet)
            continue;

          // We need to read the values for the outlet boundary condition
          if (this->boundary_conditions.type.at(
                face_integrator.boundary_id()) ==
              BoundaryConditions::BoundaryType::outlet)
            {
              face_integrator.read_dof_values_plain(newton_step);
              face_integrator.evaluate(
                EvaluationFlags::EvaluationFlags::values);
            }

          for (const auto q : face_integrator.quadrature_point_indices())
            {
              if (this->boundary_conditions.type.at(
                    face_integrator.boundary_id()) ==
                  BoundaryConditions::BoundaryType::function_weak)
                {
                  Point<dim, VectorizedArray<number>> point_batch =
                    face_integrator.quadrature_point(q);

                  Tensor<1, dim, VectorizedArray<number>> target_velocity_value;

                  target_velocity_value[0] = evaluate_function<dim, number>(
                    boundary_conditions.navier_stokes_functions
                      .at(face_integrator.boundary_id())
                      ->u,
                    point_batch);

                  target_velocity_value[1] = evaluate_function<dim, number>(
                    boundary_conditions.navier_stokes_functions
                      .at(face_integrator.boundary_id())
                      ->v,
                    point_batch);

                  if constexpr (dim == 3)
                    target_velocity_value[2] = evaluate_function<dim, number>(
                      boundary_conditions.navier_stokes_functions
                        .at(face_integrator.boundary_id())
                        ->w,
                      point_batch);

                  face_target_velocity(face - n_inner_faces, q) =
                    target_velocity_value;
                }
              else if (this->boundary_conditions.type.at(
                         face_integrator.boundary_id()) ==
                       BoundaryConditions::BoundaryType::outlet)
                {
                  face_nonlinear_previous_values[face - n_inner_faces][q] =
                    face_integrator.get_value(q);
                }
            }

          // Calculate the element size and use it to calculate the penalty
          // parameter with the beta factor given in the parameter file
          VectorizedArray<number> cell_size = 1.0;

          const auto [cell_it, _] = matrix_free.get_face_iterator(face, 0);

          cell_size = compute_cell_diameter<dim>(cell_it->measure(), fe_degree);

          const double beta =
            this->boundary_conditions.beta.at(face_integrator.boundary_id());

          effective_beta_face(face - n_inner_faces) =
            beta / std::pow(cell_size, static_cast<number>(fe_degree + 1));
        }
    }

  this->timer.leave_subsection("operator::evaluate_non_linear_term");
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::
  evaluate_time_derivative_previous_solutions(
    const VectorType &time_derivative_previous_solutions)
{
  this->timer.enter_subsection(
    "operator::evaluate_time_derivative_previous_solutions");

  const unsigned int n_cells = matrix_free.n_cell_batches();
  FECellIntegrator   integrator(matrix_free);

  time_derivatives_previous_solutions.reinit(n_cells, integrator.n_q_points);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      integrator.reinit(cell);
      integrator.read_dof_values_plain(time_derivative_previous_solutions);
      integrator.evaluate(EvaluationFlags::values);
      for (const auto q : integrator.quadrature_point_indices())
        time_derivatives_previous_solutions(cell, q) += integrator.get_value(q);
    }

  this->timer.leave_subsection(
    "operator::evaluate_time_derivative_previous_solutions");
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::evaluate_velocity_ale(
  const Mapping<dim>                             &mapping,
  const double                                    radius,
  const Point<dim>                                center_of_rotation,
  std::shared_ptr<Functions::ParsedFunction<dim>> rotor_angular_velocity)
{
  this->timer.enter_subsection("operator::evaluate_velocity_ale");

  const unsigned int n_cells =
    matrix_free.n_cell_batches() + matrix_free.n_ghost_cell_batches();

  FECellIntegrator integrator(matrix_free);
  matrix_free.update_mapping(mapping);

  // 1. Precompute values on cells:

  // Set appropriate size for tables
  velocity_ale.reinit(n_cells, integrator.n_q_points);

  // Get updated rotor angular velocity
  rotor_angular_velocity->set_time(
    this->simulation_control->get_current_time());
  const double rotor_angular_velocity_value =
    rotor_angular_velocity->value(Point<dim>());

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      integrator.reinit(cell);
      for (auto lane = 0u;
           lane < matrix_free.n_active_entries_per_cell_batch(cell);
           lane++)
        {
          // Compute radial distance from the current cell
          const auto cell_center =
            matrix_free.get_cell_iterator(cell, lane)->center();
          const auto radius_current = cell_center.distance(center_of_rotation);

          // Use prescribed rotor angular velocity only if cell is part of the
          // rotor
          double cell_rotor_angular_velocity;
          if (radius_current > radius)
            cell_rotor_angular_velocity = 0.0;
          else
            cell_rotor_angular_velocity = rotor_angular_velocity_value;

          // Compute linear velocity at quadrature points
          for (const auto q : integrator.quadrature_point_indices())
            {
              // Assumption in 3D case: rotation axis is in z
              // TODO generalize rotation axis
              const auto x = integrator.quadrature_point(q)[0][lane];
              const auto y = integrator.quadrature_point(q)[1][lane];

              this->velocity_ale[cell][q][0][lane] =
                -cell_rotor_angular_velocity * y;
              this->velocity_ale[cell][q][1][lane] =
                cell_rotor_angular_velocity * x;
            }
        }
    }

  this->timer.leave_subsection("operator::evaluate_velocity_ale");
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::update_beta_force(
  const Tensor<1, dim> &beta_force)
{
  for (unsigned int v = 0; v < VectorizedArray<number>::size(); ++v)
    {
      for (int d = 0; d < dim; ++d)
        this->beta_force[d][v] = beta_force[d];
    }
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::evaluate_residual(VectorType       &dst,
                                                         const VectorType &src)
{
  this->timer.enter_subsection("operator::evaluate_residual");

  if (enable_face_terms)
    this->matrix_free.loop(
      &NavierStokesOperatorBase::local_evaluate_residual,
      &NavierStokesOperatorBase::do_internal_face_integral_range,
      &NavierStokesOperatorBase::do_boundary_face_integral_range<true>,
      this,
      dst,
      src,
      true);
  else
    this->matrix_free.cell_loop(
      &NavierStokesOperatorBase::local_evaluate_residual, this, dst, src, true);

  this->timer.leave_subsection("operator::evaluate_residual");
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::compute_inverse_diagonal(
  VectorType &diagonal) const
{
  this->timer.enter_subsection("operator::compute_inverse_diagonal");

  std::function<void(FEFaceIntegrator &)> boundary_function;

  if ((enable_face_terms))
    boundary_function = [&](auto &integrator) {
      this->do_boundary_face_integral_local<false>(integrator);
    };

  matrix_free.initialize_dof_vector(diagonal);
  MatrixFreeTools::
    compute_diagonal<dim, -1, 0, dim + 1, number, VectorizedArray<number>>(
      matrix_free,
      diagonal,
      [&](auto &integrator) { this->do_cell_integral_local(integrator); },
      {},
      boundary_function);

  // If mortar is enabled, add diagonal entries
  if (this->enable_mortar)
    this->mortar_coupling_operator_mf->add_diagonal_entries(diagonal);

  for (const auto &i : edge_constrained_indices)
    diagonal.local_element(i) = 0.0;

  for (auto &i : diagonal)
    i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0;

  this->timer.leave_subsection("operator::compute_inverse_diagonal");
}

template <int dim, typename number>
void
NavierStokesOperatorBase<dim, number>::do_cell_integral_range(
  const MatrixFree<dim, number>               &matrix_free,
  VectorType                                  &dst,
  const VectorType                            &src,
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
void
NavierStokesOperatorBase<dim, number>::do_internal_face_integral_range(
  const MatrixFree<dim, number>               &matrix_free,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::pair<unsigned int, unsigned int> &range) const
{
  (void)matrix_free;
  (void)dst;
  (void)src;
  (void)range;

  // nothing to do
}

template <int dim, typename number>
template <bool assemble_residual>
void
NavierStokesOperatorBase<dim, number>::do_boundary_face_integral_range(
  const MatrixFree<dim, number>               &matrix_free,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::pair<unsigned int, unsigned int> &range) const
{
  FEFaceIntegrator phi(matrix_free, true, 0);

  for (auto cell = range.first; cell < range.second; ++cell)
    {
      phi.reinit(cell);

      if (assemble_residual)
        phi.read_dof_values_plain(src);
      else
        phi.read_dof_values(src);

      do_boundary_face_integral_local<assemble_residual>(phi);

      phi.distribute_local_to_global(dst);
    }
}


// Assembles face terms of boundary conditions for both the Jacobian matrix
// and the residual. There are two types implemented:
// 1. Weak Dirichlet boundary conditions using Nitsche's symmetric
// penalty method. It adds the following term:
// (v,β·(u-u_target)) - ν(v,∇δu·n) - ν(∇v·n,(u-u_target))
// 2. Outlet boundary conditions using the directional do-nothing method.
// It adds the following term: -(v,β(u·n)_·u) where (u·n)_=min(0,u·n)
template <int dim, typename number>
template <bool assemble_residual>
void
NavierStokesOperatorBase<dim, number>::do_boundary_face_integral_local(
  FEFaceIntegrator &integrator) const
{
  if (!enable_face_terms)
    return;

  // If the boundary condition is not in our list of boundary
  // conditions or the boundary condition that is set in the list is
  // not a weak function or outlet BC, there is nothing to do, so we set the
  // values to zero and return
  if (this->boundary_conditions.type.at(integrator.boundary_id()) !=
        BoundaryConditions::BoundaryType::function_weak &&
      this->boundary_conditions.type.at(integrator.boundary_id()) !=
        BoundaryConditions::BoundaryType::outlet)
    {
      const VectorizedArray<number> zero = 0.0;

      for (unsigned int i = 0; i < integrator.dofs_per_cell; ++i)
        integrator.begin_dof_values()[i] = zero;
      return;
    }

  const auto         face       = integrator.get_current_cell_index();
  const unsigned int face_index = face - matrix_free.n_inner_face_batches();

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  const auto penalty_parameter = this->effective_beta_face[face_index];

  if (this->boundary_conditions.type.at(integrator.boundary_id()) ==
      BoundaryConditions::BoundaryType::function_weak)
    {
      integrator.evaluate(EvaluationFlags::EvaluationFlags::values |
                          EvaluationFlags::EvaluationFlags::gradients);

      for (const auto q : integrator.quadrature_point_indices())
        {
          typename FEFaceIntegrator::value_type    value_result    = {};
          typename FEFaceIntegrator::gradient_type gradient_result = {};

          const auto normal_vector = integrator.normal_vector(q);

          auto       value    = integrator.get_value(q);
          const auto gradient = integrator.get_gradient(q);


          // If we are assembling the residual, substract the target velocity
          // from the velocity value.
          if constexpr (assemble_residual)
            for (int d = 0; d < dim; ++d)
              value[d] -= this->face_target_velocity[face_index][q][d];

          for (int d = 0; d < dim; ++d)
            {
              // Assemble (v,beta (u-u_target)) for the main penalization
              value_result[d] += penalty_parameter * value[d];

              // Assemble ν(v,∇δu·n)
              for (unsigned int i = 0; i < dim; ++i)
                value_result[d] -=
                  kinematic_viscosity * gradient[d][i] * normal_vector[i];
              // Assemble ν(∇v·n,(u-u_target))
              for (unsigned int i = 0; i < dim; ++i)
                gradient_result[d][i] -=
                  kinematic_viscosity * value[d] * normal_vector[i];
            }
          integrator.submit_value(value_result, q);
          integrator.submit_gradient(gradient_result, q);
        }

      integrator.integrate(EvaluationFlags::EvaluationFlags::values |
                           EvaluationFlags::EvaluationFlags::gradients);
    }

  else if (this->boundary_conditions.type.at(integrator.boundary_id()) ==
           BoundaryConditions::BoundaryType::outlet)
    {
      integrator.evaluate(EvaluationFlags::values);

      for (const auto q : integrator.quadrature_point_indices())
        {
          typename FEFaceIntegrator::value_type value_result = {};

          const auto normal_vector = integrator.normal_vector(q);

          auto value = integrator.get_value(q);

          Tensor<1, dim, VectorizedArray<number>> velocity;
          if constexpr (assemble_residual)
            {
              for (int d = 0; d < dim; ++d)
                velocity[d] = value[d];
            }
          else
            {
              // Only use the previous velocity values for the Jacobian
              for (int d = 0; d < dim; ++d)
                velocity[d] = face_nonlinear_previous_values[face][q][d];
            }

          // (u·n)_=min(0,u·n)
          VectorizedArray<number> normal_outflux = velocity * normal_vector;

          normal_outflux =
            std::min(VectorizedArray<number>(0.0), normal_outflux);

          // -(v,β(u·n)_·u)
          for (int d = 0; d < dim; ++d)
            value_result[d] -= penalty_parameter * normal_outflux * value[d];

          integrator.submit_value(value_result, q);
        }

      integrator.integrate(EvaluationFlags::EvaluationFlags::values);
    }
}


template <int dim, typename number>
IndexSet
NavierStokesOperatorBase<dim, number>::get_refinement_edges(
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

template class NavierStokesOperatorBase<2, double>;
template class NavierStokesOperatorBase<3, double>;
template class NavierStokesOperatorBase<2, float>;
template class NavierStokesOperatorBase<3, float>;

template <int dim, typename number>
NavierStokesStabilizedOperator<dim, number>::NavierStokesStabilizedOperator() =
  default;

/**
 * The expressions calculated in this cell integral are:
 * (q,∇δu) + (v,∂t δu) + (v,(u·∇)δu) + (v,(δu·∇)u) - (∇·v,δp) + ν(∇v,∇δu)
 * (Weak form Jacobian), plus three additional terms in the case of SUPG-PSPG
 * stabilization:
 * \+ (∂t δu +(u·∇)δu + (δu·∇)u + ∇δp - ν∆δu)τ·∇q (PSPG Jacobian)
 * \+ (∂t δu +(u·∇)δu + (δu·∇)u + ∇δp - ν∆δu)τu·∇v (SUPG Jacobian Part 1)
 * \+ (∂t u +(u·∇)u + ∇p - ν∆u - f )τδu·∇v (SUPG Jacobian Part 2),
 * plus two additional terms in the case of full gls stabilization:
 * \+ (∂t δu +(u·∇)δu + (δu·∇)u + ∇δp - ν∆δu)τ(−ν∆v) (GLS Jacobian)
 * \+ (∇·δu)τ'(∇·v) (LSIC Jacobian).
 */
template <int dim, typename number>
void
NavierStokesStabilizedOperator<dim, number>::do_cell_integral_local(
  FECellIntegrator &integrator) const
{
  if (this->enable_hessians_jacobian)
    integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                        EvaluationFlags::hessians);
  else
    integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

  const unsigned int cell = integrator.get_current_cell_index();

  // To identify whether the problem is transient or steady
  bool transient =
    (is_bdf(this->simulation_control->get_assembly_method())) ? true : false;

  // To identify whether mortar feature is enabled (use local variable for
  // efficiency)
  bool enable_mortar = this->enable_mortar;

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  // Vector for BDF coefficients
  const Vector<double> *bdf_coefs;
  if (transient)
    bdf_coefs = &this->simulation_control->get_bdf_coefficients();

  for (const auto q : integrator.quadrature_point_indices())
    {
      Tensor<1, dim, VectorizedArray<number>> source_value;

      // Evaluate source term function if enabled
      if (!this->buoyancy_term.empty())
        source_value = this->buoyancy_term(cell, q);
      else if (this->forcing_function)
        source_value = this->forcing_terms(cell, q);

      // Add to source term the dynamic flow control force (zero if not
      // enabled)
      source_value += this->beta_force;

      // Gather the original value/gradient
      typename FECellIntegrator::value_type    value = integrator.get_value(q);
      typename FECellIntegrator::gradient_type gradient =
        integrator.get_gradient(q);
      typename FECellIntegrator::gradient_type hessian_diagonal;

      if (this->enable_hessians_jacobian)
        hessian_diagonal = integrator.get_hessian_diagonal(q);

      // Result value/gradient we will use
      typename FECellIntegrator::value_type    value_result;
      typename FECellIntegrator::gradient_type gradient_result;
      typename FECellIntegrator::hessian_type  hessian_result;

      // Gather previous values of the velocity and the pressure
      auto previous_values   = this->nonlinear_previous_values(cell, q);
      auto previous_gradient = this->nonlinear_previous_gradient(cell, q);
      auto previous_hessian_diagonal =
        this->nonlinear_previous_hessian_diagonal(cell, q);

      Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
      if (transient)
        previous_time_derivatives =
          this->time_derivatives_previous_solutions(cell, q);

      Tensor<1, dim, VectorizedArray<number>> u_ale;
      if (enable_mortar)
        u_ale = this->velocity_ale[cell][q];

      // Get stabilization parameter
      const auto tau      = this->stabilization_parameter[cell][q];
      const auto tau_lsic = this->stabilization_parameter_lsic[cell][q];

      // Weak form Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          // ν(∇v,∇δu)
          gradient_result[i] = kinematic_viscosity * gradient[i];
          // -(∇·v,δp)
          gradient_result[i][i] += -value[dim];
          // +(q,∇δu)
          value_result[dim] += gradient[i][i];

          for (unsigned int k = 0; k < dim; ++k)
            {
              // +(v,(u·∇)δu + (δu·∇)u)
              value_result[i] += gradient[i][k] * previous_values[k] +
                                 previous_gradient[i][k] * value[k];

              // -(v, (u_ale·∇)δu)
              if (enable_mortar)
                value_result[i] -= gradient[i][k] * u_ale[k];
            }
          // +(v,∂t δu)
          if (transient)
            value_result[i] += (*bdf_coefs)[0] * value[i];
        }

      // PSPG Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              // (-ν∆δu + (u·∇)δu + (δu·∇)u)·τ∇q
              gradient_result[dim][i] +=
                tau * (-kinematic_viscosity * hessian_diagonal[i][k] +
                       gradient[i][k] * previous_values[k] +
                       previous_gradient[i][k] * value[k]);

              // -(u_ale·∇δu)·τ∇q
              if (enable_mortar)
                gradient_result[dim][i] -= tau * (gradient[i][k] * u_ale[k]);
            }
          // +(∂t δu)·τ∇q
          if (transient)
            gradient_result[dim][i] += tau * (*bdf_coefs)[0] * value[i];
        }
      // (∇δp)τ·∇q
      gradient_result[dim] += tau * gradient[dim];

      // SUPG Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              // Part 1
              for (unsigned int l = 0; l < dim; ++l)
                {
                  // +((u·∇)δu + (δu·∇)u - ν∆δu)τ(u·∇)v
                  gradient_result[i][k] +=
                    tau * previous_values[k] *
                    (gradient[i][l] * previous_values[l] +
                     previous_gradient[i][l] * value[l] -
                     kinematic_viscosity * hessian_diagonal[i][l]);

                  if (enable_mortar)
                    {
                      // -((u_ale·∇)δu)τ(u·∇)v
                      gradient_result[i][k] -=
                        tau * previous_values[k] * gradient[i][l] * u_ale[l];
                      // (-(u_ale·∇)δu)τ(-u_ale·∇)v
                      gradient_result[i][k] +=
                        tau * u_ale[k] * gradient[i][l] * u_ale[l];
                      // -((u·∇)δu + (δu·∇)u - ν∆δu)τ(u_ale·∇)v
                      gradient_result[i][k] -=
                        tau * u_ale[k] *
                        (gradient[i][l] * previous_values[l] +
                         previous_gradient[i][l] * value[l] -
                         kinematic_viscosity * hessian_diagonal[i][l]);
                    }
                }
              // +(∇δp)τ(u·∇)v
              gradient_result[i][k] +=
                tau * previous_values[k] * (gradient[dim][i]);

              // +(∂t δu)τ(u·∇)v
              if (transient)
                gradient_result[i][k] +=
                  tau * previous_values[k] * ((*bdf_coefs)[0] * value[i]);

              if (enable_mortar)
                {
                  // -(∇δp)τ(u_ale·∇)v
                  gradient_result[i][k] -= tau * u_ale[k] * (gradient[dim][i]);
                  // -(∂t δu)τ(u_ale·∇)v
                  if (transient)
                    gradient_result[i][k] -=
                      tau * u_ale[k] * ((*bdf_coefs)[0] * value[i]);
                }

              // Part 2
              for (unsigned int l = 0; l < dim; ++l)
                {
                  // +((u·∇)u - ν∆u)τ(δu·∇)v
                  gradient_result[i][k] +=
                    tau * value[k] *
                    (previous_gradient[i][l] * previous_values[l] -
                     kinematic_viscosity * previous_hessian_diagonal[i][l]);

                  // -((u_ale·∇)u)τ(δu·∇)v
                  if (enable_mortar)
                    gradient_result[i][k] -=
                      tau * value[k] * previous_gradient[i][l] * u_ale[l];
                }
              // +(∇p - f)τ(δu·∇)v
              gradient_result[i][k] +=
                tau * value[k] * (previous_gradient[dim][i] - source_value[i]);

              // +(∂t u)τ(δu·∇)v
              if (transient)
                gradient_result[i][k] += tau * value[k] *
                                         ((*bdf_coefs)[0] * previous_values[i] +
                                          previous_time_derivatives[i]);
            }
        }

      if (this->stabilization ==
          Parameters::Stabilization::NavierStokesStabilization::gls)
        {
          // GLS Jacobian
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  if (this->enable_hessians_jacobian)
                    {
                      for (unsigned int l = 0; l < dim; ++l)
                        {
                          // +((u·∇)δu + (δu·∇)u - ν∆δu)τ(−ν∆v)
                          hessian_result[i][k][k] +=
                            tau * -kinematic_viscosity *
                            (gradient[i][l] * previous_values[l] +
                             previous_gradient[i][l] * value[l] -
                             kinematic_viscosity * hessian_diagonal[i][l]);
                        }

                      // +(∇δp)τ(−ν∆v)
                      hessian_result[i][k][k] +=
                        tau * -kinematic_viscosity * (gradient[dim][i]);

                      // +(∂t δu)τ(−ν∆v)
                      if (transient)
                        hessian_result[i][k][k] += tau * -kinematic_viscosity *
                                                   ((*bdf_coefs)[0] * value[i]);
                    }

                  // LSIC term
                  // (∇·δu)τ'(∇·v)
                  gradient_result[i][i] += tau_lsic * gradient[k][k];
                }
            }
        }

      integrator.submit_gradient(gradient_result, q);
      integrator.submit_value(value_result, q);

      if (this->enable_hessians_jacobian)
        integrator.submit_hessian(hessian_result, q);
    }

  if (this->enable_hessians_jacobian)
    integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients |
                         EvaluationFlags::hessians);
  else
    integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

/**
 * The expressions calculated in this cell integral are:
 * (q, ∇·u) + (v,∂t u) + (v,(u·∇)u) - (∇·v,p) + ν(∇v,∇u) - (v,f) (Weak form),
 * plus two additional terms in the case of SUPG-PSPG stabilization:
 * \+ (∂t u +(u·∇)u + ∇p - ν∆u - f)τ∇·q (PSPG term)
 * \+ (∂t u +(u·∇)u + ∇p - ν∆u - f)τu·∇v (SUPG term),
 * plus two additional terms in the case of full gls stabilization:
 * \+ (∂t u +(u·∇)u + ∇p - ν∆u - f)τ(−ν∆v) (GLS term)
 * \+ (∇·u)τ'(∇·v) (LSIC term).
 */
template <int dim, typename number>
void
NavierStokesStabilizedOperator<dim, number>::local_evaluate_residual(
  const MatrixFree<dim, number>               &matrix_free,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::pair<unsigned int, unsigned int> &range) const
{
  FECellIntegrator integrator(matrix_free);

  const double kinematic_viscosity =
    this->properties_manager->get_rheology()->get_kinematic_viscosity();

  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      integrator.reinit(cell);
      integrator.read_dof_values_plain(src);

      if (this->enable_hessians_residual)
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients |
                            EvaluationFlags::hessians);
      else
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients);

      // To identify whether the problem is transient or steady
      bool transient =
        (is_bdf(this->simulation_control->get_assembly_method())) ? true :
                                                                    false;

      // To identify whether mortar feature is enabled (use local variable for
      // efficiency)
      bool enable_mortar = this->enable_mortar;

      // Vector for BDF coefficients
      const Vector<double> *bdf_coefs;
      if (transient)
        bdf_coefs = &this->simulation_control->get_bdf_coefficients();

      for (const auto q : integrator.quadrature_point_indices())
        {
          Tensor<1, dim, VectorizedArray<number>> source_value;

          // Evaluate source term function if enabled
          if (!this->buoyancy_term.empty())
            source_value = this->buoyancy_term(cell, q);
          else if (this->forcing_function)
            source_value = this->forcing_terms(cell, q);

          // Add to source term the dynamic flow control force (zero if not
          // enabled)
          source_value += this->beta_force;

          // Gather the original value/gradient
          typename FECellIntegrator::value_type value = integrator.get_value(q);
          typename FECellIntegrator::gradient_type gradient =
            integrator.get_gradient(q);
          typename FECellIntegrator::gradient_type hessian_diagonal;

          if (this->enable_hessians_residual)
            hessian_diagonal = integrator.get_hessian_diagonal(q);

          // Time derivatives of previous solutions
          Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
          if (transient)
            previous_time_derivatives =
              this->time_derivatives_previous_solutions(cell, q);

          Tensor<1, dim, VectorizedArray<number>> u_ale;
          if (enable_mortar)
            u_ale = this->velocity_ale[cell][q];

          // Get stabilization parameter
          const auto tau      = this->stabilization_parameter[cell][q];
          const auto tau_lsic = this->stabilization_parameter_lsic[cell][q];

          // Result value/gradient we will use
          typename FECellIntegrator::value_type    value_result;
          typename FECellIntegrator::gradient_type gradient_result;
          typename FECellIntegrator::hessian_type  hessian_result;

          // Weak form
          for (unsigned int i = 0; i < dim; ++i)
            {
              // ν(∇v,∇u)
              gradient_result[i] = kinematic_viscosity * gradient[i];
              // -(∇·v,p)
              gradient_result[i][i] += -value[dim];
              // +(v,-f)
              value_result[i] = -source_value[i];

              // +(v,∂t u)
              if (transient)
                value_result[i] +=
                  (*bdf_coefs)[0] * value[i] + previous_time_derivatives[i];


              // +(q,∇·u)
              value_result[dim] += gradient[i][i];

              for (unsigned int k = 0; k < dim; ++k)
                {
                  // +(v,(u·∇)u)
                  value_result[i] += gradient[i][k] * value[k];

                  // -(v, (u_ale·∇)u)
                  if (enable_mortar)
                    value_result[i] -= gradient[i][k] * u_ale[k];
                }
            }

          // PSPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  // (-ν∆u + (u·∇)u)·τ∇q
                  gradient_result[dim][i] +=
                    tau * (-kinematic_viscosity * hessian_diagonal[i][k] +
                           gradient[i][k] * value[k]);

                  // -((u_ale·∇)u)·τ∇q
                  if (enable_mortar)
                    gradient_result[dim][i] -= tau * gradient[i][k] * u_ale[k];
                }
              // +(-f)·τ∇q
              gradient_result[dim][i] += tau * (-source_value[i]);

              // +(∂t u)·τ∇q
              if (transient)
                gradient_result[dim][i] += tau * ((*bdf_coefs)[0] * value[i] +
                                                  previous_time_derivatives[i]);
            }
          // +(∇p)τ∇·q
          gradient_result[dim] += tau * gradient[dim];

          // SUPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  for (unsigned int l = 0; l < dim; ++l)
                    {
                      // (-ν∆u)τ(u·∇)v
                      gradient_result[i][k] += -tau * kinematic_viscosity *
                                               value[k] *
                                               hessian_diagonal[i][l];

                      // + ((u·∇)u)τ(u·∇)v
                      gradient_result[i][k] +=
                        tau * value[k] * gradient[i][l] * value[l];

                      if (enable_mortar)
                        {
                          // -((u_ale·∇)u)τ(u·∇)v
                          gradient_result[i][k] -=
                            tau * value[k] * gradient[i][l] * u_ale[l];

                          // -(-ν∆u)τ(u_ale·∇)v
                          gradient_result[i][k] += tau * kinematic_viscosity *
                                                   u_ale[k] *
                                                   hessian_diagonal[i][l];

                          // -((u·∇)u)τ(u_ale·∇)v
                          gradient_result[i][k] -=
                            tau * u_ale[k] * gradient[i][l] * value[l];

                          // -(-(u_ale·∇)u)τ(u_ale·∇)v
                          gradient_result[i][k] +=
                            tau * u_ale[k] * gradient[i][l] * u_ale[l];
                        }
                    }
                  // + (∇p - f)τ(u·∇)v
                  gradient_result[i][k] +=
                    tau * value[k] * (gradient[dim][i] - source_value[i]);

                  // + (∂t u)τ(u·∇)v
                  if (transient)
                    gradient_result[i][k] += tau * value[k] *
                                             ((*bdf_coefs)[0] * value[i] +
                                              previous_time_derivatives[i]);

                  if (enable_mortar)
                    {
                      // - (∇p - f)τ(u_ale·∇)v
                      gradient_result[i][k] -=
                        tau * u_ale[k] * (gradient[dim][i] - source_value[i]);

                      // - (∂t u)τ(u_ale·∇)v
                      if (transient)
                        gradient_result[i][k] -= tau * u_ale[k] *
                                                 ((*bdf_coefs)[0] * value[i] +
                                                  previous_time_derivatives[i]);
                    }
                }
            }

          if (this->stabilization ==
              Parameters::Stabilization::NavierStokesStabilization::gls)
            {
              // GLS term
              for (unsigned int i = 0; i < dim; ++i)
                {
                  for (unsigned int k = 0; k < dim; ++k)
                    {
                      if (this->enable_hessians_residual)
                        {
                          for (unsigned int l = 0; l < dim; ++l)
                            {
                              // (-ν∆u + (u·∇)u)τ(−ν∆v)
                              hessian_result[i][k][k] +=
                                tau * -kinematic_viscosity *
                                (-kinematic_viscosity * hessian_diagonal[i][l] +
                                 gradient[i][l] * value[l]);
                            }
                          // + (∇p - f)τ(−ν∆v)
                          hessian_result[i][k][k] +=
                            tau * -kinematic_viscosity *
                            (gradient[dim][i] - source_value[i]);

                          // + (∂t u)τ(−ν∆v)
                          if (transient)
                            hessian_result[i][k][k] +=
                              tau * -kinematic_viscosity *
                              ((*bdf_coefs)[0] * value[i] +
                               previous_time_derivatives[i]);
                        }

                      // LSIC term
                      // (∇·u)τ'(∇·v)
                      gradient_result[i][i] += tau_lsic * gradient[k][k];
                    }
                }
            }

          integrator.submit_gradient(gradient_result, q);
          integrator.submit_value(value_result, q);
          if (this->enable_hessians_residual)
            integrator.submit_hessian(hessian_result, q);
        }

      if (this->enable_hessians_residual)
        integrator.integrate_scatter(EvaluationFlags::values |
                                       EvaluationFlags::gradients |
                                       EvaluationFlags::hessians,
                                     dst);
      else
        integrator.integrate_scatter(EvaluationFlags::values |
                                       EvaluationFlags::gradients,
                                     dst);
    }
}

template class NavierStokesStabilizedOperator<2, double>;
template class NavierStokesStabilizedOperator<3, double>;
template class NavierStokesStabilizedOperator<2, float>;
template class NavierStokesStabilizedOperator<3, float>;


template <int dim, typename number>
NavierStokesNonNewtonianStabilizedOperator<dim, number>::
  NavierStokesNonNewtonianStabilizedOperator() = default;

/**
 * The expressions calculated in this cell integral are:
 * (q,∇δu) + (v,∂t δu) + (v,(u·∇)δu) + (v,(δu·∇)u) - (∇·v,δp) + ν(∇v,(∇δu +
 * ∇δuT)) + (∇v, 0.5/γ_dot (∂ν/∂γ_dot)(∇u + ∇uT)(∇δu + ∇δuT)(∇u + ∇uT)) (Weak
 * form Jacobian), plus three additional terms in the case of SUPG-PSPG
 * stabilization:
 * \+ (∂t δu +(u·∇)δu + (δu·∇)u + ∇δp - ν∆δu - (∇ν)(∇δu + ∇δuT))τ·∇q (PSPG
 * Jacobian)
 * \+ (∂t δu +(u·∇)δu + (δu·∇)u + ∇δp - ν∆δu - (∇ν)(∇δu + ∇δuT))τu·∇v (SUPG
 * Jacobian Part 1)
 * \+ (∂t u +(u·∇)u + ∇p - ν∆u - (∇ν)((∇u + ∇uT) - f )τδu·∇v (SUPG Jacobian
 * Part 2),
 */
template <int dim, typename number>
void
NavierStokesNonNewtonianStabilizedOperator<dim, number>::do_cell_integral_local(
  FECellIntegrator &integrator) const
{
  if (this->enable_hessians_jacobian)
    integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients |
                        EvaluationFlags::hessians);
  else
    integrator.evaluate(EvaluationFlags::values | EvaluationFlags::gradients);

  const unsigned int cell = integrator.get_current_cell_index();

  // To identify whether the problem is transient or steady
  bool transient =
    (is_bdf(this->simulation_control->get_assembly_method())) ? true : false;

  // Vector for BDF coefficients
  const Vector<double> *bdf_coefs;
  if (transient)
    bdf_coefs = &this->simulation_control->get_bdf_coefficients();

  for (const auto q : integrator.quadrature_point_indices())
    {
      Tensor<1, dim, VectorizedArray<number>> source_value;

      // Evaluate source term function if enabled
      if (this->forcing_function)
        source_value = this->forcing_terms(cell, q);

      // Add to source term the dynamic flow control force (zero if not
      // enabled)
      source_value += this->beta_force;

      // Gather the original value/gradient
      typename FECellIntegrator::value_type    value = integrator.get_value(q);
      typename FECellIntegrator::gradient_type gradient =
        integrator.get_gradient(q);
      typename FECellIntegrator::gradient_type hessian_diagonal;

      if (this->enable_hessians_jacobian)
        hessian_diagonal = integrator.get_hessian_diagonal(q);

      // Result value/gradient we will use
      typename FECellIntegrator::value_type    value_result;
      typename FECellIntegrator::gradient_type gradient_result;
      typename FECellIntegrator::hessian_type  hessian_result;

      // Gather previous values of the velocity and the pressure
      auto previous_values   = this->nonlinear_previous_values(cell, q);
      auto previous_gradient = this->nonlinear_previous_gradient(cell, q);
      auto previous_hessian_diagonal =
        this->nonlinear_previous_hessian_diagonal(cell, q);

      Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
      if (transient)
        previous_time_derivatives =
          this->time_derivatives_previous_solutions(cell, q);

      // Gather previous values of the shear_rate and shear_rate_magnitude
      auto previous_shear_rate = this->previous_shear_rate(cell, q);
      auto previous_shear_rate_magnitude =
        this->previous_shear_rate_magnitude(cell, q);

      previous_shear_rate_magnitude =
        std::max(previous_shear_rate_magnitude, VectorizedArray<number>(1e-3));

      // Get kinematic viscosity from rheology model
      const auto kinematic_viscosity =
        this->kinematic_viscosity_vector[cell][q];
      const auto kinematic_viscosity_gradient =
        this->kinematic_viscosity_gradient[cell][q];
      const auto grad_kinematic_viscosity_shear_rate =
        this->grad_kinematic_viscosity_shear_rate[cell][q];

      // Get stabilization parameter
      const auto tau = this->stabilization_parameter[cell][q];

      // Calculate shear rate per component
      typename FECellIntegrator::gradient_type shear_rate;
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int j = 0; j < dim; ++j)
            {
              shear_rate[i][j] = gradient[i][j] + gradient[j][i];
            }
        }

      // Precompute double contraction between shear rate and previous shear
      // rate
      VectorizedArray<number> shear_rate_previous_product =
        VectorizedArray<number>(0.);

      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              shear_rate_previous_product +=
                shear_rate[i][k] * previous_shear_rate[k][i];
            }
        }

      // Weak form Jacobian
      value_result[dim] = VectorizedArray<number>(0.);
      for (unsigned int i = 0; i < dim; ++i)
        {
          // -(∇·v,δp)
          gradient_result[i][i] = -value[dim];
          // +(q,∇δu)
          value_result[dim] += gradient[i][i];
          // ν(∇v,(∇δu + ∇δuT))
          gradient_result[i] += kinematic_viscosity * shear_rate[i];
          for (unsigned int k = 0; k < dim; ++k)
            {
              // (∇v, 0.5/γ_dot (∂ν/∂γ_dot)(∇u + ∇uT)(∇δu + ∇δuT)(∇u + ∇uT))
              gradient_result[i][k] += 0.5 / previous_shear_rate_magnitude *
                                       grad_kinematic_viscosity_shear_rate *
                                       shear_rate_previous_product *
                                       previous_shear_rate[i][k];
              // +(v,(u·∇)δu + (δu·∇)u)
              value_result[i] += gradient[i][k] * previous_values[k] +
                                 previous_gradient[i][k] * value[k];
            }
          // +(v,∂t δu)
          if (transient)
            value_result[i] += (*bdf_coefs)[0] * value[i];
        }

      // PSPG Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              // (-ν∆δu - (∇ν)(∇δu + ∇δuT) + (u·∇)δu + (δu·∇)u)·τ∇q
              gradient_result[dim][i] +=
                tau * (-kinematic_viscosity * hessian_diagonal[i][k] -
                       kinematic_viscosity_gradient[k] * shear_rate[k][i] +
                       gradient[i][k] * previous_values[k] +
                       previous_gradient[i][k] * value[k]);
            }
          // +(∂t δu)·τ∇q
          if (transient)
            gradient_result[dim][i] += tau * (*bdf_coefs)[0] * value[i];
        }
      // (∇δp)τ·∇q
      gradient_result[dim] += tau * gradient[dim];

      // SUPG Jacobian
      for (unsigned int i = 0; i < dim; ++i)
        {
          for (unsigned int k = 0; k < dim; ++k)
            {
              // Part 1
              for (unsigned int l = 0; l < dim; ++l)
                {
                  // +((u·∇)δu + (δu·∇)u - ν∆δu  - (∇ν)(∇δu + ∇δuT))τ(u·∇)v
                  gradient_result[i][k] +=
                    tau * previous_values[k] *
                    (gradient[i][l] * previous_values[l] +
                     previous_gradient[i][l] * value[l] -
                     kinematic_viscosity * hessian_diagonal[i][l] -
                     kinematic_viscosity_gradient[i] * shear_rate[i][l]);
                }
              // +(∇δp)τ(u·∇)v
              gradient_result[i][k] +=
                tau * previous_values[k] * (gradient[dim][i]);

              // +(∂t δu)τ(u·∇)v
              if (transient)
                gradient_result[i][k] +=
                  tau * previous_values[k] * ((*bdf_coefs)[0] * value[i]);


              // Part 2
              for (unsigned int l = 0; l < dim; ++l)
                {
                  // +((u·∇)u - ν∆u - (∇ν)((∇u + ∇uT))τ(δu·∇)v
                  gradient_result[i][k] +=
                    tau * value[k] *
                    (previous_gradient[i][l] * previous_values[l] -
                     kinematic_viscosity * previous_hessian_diagonal[i][l] -
                     kinematic_viscosity_gradient[i] *
                       previous_shear_rate[i][l]);
                }
              // +(∇p - f)τ(δu·∇)v
              gradient_result[i][k] +=
                tau * value[k] * (previous_gradient[dim][i] - source_value[i]);

              // +(∂t u)τ(δu·∇)v
              if (transient)
                gradient_result[i][k] += tau * value[k] *
                                         ((*bdf_coefs)[0] * previous_values[i] +
                                          previous_time_derivatives[i]);
            }
        }

      integrator.submit_gradient(gradient_result, q);
      integrator.submit_value(value_result, q);

      if (this->enable_hessians_jacobian)
        integrator.submit_hessian(hessian_result, q);
    }

  if (this->enable_hessians_jacobian)
    integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients |
                         EvaluationFlags::hessians);
  else
    integrator.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

/**
 * The expressions calculated in this cell integral are:
 * (q, ∇·u) + (v,∂t u) + (v,(u·∇)u) - (∇·v,p) + ν(∇v,(∇u + ∇uT)) - (v,f) (Weak
 * form), plus two additional terms in the case of SUPG-PSPG stabilization:
 * \+ (∂t u +(u·∇)u + ∇p -ν∆u - (∇ν)((∇u + ∇uT) - f)τ∇·q (PSPG term)
 * \+ (∂t u +(u·∇)u + ∇p -ν∆u - (∇ν)((∇u + ∇uT) - f)τu·∇v (SUPG term),
 */
template <int dim, typename number>
void
NavierStokesNonNewtonianStabilizedOperator<dim, number>::
  local_evaluate_residual(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const
{
  FECellIntegrator integrator(matrix_free);

  for (unsigned int cell = range.first; cell < range.second; ++cell)
    {
      integrator.reinit(cell);
      integrator.read_dof_values_plain(src);

      if (this->enable_hessians_residual)
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients |
                            EvaluationFlags::hessians);
      else
        integrator.evaluate(EvaluationFlags::values |
                            EvaluationFlags::gradients);

      // To identify whether the problem is transient or steady
      bool transient =
        (is_bdf(this->simulation_control->get_assembly_method())) ? true :
                                                                    false;

      // Vector for BDF coefficients
      const Vector<double> *bdf_coefs;
      if (transient)
        bdf_coefs = &this->simulation_control->get_bdf_coefficients();

      for (const auto q : integrator.quadrature_point_indices())
        {
          Tensor<1, dim, VectorizedArray<number>> source_value;

          // Evaluate source term function if enabled
          if (this->forcing_function)
            source_value = this->forcing_terms(cell, q);

          // Add to source term the dynamic flow control force (zero if not
          // enabled)
          source_value += this->beta_force;

          // Gather the original value/gradient
          typename FECellIntegrator::value_type value = integrator.get_value(q);
          typename FECellIntegrator::gradient_type gradient =
            integrator.get_gradient(q);
          typename FECellIntegrator::gradient_type hessian_diagonal;

          if (this->enable_hessians_residual)
            hessian_diagonal = integrator.get_hessian_diagonal(q);

          // Time derivatives of previous solutions
          Tensor<1, dim + 1, VectorizedArray<number>> previous_time_derivatives;
          if (transient)
            previous_time_derivatives =
              this->time_derivatives_previous_solutions(cell, q);

          // Get stabilization parameter
          const auto tau = this->stabilization_parameter[cell][q];

          // Get kinematic viscosity from rheology model
          const auto kinematic_viscosity =
            this->kinematic_viscosity_vector[cell][q];

          // Get kinematic viscosity gradient
          const auto kinematic_viscosity_gradient =
            this->kinematic_viscosity_gradient[cell][q];

          // Result value/gradient we will use
          typename FECellIntegrator::value_type    value_result;
          typename FECellIntegrator::gradient_type gradient_result;
          typename FECellIntegrator::hessian_type  hessian_result;

          // Additional values/gradients needed
          typename FECellIntegrator::gradient_type shear_rate;

          // Calculate shear rate per component
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int j = 0; j < dim; ++j)
                {
                  shear_rate[i][j] = gradient[i][j] + gradient[j][i];
                }
            }

          // Weak form
          for (unsigned int i = 0; i < dim; ++i)
            {
              // ν(∇v,(∇u + ∇uT))
              gradient_result[i] = kinematic_viscosity * shear_rate[i];
              // -(∇·v,p)
              gradient_result[i][i] += -value[dim];
              // +(v,-f)
              value_result[i] = -source_value[i];

              // +(v,∂t u)
              if (transient)
                value_result[i] +=
                  (*bdf_coefs)[0] * value[i] + previous_time_derivatives[i];


              // +(q,∇·u)
              value_result[dim] += gradient[i][i];

              for (unsigned int k = 0; k < dim; ++k)
                {
                  // +(v,(u·∇)u)
                  value_result[i] += gradient[i][k] * value[k];
                }
            }

          // PSPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  // (-ν∆u - (∇ν)((∇u + ∇uT)) + (u·∇)u)·τ∇q
                  gradient_result[dim][i] +=
                    tau * (-kinematic_viscosity * hessian_diagonal[i][k] -
                           kinematic_viscosity_gradient[k] * shear_rate[k][i] +
                           gradient[i][k] * value[k]);
                }
              // +(-f)·τ∇q
              gradient_result[dim][i] += tau * (-source_value[i]);

              // +(∂t u)·τ∇q
              if (transient)
                gradient_result[dim][i] += tau * ((*bdf_coefs)[0] * value[i] +
                                                  previous_time_derivatives[i]);
            }
          // +(∇p)τ∇·q
          gradient_result[dim] += tau * gradient[dim];

          // SUPG term
          for (unsigned int i = 0; i < dim; ++i)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  for (unsigned int l = 0; l < dim; ++l)
                    {
                      // (-ν∆u - (∇ν)((∇u + ∇uT)))τ(u·∇)v
                      gradient_result[i][k] +=
                        tau * value[k] *
                        (-kinematic_viscosity * hessian_diagonal[i][l] -
                         kinematic_viscosity_gradient[l] * shear_rate[l][i]);

                      // + ((u·∇)u)τ(u·∇)v
                      gradient_result[i][k] +=
                        tau * value[k] * gradient[i][l] * value[l];
                    }
                  // + (∇p - f)τ(u·∇)v
                  gradient_result[i][k] +=
                    tau * value[k] * (gradient[dim][i] - source_value[i]);

                  // + (∂t u)τ(u·∇)v
                  if (transient)
                    gradient_result[i][k] += tau * value[k] *
                                             ((*bdf_coefs)[0] * value[i] +
                                              previous_time_derivatives[i]);
                }
            }

          integrator.submit_gradient(gradient_result, q);
          integrator.submit_value(value_result, q);
          if (this->enable_hessians_residual)
            integrator.submit_hessian(hessian_result, q);
        }

      if (this->enable_hessians_residual)
        integrator.integrate_scatter(EvaluationFlags::values |
                                       EvaluationFlags::gradients |
                                       EvaluationFlags::hessians,
                                     dst);
      else
        integrator.integrate_scatter(EvaluationFlags::values |
                                       EvaluationFlags::gradients,
                                     dst);
    }
}

template class NavierStokesNonNewtonianStabilizedOperator<2, double>;
template class NavierStokesNonNewtonianStabilizedOperator<3, double>;
template class NavierStokesNonNewtonianStabilizedOperator<2, float>;
template class NavierStokesNonNewtonianStabilizedOperator<3, float>;
