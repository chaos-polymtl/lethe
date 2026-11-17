// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/ib_scalar_constraints.h>
#include <core/ib_stencil.h>
#include <core/lethe_grid_tools.h>
#include <core/tracer_reaction_model.h>

#include <deal.II/base/mpi.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>

namespace
{
  /**
   * @brief Return the derivatives of the shape functions of a cell at a point,
   * along a direction.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] mapping Mapping of the cell.
   *
   * @param[in] cell Cell whose shape functions are differentiated.
   *
   * @param[in] point Point at which the derivatives are evaluated. It can lie
   * outside of the cell, in which case the polynomials of the cell are
   * extrapolated.
   *
   * @param[in] direction Unit vector along which the shape functions are
   * differentiated.
   *
   * @return Derivatives of the shape functions, in the order of the degrees of
   * freedom of the cell.
   */
  template <int dim>
  std::vector<double>
  compute_shape_function_directional_derivatives(
    const Mapping<dim>                                   &mapping,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    const Point<dim>                                     &point,
    const Tensor<1, dim>                                 &direction)
  {
    const FiniteElement<dim> &fe = cell->get_fe();
    FEValues<dim>             fe_values(mapping,
                            fe,
                            Quadrature<dim>(
                              mapping.transform_real_to_unit_cell(cell, point)),
                            update_gradients);
    fe_values.reinit(cell);

    std::vector<double> derivatives(fe.n_dofs_per_cell());
    for (unsigned int m = 0; m < fe.n_dofs_per_cell(); ++m)
      derivatives[m] = fe_values.shape_grad(m, 0) * direction;
    return derivatives;
  }

  /**
   * @brief Return the rate of a surface reaction of order n,
   * \f$ \alpha C_\text{eff}^{n-1} c \f$ with
   * \f$ C_\text{eff} = \sqrt{c^2 + \epsilon^2} \f$, and its derivative with
   * respect to the concentration \f$ c \f$,
   * \f$ \alpha C_\text{eff}^{n-1} \left(1 + (n-1) (c / C_\text{eff})^2
   * \right) \f$. The effective concentration is the one of the volumetric
   * reaction (see calculate_effective_concentration()).
   *
   * @param[in] concentration Concentration \f$ c \f$ on the surface.
   *
   * @param[in] reaction_constant Reaction constant \f$ \alpha \f$.
   *
   * @param[in] reaction_order Reaction order \f$ n \f$. If it is lower than 1,
   * the reaction threshold must be positive.
   *
   * @param[in] reaction_threshold Reaction threshold \f$ \epsilon \f$.
   *
   * @return Pair containing the reaction rate and its derivative.
   */
  std::pair<double, double>
  compute_surface_reaction_rate(const double concentration,
                                const double reaction_constant,
                                const double reaction_order,
                                const double reaction_threshold)
  {
    const double effective_concentration =
      calculate_effective_concentration(concentration, reaction_threshold);
    const double prefactor =
      reaction_constant *
      std::pow(effective_concentration, reaction_order - 1.);

    // Derivative of the effective concentration, c / C_eff. Its magnitude does
    // not exceed 1, and it is zero if the effective concentration is zero.
    const double effective_concentration_derivative =
      calculate_d_effective_concentration_dC(concentration, reaction_threshold);

    return {prefactor * concentration,
            prefactor *
              (1. + (reaction_order - 1.) * effective_concentration_derivative *
                      effective_concentration_derivative)};
  }
} // namespace

template <int dim>
IBScalarConstraints<dim>::IBScalarConstraints(
  const unsigned int                                p_stencil_degree,
  const double                                      p_length_ratio,
  const bool                                        p_enable_extrapolation,
  const Parameters::IBTracerBoundaryCondition<dim> &p_boundary_condition)
  : stencil_degree(p_stencil_degree)
  , length_ratio(p_length_ratio)
  , enable_extrapolation(p_enable_extrapolation)
  , boundary_condition(p_boundary_condition)
  , imposes_flux(
      p_boundary_condition.type == Parameters::IBTracerBoundaryType::neumann ||
      p_boundary_condition.type == Parameters::IBTracerBoundaryType::robin)
  , diffusivity(0.)
  , inverse_time_step(0.)
  , n_interface_dofs(0)
  , n_dofs_without_stencil(0)
{
  AssertThrow(stencil_degree > 0,
              ExcMessage("The degree of the immersed boundary stencil is " +
                         std::to_string(stencil_degree) +
                         ". It must be at least 1. Increase the parameter "
                         "<stencil degree> of the subsection <extrapolation "
                         "function> of the subsection <particles>."));
  AssertThrow(length_ratio > 0.,
              ExcMessage("The length ratio of the immersed boundary stencil "
                         "is " +
                         std::to_string(length_ratio) +
                         ". It must be positive. Change the parameter <length "
                         "ratio> of the subsection <extrapolation function> of "
                         "the subsection <particles>."));
  AssertThrow(boundary_condition.type != Parameters::IBTracerBoundaryType::none,
              ExcMessage("The immersed boundary constraints of a scalar field "
                         "require a boundary condition whose type is not "
                         "none."));
}

template <int dim>
void
IBScalarConstraints<dim>::update(const DoFHandler<dim> &dof_handler,
                                 const DoFHandler<dim> &geometry_dof_handler,
                                 const Mapping<dim>    &mapping,
                                 const std::vector<IBParticle<dim>> &particles,
                                 const AffineConstraints<double> &constraints,
                                 const double                     p_diffusivity,
                                 const double p_inverse_time_step)
{
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  AssertThrow(fe.n_components() == 1 && fe.has_support_points() &&
                fe.conforms(FiniteElementData<dim>::H1),
              ExcMessage("The boundary conditions on the particles of the "
                         "sharp-interface immersed boundary method require a "
                         "scalar and continuous finite element with support "
                         "points (FE_Q)."));
  AssertThrow(p_diffusivity >= 0. && p_inverse_time_step >= 0.,
              ExcMessage("The diffusivity and the inverse of the time step "
                         "used to scale the immersed boundary equations must "
                         "not be negative."));
  AssertThrow(!imposes_flux || p_diffusivity > 0.,
              ExcMessage("A flux is imposed on the surface of the particles "
                         "(neumann or robin type), but the diffusivity is " +
                         std::to_string(p_diffusivity) +
                         ". The diffusivity must be positive to impose a "
                         "flux."));

  diffusivity       = p_diffusivity;
  inverse_time_step = p_inverse_time_step;

  const MPI_Comm  mpi_communicator   = dof_handler.get_mpi_communicator();
  const IndexSet &locally_owned_dofs = dof_handler.locally_owned_dofs();
  locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

  const unsigned int            n_dofs_per_cell = fe.n_dofs_per_cell();
  const types::global_dof_index n_relevant_dofs =
    locally_relevant_dofs.n_elements();
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // The real support points of the degrees of freedom of a cell are the
  // quadrature points of a quadrature built on the unit support points.
  FEValues<dim> fe_values_support_points(mapping,
                                         fe,
                                         Quadrature<dim>(
                                           fe.get_unit_support_points()),
                                         update_quadrature_points);

  // The shapes of the particles use the cells of the geometry DoFHandler as
  // cell guesses, since their precalculations are done on this DoFHandler.
  const auto get_geometry_cell =
    [&geometry_dof_handler](
      const typename DoFHandler<dim>::active_cell_iterator &cell) {
      return typename DoFHandler<dim>::active_cell_iterator(
        &cell->get_triangulation(),
        cell->level(),
        cell->index(),
        &geometry_dof_handler);
    };

  // Signed distance to the closest particle, and index of this particle, at
  // the support point of every locally relevant degree of freedom. A negative
  // distance means that the support point is inside a particle.
  std::vector<double>       dof_signed_distance(n_relevant_dofs,
                                          std::numeric_limits<double>::max());
  std::vector<unsigned int> dof_closest_particle(n_relevant_dofs,
                                                 numbers::invalid_unsigned_int);
  std::vector<bool> dof_signed_distance_is_computed(n_relevant_dofs, false);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() && !cell->is_ghost())
        continue;

      fe_values_support_points.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      const auto geometry_cell = get_geometry_cell(cell);

      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          const types::global_dof_index relevant_index =
            locally_relevant_dofs.index_within_set(local_dof_indices[i]);
          if (dof_signed_distance_is_computed[relevant_index])
            continue;
          dof_signed_distance_is_computed[relevant_index] = true;

          const Point<dim> &support_point =
            fe_values_support_points.quadrature_point(i);
          for (unsigned int p = 0; p < particles.size(); ++p)
            {
              const double signed_distance =
                particles[p].get_levelset(support_point, geometry_cell);
              if (signed_distance < dof_signed_distance[relevant_index])
                {
                  dof_signed_distance[relevant_index]  = signed_distance;
                  dof_closest_particle[relevant_index] = p;
                }
            }
        }
    }

  // Classify the locally owned and ghost cells according to the number of
  // their support points that are inside the particles.
  cell_categories.assign(dof_handler.get_triangulation().n_active_cells(),
                         CellCategory::fluid);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() && !cell->is_ghost())
        continue;

      cell->get_dof_indices(local_dof_indices);
      const unsigned int n_dofs_inside =
        static_cast<unsigned int>(std::ranges::count_if(
          local_dof_indices, [&](const types::global_dof_index dof_index) {
            return dof_signed_distance[locally_relevant_dofs.index_within_set(
                     dof_index)] <= 0.;
          }));

      CellCategory category = CellCategory::cut;
      if (n_dofs_inside == 0)
        category = CellCategory::fluid;
      else if (n_dofs_inside == n_dofs_per_cell)
        category = CellCategory::inside;
      cell_categories[cell->active_cell_index()] = category;
    }

  // With elements of degree 1, a fluid cell whose degrees of freedom all
  // belong to cut cells or are constrained has an equation that only involves
  // constrained degrees of freedom. As in the fluid dynamics solver, such a
  // cell is treated as a cut cell. A degree of freedom can belong to a cut cell
  // that is not visible on this process, hence the communication.
  if (fe.degree == 1)
    {
      LinearAlgebra::distributed::Vector<double> dof_is_in_cut_cell(
        locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (!cell->is_locally_owned() ||
              cell_categories[cell->active_cell_index()] != CellCategory::cut)
            continue;
          cell->get_dof_indices(local_dof_indices);
          for (const auto dof_index : local_dof_indices)
            dof_is_in_cut_cell(dof_index) = 1.;
        }
      dof_is_in_cut_cell.compress(VectorOperation::add);
      dof_is_in_cut_cell.update_ghost_values();

      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if ((!cell->is_locally_owned() && !cell->is_ghost()) ||
              cell_categories[cell->active_cell_index()] != CellCategory::fluid)
            continue;
          cell->get_dof_indices(local_dof_indices);
          const bool cell_is_overconstrained = std::ranges::all_of(
            local_dof_indices, [&](const types::global_dof_index dof_index) {
              return std::as_const(dof_is_in_cut_cell)(dof_index) > 0. ||
                     constraints.is_constrained(dof_index);
            });
          if (cell_is_overconstrained)
            cell_categories[cell->active_cell_index()] = CellCategory::cut;
        }
    }

  // Classify the locally relevant degrees of freedom according to the
  // categories of the cells they belong to. These cells can be owned by other
  // processes, hence the communication.
  LinearAlgebra::distributed::Vector<double> dof_is_in_fluid_cell(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  LinearAlgebra::distributed::Vector<double> dof_is_in_non_fluid_cell(
    locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;
      cell->get_dof_indices(local_dof_indices);
      auto &dof_flags =
        (cell_categories[cell->active_cell_index()] == CellCategory::fluid) ?
          dof_is_in_fluid_cell :
          dof_is_in_non_fluid_cell;
      for (const auto dof_index : local_dof_indices)
        dof_flags(dof_index) = 1.;
    }
  dof_is_in_fluid_cell.compress(VectorOperation::add);
  dof_is_in_fluid_cell.update_ghost_values();
  dof_is_in_non_fluid_cell.compress(VectorOperation::add);
  dof_is_in_non_fluid_cell.update_ghost_values();

  // The locally relevant degrees of freedom are iterated in increasing order,
  // which is the order of their position in the index set.
  dof_categories.assign(n_relevant_dofs, DoFCategory::fluid);
  types::global_dof_index dof_position = 0;
  for (const types::global_dof_index dof_index : locally_relevant_dofs)
    {
      const bool is_in_fluid_cell =
        std::as_const(dof_is_in_fluid_cell)(dof_index) > 0.;
      const bool is_in_non_fluid_cell =
        std::as_const(dof_is_in_non_fluid_cell)(dof_index) > 0.;
      if (is_in_non_fluid_cell)
        dof_categories[dof_position] =
          is_in_fluid_cell ? DoFCategory::interface : DoFCategory::extension;
      ++dof_position;
    }

  // Build the equation of every locally owned interface degree of freedom that
  // is not constrained. The stencil cell is a fluid cell that contains the
  // degree of freedom, so that the sparsity pattern already couples the degree
  // of freedom with the degrees of freedom of the stencil cell.
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cells;
  LetheGridTools::vertices_cell_mapping(dof_handler, vertices_to_cells);

  IBStencil<dim>            stencil;
  const std::vector<double> stencil_coefficients =
    stencil.coefficients(stencil_degree, length_ratio);
  const unsigned int n_stencil_points = stencil_coefficients.size();
  const unsigned int n_interpolation_points =
    stencil.number_of_interpolation_support_points(stencil_degree);

  // A stencil with a single coefficient evaluates the polynomial of the
  // stencil cell directly on the surface (see IBStencil). The normal
  // derivative is then the derivative of this polynomial. Otherwise, it is the
  // derivative of the polynomial of the 1D stencil.
  const bool stencil_uses_cell_directly = (n_stencil_points == 1);
  const std::vector<double> stencil_derivative_coefficients =
    (imposes_flux && !stencil_uses_cell_directly) ?
      stencil.derivative_coefficients(stencil_degree, length_ratio) :
      std::vector<double>();

  // Relative distance under which a degree of freedom is considered to be on
  // the surface of the particle, in which case the boundary condition is
  // imposed on the degree of freedom, as in the fluid dynamics solver.
  constexpr double relative_distance_tolerance = 1e-12;

  std::vector<bool>                    row_is_built(n_relevant_dofs, false);
  std::vector<types::global_dof_index> stencil_dof_indices(n_dofs_per_cell);
  std::vector<Point<dim>>              unit_stencil_points(n_stencil_points);
  types::global_dof_index              n_local_dofs_without_stencil = 0;
  interface_rows.clear();

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if ((!cell->is_locally_owned() && !cell->is_ghost()) ||
          cell_categories[cell->active_cell_index()] != CellCategory::cut)
        continue;

      fe_values_support_points.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      const auto   geometry_cell = get_geometry_cell(cell);
      const double cell_size     = cell->diameter();

      // The equation imposing the flux is divided by the diffusivity and
      // multiplied by the cell size to have the dimension of the solution, as
      // the equation imposing the value
      double scaling = cell->measure() * get_reference_rate(cell);
      if (imposes_flux)
        scaling *= cell_size / diffusivity;

      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          const types::global_dof_index dof_index = local_dof_indices[i];
          if (!locally_owned_dofs.is_element(dof_index) ||
              constraints.is_constrained(dof_index))
            continue;

          const types::global_dof_index dof_relevant_index =
            locally_relevant_dofs.index_within_set(dof_index);
          if (dof_categories[dof_relevant_index] != DoFCategory::interface ||
              row_is_built[dof_relevant_index])
            continue;
          row_is_built[dof_relevant_index] = true;

          const Point<dim> &dof_point =
            fe_values_support_points.quadrature_point(i);
          const IBParticle<dim> &particle =
            particles[dof_closest_particle[dof_relevant_index]];

          InterfaceRow row;
          row.dof_index = dof_index;
          row.scaling   = scaling;
          particle.closest_surface_point(dof_point,
                                         row.boundary_point,
                                         geometry_cell);

          const double distance_to_surface =
            dof_point.distance(row.boundary_point);
          const bool dof_is_on_surface =
            distance_to_surface <= relative_distance_tolerance * cell_size;

          // Normal to the surface pointing out of the particle, only used to
          // impose a flux. On the surface, it is the direction of the
          // gradient of the signed distance.
          Tensor<1, dim> normal;
          if (imposes_flux)
            {
              normal =
                dof_is_on_surface ?
                  particle.shape->gradient_with_cell_guess(dof_point,
                                                           geometry_cell) :
                  (dof_point - row.boundary_point) / distance_to_surface;
              const double normal_norm = normal.norm();
              if (normal_norm > 0.)
                normal /= normal_norm;
            }

          // Find the stencil cell: a fluid cell sharing a vertex with the cut
          // cell, containing the point used for the cell detection and the
          // degree of freedom. The farthest interpolation point must not be in
          // the cut cell, otherwise the stencil would extrapolate the solution
          // of the cut cell.
          bool stencil_cell_is_found = false;
          if (enable_extrapolation && !dof_is_on_surface)
            {
              const auto [stencil_point, interpolation_points] =
                stencil.support_points_for_interpolation(stencil_degree,
                                                         length_ratio,
                                                         particle,
                                                         dof_point,
                                                         geometry_cell);
              const Point<dim> detection_point =
                stencil.point_for_cell_detection(particle,
                                                 dof_point,
                                                 geometry_cell);

              if (!cell->point_inside(
                    interpolation_points[n_interpolation_points - 1]))
                {
                  for (const auto &candidate_cell :
                       LetheGridTools::find_cells_around_cell<dim>(
                         vertices_to_cells, cell))
                    {
                      if (cell_categories[candidate_cell
                                            ->active_cell_index()] !=
                            CellCategory::fluid ||
                          !candidate_cell->point_inside(detection_point))
                        continue;

                      candidate_cell->get_dof_indices(stencil_dof_indices);
                      if (std::ranges::find(stencil_dof_indices, dof_index) ==
                          stencil_dof_indices.end())
                        continue;

                      // The stencil points are expressed in the reference cell
                      // of the stencil cell. They can lie outside of it, in
                      // which case the polynomial of the stencil cell is
                      // extrapolated.
                      unit_stencil_points[0] =
                        mapping.transform_real_to_unit_cell(candidate_cell,
                                                            stencil_point);
                      for (unsigned int k = 1; k < n_stencil_points; ++k)
                        unit_stencil_points[k] =
                          mapping.transform_real_to_unit_cell(
                            candidate_cell, interpolation_points[k - 1]);

                      row.column_indices = stencil_dof_indices;
                      row.value_coefficients.assign(n_dofs_per_cell, 0.);
                      for (unsigned int m = 0; m < n_dofs_per_cell; ++m)
                        for (unsigned int k = 0; k < n_stencil_points; ++k)
                          row.value_coefficients[m] +=
                            stencil_coefficients[k] *
                            fe.shape_value(m, unit_stencil_points[k]);

                      if (imposes_flux && stencil_uses_cell_directly)
                        row.normal_derivative_coefficients =
                          compute_shape_function_directional_derivatives(
                            mapping,
                            candidate_cell,
                            row.boundary_point,
                            normal);
                      else if (imposes_flux)
                        {
                          // The coordinate of the 1D stencil increases towards
                          // the particle, and the distance between the degree
                          // of freedom and the surface is length_ratio on the
                          // 1D stencil.
                          const double derivative_factor =
                            -length_ratio / distance_to_surface;
                          row.normal_derivative_coefficients.assign(
                            n_dofs_per_cell, 0.);
                          for (unsigned int m = 0; m < n_dofs_per_cell; ++m)
                            for (unsigned int k = 0; k < n_stencil_points; ++k)
                              row.normal_derivative_coefficients[m] +=
                                derivative_factor *
                                stencil_derivative_coefficients[k] *
                                fe.shape_value(m, unit_stencil_points[k]);
                        }

                      stencil_cell_is_found = true;
                      break;
                    }
                }
            }

          if (!stencil_cell_is_found)
            {
              ++n_local_dofs_without_stencil;
              if (!imposes_flux)
                {
                  // The value is imposed on the degree of freedom
                  row.column_indices     = {dof_index};
                  row.value_coefficients = {1.};
                }
              else
                {
                  // The flux is computed with the value and the gradient of
                  // the solution at the degree of freedom, in a fluid cell that
                  // contains it. Such a cell exists, since the degree of
                  // freedom is an interface degree of freedom. If it is not
                  // found or if the normal is undefined, the equation of the
                  // fluid cells is kept, which imposes a zero flux.
                  if (normal.norm() <= 0.)
                    continue;

                  bool fluid_cell_is_found = false;
                  for (const auto &fluid_cell :
                       LetheGridTools::find_cells_around_cell<dim>(
                         vertices_to_cells, cell))
                    {
                      if (cell_categories[fluid_cell->active_cell_index()] !=
                          CellCategory::fluid)
                        continue;

                      fluid_cell->get_dof_indices(stencil_dof_indices);
                      const auto dof_iterator =
                        std::ranges::find(stencil_dof_indices, dof_index);
                      if (dof_iterator == stencil_dof_indices.end())
                        continue;

                      const auto local_dof_position = static_cast<unsigned int>(
                        std::distance(stencil_dof_indices.begin(),
                                      dof_iterator));
                      row.column_indices = stencil_dof_indices;
                      row.value_coefficients.assign(n_dofs_per_cell, 0.);
                      row.value_coefficients[local_dof_position] = 1.;
                      row.normal_derivative_coefficients =
                        compute_shape_function_directional_derivatives(
                          mapping, fluid_cell, dof_point, normal);
                      fluid_cell_is_found = true;
                      break;
                    }
                  if (!fluid_cell_is_found)
                    continue;
                }
            }

          interface_rows.push_back(std::move(row));
        }
    }

  n_interface_dofs = Utilities::MPI::sum(static_cast<types::global_dof_index>(
                                           interface_rows.size()),
                                         mpi_communicator);
  n_dofs_without_stencil =
    Utilities::MPI::sum(n_local_dofs_without_stencil, mpi_communicator);
}

template <int dim>
void
IBScalarConstraints<dim>::assemble_matrix(
  const DoFHandler<dim>           &dof_handler,
  const Mapping<dim>              &mapping,
  const Quadrature<dim>           &quadrature,
  const AffineConstraints<double> &constraints,
  const GlobalVectorType          &evaluation_point,
  TrilinosWrappers::SparseMatrix  &system_matrix) const
{
  const FiniteElement<dim> &fe              = dof_handler.get_fe();
  const unsigned int        n_dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int        n_q_points      = quadrature.size();

  FEValues<dim>      fe_values(mapping,
                          fe,
                          quadrature,
                          update_gradients | update_JxW_values);
  FullMatrix<double> local_matrix(n_dofs_per_cell, n_dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // Laplace equation of the extension degrees of freedom, assembled on the cut
  // cells and on the cells inside the particles. The rows of the other degrees
  // of freedom are left empty, so that the equations of the fluid and interface
  // degrees of freedom are not modified.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() || is_fluid_cell(cell))
        continue;

      fe_values.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
      const double extension_diffusivity =
        get_reference_rate(cell) * cell->diameter() * cell->diameter();

      local_matrix = 0;
      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          if (get_dof_category(local_dof_indices[i]) != DoFCategory::extension)
            continue;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double          JxW        = fe_values.JxW(q);
              const Tensor<1, dim> &grad_phi_i = fe_values.shape_grad(i, q);
              for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                {
                  // term: D_ext * \nabla \phi_i \cdot \nabla \phi_j
                  local_matrix(i, j) += extension_diffusivity * grad_phi_i *
                                        fe_values.shape_grad(j, q) * JxW;
                }
            }
        }
      constraints.distribute_local_to_global(local_matrix,
                                             local_dof_indices,
                                             system_matrix);
    }
  system_matrix.compress(VectorOperation::add);

  // Replace the equation of the interface degrees of freedom by their stencil.
  FullMatrix<double>                   row_matrix;
  std::vector<types::global_dof_index> row_indices(1);
  for (const auto &row : interface_rows)
    {
      system_matrix.clear_row(row.dof_index);

      row_indices[0] = row.dof_index;
      row_matrix.reinit(1, row.column_indices.size());
      if (!imposes_flux)
        {
          // term: s * v_m
          for (unsigned int m = 0; m < row.column_indices.size(); ++m)
            row_matrix(0, m) = row.scaling * row.value_coefficients[m];
        }
      else
        {
          // Derivative of the imposed flux with respect to the value on the
          // surface, which is zero for the neumann type
          double flux_derivative = 0.;
          if (boundary_condition.type ==
              Parameters::IBTracerBoundaryType::robin)
            {
              double value_at_boundary = 0.;
              for (unsigned int m = 0; m < row.column_indices.size(); ++m)
                value_at_boundary += row.value_coefficients[m] *
                                     evaluation_point(row.column_indices[m]);
              flux_derivative = compute_surface_reaction_rate(
                                  value_at_boundary,
                                  boundary_condition.reaction_constant,
                                  boundary_condition.reaction_order,
                                  boundary_condition.reaction_threshold)
                                  .second;
            }

          // term: s * (dJ/du * v_m - D * g_m)
          for (unsigned int m = 0; m < row.column_indices.size(); ++m)
            row_matrix(0, m) =
              row.scaling *
              (flux_derivative * row.value_coefficients[m] -
               diffusivity * row.normal_derivative_coefficients[m]);
        }

      // The constrained columns are distributed to their master degrees of
      // freedom, which share a cell with the interface degree of freedom.
      constraints.distribute_local_to_global(row_matrix,
                                             row_indices,
                                             row.column_indices,
                                             system_matrix);
    }
  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
IBScalarConstraints<dim>::assemble_rhs(
  const DoFHandler<dim>           &dof_handler,
  const Mapping<dim>              &mapping,
  const Quadrature<dim>           &quadrature,
  const AffineConstraints<double> &constraints,
  const GlobalVectorType          &evaluation_point,
  GlobalVectorType                &system_rhs) const
{
  const FiniteElement<dim> &fe              = dof_handler.get_fe();
  const unsigned int        n_dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int        n_q_points      = quadrature.size();

  FEValues<dim>                        fe_values(mapping,
                          fe,
                          quadrature,
                          update_gradients | update_JxW_values);
  Vector<double>                       local_rhs(n_dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);
  std::vector<Tensor<1, dim>>          solution_gradients(n_q_points);

  // Residual of the Laplace equation of the extension degrees of freedom.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned() || is_fluid_cell(cell))
        continue;

      fe_values.reinit(cell);
      fe_values.get_function_gradients(evaluation_point, solution_gradients);
      cell->get_dof_indices(local_dof_indices);
      const double extension_diffusivity =
        get_reference_rate(cell) * cell->diameter() * cell->diameter();

      local_rhs = 0;
      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          if (get_dof_category(local_dof_indices[i]) != DoFCategory::extension)
            continue;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // term: -D_ext * \nabla \phi_i \cdot \nabla u
              local_rhs(i) -= extension_diffusivity *
                              fe_values.shape_grad(i, q) *
                              solution_gradients[q] * fe_values.JxW(q);
            }
        }
      constraints.distribute_local_to_global(local_rhs,
                                             local_dof_indices,
                                             system_rhs);
    }
  system_rhs.compress(VectorOperation::add);

  // Residual of the stencil equation of the interface degrees of freedom.
  for (const auto &row : interface_rows)
    {
      double value_at_boundary = 0.;
      for (unsigned int m = 0; m < row.column_indices.size(); ++m)
        value_at_boundary +=
          row.value_coefficients[m] * evaluation_point(row.column_indices[m]);

      double residual = 0.;
      if (!imposes_flux)
        {
          // u_Gamma - u_D
          residual = value_at_boundary -
                     boundary_condition.dirichlet->value(row.boundary_point);
        }
      else
        {
          double normal_derivative_at_boundary = 0.;
          for (unsigned int m = 0; m < row.column_indices.size(); ++m)
            normal_derivative_at_boundary +=
              row.normal_derivative_coefficients[m] *
              evaluation_point(row.column_indices[m]);

          // Imposed flux J into the particle: q for the neumann type, and
          // alpha * C_eff^(n-1) * u_Gamma + q for the robin type
          double imposed_flux =
            boundary_condition.flux->value(row.boundary_point);
          if (boundary_condition.type ==
              Parameters::IBTracerBoundaryType::robin)
            imposed_flux += compute_surface_reaction_rate(
                              value_at_boundary,
                              boundary_condition.reaction_constant,
                              boundary_condition.reaction_order,
                              boundary_condition.reaction_threshold)
                              .first;

          // J - D * du/dn. With this sign, the diagonal coefficient is
          // usually positive, as for the equations assembled on the fluid
          // cells, since the normal derivative on the surface decreases when
          // the value of the degree of freedom increases.
          residual = imposed_flux - diffusivity * normal_derivative_at_boundary;
        }

      system_rhs(row.dof_index) = -row.scaling * residual;
    }
  system_rhs.compress(VectorOperation::insert);
}

template <int dim>
typename IBScalarConstraints<dim>::DoFCategory
IBScalarConstraints<dim>::get_dof_category(
  const types::global_dof_index dof_index) const
{
  const types::global_dof_index relevant_index =
    locally_relevant_dofs.index_within_set(dof_index);
  AssertIndexRange(relevant_index, dof_categories.size());
  return dof_categories[relevant_index];
}

template <int dim>
double
IBScalarConstraints<dim>::get_reference_rate(
  const typename DoFHandler<dim>::active_cell_iterator &cell) const
{
  const double cell_size                 = cell->diameter();
  const double inverse_squared_cell_size = 1. / (cell_size * cell_size);
  const double reference_rate =
    inverse_time_step + diffusivity * inverse_squared_cell_size;

  // In a steady simulation without diffusivity, a unit diffusivity is used.
  return (reference_rate > 0.) ? reference_rate : inverse_squared_cell_size;
}

template class IBScalarConstraints<2>;
template class IBScalarConstraints<3>;
