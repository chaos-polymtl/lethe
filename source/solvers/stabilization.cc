// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/lethe_grid_tools.h>

#include <solvers/stabilization.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <vector>

template <int dim>
void
moe_scalar_limiter(const DoFHandler<dim> &dof_handler,
                   const Mapping<dim>    &mapping,
                   const Quadrature<dim> &cell_quadrature,
                   GlobalVectorType      &locally_relevant_vector,
                   GlobalVectorType      &locally_owned_vector)
{
  // The limiter requires a vector in which the value can be modified.
  // For this we use the local_evaluation_point. We thus begin by making sure
  // that the local_evaluation_point is set at the present value of the solution
  // to limit.
  locally_owned_vector = locally_relevant_vector;

  // We define the cutoff function from the Moe paper using a lambda function.
  // The cutoff factor is larger than one so that the cutoff function is smaller
  // than one when a cell extremum is equal to its bound. According to the Moe
  // paper, this is necessary to suppress the oscillations in more than one
  // dimension.
  const double cutoff_factor = 1.1;
  auto         phi_limit     = [cutoff_factor](const double y) {
    return std::min(1., y / cutoff_factor);
  };

  // A variation of the solution within a cell that is below this fraction of
  // the magnitude of the solution is round-off noise. It is not limited since
  // the ratio between the bounds and such a variation is meaningless.
  const double relative_variation_tolerance =
    100. * std::numeric_limits<double>::epsilon();

  // The strategy for limiting requires looping over the cells twice.
  // First loop over all active cells (local and ghost):
  // 1. Calculate the max and the min value of the field we wish to limit
  // Loop over every local cells:
  // 2. Calculate approximate upper and lower bounds using the neighbors
  // 3. Calculate the value of the theta limiting parameter
  // 4. Rescale the nodal values of the solution using the average value of the
  // solution within the element and the theta limiter.

  // We need the vertices-to-cell map to have access to the neighbors rapidly
  // We obtain this map using the LetheGridTools functionnalities.
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;
  LetheGridTools::vertices_cell_mapping(dof_handler, vertices_to_cell);

  // The cells that share a vertex through a periodic boundary are also
  // neighbors. They are stored in a different map.
  const bool has_periodic_faces =
    !dof_handler.get_triangulation().get_periodic_face_map().empty();
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_periodic_cell;
  if (has_periodic_faces)
    LetheGridTools::vertices_cell_mapping_with_periodic_boundaries(
      dof_handler, vertices_to_periodic_cell);

  // The values of the degrees of freedom are used as the values of the solution
  // at the support points. This requires a finite element with support points.
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  Assert(fe.has_support_points(),
         ExcMessage(
           "The Moe limiter requires a finite element with support points."));

  // The mean is obtained by integrating the solution over the cell. The
  // solution at the quadrature points is also used to sample the extrema.
  FEValues<dim>       fe_values(mapping,
                          fe,
                          cell_quadrature,
                          update_values | update_JxW_values);
  const unsigned int  n_q_points = cell_quadrature.size();
  std::vector<double> values_q(n_q_points);

  // Step 1. loop over every cell and calculate the max, the min and the mean.
  // These values are stored in std::map for rapid access.
  std::map<typename ::dealii::types::global_cell_index, double>
    max_value_per_cells;
  std::map<typename ::dealii::types::global_cell_index, double>
    min_value_per_cells;
  std::map<typename ::dealii::types::global_cell_index, double>
    mean_value_per_cells;

  std::vector<types::global_dof_index> local_dof_indices(
    dof_handler.get_fe().n_dofs_per_cell());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // We need to loop over cells that are either locally_owned or ghost to
      // gather all the neighbors.
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          double max_value = std::numeric_limits<double>::lowest();
          double min_value = std::numeric_limits<double>::max();
          for (const auto &i_dof : local_dof_indices)
            {
              // Get the max and the min value of the solution at the support
              // points
              const double dof_value = locally_relevant_vector(i_dof);
              max_value              = std::max(max_value, dof_value);
              min_value              = std::min(min_value, dof_value);
            }

          // Calculate the mean solution \bar{w} by integrating the solution
          // over the cell. The arithmetic average of the nodal values is not
          // the mean of the solution for degrees above one, nor on non-affine
          // cells, and rescaling around it does not conserve the integral of
          // the solution.
          fe_values.reinit(cell);
          fe_values.get_function_values(locally_relevant_vector, values_q);
          double cell_integral = 0;
          double cell_measure  = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // The extrema of the solution of degree above one are not
              // necessarily located at the support points.
              max_value = std::max(max_value, values_q[q]);
              min_value = std::min(min_value, values_q[q]);

              cell_integral += values_q[q] * fe_values.JxW(q);
              cell_measure += fe_values.JxW(q);
            }
          const double mean_value = cell_integral / cell_measure;
          min_value_per_cells[cell->global_active_cell_index()]  = min_value;
          max_value_per_cells[cell->global_active_cell_index()]  = max_value;
          mean_value_per_cells[cell->global_active_cell_index()] = mean_value;
        }
    }

  // Loop over every cell and do:
  // 2. Calculate approximate upper and lower bounds using the neighbhors
  // 3. Calculate the theta limiter
  // 4. Rescale the nodal values of the solution

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices(local_dof_indices);

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);

          // Add the neighbors through the periodic boundaries. Only the cells
          // for which the solution is available are kept.
          if (has_periodic_faces)
            for (const auto &periodic_neighbor :
                 LetheGridTools::find_cells_around_cell<dim>(
                   vertices_to_periodic_cell, cell))
              if (periodic_neighbor->is_locally_owned() ||
                  periodic_neighbor->is_ghost())
                active_neighbors.push_back(periodic_neighbor);

          // 2. Calculate approximate upper and lower bounds using the
          // neighbhors

          // Here we assume that alpha = 0 (see the Moe paper). The definition
          // of alpha in the original paper is flaky in my opinion, let's use
          // the more robust version as a starting point.
          double upper_bound =
            mean_value_per_cells.at(cell->global_active_cell_index());
          double lower_bound =
            mean_value_per_cells.at(cell->global_active_cell_index());


          for (const auto &neighbor : active_neighbors)
            {
              // The shock capture mechanism theta calculation should not
              // consider the cell itself.
              if (cell->global_active_cell_index() ==
                  neighbor->global_active_cell_index())
                continue;

              // The original paper uses the max and the min of the neighbor
              // solutions. I have also tested with the mean values of the
              // neighbor. It seemed a bit more robust, but I'd rather follow
              // the paper for now. Let's just keep in mind that this can be a
              // fallback.

              upper_bound = std::max(upper_bound,
                                     max_value_per_cells.at(
                                       neighbor->global_active_cell_index()));

              lower_bound = std::min(lower_bound,
                                     min_value_per_cells.at(
                                       neighbor->global_active_cell_index()));
            }

          // 3. Calculate the value of the theta limiting using the max, min and
          // mean values as well as the bounds.
          const double max_value =
            max_value_per_cells.at(cell->global_active_cell_index());
          const double min_value =
            min_value_per_cells.at(cell->global_active_cell_index());
          const double mean_value =
            mean_value_per_cells.at(cell->global_active_cell_index());

          // The distances between the mean and the extrema of the cell, and
          // between the mean and the bounds, are all positive. The ratio is
          // only calculated when the variation within the cell is above
          // round-off. Otherwise there is nothing to limit.
          const double max_variation = max_value - mean_value;
          const double min_variation = mean_value - min_value;
          const double variation_tolerance =
            relative_variation_tolerance *
            std::max(std::abs(max_value), std::abs(min_value));

          double theta = 1.;
          if (max_variation > variation_tolerance)
            theta =
              std::min(theta,
                       phi_limit((upper_bound - mean_value) / max_variation));
          if (min_variation > variation_tolerance)
            theta =
              std::min(theta,
                       phi_limit((mean_value - lower_bound) / min_variation));

          // The rescaling must never invert the solution around its mean
          theta = std::max(theta, 0.);

          // 4. Rescale the solution within the element
          for (const auto &i : local_dof_indices)
            {
              // Get the value of the solution at the DOFs and rescale it
              // using the mean value and the theta limiter.
              const double dof_value = locally_relevant_vector(i);
              locally_owned_vector(i) =
                mean_value + theta * (dof_value - mean_value);
            }
        }
    }

  // Reset the present solution and the evaluation point to the new updated
  // solution.
  locally_relevant_vector = locally_owned_vector;
}


template void
moe_scalar_limiter(const DoFHandler<2> &dof_handler,
                   const Mapping<2>    &mapping,
                   const Quadrature<2> &cell_quadrature,
                   GlobalVectorType    &locally_relevant_vector,
                   GlobalVectorType    &locally_owned_vector);

template void
moe_scalar_limiter(const DoFHandler<3> &dof_handler,
                   const Mapping<3>    &mapping,
                   const Quadrature<3> &cell_quadrature,
                   GlobalVectorType    &locally_relevant_vector,
                   GlobalVectorType    &locally_owned_vector);


template <int dim>
void
kuzmin_scalar_limiter(const DoFHandler<dim> &dof_handler,
                      const Mapping<dim>    &mapping,
                      const Quadrature<dim> &cell_quadrature,
                      GlobalVectorType      &locally_relevant_vector,
                      GlobalVectorType      &locally_owned_vector)
{
  // The limiter requires a vector in which the value can be modified. We begin
  // by making sure that it is set at the present value of the solution to
  // limit.
  locally_owned_vector = locally_relevant_vector;

  // The values of the degrees of freedom are used as the values of the solution
  // at the support points. This requires a finite element with support points.
  const FiniteElement<dim> &fe = dof_handler.get_fe();
  Assert(
    fe.has_support_points(),
    ExcMessage(
      "The Kuzmin limiter requires a finite element with support points."));
  Assert(
    fe.reference_cell().is_hyper_cube(),
    ExcMessage(
      "The Kuzmin limiter is only implemented for quadrilateral and hexahedral cells."));

  // A piecewise constant solution has nothing to limit.
  if (fe.degree == 0)
    return;

  // The hierarchical limiting requires second derivatives. For a degree of
  // one, a single factor is calculated from the values at the vertices.
  const bool is_hierarchical = fe.degree >= 2;

  // A variation that is below this fraction of the magnitude of the quantities
  // it is compared to is round-off noise. It is not limited since the ratio
  // between the bounds and such a variation is meaningless.
  const double relative_variation_tolerance =
    100. * std::numeric_limits<double>::epsilon();

  // Largest factor between zero and one by which the variation of a quantity
  // between the centroid of a cell and one of its vertices can be multiplied
  // for the value at the vertex to remain within the bounds of this vertex.
  // The bounds always include the value at the centroid.
  const auto calculate_correction_factor =
    [relative_variation_tolerance](const double centroid_value,
                                   const double vertex_value,
                                   const double lower_bound,
                                   const double upper_bound) {
      const double variation = vertex_value - centroid_value;
      const double variation_tolerance =
        relative_variation_tolerance * std::max({std::abs(centroid_value),
                                                 std::abs(vertex_value),
                                                 std::abs(lower_bound),
                                                 std::abs(upper_bound)});
      if (variation > variation_tolerance)
        return std::clamp((upper_bound - centroid_value) / variation, 0., 1.);
      if (variation < -variation_tolerance)
        return std::clamp((lower_bound - centroid_value) / variation, 0., 1.);
      return 1.;
    };

  // The strategy for limiting requires looping over the cells twice.
  // First loop over all active cells (local and ghost):
  // 1. Calculate the mean of the solution and, for the hierarchical limiting,
  // its gradient and second derivatives.
  // 2. Calculate the bounds at every vertex using the cells that share it.
  // Loop over every local cells:
  // 3. Calculate the correction factors using the bounds at the vertices.
  // 4. Rescale the linear part and the remainder of the solution.

  const auto        &triangulation   = dof_handler.get_triangulation();
  const unsigned int n_active_cells  = triangulation.n_active_cells();
  const unsigned int n_vertices      = triangulation.n_vertices();
  const unsigned int n_dofs_per_cell = fe.n_dofs_per_cell();
  const unsigned int n_q_points      = cell_quadrature.size();

  FEValues<dim>       fe_values(mapping,
                          fe,
                          cell_quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  std::vector<double> values_q(n_q_points);

  // Values stored for every cell and accessed with the active cell index.
  std::vector<double>                  mean_per_cell(n_active_cells, 0.);
  std::vector<Point<dim>>              centroid_per_cell(n_active_cells);
  std::vector<Tensor<1, dim>>          gradient_per_cell(n_active_cells);
  std::vector<SymmetricTensor<2, dim>> hessian_per_cell(
    is_hierarchical ? n_active_cells : 0);

  // Bounds stored for every vertex and accessed with the vertex index.
  const double   largest_value = std::numeric_limits<double>::max();
  Tensor<1, dim> largest_gradient;
  for (unsigned int d = 0; d < dim; ++d)
    largest_gradient[d] = largest_value;

  std::vector<double>         min_mean_per_vertex(n_vertices, largest_value);
  std::vector<double>         max_mean_per_vertex(n_vertices, -largest_value);
  std::vector<Tensor<1, dim>> min_gradient_per_vertex(n_vertices,
                                                      largest_gradient);
  std::vector<Tensor<1, dim>> max_gradient_per_vertex(n_vertices,
                                                      -largest_gradient);
  std::vector<bool>           vertex_is_at_boundary(n_vertices, false);

  // Extend the bounds of a vertex with the mean and the gradient of a cell.
  const auto update_vertex_bounds = [&](const unsigned int vertex_index,
                                        const unsigned int cell_index) {
    min_mean_per_vertex[vertex_index] =
      std::min(min_mean_per_vertex[vertex_index], mean_per_cell[cell_index]);
    max_mean_per_vertex[vertex_index] =
      std::max(max_mean_per_vertex[vertex_index], mean_per_cell[cell_index]);
    for (unsigned int d = 0; d < dim; ++d)
      {
        min_gradient_per_vertex[vertex_index][d] =
          std::min(min_gradient_per_vertex[vertex_index][d],
                   gradient_per_cell[cell_index][d]);
        max_gradient_per_vertex[vertex_index][d] =
          std::max(max_gradient_per_vertex[vertex_index][d],
                   gradient_per_cell[cell_index][d]);
      }
  };

  // The gradient and the second derivatives are those of the L2 projection of
  // the solution onto the polynomials of degree two. The basis is made of 1,
  // the coordinates r_i relative to the centroid, and the products r_i r_j / 2
  // (i = j) and r_i r_j (i < j), so that the coefficients are the mean value of
  // the projection at the centroid, its gradient at the centroid and its second
  // derivatives. The coordinates are divided by the size of the cell to obtain
  // a well-conditioned system.
  const unsigned int  n_basis_functions = 1 + dim + (dim * (dim + 1)) / 2;
  FullMatrix<double>  gram_matrix(n_basis_functions, n_basis_functions);
  Vector<double>      projection_rhs(n_basis_functions);
  Vector<double>      projection_coefficients(n_basis_functions);
  std::vector<double> basis_values(n_basis_functions);
  Assert(
    !is_hierarchical || n_q_points >= n_basis_functions,
    ExcMessage(
      "The quadrature does not have enough points to project the solution onto the polynomials of degree two."));

  // Step 1. Loop over every cell and calculate the mean, the gradient and the
  // second derivatives.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // We need to loop over cells that are either locally_owned or ghost to
      // gather all the cells that share a vertex with a locally owned cell.
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          const unsigned int cell_index = cell->active_cell_index();

          fe_values.reinit(cell);
          fe_values.get_function_values(locally_relevant_vector, values_q);

          // Mean of the solution and centroid of the cell
          double         cell_integral = 0;
          double         cell_measure  = 0;
          Tensor<1, dim> first_moment;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              const double JxW = fe_values.JxW(q);
              cell_integral += values_q[q] * JxW;
              cell_measure += JxW;
              first_moment += JxW * fe_values.quadrature_point(q);
            }
          mean_per_cell[cell_index] = cell_integral / cell_measure;
          centroid_per_cell[cell_index] =
            Point<dim>(first_moment / cell_measure);

          if (is_hierarchical)
            {
              const double cell_size = cell->diameter();

              gram_matrix    = 0;
              projection_rhs = 0;
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  const double         JxW = fe_values.JxW(q);
                  const Tensor<1, dim> relative_position =
                    (fe_values.quadrature_point(q) -
                     centroid_per_cell[cell_index]) /
                    cell_size;

                  basis_values[0] = 1.;
                  for (unsigned int i = 0; i < dim; ++i)
                    basis_values[1 + i] = relative_position[i];
                  unsigned int k = 1 + dim;
                  for (unsigned int i = 0; i < dim; ++i)
                    for (unsigned int j = i; j < dim; ++j)
                      basis_values[k++] = (i == j ? 0.5 : 1.) *
                                          relative_position[i] *
                                          relative_position[j];

                  for (unsigned int a = 0; a < n_basis_functions; ++a)
                    {
                      projection_rhs(a) += values_q[q] * basis_values[a] * JxW;
                      for (unsigned int b = 0; b < n_basis_functions; ++b)
                        gram_matrix(a, b) +=
                          basis_values[a] * basis_values[b] * JxW;
                    }
                }
              gram_matrix.gauss_jordan();
              gram_matrix.vmult(projection_coefficients, projection_rhs);

              // The coefficients are those of the scaled coordinates
              for (unsigned int i = 0; i < dim; ++i)
                gradient_per_cell[cell_index][i] =
                  projection_coefficients(1 + i) / cell_size;
              unsigned int k = 1 + dim;
              for (unsigned int i = 0; i < dim; ++i)
                for (unsigned int j = i; j < dim; ++j)
                  hessian_per_cell[cell_index][i][j] =
                    projection_coefficients(k++) / (cell_size * cell_size);
            }

          // Step 2. The cell contributes to the bounds of its vertices.
          for (unsigned int v = 0; v < cell->n_vertices(); ++v)
            update_vertex_bounds(cell->vertex_index(v), cell_index);

          // Identify the vertices on the boundaries of the domain. The
          // periodic boundaries are not boundaries of the domain.
          for (const auto f : cell->face_indices())
            if (cell->at_boundary(f) && !cell->has_periodic_neighbor(f))
              for (unsigned int v = 0; v < cell->face(f)->n_vertices(); ++v)
                vertex_is_at_boundary[cell->face(f)->vertex_index(v)] = true;
        }
    }

  // A hanging vertex is located on a face of a coarser cell without being one
  // of its vertices. This coarser cell is added to the bounds of all the
  // vertices of the face it shares with a finer cell.
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          for (const auto f : cell->face_indices())
            {
              if (cell->at_boundary(f) || !cell->neighbor_is_coarser(f))
                continue;

              const auto neighbor = cell->neighbor(f);
              if (neighbor->is_locally_owned() || neighbor->is_ghost())
                for (unsigned int v = 0; v < cell->face(f)->n_vertices(); ++v)
                  update_vertex_bounds(cell->face(f)->vertex_index(v),
                                       neighbor->active_cell_index());
            }
        }
    }

  // The vertices that coincide through a periodic boundary share the same
  // bounds.
  if (!triangulation.get_periodic_face_map().empty())
    {
      std::map<unsigned int, std::vector<unsigned int>>
                                           coinciding_vertex_groups;
      std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;
      GridTools::collect_coinciding_vertices(triangulation,
                                             coinciding_vertex_groups,
                                             vertex_to_coinciding_vertex_group);

      for (const auto &[group, coinciding_vertices] : coinciding_vertex_groups)
        {
          double         min_mean     = largest_value;
          double         max_mean     = -largest_value;
          Tensor<1, dim> min_gradient = largest_gradient;
          Tensor<1, dim> max_gradient = -largest_gradient;
          for (const auto &vertex_index : coinciding_vertices)
            {
              min_mean = std::min(min_mean, min_mean_per_vertex[vertex_index]);
              max_mean = std::max(max_mean, max_mean_per_vertex[vertex_index]);
              for (unsigned int d = 0; d < dim; ++d)
                {
                  min_gradient[d] =
                    std::min(min_gradient[d],
                             min_gradient_per_vertex[vertex_index][d]);
                  max_gradient[d] =
                    std::max(max_gradient[d],
                             max_gradient_per_vertex[vertex_index][d]);
                }
            }
          for (const auto &vertex_index : coinciding_vertices)
            {
              min_mean_per_vertex[vertex_index]     = min_mean;
              max_mean_per_vertex[vertex_index]     = max_mean;
              min_gradient_per_vertex[vertex_index] = min_gradient;
              max_gradient_per_vertex[vertex_index] = max_gradient;
            }
        }
    }

  // For a degree of one, the values at the vertices are those of the degrees
  // of freedom whose support points are the vertices of the reference cell.
  const double                   support_point_tolerance = 1e-10;
  const std::vector<Point<dim>> &unit_support_points =
    fe.get_unit_support_points();
  std::vector<unsigned int> vertex_to_dof(GeometryInfo<dim>::vertices_per_cell,
                                          numbers::invalid_unsigned_int);
  if (!is_hierarchical)
    for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
          if (unit_support_points[i].distance(
                GeometryInfo<dim>::unit_cell_vertex(v)) <
              support_point_tolerance)
            vertex_to_dof[v] = i;
        Assert(
          vertex_to_dof[v] != numbers::invalid_unsigned_int,
          ExcMessage(
            "The finite element does not have a support point at every vertex."));
      }

  // The linear part of the solution is evaluated at the support points.
  FEValues<dim> fe_values_support_points(mapping,
                                         fe,
                                         Quadrature<dim>(unit_support_points),
                                         update_quadrature_points);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);
  std::vector<double>                  dof_values(n_dofs_per_cell);
  std::vector<double>                  limited_values(n_dofs_per_cell);

  // Loop over every cell and do:
  // 3. Calculate the correction factors
  // 4. Rescale the linear part and the remainder of the solution
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int    cell_index = cell->active_cell_index();
          const double          mean_value = mean_per_cell[cell_index];
          const Point<dim>     &centroid   = centroid_per_cell[cell_index];
          const Tensor<1, dim> &gradient   = gradient_per_cell[cell_index];

          cell->get_dof_indices(local_dof_indices);
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
            dof_values[i] = locally_relevant_vector(local_dof_indices[i]);

          // 3. Calculate the correction factor of the linear part (alpha_1)
          // and of the remainder (alpha_2).
          double alpha_1 = 1.;
          double alpha_2 = 1.;
          for (unsigned int v = 0; v < cell->n_vertices(); ++v)
            {
              const unsigned int vertex_index = cell->vertex_index(v);
              if (vertex_is_at_boundary[vertex_index])
                continue;

              if (is_hierarchical)
                {
                  const Tensor<1, dim> vertex_position =
                    cell->vertex(v) - centroid;

                  // Linear part of the solution at the vertex
                  alpha_1 = std::min(alpha_1,
                                     calculate_correction_factor(
                                       mean_value,
                                       mean_value + gradient * vertex_position,
                                       min_mean_per_vertex[vertex_index],
                                       max_mean_per_vertex[vertex_index]));

                  // Linear reconstruction of the gradient at the vertex
                  const Tensor<1, dim> vertex_gradient =
                    gradient + hessian_per_cell[cell_index] * vertex_position;
                  for (unsigned int d = 0; d < dim; ++d)
                    alpha_2 =
                      std::min(alpha_2,
                               calculate_correction_factor(
                                 gradient[d],
                                 vertex_gradient[d],
                                 min_gradient_per_vertex[vertex_index][d],
                                 max_gradient_per_vertex[vertex_index][d]));
                }
              else
                {
                  // Solution at the vertex
                  alpha_2 = std::min(alpha_2,
                                     calculate_correction_factor(
                                       mean_value,
                                       dof_values[vertex_to_dof[v]],
                                       min_mean_per_vertex[vertex_index],
                                       max_mean_per_vertex[vertex_index]));
                }
            }

          // The linear part is not limited more than the remainder. This is
          // what preserves the smooth extrema. For a degree of one, the single
          // factor is applied to the whole deviation from the mean.
          alpha_1 = is_hierarchical ? std::max(alpha_1, alpha_2) : alpha_2;

          // Nothing to do if the solution is not limited within the cell
          if (alpha_1 >= 1. && alpha_2 >= 1.)
            continue;

          // 4. Rescale the linear part and the remainder of the solution.
          if (is_hierarchical)
            fe_values_support_points.reinit(cell);
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
            {
              const double linear_part =
                is_hierarchical ?
                  gradient *
                    (fe_values_support_points.quadrature_point(i) - centroid) :
                  0.;
              limited_values[i] =
                mean_value + alpha_1 * linear_part +
                alpha_2 * (dof_values[i] - mean_value - linear_part);
            }

          // The linear part has a mean of zero when the mapping can be
          // represented by the finite element. To conserve the integral of the
          // solution in every case, the limited solution is shifted to recover
          // the mean of the solution.
          fe_values.reinit(cell);
          double limited_integral = 0;
          double cell_measure     = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              double limited_value_q = 0;
              for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                limited_value_q +=
                  limited_values[i] * fe_values.shape_value(i, q);
              limited_integral += limited_value_q * fe_values.JxW(q);
              cell_measure += fe_values.JxW(q);
            }
          const double mean_correction =
            mean_value - limited_integral / cell_measure;

          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
            locally_owned_vector(local_dof_indices[i]) =
              limited_values[i] + mean_correction;
        }
    }

  // Reset the present solution and the evaluation point to the new updated
  // solution.
  locally_relevant_vector = locally_owned_vector;
}


template void
kuzmin_scalar_limiter(const DoFHandler<2> &dof_handler,
                      const Mapping<2>    &mapping,
                      const Quadrature<2> &cell_quadrature,
                      GlobalVectorType    &locally_relevant_vector,
                      GlobalVectorType    &locally_owned_vector);

template void
kuzmin_scalar_limiter(const DoFHandler<3> &dof_handler,
                      const Mapping<3>    &mapping,
                      const Quadrature<3> &cell_quadrature,
                      GlobalVectorType    &locally_relevant_vector,
                      GlobalVectorType    &locally_owned_vector);
