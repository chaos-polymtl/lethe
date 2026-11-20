// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/lethe_grid_tools.h>

#include <fem-dem/sharp_ib_stencil_search.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/fe/mapping_q1.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>

template <int dim>
void
SharpIBStencilSearch<dim>::clear()
{
  results.clear();
  index.clear();
  cells.clear();
  index_ready = false;
}

template <int dim>
void
SharpIBStencilSearch<dim>::reinit(const dealii::DoFHandler<dim> &handler)
{
  mesh_change_connection.disconnect();
  mesh_movement_connection.disconnect();
  repartition_connection.disconnect();
  clear();
  dof_handler = &handler;
  mesh_change_connection =
    handler.get_triangulation().signals.any_change.connect(
      [this]() { clear(); });
  mesh_movement_connection =
    handler.get_triangulation().signals.mesh_movement.connect(
      [this]() { clear(); });
  repartition_connection =
    handler.get_triangulation().signals.post_distributed_repartition.connect(
      [this]() { clear(); });
}

template <int dim>
void
SharpIBStencilSearch<dim>::build_index()
{
  std::vector<BoxEntry> boxes;
  for (const auto &cell : dof_handler->active_cell_iterators())
    if (!cell->is_artificial())
      {
        auto bounds = cell->bounding_box().get_boundary_points();
        // Include points just outside a face that satisfy the reference-cell
        // tolerance. Q1 coordinate derivatives are bounded by vertex spans.
        for (int d = 0; d < dim; ++d)
          {
            const double padding =
              dim * 1e-12 * cell->diameter() +
              16 * std::numeric_limits<double>::epsilon() *
                std::max(
                  {1., std::abs(bounds.first[d]), std::abs(bounds.second[d])});
            bounds.first[d] -= padding;
            bounds.second[d] += padding;
          }
        boxes.emplace_back(dealii::BoundingBox<dim>(bounds), cells.size());
        cells.push_back(cell);
      }
  index       = dealii::pack_rtree(boxes);
  index_ready = true;
}

template <int dim>
typename SharpIBStencilSearch<dim>::Cell
SharpIBStencilSearch<dim>::find(
  const dealii::Point<dim>               &point,
  const Cell                             &seed,
  std::map<unsigned int, std::set<Cell>> &vertices_to_cell)
{
  AssertThrow(dof_handler != nullptr,
              dealii::ExcMessage("Stencil search is not initialized"));
  std::array<double, dim> coordinates;
  for (int d = 0; d < dim; ++d)
    {
      AssertThrow(std::isfinite(point[d]),
                  dealii::ExcMessage("Non-finite stencil query"));
      coordinates[d] = point[d];
    }
  const auto key    = std::make_pair(seed, coordinates);
  const auto cached = results.find(key);
  if (cached != results.end())
    {
      return cached->second;
    }
  const auto neighbors =
    LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, seed);
  for (const auto &cell : neighbors)
    if (!cell->is_artificial() && cell->point_inside(point))
      return results.emplace(key, cell).first->second;

  if (!index_ready)
    build_index();
  std::vector<BoxEntry> candidates;
  index.query(boost::geometry::index::intersects(point),
              std::back_inserter(candidates));
  // R-tree traversal order is not a geometric tie-breaker. Use mesh iteration
  // order, as the original level scan does on an unrefined mesh.
  std::sort(candidates.begin(),
            candidates.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });
  dealii::MappingQ1<dim> mapping;
  for (const auto &candidate : candidates)
    {
      const auto &cell = cells[candidate.second];
      try
        {
          const auto unit_point =
            mapping.transform_real_to_unit_cell(cell, point);
          if (dealii::GeometryInfo<dim>::distance_to_unit_cell(unit_point) <
              1e-12)
            return results.emplace(key, cell).first->second;
        }
      catch (
        const typename dealii::MappingQGeneric<dim>::ExcTransformationFailed &)
        {}
    }
  return results.emplace(key, dof_handler->end()).first->second;
}

template class SharpIBStencilSearch<2>;
template class SharpIBStencilSearch<3>;
