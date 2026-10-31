// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_sharp_ib_stencil_search_h
#define lethe_sharp_ib_stencil_search_h

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/rtree.h>

#include <boost/signals2/connection.hpp>

#include <array>
#include <map>
#include <set>
#include <utility>
#include <vector>

/**
 * @brief Cache geometric stencil searches on the locally accessible Q1 mesh.
 *
 * Owned and ghost cells are searchable; artificial cells are excluded. Exact
 * query coordinates and the starting cell identify reusable results, including
 * misses. No particle, interpolation, or matrix data are cached. Mesh signals
 * invalidate the index; callers must also reinitialize after DoF
 * redistribution.
 *
 * @tparam dim Spatial dimension.
 */
template <int dim>
class SharpIBStencilSearch
{
public:
  using Cell = typename dealii::DoFHandler<dim>::active_cell_iterator;

  /**
   * @brief Clear previous results and attach to the current DoFHandler.
   * @param dof_handler Mesh and DoFs to search. Must outlive this object.
   */
  void
  reinit(const dealii::DoFHandler<dim> &dof_handler);

  /**
   * @brief Find a containing cell, or return DoFHandler::end() on a miss.
   *
   * Nearby cells retain the existing search preference. The fallback uses a
   * bounding-box index followed by exact inverse-map containment checks.
   *
   * @param point Query in the fundamental periodic domain.
   * @param seed Starting cell used for the nearby search and cache key.
   * @param vertices_to_cell Locally accessible cells incident on each vertex.
   */
  Cell
  find(const dealii::Point<dim>               &point,
       const Cell                             &seed,
       std::map<unsigned int, std::set<Cell>> &vertices_to_cell);

private:
  /** @brief Invalidate all mesh-dependent state. */
  void
  clear();

  /** @brief Build the index once, on the first fallback search. */
  void
  build_index();

  using BoxEntry = std::pair<dealii::BoundingBox<dim>, unsigned int>;
  const dealii::DoFHandler<dim> *dof_handler = nullptr;
  dealii::RTree<BoxEntry>        index;
  std::vector<Cell>              cells;
  bool                           index_ready = false;
  std::map<std::pair<Cell, std::array<double, dim>>, Cell> results;
  // Declared last so the callback is disconnected before other members die.
  boost::signals2::scoped_connection mesh_change_connection;
  boost::signals2::scoped_connection mesh_movement_connection;
  boost::signals2::scoped_connection repartition_connection;
};

#endif
