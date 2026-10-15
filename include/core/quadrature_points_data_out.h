// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_quadrature_points_data_out_h
#define lethe_quadrature_points_data_out_h

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <functional>
#include <string>
#include <vector>


/**
 * @brief Build point patches located at quadrature points for visualization,
 * for example to display the QCM filter spheres centered at each quadrature
 * point. Each quadrature point stores its filter radius, filter diameter and
 * the global active index of the cell it belongs to.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class QuadraturePointsDataOut : public dealii::DataOutInterface<0, dim>
{
public:
  /**
   * @brief Build one point patch per quadrature point of every locally owned
   * cell. Each patch stores the filter radius, the filter diameter and the
   * global active index of the cell that owns the quadrature point.
   *
   * @param[in] dof_handler DoFHandler whose active cells are traversed.
   * @param[in] mapping Mapping used to compute the quadrature point locations.
   * @param[in] quadrature Quadrature rule defining the points in each cell.
   * @param[in] filter_radius Function returning the filter radius of a cell.
   */
  void
  build_patches(
    const dealii::DoFHandler<dim> &dof_handler,
    const dealii::Mapping<dim>    &mapping,
    const dealii::Quadrature<dim> &quadrature,
    const std::function<
      double(const typename dealii::DoFHandler<dim>::active_cell_iterator &)>
      &filter_radius);

private:
  /**
   * @brief Implementation of the corresponding function of the base class.
   */
  const std::vector<dealii::DataOutBase::Patch<0, dim>> &
  get_patches() const override
  {
    return patches;
  }

  /**
   * @brief Implementation of the corresponding function of the base class.
   */
  std::vector<std::string>
  get_dataset_names() const override
  {
    return dataset_names;
  }

  /// Point patches filled by build_patches() and written by the base class.
  std::vector<dealii::DataOutBase::Patch<0, dim>> patches;

  /// Names of the data fields stored in each patch, in row order.
  const std::vector<std::string> dataset_names = {"filter_radius",
                                                  "filter_diameter",
                                                  "cell_index"};
};

#endif
