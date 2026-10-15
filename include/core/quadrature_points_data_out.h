// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_quadrature_points_data_out_h
#define lethe_quadrature_points_data_out_h

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>

#include <string>
#include <vector>

using namespace dealii;


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
   * @brief Build point patches located at quadrature points for visualization.
   *
   * @param[in] location Position of the quadrature point in the real cell.
   * @param[in] filter_radius Radius of the QCM filter centered at the
   * quadrature point.
   * @param[in] cell_index Global active index of the cell that owns the
   * quadrature point.
   */
  void
  build_patches(const dealii::Point<dim>              &location,
                const double                           filter_radius,
                const dealii::types::global_cell_index cell_index);

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
