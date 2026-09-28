// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_immersed_solid_classifier_h
#define lethe_immersed_solid_classifier_h

#include <core/shape.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>

#include <deal.II/numerics/rtree.h>

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

/**
 * @brief Phase of a region of the domain with respect to the immersed solids.
 */
enum class CellPhase : std::uint8_t
{
  /// The region lies entirely in the fluid.
  fluid,
  /// The region lies entirely inside one immersed solid.
  solid,
  /// The region may be intersected by the surface of an immersed solid.
  cut
};

/**
 * @brief Interface classifying regions of the domain as fluid, solid or cut
 * with respect to immersed solids, and evaluating the fluid indicator
 * function (1 in the fluid, 0 in the solids).
 *
 * The classification is conservative: a region is only classified as fluid
 * (or solid) if it is guaranteed to lie entirely in the fluid (or in a solid).
 * The fluid indicator of a cut region must then be evaluated point by point
 * with fill_fluid_indicator(), which only tests the solids that may cut the
 * region. This interface is independent of the solver that owns the solids so
 * that the consumers of the fluid indicator (e.g. filters) are not coupled to
 * the internal representation of the immersed boundaries of a solver.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class ImmersedSolidClassifier
{
public:
  /**
   * @brief Default destructor.
   */
  virtual ~ImmersedSolidClassifier() = default;

  /**
   * @brief Classify a ball with respect to the immersed solids.
   *
   * @param[in] center Center of the ball.
   *
   * @param[in] radius Radius of the ball.
   *
   * @param[out] solids Indices of the solids relevant to the ball. It is empty
   * if the ball is in the fluid, it contains the index of the solid containing
   * the ball if the ball is in a solid, and it contains the indices of the
   * solids that may cut the ball if the ball is cut. The vector is cleared
   * first; its capacity is reused.
   *
   * @return Phase of the ball.
   */
  virtual CellPhase
  classify_ball(const dealii::Point<dim>  &center,
                const double               radius,
                std::vector<unsigned int> &solids) const = 0;

  /**
   * @brief Evaluate the fluid indicator at points of a cut region.
   *
   * @param[in] cutting_solids Indices of the solids that may cut the region,
   * as returned by classify_ball() or classify_cell().
   *
   * @param[in] points Points at which the fluid indicator is evaluated. They
   * must lie in the region that was classified.
   *
   * @param[out] fluid_indicator Fluid indicator at the points: 1 if the point
   * lies outside of all the cutting solids, 0 otherwise. It is resized to the
   * number of points.
   */
  virtual void
  fill_fluid_indicator(const std::vector<unsigned int>       &cutting_solids,
                       const std::vector<dealii::Point<dim>> &points,
                       std::vector<double> &fluid_indicator) const = 0;

  /**
   * @brief Classify a cell with respect to the immersed solids. The cell is
   * enclosed in the ball circumscribing the bounding box of the cell given by
   * the mapping, slightly enlarged to account for curved cells, and this ball
   * is classified with classify_ball().
   *
   * @param[in] cell Cell to classify.
   *
   * @param[in] mapping Mapping describing the geometry of the cell.
   *
   * @param[out] solids Indices of the solids relevant to the cell. See
   * classify_ball().
   *
   * @return Phase of the cell.
   */
  CellPhase
  classify_cell(const typename dealii::Triangulation<dim>::cell_iterator &cell,
                const dealii::Mapping<dim> &mapping,
                std::vector<unsigned int>  &solids) const;
};

/**
 * @brief Classification of the domain with respect to immersed solids
 * described by Shape objects, e.g. the particles of the sharp-edge immersed
 * boundary solver.
 *
 * The value of a shape is its signed distance function, negative inside the
 * solid. The classification of a ball of center \f$\mathbf{c}\f$ and radius
 * \f$r\f$ assumes that the signed distance function \f$\phi\f$ of each shape
 * is 1-Lipschitz, which holds for exact signed distances:
 * - if \f$\phi(\mathbf{c}) < -r\f$ for one solid, the ball lies in that solid;
 * - if \f$\phi(\mathbf{c}) > r\f$ for all solids, the ball lies in the fluid;
 * - otherwise, the ball is cut by the solids for which
 *   \f$|\phi(\mathbf{c})| \le r\f$.
 *
 * A point is in the fluid if the signed distance of every solid is strictly
 * positive, consistently with the sharp-edge immersed boundary solver, which
 * considers that a point with a non-positive signed distance is inside a
 * particle.
 *
 * The candidate solids of a ball are found with an R-tree of the bounding
 * boxes of the solids, so that the cost of a classification does not scale
 * with the number of solids. A bounding box is currently only built for
 * spheres. The other shapes are always considered as candidates, which is
 * correct but slower.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class ShapesImmersedSolidClassifier : public ImmersedSolidClassifier<dim>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] solid_shapes Shapes of the immersed solids, in their current
   * position and orientation. Shapes that are not defined by a signed distance
   * function (e.g. IGES CAD files) are not supported.
   */
  explicit ShapesImmersedSolidClassifier(
    const std::vector<std::shared_ptr<Shape<dim>>> &solid_shapes);

  /**
   * @copydoc ImmersedSolidClassifier::classify_ball
   */
  CellPhase
  classify_ball(const dealii::Point<dim>  &center,
                const double               radius,
                std::vector<unsigned int> &solids) const override;

  /**
   * @copydoc ImmersedSolidClassifier::fill_fluid_indicator
   */
  void
  fill_fluid_indicator(const std::vector<unsigned int>       &cutting_solids,
                       const std::vector<dealii::Point<dim>> &points,
                       std::vector<double> &fluid_indicator) const override;

  /**
   * @brief Return the number of immersed solids.
   *
   * @return Number of solids.
   */
  unsigned int
  n_solids() const
  {
    return shapes.size();
  }

  /**
   * @brief Return the number of immersed solids without bounding box, which
   * are candidates for the classification of every region.
   *
   * @return Number of solids without bounding box.
   */
  unsigned int
  n_unbounded_solids() const
  {
    return unbounded_solids.size();
  }

private:
  /// Shapes of the immersed solids.
  std::vector<std::shared_ptr<Shape<dim>>> shapes;

  /// Indices of the solids without bounding box.
  std::vector<unsigned int> unbounded_solids;

  /// R-tree of the bounding boxes of the solids that have one, paired with the
  /// index of the solid.
  dealii::RTree<std::pair<dealii::BoundingBox<dim>, unsigned int>>
    bounded_solids_tree;
};

#endif
