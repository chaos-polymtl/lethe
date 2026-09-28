// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the classification of the cells of a mesh with respect to
 * immersed spheres and the evaluation of the fluid indicator in the cut cells
 * by ShapesImmersedSolidClassifier.
 *
 * Two configurations are tested in 2D and 3D: a single sphere, and two spheres
 * touching at a point. For each configuration, the test checks that:
 * - the cells classified as fluid (or solid) have all their vertices outside
 *   (or inside) of the spheres;
 * - the classification is identical when the spheres are handled through the
 *   search path used for shapes without bounding box.
 *
 * The test then integrates the volume of a sphere, exactly in the solid cells
 * and with an iterated Gauss quadrature in the cut cells, for an increasing
 * number of subdivisions of the quadrature. Since the quadrature error of the
 * discontinuous fluid indicator oscillates with the position of the sphere
 * relative to the quadrature points, the error is averaged over several
 * positions of the sphere to exhibit the convergence rate.
 */

// Deal.II includes
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/immersed_solid_classifier.h>
#include <core/shape.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

namespace
{
  /**
   * @brief Sphere whose signed distance is evaluated directly. Since it does
   * not derive from Sphere, ShapesImmersedSolidClassifier treats it as a shape
   * without bounding box.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class UnboundedSphere : public Shape<dim>
  {
  public:
    /**
     * @brief Constructor.
     *
     * @param[in] sphere_radius Radius of the sphere.
     *
     * @param[in] sphere_center Center of the sphere.
     */
    UnboundedSphere(const double sphere_radius, const Point<dim> &sphere_center)
      : Shape<dim>(sphere_radius, sphere_center, Tensor<1, 3>())
      , sphere_radius(sphere_radius)
      , sphere_center(sphere_center)
    {}

    /**
     * @brief Return the signed distance to the sphere.
     *
     * @param[in] evaluation_point Point at which the distance is evaluated.
     *
     * @return Signed distance, negative inside the sphere.
     */
    double
    value(const Point<dim> &evaluation_point,
          const unsigned int /*component*/ = 0) const override
    {
      return evaluation_point.distance(sphere_center) - sphere_radius;
    }

    /**
     * @brief Return a copy of the shape.
     *
     * @return Copy of the shape.
     */
    std::shared_ptr<Shape<dim>>
    static_copy() const override
    {
      return std::make_shared<UnboundedSphere<dim>>(sphere_radius,
                                                    sphere_center);
    }

  private:
    /// Radius of the sphere.
    const double sphere_radius;

    /// Center of the sphere.
    const Point<dim> sphere_center;
  };

  /**
   * @brief Format a boolean as true or false.
   *
   * @param[in] value Boolean to format.
   *
   * @return Formatted boolean.
   */
  std::string
  format(const bool value)
  {
    return value ? "true" : "false";
  }

  /**
   * @brief Format a number in scientific notation with three digits.
   *
   * @param[in] value Number to format.
   *
   * @return Formatted number.
   */
  std::string
  format(const double value)
  {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(3) << value;
    return stream.str();
  }

  /**
   * @brief Format a convergence rate in fixed notation with two digits.
   *
   * @param[in] rate Convergence rate to format.
   *
   * @return Formatted rate.
   */
  std::string
  format_rate(const double rate)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << rate;
    return stream.str();
  }

  /**
   * @brief Integrate the volume of the immersed solids: the solid cells are
   * integrated exactly and the cut cells with an iterated Gauss quadrature.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] triangulation Triangulation of the domain.
   *
   * @param[in] mapping Mapping of the domain.
   *
   * @param[in] classifier Classifier of the immersed solids.
   *
   * @param[in] n_subdivisions Number of subdivisions per direction of the
   * quadrature of the cut cells.
   *
   * @return Volume of the immersed solids.
   */
  template <int dim>
  double
  integrate_solid_volume(const Triangulation<dim>           &triangulation,
                         const Mapping<dim>                 &mapping,
                         const ImmersedSolidClassifier<dim> &classifier,
                         const unsigned int                  n_subdivisions)
  {
    const FE_Q<dim>      fe(1);
    const QIterated<dim> cut_cell_quadrature(QGauss<1>(2), n_subdivisions);
    FEValues<dim>        fe_values(mapping,
                            fe,
                            cut_cell_quadrature,
                            update_quadrature_points | update_JxW_values);

    std::vector<unsigned int> solids;
    std::vector<double>       fluid_indicator;
    double                    solid_volume = 0.;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        const CellPhase phase = classifier.classify_cell(cell, mapping, solids);
        if (phase == CellPhase::solid)
          solid_volume += cell->measure();
        else if (phase == CellPhase::cut)
          {
            fe_values.reinit(cell);
            classifier.fill_fluid_indicator(solids,
                                            fe_values.get_quadrature_points(),
                                            fluid_indicator);
            for (const unsigned int q : fe_values.quadrature_point_indices())
              solid_volume += (1. - fluid_indicator[q]) * fe_values.JxW(q);
          }
      }
    return solid_volume;
  }

  /**
   * @brief Return the volume (area in 2D) of a sphere.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] radius Radius of the sphere.
   *
   * @return Volume of the sphere.
   */
  template <int dim>
  double
  sphere_volume(const double radius)
  {
    return (dim == 2) ? std::numbers::pi * radius * radius :
                        4. / 3. * std::numbers::pi * radius * radius * radius;
  }

  /**
   * @brief Classify the cells of a uniform mesh of the unit hypercube with
   * respect to a set of spheres and verify the classification.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] label Description of the configuration printed in the log.
   *
   * @param[in] spheres Radius and center of each sphere.
   */
  template <int dim>
  void
  test_classification(const std::string                                &label,
                      const std::vector<std::pair<double, Point<dim>>> &spheres)
  {
    deallog << label << ", dim = " << dim << std::endl;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, 0., 1.);
    triangulation.refine_global(4);
    const MappingQ1<dim> mapping;

    std::vector<std::shared_ptr<Shape<dim>>> shapes;
    std::vector<std::shared_ptr<Shape<dim>>> unbounded_shapes;
    for (const auto &[radius, center] : spheres)
      {
        shapes.push_back(
          std::make_shared<Sphere<dim>>(radius, center, Tensor<1, 3>()));
        unbounded_shapes.push_back(
          std::make_shared<UnboundedSphere<dim>>(radius, center));
      }

    const ShapesImmersedSolidClassifier<dim> classifier(shapes);
    const ShapesImmersedSolidClassifier<dim> unbounded_classifier(
      unbounded_shapes);

    deallog << "  unbounded solids of the spheres        : "
            << classifier.n_unbounded_solids() << std::endl;
    deallog << "  unbounded solids of the unbounded ones : "
            << unbounded_classifier.n_unbounded_solids() << std::endl;

    // Signed distance to the union of the spheres, used to verify the
    // classification at the vertices of the cells.
    const auto union_signed_distance = [&shapes](const Point<dim> &point) {
      double signed_distance = std::numeric_limits<double>::max();
      for (const auto &shape : shapes)
        signed_distance = std::min(signed_distance, shape->value(point));
      return signed_distance;
    };

    unsigned int              n_fluid_cells = 0;
    unsigned int              n_solid_cells = 0;
    unsigned int              n_cut_cells   = 0;
    bool                      consistent    = true;
    bool                      identical     = true;
    std::vector<unsigned int> solids;
    std::vector<unsigned int> unbounded_solids;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        const CellPhase phase = classifier.classify_cell(cell, mapping, solids);
        const CellPhase unbounded_phase =
          unbounded_classifier.classify_cell(cell, mapping, unbounded_solids);
        identical = identical && (phase == unbounded_phase) &&
                    (solids == unbounded_solids);

        for (const unsigned int v : cell->vertex_indices())
          {
            const double signed_distance =
              union_signed_distance(cell->vertex(v));
            if (phase == CellPhase::fluid && signed_distance <= 0.)
              consistent = false;
            if (phase == CellPhase::solid && signed_distance >= 0.)
              consistent = false;
          }

        if (phase == CellPhase::fluid)
          ++n_fluid_cells;
        else if (phase == CellPhase::solid)
          ++n_solid_cells;
        else
          ++n_cut_cells;
      }

    deallog << "  fluid cells                            : " << n_fluid_cells
            << std::endl;
    deallog << "  solid cells                            : " << n_solid_cells
            << std::endl;
    deallog << "  cut cells                              : " << n_cut_cells
            << std::endl;
    deallog << "  vertices consistent with the phases    : "
            << format(consistent) << std::endl;
    deallog << "  identical without bounding boxes       : "
            << format(identical) << std::endl;
  }

  /**
   * @brief Integrate the volume of a sphere for an increasing number of
   * subdivisions of the quadrature of the cut cells, and report the relative
   * error averaged over several positions of the sphere as well as the
   * observed convergence rate.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  void
  test_volume_convergence()
  {
    deallog << "Convergence of the solid volume, dim = " << dim << std::endl;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, 0., 1.);
    triangulation.refine_global(4);
    const MappingQ1<dim> mapping;

    constexpr double       radius      = 0.3;
    constexpr double       cell_size   = 1. / 16.;
    constexpr unsigned int n_positions = 8;
    // The positions of the sphere are spread within one cell around the
    // center of the domain with a golden-ratio sequence.
    const double golden_ratio_conjugate = (std::sqrt(5.) - 1.) / 2.;

    const std::vector<unsigned int> subdivisions = {1, 2, 4, 8};
    std::vector<double>             mean_errors(subdivisions.size(), 0.);
    for (unsigned int k = 0; k < n_positions; ++k)
      {
        Point<dim> center;
        for (unsigned int d = 0; d < dim; ++d)
          center[d] =
            0.5 +
            (std::fmod((k + 1) * golden_ratio_conjugate * (d + 1), 1.) - 0.5) *
              cell_size;

        const std::vector<std::shared_ptr<Shape<dim>>> shapes = {
          std::make_shared<Sphere<dim>>(radius, center, Tensor<1, 3>())};
        const ShapesImmersedSolidClassifier<dim> classifier(shapes);

        for (unsigned int i = 0; i < subdivisions.size(); ++i)
          {
            const double volume = integrate_solid_volume<dim>(triangulation,
                                                              mapping,
                                                              classifier,
                                                              subdivisions[i]);
            mean_errors[i] += std::abs(volume - sphere_volume<dim>(radius)) /
                              sphere_volume<dim>(radius) / n_positions;
          }
      }

    for (unsigned int i = 0; i < subdivisions.size(); ++i)
      {
        deallog << "  subdivisions " << subdivisions[i]
                << ", mean relative error : " << format(mean_errors[i]);
        if (i > 0)
          deallog << ", rate : "
                  << format_rate(
                       std::log2(mean_errors[i - 1] / mean_errors[i]));
        deallog << std::endl;
      }
  }
} // namespace

int
main()
{
  try
    {
      initlog();

      test_classification<2>("Single sphere", {{0.3, Point<2>(0.5, 0.5)}});
      test_classification<2>("Two touching spheres",
                             {{0.2, Point<2>(0.3, 0.5)},
                              {0.2, Point<2>(0.7, 0.5)}});
      test_volume_convergence<2>();
      test_classification<3>("Single sphere", {{0.3, Point<3>(0.5, 0.5, 0.5)}});
      test_classification<3>("Two touching spheres",
                             {{0.2, Point<3>(0.3, 0.5, 0.5)},
                              {0.2, Point<3>(0.7, 0.5, 0.5)}});
      test_volume_convergence<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
