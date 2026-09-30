// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Anderson-Jackson filter of a uniform velocity-pressure field
 * around an immersed sphere, in a domain that is periodic in every direction.
 *
 * A sphere of radius a = 0.2 is centered at a vertex in the middle of the unit
 * square (cube). With a top-hat kernel of radius R = 0.35, the solid volume
 * fraction at the center of the sphere is the volume of the sphere divided by
 * the volume of the kernel, a^d / R^d, and at the distance 0.25 from the
 * center it is the volume of the lens-shaped intersection of the sphere and
 * the kernel divided by the volume of the kernel. For an increasing number of
 * subdivisions of the quadrature of the cut cells, the test prints:
 * - the relative error of the volume of the sphere integrated by the source
 *   quadratures;
 * - whether the solid volume fraction at the center times the volume of the
 *   kernel equals this integrated volume, which checks that every solid
 *   source point reached the filter center;
 * - the solid volume fraction at the center and off the center.
 *
 * The quadrature error of the discontinuous fluid indicator does not decrease
 * monotonically for a sphere centered at a vertex (see the test
 * immersed_solid_classifier_01 for averaged convergence rates). The test also
 * compares the solid volume fraction at the center with a gaussian kernel to
 * its closed form, and checks that the phase averages are left undefined, and
 * zero, at the filter centers whose kernel lies entirely in a large sphere.
 * The phase-averaged velocity and pressure must recover the uniform field
 * wherever they are defined.
 *
 * The output is identical for any number of processes.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <core/immersed_solid_classifier.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/shape.h>

#include <fem-dem/anderson_jackson_filter.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

// Tests (with common definitions)
#include <../tests/fem-dem/anderson_jackson_filter_test_utilities.h>

#include <../tests/tests.h>

namespace
{
  /// Radius of the immersed sphere.
  constexpr double sphere_radius = 0.2;

  /// Distance between the center of the sphere and the off-center filter
  /// center.
  constexpr double off_center_distance = 0.25;

  /**
   * @brief Return the volume (area in 2D) of a ball.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] radius Radius of the ball.
   *
   * @return Volume of the ball.
   */
  template <int dim>
  double
  ball_volume(const double radius)
  {
    return (dim == 2) ? std::numbers::pi * radius * radius :
                        4. / 3. * std::numbers::pi * radius * radius * radius;
  }

  /**
   * @brief Return the volume (area in 2D) of the intersection of two balls
   * that partially overlap.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] r1 Radius of the first ball.
   *
   * @param[in] r2 Radius of the second ball.
   *
   * @param[in] d Distance between the centers of the balls.
   *
   * @return Volume of the intersection.
   */
  template <int dim>
  double
  intersection_volume(const double r1, const double r2, const double d)
  {
    if constexpr (dim == 2)
      return r1 * r1 * std::acos((d * d + r1 * r1 - r2 * r2) / (2. * d * r1)) +
             r2 * r2 * std::acos((d * d + r2 * r2 - r1 * r1) / (2. * d * r2)) -
             0.5 * std::sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) *
                             (d + r1 + r2));
    else
      return std::numbers::pi * (r1 + r2 - d) * (r1 + r2 - d) *
             (d * d + 2. * d * r2 - 3. * r2 * r2 + 2. * d * r1 + 6. * r1 * r2 -
              3. * r1 * r1) /
             (12. * d);
  }

  /**
   * @brief Return the mass of a centered gaussian of unit standard deviation
   * contained in the ball of radius X.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] X Radius of the ball in number of standard deviations.
   *
   * @return Mass contained in the ball.
   */
  template <int dim>
  double
  gaussian_mass_in_ball(const double X)
  {
    if constexpr (dim == 2)
      return -std::expm1(-0.5 * X * X);
    else
      return std::erf(X / std::numbers::sqrt2) -
             std::sqrt(2. / std::numbers::pi) * X * std::exp(-0.5 * X * X);
  }

  /**
   * @brief Return the point at the middle of the unit square (cube).
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @return Middle of the domain.
   */
  template <int dim>
  Point<dim>
  middle_of_domain()
  {
    Point<dim> middle;
    for (unsigned int d = 0; d < dim; ++d)
      middle[d] = 0.5;
    return middle;
  }

  /**
   * @brief Apply the filter to the uniform field around a sphere centered at
   * the middle of the domain.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] problem Uniform field in the periodic unit square (cube).
   *
   * @param[in] radius Radius of the sphere.
   *
   * @param[in,out] filter Filter to apply.
   */
  template <int dim>
  void
  filter_around_sphere(const UniformFluidProblem<dim> &problem,
                       const double                    radius,
                       AndersonJacksonFilter<dim>     &filter)
  {
    const std::vector<std::shared_ptr<Shape<dim>>> shapes = {
      std::make_shared<Sphere<dim>>(radius,
                                    middle_of_domain<dim>(),
                                    Tensor<1, 3>())};
    const ShapesImmersedSolidClassifier<dim> classifier(shapes);

    filter.apply(problem.dof_handler,
                 problem.mapping,
                 problem.solution,
                 classifier,
                 problem.periodic_boundaries);
  }

  /**
   * @brief Filter the field around a sphere with a top-hat kernel, for an
   * increasing number of subdivisions of the quadrature of the cut cells.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] n_refinements Number of global refinements of the unit square
   * (cube).
   *
   * @param[in] subdivisions Numbers of subdivisions of the quadrature of the
   * cut cells.
   */
  template <int dim>
  void
  test_top_hat(const unsigned int               n_refinements,
               const std::vector<unsigned int> &subdivisions)
  {
    constexpr double kernel_radius = 0.35;

    deallog << "Sphere and top-hat kernel, dim = " << dim << std::endl;

    const UniformFluidProblem<dim> problem(n_refinements, 1, true);
    const Point<dim>               center     = middle_of_domain<dim>();
    Point<dim>                     off_center = center;
    off_center[0] += off_center_distance;

    const double solid_volume  = ball_volume<dim>(sphere_radius);
    const double kernel_volume = ball_volume<dim>(kernel_radius);
    deallog << "  exact solid fraction at the center      : "
            << format(solid_volume / kernel_volume) << std::endl;
    deallog << "  exact solid fraction off the center     : "
            << format(intersection_volume<dim>(sphere_radius,
                                               kernel_radius,
                                               off_center_distance) /
                      kernel_volume)
            << std::endl;

    for (const unsigned int n_subdivisions : subdivisions)
      {
        Parameters::AndersonJacksonFilter filter_parameters;
        filter_parameters.kernel_type  = Parameters::FilterKernelType::top_hat;
        filter_parameters.filter_width = kernel_radius;
        filter_parameters.cut_cell_subdivisions = n_subdivisions;

        AndersonJacksonFilter<dim> filter(filter_parameters,
                                          Parameters::Timer(),
                                          MPI_COMM_WORLD);
        filter_around_sphere(problem, sphere_radius, filter);

        const auto &statistics = filter.get_statistics();
        const auto &solid_fraction =
          filter.get_filtered_fields().solid_volume_fraction;
        const double center_fraction = value_at_filter_center(
          filter.get_dof_handler(), problem.mapping, solid_fraction, center);
        const double off_center_fraction =
          value_at_filter_center(filter.get_dof_handler(),
                                 problem.mapping,
                                 solid_fraction,
                                 off_center);
        const auto [velocity_error, pressure_error] =
          maximum_phase_average_errors(filter);

        deallog << "  cut cell subdivisions " << n_subdivisions << std::endl;
        deallog << "    fluid, solid and cut cells            : "
                << statistics.n_fluid_cells << " " << statistics.n_solid_cells
                << " " << statistics.n_cut_cells << std::endl;
        deallog << "    relative error of the solid volume    : "
                << format(std::abs(statistics.source_solid_volume -
                                   solid_volume) /
                            solid_volume,
                          3)
                << std::endl;
        deallog << "    center fraction matches solid volume  : "
                << format(std::abs(center_fraction * kernel_volume -
                                   statistics.source_solid_volume) <
                          1e-12 * solid_volume)
                << std::endl;
        deallog << "    solid fraction at the center          : "
                << format(center_fraction) << std::endl;
        deallog << "    solid fraction off the center         : "
                << format(off_center_fraction) << std::endl;
        deallog << "    velocity and pressure are recovered   : "
                << format(velocity_error < 1e-12 && pressure_error < 1e-12)
                << std::endl;
        deallog << "    undefined phase averages              : "
                << statistics.n_undefined_averages << std::endl;
      }
  }

  /**
   * @brief Filter the field around a sphere with a gaussian kernel and compare
   * the solid volume fraction at the center of the sphere with its closed
   * form.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] n_refinements Number of global refinements of the unit square
   * (cube).
   */
  template <int dim>
  void
  test_gaussian(const unsigned int n_refinements)
  {
    constexpr double standard_deviation = 0.1;
    constexpr double cutoff             = 3.5;

    deallog << "Sphere and gaussian kernel, dim = " << dim << std::endl;

    const UniformFluidProblem<dim> problem(n_refinements, 1, true);

    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.kernel_type     = Parameters::FilterKernelType::gaussian;
    filter_parameters.filter_width    = standard_deviation;
    filter_parameters.gaussian_cutoff = cutoff;
    AndersonJacksonFilter<dim> filter(filter_parameters,
                                      Parameters::Timer(),
                                      MPI_COMM_WORLD);
    filter_around_sphere(problem, sphere_radius, filter);

    const double center_fraction =
      value_at_filter_center(filter.get_dof_handler(),
                             problem.mapping,
                             filter.get_filtered_fields().solid_volume_fraction,
                             middle_of_domain<dim>());
    const auto [velocity_error, pressure_error] =
      maximum_phase_average_errors(filter);

    // Mass of the truncated gaussian contained in the sphere.
    deallog << "  exact solid fraction at the center      : "
            << format(gaussian_mass_in_ball<dim>(sphere_radius /
                                                 standard_deviation) /
                      gaussian_mass_in_ball<dim>(cutoff))
            << std::endl;
    deallog << "  solid fraction at the center            : "
            << format(center_fraction) << std::endl;
    deallog << "  velocity and pressure are recovered     : "
            << format(velocity_error < 1e-12 && pressure_error < 1e-12)
            << std::endl;
    deallog << "  undefined phase averages                : "
            << filter.get_statistics().n_undefined_averages << std::endl;
  }

  /**
   * @brief Filter the field around a sphere that is larger than the kernel,
   * so that the kernel of the filter centers near the center of the sphere
   * contains no fluid.
   */
  void
  test_undefined_averages()
  {
    constexpr unsigned int dim                 = 2;
    constexpr double       large_sphere_radius = 0.4;

    deallog << "Large sphere and small top-hat kernel, dim = " << dim
            << std::endl;

    const UniformFluidProblem<dim> problem(4, 1, true);

    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.kernel_type  = Parameters::FilterKernelType::top_hat;
    filter_parameters.filter_width = 0.1;
    AndersonJacksonFilter<dim> filter(filter_parameters,
                                      Parameters::Timer(),
                                      MPI_COMM_WORLD);
    filter_around_sphere(problem, large_sphere_radius, filter);

    const auto &statistics = filter.get_statistics();
    const auto [velocity_error, pressure_error] =
      maximum_phase_average_errors(filter);

    deallog << "  fluid, solid and cut cells              : "
            << statistics.n_fluid_cells << " " << statistics.n_solid_cells
            << " " << statistics.n_cut_cells << std::endl;
    deallog << "  undefined phase averages                : "
            << statistics.n_undefined_averages << std::endl;
    deallog << "  undefined phase averages are zero       : "
            << format(undefined_averages_are_zero(filter)) << std::endl;
    deallog << "  velocity and pressure are recovered     : "
            << format(velocity_error < 1e-12 && pressure_error < 1e-12)
            << std::endl;
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      test_top_hat<2>(4, {1, 2, 4, 8});
      test_top_hat<3>(3, {1, 2, 4});
      test_gaussian<2>(4);
      test_gaussian<3>(3);
      test_undefined_averages();
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
