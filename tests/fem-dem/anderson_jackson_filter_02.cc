// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Anderson-Jackson filter near the walls of a non-periodic
 * domain, for a uniform velocity-pressure field without immersed solid.
 *
 * The kernel is truncated by the walls. On a uniform mesh of Q1 elements, the
 * source points seen by a filter center located on a face, an edge or a
 * corner of the unit square (cube) are the half, the quarter or the eighth of
 * the symmetric set of points seen by a filter center at the middle of the
 * domain. The kernel mass on the faces, the edges and the corners must thus be
 * exactly 1/2, 1/4 and 1/8 of the kernel mass at the middle of the domain,
 * which checks that the kernel is symmetric and that no source point is
 * missed or counted twice. When the volume fractions are renormalized by the
 * kernel mass, the fluid volume fraction must be one everywhere. The
 * phase-averaged velocity and pressure must recover the uniform field.
 *
 * The output is identical for any number of processes.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

// Lethe
#include <core/immersed_solid_classifier.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/shape.h>

#include <fem-dem/anderson_jackson_filter.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// Tests (with common definitions)
#include <../tests/fem-dem/anderson_jackson_filter_test_utilities.h>

#include <../tests/tests.h>

namespace
{
  /**
   * @brief Filter a uniform field in the unit square (cube) with walls and
   * print the kernel mass on the boundaries relative to the kernel mass at the
   * middle of the domain.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] label Description of the case printed in the log.
   *
   * @param[in] n_refinements Number of global refinements of the unit square
   * (cube).
   *
   * @param[in] kernel_type Kernel of the filter.
   *
   * @param[in] filter_width Standard deviation of the gaussian kernel or
   * radius of the top-hat kernel. The support of the kernel centered at the
   * middle of the domain must not reach the walls.
   *
   * @param[in] gaussian_cutoff Truncation radius of the gaussian kernel in
   * number of standard deviations.
   */
  template <int dim>
  void
  test(const std::string                 &label,
       const unsigned int                 n_refinements,
       const Parameters::FilterKernelType kernel_type,
       const double                       filter_width,
       const double                       gaussian_cutoff)
  {
    deallog << label << ", dim = " << dim << std::endl;

    const UniformFluidProblem<dim> problem(n_refinements, 1, false);

    const std::vector<std::shared_ptr<Shape<dim>>> no_solids;
    const ShapesImmersedSolidClassifier<dim>       classifier(no_solids);

    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.kernel_type     = kernel_type;
    filter_parameters.filter_width    = filter_width;
    filter_parameters.gaussian_cutoff = gaussian_cutoff;

    // Kernel truncated by the walls.
    {
      AndersonJacksonFilter<dim> filter(filter_parameters,
                                        Parameters::Timer(),
                                        MPI_COMM_WORLD);
      filter.apply(problem.dof_handler,
                   problem.mapping,
                   problem.solution,
                   classifier,
                   problem.periodic_boundaries);

      const auto &kernel_mass = filter.get_filtered_fields().kernel_mass;
      const auto  mass_at     = [&](const Point<dim> &point) {
        return value_at_filter_center(filter.get_dof_handler(),
                                      problem.mapping,
                                      kernel_mass,
                                      point);
      };

      Point<dim> middle;
      for (unsigned int d = 0; d < dim; ++d)
        middle[d] = 0.5;
      Point<dim> face_center = middle;
      face_center[0]         = 0.;
      Point<dim> edge_center = face_center;
      edge_center[1]         = 0.;
      const Point<dim> corner;

      const double middle_mass = mass_at(middle);
      deallog << "  kernel mass at the middle            : "
              << format(middle_mass) << std::endl;
      deallog << "  face to middle kernel mass ratio     : "
              << format(mass_at(face_center) / middle_mass) << std::endl;
      if (dim == 3)
        deallog << "  edge to middle kernel mass ratio     : "
                << format(mass_at(edge_center) / middle_mass) << std::endl;
      deallog << "  corner to middle kernel mass ratio   : "
              << format(mass_at(corner) / middle_mass) << std::endl;

      const auto [velocity_error, pressure_error] =
        maximum_phase_average_errors(filter);
      deallog << "  velocity is recovered                : "
              << format(velocity_error < 1e-12) << std::endl;
      deallog << "  pressure is recovered                : "
              << format(pressure_error < 1e-12) << std::endl;
      deallog << "  undefined phase averages             : "
              << filter.get_statistics().n_undefined_averages << std::endl;
    }

    // Kernel renormalized by its mass.
    {
      filter_parameters.normalize_at_domain_boundaries = true;
      AndersonJacksonFilter<dim> filter(filter_parameters,
                                        Parameters::Timer(),
                                        MPI_COMM_WORLD);
      filter.apply(problem.dof_handler,
                   problem.mapping,
                   problem.solution,
                   classifier,
                   problem.periodic_boundaries);

      const auto &fields         = filter.get_filtered_fields();
      double      fluid_error    = 0.;
      double      solid_fraction = 0.;
      for (const types::global_dof_index dof :
           filter.get_dof_handler().locally_owned_dofs())
        {
          fluid_error =
            std::max(fluid_error,
                     std::abs(fields.fluid_volume_fraction(dof) - 1.));
          solid_fraction =
            std::max(solid_fraction,
                     std::abs(fields.solid_volume_fraction(dof)));
        }
      fluid_error    = Utilities::MPI::max(fluid_error, MPI_COMM_WORLD);
      solid_fraction = Utilities::MPI::max(solid_fraction, MPI_COMM_WORLD);

      deallog << "  renormalized fluid fraction is one   : "
              << format(fluid_error < 1e-14) << std::endl;
      deallog << "  renormalized solid fraction is zero  : "
              << format(solid_fraction < 1e-14) << std::endl;
    }
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      test<2>(
        "Gaussian kernel", 4, Parameters::FilterKernelType::gaussian, 0.1, 3.);
      test<2>(
        "Top-hat kernel", 4, Parameters::FilterKernelType::top_hat, 0.3, 3.);
      test<3>(
        "Gaussian kernel", 3, Parameters::FilterKernelType::gaussian, 0.1, 3.5);
      test<3>(
        "Top-hat kernel", 3, Parameters::FilterKernelType::top_hat, 0.3, 3.);
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
