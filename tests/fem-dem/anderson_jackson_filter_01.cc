// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Anderson-Jackson filter on a uniform velocity-pressure
 * field in a domain that is periodic in every direction and contains no
 * immersed solid.
 *
 * On a uniform periodic mesh of Q1 elements, every filter center sees the
 * same arrangement of source points, so the kernel mass must be identical at
 * every filter center. A contribution missed across the boundaries of the
 * subdomains of the processes or across the periodic boundaries breaks this
 * uniformity. The kernel mass differs from one by the quadrature error of the
 * truncated kernel. Without solids, the fluid volume fraction equals the
 * kernel mass and the phase-averaged velocity and pressure recover the
 * uniform field exactly. A case with Q2 elements checks the filter with
 * higher-order elements, whose filter centers are not all equivalent, and a
 * case with the top-hat kernel checks the second kernel.
 *
 * The output is identical for any number of processes.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>

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
   * @brief Filter a uniform field in a periodic domain and print the
   * statistics of the filtered fields.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] label Description of the case printed in the log.
   *
   * @param[in] velocity_degree Degree of the velocity and pressure elements.
   * The filtered fields use the same degree.
   *
   * @param[in] n_refinements Number of global refinements of the unit square
   * (cube).
   *
   * @param[in] kernel_type Kernel of the filter.
   *
   * @param[in] filter_width Standard deviation of the gaussian kernel or
   * radius of the top-hat kernel.
   *
   * @param[in] gaussian_cutoff Truncation radius of the gaussian kernel in
   * number of standard deviations.
   */
  template <int dim>
  void
  test(const std::string                 &label,
       const unsigned int                 velocity_degree,
       const unsigned int                 n_refinements,
       const Parameters::FilterKernelType kernel_type,
       const double                       filter_width,
       const double                       gaussian_cutoff)
  {
    deallog << label << ", dim = " << dim
            << ", velocity degree = " << velocity_degree << std::endl;

    const UniformFluidProblem<dim> problem(n_refinements,
                                           velocity_degree,
                                           true);

    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.kernel_type     = kernel_type;
    filter_parameters.filter_width    = filter_width;
    filter_parameters.gaussian_cutoff = gaussian_cutoff;

    const std::vector<std::shared_ptr<Shape<dim>>> no_solids;
    const ShapesImmersedSolidClassifier<dim>       classifier(no_solids);

    AndersonJacksonFilter<dim> filter(filter_parameters,
                                      Parameters::Timer(),
                                      MPI_COMM_WORLD);
    filter.apply(problem.dof_handler,
                 problem.mapping,
                 problem.solution,
                 classifier,
                 problem.periodic_boundaries);

    const auto &statistics = filter.get_statistics();
    const auto &fields     = filter.get_filtered_fields();

    // Without solids, the fluid volume fraction equals the kernel mass and the
    // solid volume fraction is zero.
    double fraction_difference = 0.;
    double solid_fraction      = 0.;
    for (const types::global_dof_index dof :
         filter.get_dof_handler().locally_owned_dofs())
      {
        fraction_difference =
          std::max(fraction_difference,
                   std::abs(fields.fluid_volume_fraction(dof) -
                            fields.kernel_mass(dof)));
        solid_fraction =
          std::max(solid_fraction, std::abs(fields.solid_volume_fraction(dof)));
      }
    fraction_difference =
      Utilities::MPI::max(fraction_difference, MPI_COMM_WORLD);
    solid_fraction = Utilities::MPI::max(solid_fraction, MPI_COMM_WORLD);

    const auto [velocity_error, pressure_error] =
      maximum_phase_average_errors(filter);

    deallog << "  number of filter centers             : "
            << statistics.n_filter_centers << std::endl;
    deallog << "  number of fluid, solid and cut cells : "
            << statistics.n_fluid_cells << " " << statistics.n_solid_cells
            << " " << statistics.n_cut_cells << std::endl;
    deallog << "  minimum kernel mass                  : "
            << format(statistics.minimum_kernel_mass) << std::endl;
    deallog << "  maximum kernel mass                  : "
            << format(statistics.maximum_kernel_mass) << std::endl;
    // The filter centers of Q1 elements are all equivalent.
    if (velocity_degree == 1)
      deallog << "  kernel mass is uniform               : "
              << format(statistics.maximum_kernel_mass -
                          statistics.minimum_kernel_mass <
                        1e-12)
              << std::endl;
    deallog << "  fluid fraction equals kernel mass    : "
            << format(fraction_difference < 1e-14) << std::endl;
    deallog << "  solid fraction is zero               : "
            << format(solid_fraction < 1e-14) << std::endl;
    deallog << "  velocity is recovered                : "
            << format(velocity_error < 1e-12) << std::endl;
    deallog << "  pressure is recovered                : "
            << format(pressure_error < 1e-12) << std::endl;
    deallog << "  undefined phase averages             : "
            << statistics.n_undefined_averages << std::endl;
    deallog << "  volume of the source cells is one    : "
            << format(std::abs(statistics.source_volume - 1.) < 1e-12)
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

      test<2>("Gaussian kernel",
              1,
              4,
              Parameters::FilterKernelType::gaussian,
              0.1,
              3.);
      test<2>("Gaussian kernel",
              2,
              3,
              Parameters::FilterKernelType::gaussian,
              0.1,
              3.);
      test<3>("Gaussian kernel",
              1,
              3,
              Parameters::FilterKernelType::gaussian,
              0.1,
              3.5);
      test<3>(
        "Top-hat kernel", 1, 3, Parameters::FilterKernelType::top_hat, 0.3, 3.);
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
