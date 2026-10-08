// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Anderson-Jackson filter of a Couette flow between a
 * stationary wall and a moving wall, with and without the extension of the
 * velocity beyond the walls.
 *
 * The unit square is periodic in x. The velocity is u = (y, 0), which vanishes
 * on the wall at y = 0 and equals the velocity (1, 0) of the wall at y = 1.
 * The filtered velocity is printed at filter centers located at increasing
 * distances from the stationary wall. Without extension, the kernel is
 * truncated by the wall and the filtered velocity is the velocity at the
 * centroid of the truncated kernel, which lies inside the domain. With
 * extension, the part of the kernel beyond the wall is filled with fluid
 * moving at the velocity of the wall, which halves the deviation from the
 * exact profile at the wall. Beyond the support radius of the kernel (0.3),
 * both are exact. The test also checks that the filtered profile is
 * antisymmetric about the middle of the channel, like the exact one, that it
 * is exact at the middle, and that the filtered vertical velocity is zero.
 *
 * The fluid volume fraction is printed along with the velocity. Without
 * extension, it is the mass of the truncated kernel, which is one half on the
 * wall. With extension, the fluid beyond the wall is counted in the fluid
 * volume fraction, which is one everywhere, so that the fluid and solid
 * volume fractions sum to one and the product of the fluid volume fraction
 * and of the filtered velocity remains the flux of fluid seen by the kernel.
 *
 * The output is identical for any number of processes.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/grids.h>
#include <core/immersed_solid_classifier.h>
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>
#include <core/periodic_boundary.h>
#include <core/shape.h>
#include <core/vector.h>

#include <fem-dem/anderson_jackson_filter.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

// Tests (with common definitions)
#include <../tests/fem-dem/anderson_jackson_filter_test_utilities.h>

#include <../tests/tests.h>

namespace
{
  /**
   * @brief Couette velocity u = (y, 0) and zero pressure.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  class CouetteFlow : public Function<dim>
  {
  public:
    /**
     * @brief Constructor.
     */
    CouetteFlow()
      : Function<dim>(dim + 1)
    {}

    /**
     * @brief Return a component of the velocity-pressure field.
     *
     * @param[in] point Evaluation point.
     *
     * @param[in] component Component of the field.
     *
     * @return Value of the component.
     */
    double
    value(const Point<dim> &point, const unsigned int component) const override
    {
      return (component == 0) ? point[1] : 0.;
    }
  };

  /**
   * @brief Filter the Couette flow and print the filtered velocity near the
   * stationary wall.
   *
   * @param[in] extend_velocity_beyond_walls Fill the part of the kernel beyond
   * the walls with fluid moving at the velocity of the walls.
   */
  void
  test(const bool extend_velocity_beyond_walls)
  {
    constexpr int          dim           = 2;
    constexpr unsigned int n_refinements = 4;
    constexpr unsigned int n_cells       = 16;

    deallog << "Couette flow, velocity beyond the walls: "
            << (extend_velocity_beyond_walls ? "velocity of the wall" :
                                               "kernel truncated")
            << std::endl;

    // Unit square, periodic in x, with walls at y = 0 (boundary 2) and y = 1
    // (boundary 3).
    parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(triangulation, 0., 1., true);
    Parameters::PeriodicBoundaries periodic_boundaries;
    periodic_boundaries[0] = {1, 0};
    setup_periodic_boundary_conditions(triangulation, periodic_boundaries);
    triangulation.refine_global(n_refinements);

    const FESystem<dim> fe(FE_Q<dim>(1), dim, FE_Q<dim>(1), 1);
    DoFHandler<dim>     dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);
    const MappingQ<dim> mapping(1);

    const IndexSet   locally_owned_dofs = dof_handler.locally_owned_dofs();
    GlobalVectorType owned_solution(locally_owned_dofs, MPI_COMM_WORLD);
    VectorTools::interpolate(mapping,
                             dof_handler,
                             CouetteFlow<dim>(),
                             owned_solution);
    GlobalVectorType solution;
    solution.reinit(locally_owned_dofs,
                    DoFTools::extract_locally_relevant_dofs(dof_handler),
                    MPI_COMM_WORLD);
    solution = owned_solution;

    AndersonJacksonFilter<dim>::WallVelocities wall_velocities;
    wall_velocities[2] = std::make_shared<Functions::ZeroFunction<dim>>(dim);
    wall_velocities[3] = std::make_shared<Functions::ConstantFunction<dim>>(
      std::vector<double>{1., 0.});

    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.kernel_type     = Parameters::FilterKernelType::gaussian;
    filter_parameters.filter_width    = 0.1;
    filter_parameters.gaussian_cutoff = 3.;
    filter_parameters.extend_velocity_beyond_walls =
      extend_velocity_beyond_walls;

    const std::vector<std::shared_ptr<Shape<dim>>> no_solids;
    const ShapesImmersedSolidClassifier<dim>       classifier(no_solids);

    AndersonJacksonFilter<dim> filter(filter_parameters,
                                      Parameters::Timer(),
                                      MPI_COMM_WORLD);
    filter.apply(dof_handler,
                 mapping,
                 solution,
                 classifier,
                 periodic_boundaries,
                 wall_velocities);

    const auto &fields     = filter.get_filtered_fields();
    const auto  filtered_u = [&](const double y) {
      return value_at_filter_center(filter.get_dof_handler(),
                                    mapping,
                                    fields.velocity[0],
                                    Point<dim>(0.5, y));
    };
    const auto fluid_fraction = [&](const double y) {
      return value_at_filter_center(filter.get_dof_handler(),
                                    mapping,
                                    fields.fluid_volume_fraction,
                                    Point<dim>(0.5, y));
    };

    // The exact profile is antisymmetric about the middle of the channel:
    // u(1 - y) = 1 - u(y).
    bool antisymmetric = true;
    for (unsigned int k = 0; k <= 5; ++k)
      {
        const double y = static_cast<double>(k) / n_cells;
        const double u = filtered_u(y);
        deallog << "  y = " << format(y) << ", filtered u = " << format(u)
                << ", exact u = " << format(y)
                << ", fluid fraction = " << format(fluid_fraction(y))
                << std::endl;
        antisymmetric =
          antisymmetric && std::abs(u + filtered_u(1. - y) - 1.) < 1e-12;
      }

    double maximum_v      = 0.;
    double fraction_error = 0.;
    for (const types::global_dof_index dof :
         filter.get_dof_handler().locally_owned_dofs())
      {
        maximum_v = std::max(maximum_v, std::abs(fields.velocity[1](dof)));
        fraction_error =
          std::max(fraction_error,
                   std::abs(fields.fluid_volume_fraction(dof) +
                            fields.solid_volume_fraction(dof) - 1.));
      }
    maximum_v      = Utilities::MPI::max(maximum_v, MPI_COMM_WORLD);
    fraction_error = Utilities::MPI::max(fraction_error, MPI_COMM_WORLD);

    deallog << "  profile antisymmetric about the middle : "
            << format(antisymmetric) << std::endl;
    deallog << "  exact at the middle                    : "
            << format(std::abs(filtered_u(0.5) - 0.5) < 1e-12) << std::endl;
    deallog << "  filtered v is zero                     : "
            << format(maximum_v < 1e-14) << std::endl;
    deallog << "  volume fractions sum to one            : "
            << format(fraction_error < 1e-12) << std::endl;
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      test(false);
      test(true);
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
