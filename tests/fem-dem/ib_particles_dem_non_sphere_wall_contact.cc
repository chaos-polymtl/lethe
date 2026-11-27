// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the contact of a non-spherical particle with a wall in the sharp
 * immersed boundary DEM.
 *
 * An ellipsoid with equal radii is a sphere, but it uses the generic detection
 * and search of the contact. When it overlaps a wall, the force must be the
 * Hertz force along the normal of the wall, for a boundary of the domain and
 * for a floating wall, also when the point which defines the floating wall is
 * far from the particle. There is no force when it does not touch the wall.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/ib_particle.h>
#include <core/parameters.h>
#include <core/parameters_lagrangian.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/fem-dem/ib_particles_dem_test_utilities.h>

#include <../tests/tests.h>

#include <cmath>
#include <memory>
#include <vector>

using namespace dealii;

constexpr double dt_dem = 1e-5;

// The DEM uses 0.66665 in place of 2/3 in the normal spring constant.
constexpr double hertz_tolerance = 1e-4;

template <int dim>
void
test()
{
  auto parameters =
    make_ib_particles_parameters<dim>((dim == 2) ? "0; 0" : "0; 0; 0");
  // The effective radius of the ellipsoid is then used in the Hertz model
  parameters->approximate_radius_for_contact = true;

  const auto hertz_force = [&](const IBParticle<dim> &particle,
                               const double           overlap) {
    const double effective_youngs_modulus =
      1. /
      ((1. - particle.poisson_ratio * particle.poisson_ratio) /
         particle.youngs_modulus +
       (1. - parameters->wall_poisson_ratio * parameters->wall_poisson_ratio) /
         parameters->wall_youngs_modulus);
    // The DEM gives the wall 100 times the radius of the particle
    const double effective_radius = particle.radius * 100. / 101.;
    return 4. / 3. * effective_youngs_modulus * std::sqrt(effective_radius) *
           std::pow(overlap, 1.5);
  };

  deallog << "Dimension " << dim << std::endl;

  // Boundary x = 0 of the unit box
  {
    constexpr double radius  = 0.1;
    constexpr double overlap = 1e-3;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, 0., 1., true);
    triangulation.refine_global(3);
    const FESystem<dim> fe(FE_Q<dim>(1), dim + 1);
    DoFHandler<dim>     dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    const auto compute_force = [&](const double signed_overlap) {
      Point<dim> position;
      for (int d = 0; d < dim; ++d)
        position[d] = 0.5;
      position[0] = radius - signed_overlap;
      std::vector<IBParticle<dim>> particles{
        make_ib_ellipsoid<dim>(0, radius, position)};

      IBParticlesDEM<dim> ib_dem;
      ib_dem.initialize(
        parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
      ib_dem.update_particles(particles, 0.);
      ib_dem.update_particles_boundary_contact(particles,
                                               dof_handler,
                                               QGauss<dim - 1>(2),
                                               MappingQ1<dim>());
      std::vector<Tensor<1, 3>> force(particles.size());
      std::vector<Tensor<1, 3>> torque(particles.size());
      ib_dem.calculate_pw_contact_force(dt_dem, force, torque);
      return force[0];
    };

    Tensor<1, 3> expected;
    expected[0] =
      hertz_force(make_ib_ellipsoid<dim>(0, radius, Point<dim>()), overlap);

    log_check("Force of a boundary is the Hertz force",
              (compute_force(overlap) - expected).norm() <
                hertz_tolerance * expected.norm());
    log_check("No force when the boundary is not touched",
              compute_force(-overlap).norm() <
                hertz_tolerance * expected.norm());
  }

  // Oblique floating wall through the origin
  {
    constexpr double radius  = 1.;
    constexpr double overlap = 1e-3;

    const Tensor<1, dim> normal = get_pair_direction<dim>();

    auto floating_walls =
      std::make_shared<Parameters::Lagrangian::FloatingWalls<dim>>();
    floating_walls->points_on_walls.push_back(Point<dim>());
    floating_walls->floating_walls_normal_vectors.push_back(normal);
    floating_walls->time_start.push_back(0.);
    floating_walls->time_end.push_back(2.);

    const auto compute_force = [&](const double signed_overlap,
                                   const double tangential_distance) {
      const Point<dim> position = Point<dim>() +
                                  (radius - signed_overlap) * normal +
                                  tangential_distance * get_pair_tangent<dim>();
      const std::vector<IBParticle<dim>> particles{
        make_ib_ellipsoid<dim>(0, radius, position)};

      IBParticlesDEM<dim> ib_dem;
      ib_dem.initialize(
        parameters, floating_walls, nullptr, MPI_COMM_WORLD, particles);
      // The wall is active at this time
      ib_dem.update_particles(particles, 1.);
      std::vector<Tensor<1, 3>> force(particles.size());
      std::vector<Tensor<1, 3>> torque(particles.size());
      ib_dem.calculate_pw_contact_force(dt_dem, force, torque);
      return force[0];
    };

    const Tensor<1, 3> expected =
      hertz_force(make_ib_ellipsoid<dim>(0, radius, Point<dim>()), overlap) *
      tensor_nd_to_3d(normal);

    log_check("Force of a floating wall is the Hertz force",
              (compute_force(overlap, 0.) - expected).norm() <
                hertz_tolerance * expected.norm());
    log_check("Floating wall is found far from the point which defines it",
              (compute_force(overlap, 50. * radius) - expected).norm() <
                hertz_tolerance * expected.norm());
    log_check("No force when the floating wall is not touched",
              compute_force(-overlap, 50. * radius).norm() <
                hertz_tolerance * expected.norm());
  }
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();
      test<2>();
      test<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
