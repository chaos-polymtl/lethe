// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test that the lubrication force of the sharp immersed boundary DEM on
 * a particle is the sum of the contributions of all its neighbors.
 *
 * Particle-particle: three spheres are aligned and only the middle one moves.
 * The force on each sphere must be the sum of the forces obtained with the
 * two pairs of neighboring spheres taken separately. The two outer spheres are
 * contact candidates of each other but are outside of the lubrication range.
 *
 * Particle-wall: a sphere moves toward the corner of a box. The force must be
 * the sum of the forces obtained with the sphere close to each of the two
 * walls separately.
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

#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/fem-dem/ib_particles_dem_test_utilities.h>

#include <../tests/tests.h>

#include <memory>
#include <vector>

using namespace dealii;

// Parameters of the lubrication model.
constexpr double h_max  = 0.05;
constexpr double h_min  = 0.001;
constexpr double mu     = 0.01;
constexpr double dt_dem = 1e-5;
constexpr double radius = 0.1;

template <int dim>
std::vector<Tensor<1, 3>>
particle_particle_lubrication(
  const std::shared_ptr<Parameters::IBParticles<dim>> &parameters,
  const std::vector<IBParticle<dim>>                  &particles)
{
  IBParticlesDEM<dim> ib_dem;
  ib_dem.initialize(parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
  ib_dem.update_contact_candidates();

  std::vector<Tensor<1, 3>> force(particles.size());
  std::vector<Tensor<1, 3>> torque(particles.size());
  ib_dem.calculate_pp_lubrication_force(
    dt_dem, h_max, h_min, mu, force, torque);
  return force;
}

template <int dim>
Tensor<1, 3>
particle_wall_lubrication(
  const std::shared_ptr<Parameters::IBParticles<dim>> &parameters,
  const DoFHandler<dim>                               &dof_handler,
  const Point<dim>                                    &position,
  const Tensor<1, 3>                                  &velocity)
{
  std::vector<IBParticle<dim>> particles{
    make_ib_sphere<dim>(0, radius, position)};
  particles[0].velocity = velocity;
  particles[0].initialize_previous_solution();

  IBParticlesDEM<dim> ib_dem;
  ib_dem.initialize(parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
  ib_dem.update_particles_boundary_contact(particles,
                                           dof_handler,
                                           QGauss<dim - 1>(2),
                                           MappingQ1<dim>());

  std::vector<Tensor<1, 3>> force(particles.size());
  std::vector<Tensor<1, 3>> torque(particles.size());
  ib_dem.calculate_pw_lubrication_force(
    dt_dem, h_max, h_min, mu, force, torque);
  return force[0];
}

template <int dim>
void
test()
{
  const double tolerance = 1e-12;
  auto         parameters =
    make_ib_particles_parameters<dim>((dim == 2) ? "0; 0" : "0; 0; 0");

  auto to_string = [](const bool value) { return value ? "true" : "false"; };

  deallog << "Dimension " << dim << std::endl;

  // Particle-particle lubrication. The gaps between the spheres are 0.01 and
  // 0.02, the gap between the two outer spheres is larger than h_max.
  {
    Point<dim> position_a;
    Point<dim> position_b;
    position_b[0] = 2. * radius + 0.01;
    Point<dim> position_c;
    position_c[0] = position_b[0] + 2. * radius + 0.02;

    Tensor<1, 3> velocity_b;
    velocity_b[0] = 1.;
    velocity_b[1] = 0.5;

    auto make_particle = [&](const unsigned int  id,
                             const Point<dim>   &position,
                             const Tensor<1, 3> &velocity) {
      IBParticle<dim> particle = make_ib_sphere<dim>(id, radius, position);
      particle.velocity        = velocity;
      particle.initialize_previous_solution();
      return particle;
    };

    const auto force_abc = particle_particle_lubrication<dim>(
      parameters,
      {make_particle(0, position_a, Tensor<1, 3>()),
       make_particle(1, position_b, velocity_b),
       make_particle(2, position_c, Tensor<1, 3>())});
    const auto force_ab = particle_particle_lubrication<dim>(
      parameters,
      {make_particle(0, position_a, Tensor<1, 3>()),
       make_particle(1, position_b, velocity_b)});
    const auto force_bc = particle_particle_lubrication<dim>(
      parameters,
      {make_particle(0, position_b, velocity_b),
       make_particle(1, position_c, Tensor<1, 3>())});

    const double scale = force_ab[1].norm() + force_bc[0].norm();
    const bool   sum_is_right =
      (force_abc[0] - force_ab[0]).norm() < tolerance * scale &&
      (force_abc[1] - (force_ab[1] + force_bc[0])).norm() < tolerance * scale &&
      (force_abc[2] - force_bc[1]).norm() < tolerance * scale;

    deallog << "Pair lubrication forces are non-zero: "
            << to_string(force_ab[1].norm() > 0. && force_bc[0].norm() > 0.)
            << std::endl;
    deallog << "Lubrication force on each particle is the sum over its "
               "neighbors: "
            << to_string(sum_is_right) << std::endl;
  }

  // Particle-wall lubrication in the unit box, whose walls x = 0 and y = 0
  // have the boundary ids 0 and 2. The gaps to these two walls are 0.01 and
  // 0.02, the other walls are outside of the contact search radius.
  {
    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, 0., 1., true);
    triangulation.refine_global(3);
    const FESystem<dim> fe(FE_Q<dim>(1), dim + 1);
    DoFHandler<dim>     dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    Tensor<1, 3> velocity;
    velocity[0] = -0.3;
    velocity[1] = -0.2;

    Point<dim> corner_position;
    for (int d = 0; d < dim; ++d)
      corner_position[d] = 0.5;
    corner_position[0] = radius + 0.01;
    corner_position[1] = radius + 0.02;

    Point<dim> x_wall_position = corner_position;
    x_wall_position[1]         = 0.5;
    Point<dim> y_wall_position = corner_position;
    y_wall_position[0]         = 0.5;

    const Tensor<1, 3> force_corner = particle_wall_lubrication<dim>(
      parameters, dof_handler, corner_position, velocity);
    const Tensor<1, 3> force_x_wall = particle_wall_lubrication<dim>(
      parameters, dof_handler, x_wall_position, velocity);
    const Tensor<1, 3> force_y_wall = particle_wall_lubrication<dim>(
      parameters, dof_handler, y_wall_position, velocity);

    const double scale = force_x_wall.norm() + force_y_wall.norm();

    deallog << "Wall lubrication forces are non-zero: "
            << to_string(force_x_wall.norm() > 0. && force_y_wall.norm() > 0.)
            << std::endl;
    deallog << "Lubrication force near a corner is the sum over the two "
               "walls: "
            << to_string((force_corner - (force_x_wall + force_y_wall)).norm() <
                         tolerance * scale)
            << std::endl;
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
