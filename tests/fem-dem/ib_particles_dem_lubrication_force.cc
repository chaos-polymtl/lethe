// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the lubrication force of the sharp immersed boundary DEM against
 * its closed form.
 *
 * Particle-particle: two spheres of different radii move with arbitrary
 * velocities and their line of centers is not aligned with an axis. The force
 * on the first sphere must be
 * 3/2 pi mu (d_1 d_2 / (d_1 + d_2))^2 (1/h - 1/h_max) ((v_2 - v_1) . e) e,
 * where e is the unit vector from the second sphere to the first one, and the
 * force on the second sphere must be its opposite.
 *
 * Particle-wall: a sphere moves toward a flat wall. The diameter of the wall
 * is infinite, so the force must be
 * 6 pi mu R^2 (1/h - 1/h_max) ((v_wall - v) . e) e,
 * where e is the unit normal from the wall to the sphere.
 *
 * The distance between the centers and the distance from the center to the
 * wall differ from 1, so a force that depends on the length unit is detected.
 *
 * For both, the gap is limited to h_min and the force is zero beyond h_max.
 * The velocity of a moving wall is taken into account and an inverted range is
 * rejected.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
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

#include <algorithm>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

using namespace dealii;

// Parameters of the lubrication model.
constexpr double h_max  = 0.05;
constexpr double h_min  = 0.001;
constexpr double mu     = 0.01;
constexpr double dt_dem = 1e-5;
constexpr double gap    = 0.02;

/**
 * @brief Closed form of the lubrication force between two bodies.
 *
 * @param[in] effective_diameter Effective diameter d_1 d_2 / (d_1 + d_2) of the
 * two bodies.
 * @param[in] unit_vector Unit vector from the second body to the first one.
 * @param[in] relative_velocity Velocity of the second body minus the velocity
 * of the first one.
 * @param[in] surface_gap Gap between the two surfaces.
 *
 * @return The lubrication force on the first body.
 */
Tensor<1, 3>
expected_lubrication_force(const double        effective_diameter,
                           const Tensor<1, 3> &unit_vector,
                           const Tensor<1, 3> &relative_velocity,
                           const double        surface_gap = gap)
{
  if (surface_gap >= h_max)
    return Tensor<1, 3>();
  return 1.5 * std::numbers::pi * mu * effective_diameter * effective_diameter *
         (1. / std::max(surface_gap, h_min) - 1. / h_max) *
         (relative_velocity * unit_vector) * unit_vector;
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

  // Particle-particle lubrication. The distance between the centers is
  // radius_one + radius_two + gap = 0.37.
  {
    const double radius_one = 0.1;
    const double radius_two = 0.25;

    // Unit vector from the first sphere to the second one.
    Tensor<1, dim> direction;
    if constexpr (dim == 2)
      {
        direction[0] = 0.6;
        direction[1] = 0.8;
      }
    else
      {
        direction[0] = 2. / 3.;
        direction[1] = 1. / 3.;
        direction[2] = 2. / 3.;
      }

    Point<dim> position_one;
    for (int d = 0; d < dim; ++d)
      position_one[d] = 0.1 * (d + 1);
    const Point<dim> position_two =
      position_one + (radius_one + radius_two + gap) * direction;

    Tensor<1, 3> velocity_one;
    velocity_one[0] = 0.3;
    velocity_one[1] = 0.1;
    Tensor<1, 3> velocity_two;
    velocity_two[0] = -0.2;
    velocity_two[1] = 0.4;
    if constexpr (dim == 3)
      {
        velocity_one[2] = -0.2;
        velocity_two[2] = 0.1;
      }

    std::vector<IBParticle<dim>> particles{
      make_ib_sphere<dim>(0, radius_one, position_one),
      make_ib_sphere<dim>(1, radius_two, position_two)};
    particles[0].velocity = velocity_one;
    particles[1].velocity = velocity_two;
    for (auto &particle : particles)
      particle.initialize_previous_solution();

    IBParticlesDEM<dim> ib_dem;
    ib_dem.initialize(parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
    ib_dem.update_contact_candidates();

    std::vector<Tensor<1, 3>> force(particles.size());
    std::vector<Tensor<1, 3>> torque(particles.size());
    ib_dem.calculate_pp_lubrication_force(
      dt_dem, h_max, h_min, mu, force, torque);

    const double       diameter_one = 2. * radius_one;
    const double       diameter_two = 2. * radius_two;
    const Tensor<1, 3> expected_force =
      expected_lubrication_force(diameter_one * diameter_two /
                                   (diameter_one + diameter_two),
                                 -tensor_nd_to_3d(direction),
                                 velocity_two - velocity_one);

    deallog << "Particle-particle lubrication force matches the closed form: "
            << to_string((force[0] - expected_force).norm() <
                         tolerance * expected_force.norm())
            << std::endl;
    deallog << "Particle-particle lubrication forces are opposite: "
            << to_string((force[0] + force[1]).norm() <
                         tolerance * expected_force.norm())
            << std::endl;
  }

  // Particle-wall lubrication in the unit box. The sphere is at a distance
  // radius + gap = 0.12 of the wall x = 0, the other walls are outside of the
  // lubrication range.
  {
    const double radius = 0.1;

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, 0., 1., true);
    triangulation.refine_global(3);
    const FESystem<dim> fe(FE_Q<dim>(1), dim + 1);
    DoFHandler<dim>     dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    Point<dim> position;
    for (int d = 0; d < dim; ++d)
      position[d] = 0.5;
    position[0] = radius + gap;

    Tensor<1, 3> velocity;
    velocity[0] = -0.3;
    velocity[1] = 0.2;
    if constexpr (dim == 3)
      velocity[2] = 0.1;

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

    // The wall is at rest and its normal toward the sphere is along x.
    Tensor<1, 3> wall_to_particle;
    wall_to_particle[0] = 1.;
    const Tensor<1, 3> expected_force =
      expected_lubrication_force(2. * radius, wall_to_particle, -velocity);

    deallog << "Particle-wall lubrication force matches the closed form: "
            << to_string((force[0] - expected_force).norm() <
                         tolerance * expected_force.norm())
            << std::endl;
  }

  // Variations of the gap and of the wall motion.
  {
    const double         radius_one = 0.1;
    const double         radius_two = 0.25;
    const Point<dim>     origin     = get_pair_origin<dim>();
    const Tensor<1, dim> direction  = get_pair_direction<dim>();

    Tensor<1, 3> velocity;
    velocity[0] = -0.3;
    velocity[1] = 0.2;

    const auto pair_is_right = [&](const double surface_gap) {
      std::vector<IBParticle<dim>> particles{
        make_ib_sphere<dim>(0, radius_one, origin),
        make_ib_sphere<dim>(1,
                            radius_two,
                            origin + (radius_one + radius_two + surface_gap) *
                                       direction)};
      particles[0].velocity = velocity;

      IBParticlesDEM<dim> ib_dem;
      ib_dem.initialize(
        parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
      ib_dem.update_contact_candidates();
      std::vector<Tensor<1, 3>> force(particles.size());
      std::vector<Tensor<1, 3>> torque(particles.size());
      ib_dem.calculate_pp_lubrication_force(
        dt_dem, h_max, h_min, mu, force, torque);

      const Tensor<1, 3> expected =
        expected_lubrication_force(2. * radius_one * radius_two /
                                     (radius_one + radius_two),
                                   -tensor_nd_to_3d(direction),
                                   -velocity,
                                   surface_gap);
      return (force[0] - expected).norm() <= tolerance * expected.norm();
    };

    Triangulation<dim> triangulation;
    GridGenerator::hyper_cube(triangulation, 0., 1., true);
    triangulation.refine_global(3);
    const FESystem<dim> fe(FE_Q<dim>(1), dim + 1);
    DoFHandler<dim>     dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);

    const auto wall_is_right = [&](const double        surface_gap,
                                   const Tensor<1, 3> &wall_velocity) {
      auto boundary_conditions =
        std::make_shared<Parameters::Lagrangian::BCDEM>();
      // The wall x = 0 has the boundary id 0
      boundary_conditions->boundary_translational_velocity[0] = wall_velocity;

      Point<dim> position;
      for (int d = 0; d < dim; ++d)
        position[d] = 0.5;
      position[0] = radius_one + surface_gap;
      std::vector<IBParticle<dim>> particles{
        make_ib_sphere<dim>(0, radius_one, position)};
      particles[0].velocity = velocity;

      IBParticlesDEM<dim> ib_dem;
      ib_dem.initialize(
        parameters, nullptr, boundary_conditions, MPI_COMM_WORLD, particles);
      ib_dem.update_particles_boundary_contact(particles,
                                               dof_handler,
                                               QGauss<dim - 1>(2),
                                               MappingQ1<dim>());
      std::vector<Tensor<1, 3>> force(particles.size());
      std::vector<Tensor<1, 3>> torque(particles.size());
      ib_dem.calculate_pw_lubrication_force(
        dt_dem, h_max, h_min, mu, force, torque);

      Tensor<1, 3> wall_to_particle;
      wall_to_particle[0] = 1.;
      const Tensor<1, 3> expected =
        expected_lubrication_force(2. * radius_one,
                                   wall_to_particle,
                                   wall_velocity - velocity,
                                   surface_gap);
      return (force[0] - expected).norm() <= tolerance * expected.norm();
    };

    const Tensor<1, 3> wall_at_rest;
    Tensor<1, 3>       moving_wall;
    moving_wall[0] = 0.15;
    moving_wall[1] = 0.4;

    const double gap_below_minimum = 0.5 * h_min;
    const double gap_above_maximum = 1.4 * h_max;

    log_check("Particle-particle gap is limited to h_min",
              pair_is_right(gap_below_minimum));
    log_check("Particle-wall gap is limited to h_min",
              wall_is_right(gap_below_minimum, wall_at_rest));
    log_check("Particle-particle force is zero beyond h_max",
              pair_is_right(gap_above_maximum));
    log_check("Particle-wall force is zero beyond h_max",
              wall_is_right(gap_above_maximum, wall_at_rest));
    log_check("Particle-wall force accounts for the velocity of the wall",
              wall_is_right(gap, moving_wall));
  }

  // Inverted range.
  {
    Parameters::IBParticles<dim> parsed_parameters;
    ParameterHandler             prm;
    parsed_parameters.declare_parameters(prm);
    prm.enter_subsection("particles");
    prm.enter_subsection("DEM");
    prm.set("enable lubrication force", "true");
    prm.set("lubrication range max", "0.05");
    prm.set("lubrication range min", "0.2");
    prm.leave_subsection();
    prm.leave_subsection();

    bool range_is_rejected = false;
    try
      {
        parsed_parameters.parse_parameters(prm);
      }
    catch (const ExceptionBase &exception)
      {
        range_is_rejected =
          std::string(exception.what()).find("lubrication range min") !=
          std::string::npos;
      }
    log_check("Inverted lubrication range is rejected", range_is_rejected);
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
