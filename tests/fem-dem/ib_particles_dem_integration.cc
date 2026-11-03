// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the time integration of the DEM sub-stepping of the sharp
 * immersed boundary solver, without fluid and without contact.
 *
 * The first particle is integrated under a constant gravity and a constant
 * torque around z. The RK4 scheme is exact for a constant acceleration, so its
 * position, velocity, angular velocity and orientation at the end of the CFD
 * time step must match the analytical solution to round-off.
 *
 * The second particle is loaded from a file and is not integrated. It moves
 * with its constant velocity and angular velocity, so its position and
 * orientation must have advanced by one CFD time step.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <core/ib_particle.h>
#include <core/parameters.h>

#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/fem-dem/ib_particles_dem_test_utilities.h>

#include <../tests/tests.h>

#include <vector>

using namespace dealii;

template <int dim>
void
test()
{
  const double dt        = 0.01;
  const double tolerance = 1e-12;

  Tensor<1, 3> gravity;
  gravity[0] = 1.5;
  gravity[1] = -9.81;
  if constexpr (dim == 3)
    gravity[2] = 0.5;

  auto parameters = make_ib_particles_parameters<dim>(
    (dim == 2) ? "1.5; -9.81" : "1.5; -9.81; 0.5");
  parameters->coupling_frequency       = 4;
  parameters->load_particles_from_file = true;

  // Integrated particle.
  Point<dim> initial_position;
  initial_position[0] = 0.1;
  initial_position[1] = 0.2;
  if constexpr (dim == 3)
    initial_position[2] = 0.3;
  IBParticle<dim> integrated = make_ib_sphere<dim>(0, 0.1, initial_position);

  Tensor<1, 3> initial_velocity;
  initial_velocity[0] = 1.;
  initial_velocity[1] = 2.;
  if constexpr (dim == 3)
    initial_velocity[2] = -0.5;
  const double initial_omega       = 1.;
  const double initial_orientation = 0.2;
  const double inertia             = 0.25;
  const double torque              = 0.5;

  integrated.velocity       = initial_velocity;
  integrated.omega[2]       = initial_omega;
  integrated.orientation[2] = initial_orientation;
  integrated.inertia        = 0;
  for (unsigned int d = 0; d < 3; ++d)
    integrated.inertia[d][d] = inertia;
  integrated.fluid_torque[2] = torque;
  integrated.set_orientation(integrated.orientation);
  integrated.initialize_previous_solution();

  // Prescribed particle, far from the integrated one.
  Point<dim> prescribed_initial_position;
  prescribed_initial_position[0] = 5.;
  IBParticle<dim> prescribed =
    make_ib_sphere<dim>(1, 0.1, prescribed_initial_position);
  prescribed.integrate_motion = false;
  prescribed.velocity[0]      = 0.5;
  prescribed.velocity[1]      = -0.25;
  if constexpr (dim == 3)
    prescribed.velocity[2] = 0.75;
  prescribed.omega[2] = 3.;
  if constexpr (dim == 3)
    {
      prescribed.omega[0] = 1.;
      prescribed.omega[1] = -2.;
    }
  prescribed.initialize_previous_solution();

  std::vector<IBParticle<dim>> particles{integrated, prescribed};

  IBParticlesDEM<dim> ib_dem;
  ib_dem.initialize(parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
  ib_dem.update_particles(particles, 0.);
  ib_dem.update_contact_candidates();
  ib_dem.integrate_particles_motion(dt, 0., 0., 0., 0.);

  const IBParticle<dim> &integrated_result = ib_dem.dem_particles[0];
  const IBParticle<dim> &prescribed_result = ib_dem.dem_particles[1];

  // Analytical solution of the integrated particle.
  const double angular_acceleration = torque / inertia;
  Point<dim>   expected_position;
  for (int d = 0; d < dim; ++d)
    expected_position[d] = initial_position[d] + initial_velocity[d] * dt +
                           0.5 * gravity[d] * dt * dt;
  const Tensor<1, 3> expected_velocity = initial_velocity + gravity * dt;
  Tensor<1, 3>       expected_omega;
  expected_omega[2] = initial_omega + angular_acceleration * dt;
  Tensor<1, 3> expected_orientation;
  expected_orientation[2] = initial_orientation + initial_omega * dt +
                            0.5 * angular_acceleration * dt * dt;

  // Constant motion of the prescribed particle.
  Point<dim> expected_prescribed_position;
  for (int d = 0; d < dim; ++d)
    expected_prescribed_position[d] =
      prescribed_initial_position[d] + prescribed.velocity[d] * dt;
  const Tensor<1, 3> expected_prescribed_orientation = prescribed.omega * dt;

  auto to_string = [](const bool value) { return value ? "true" : "false"; };

  deallog << "Dimension " << dim << std::endl;
  deallog << "Integrated particle position is exact: "
          << to_string(integrated_result.position.distance(expected_position) <
                       tolerance)
          << std::endl;
  deallog << "Integrated particle velocity is exact: "
          << to_string((integrated_result.velocity - expected_velocity).norm() <
                       tolerance)
          << std::endl;
  deallog << "Integrated particle angular velocity is exact: "
          << to_string((integrated_result.omega - expected_omega).norm() <
                       tolerance)
          << std::endl;
  deallog << "Integrated particle orientation is exact: "
          << to_string(
               (integrated_result.orientation - expected_orientation).norm() <
               tolerance)
          << std::endl;
  deallog << "Prescribed particle position is exact: "
          << to_string(prescribed_result.position.distance(
                         expected_prescribed_position) < tolerance)
          << std::endl;
  deallog << "Prescribed particle orientation is exact: "
          << to_string(
               (prescribed_result.orientation - expected_prescribed_orientation)
                 .norm() < tolerance)
          << std::endl;
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
