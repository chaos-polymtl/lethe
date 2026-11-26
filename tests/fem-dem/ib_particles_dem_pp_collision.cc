// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test a head-on collision of two spheres integrated by the sharp
 * immersed boundary DEM.
 *
 * With the damping of the Hertz model, the ratio of the relative velocities
 * after and before a collision is the restitution coefficient, whatever the
 * masses and the velocity. The linear momentum must be conserved.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <core/ib_particle.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <fem-dem/ib_particles_dem.h>

// Tests
#include <../tests/fem-dem/ib_particles_dem_test_utilities.h>

#include <../tests/tests.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

using namespace dealii;

template <int dim>
void
test()
{
  constexpr double radius_one = 0.1;
  constexpr double radius_two = 0.25;
  constexpr double gap        = 4e-5;
  // The relative velocity is 1
  constexpr double speed_one = 0.7;
  constexpr double speed_two = -0.3;

  // Long enough for the collision to end. The particles are heavier in 2D,
  // where the mass is proportional to the area.
  const double     duration              = (dim == 2) ? 4e-3 : 2e-3;
  constexpr double restitution_tolerance = 1e-3;
  constexpr double momentum_tolerance    = 1e-10;

  auto parameters =
    make_ib_particles_parameters<dim>((dim == 2) ? "0; 0" : "0; 0; 0");
  parameters->coupling_frequency = 2000;

  const Point<dim>     origin    = get_pair_origin<dim>();
  const Tensor<1, dim> direction = get_pair_direction<dim>();
  const Tensor<1, 3>   normal    = tensor_nd_to_3d(direction);

  deallog << "Dimension " << dim << std::endl;

  const std::vector<std::pair<std::string, double>> cases{{"1", 1.},
                                                          {"0.9", 0.9},
                                                          {"0.5", 0.5}};
  for (const auto &[name, restitution] : cases)
    {
      std::vector<IBParticle<dim>> particles{
        make_ib_sphere<dim>(0, radius_one, origin),
        make_ib_sphere<dim>(
          1, radius_two, origin + (radius_one + radius_two + gap) * direction)};
      particles[0].velocity = speed_one * normal;
      particles[1].velocity = speed_two * normal;
      for (auto &particle : particles)
        {
          particle.restitution_coefficient = restitution;
          particle.initialize_previous_solution();
        }
      const Tensor<1, 3> initial_momentum =
        particles[0].mass * particles[0].velocity +
        particles[1].mass * particles[1].velocity;

      IBParticlesDEM<dim> ib_dem;
      ib_dem.initialize(
        parameters, nullptr, nullptr, MPI_COMM_WORLD, particles);
      ib_dem.update_particles(particles, 0.);
      ib_dem.update_contact_candidates();
      ib_dem.integrate_particles_motion(duration, 0., 0., 0., 0.);

      const auto  &result = ib_dem.dem_particles;
      const double rebound_speed =
        (result[1].velocity - result[0].velocity) * normal;
      const Tensor<1, 3> momentum = result[0].mass * result[0].velocity +
                                    result[1].mass * result[1].velocity;

      log_check("Rebound is the restitution coefficient " + name,
                std::abs(rebound_speed - restitution) < restitution_tolerance);
      log_check("Momentum is conserved with the restitution coefficient " +
                  name,
                (momentum - initial_momentum).norm() <
                  momentum_tolerance * initial_momentum.norm());
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
