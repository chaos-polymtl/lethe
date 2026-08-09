// SPDX-FileCopyrightText: Copyright (c) 2020-2022, 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the performance of the velocity verlet integrator
 * class. A single particle is subjected to a constant body force, a constant
 * force and a constant torque, and the complete sequence of the staggered
 * scheme is carried out: an opening half step, a series of regular steps and a
 * closing half step. For a constant acceleration, the scheme is exact, so the
 * position, the velocity and the angular velocity must match the analytical
 * solution to machine precision.
 */

// Deal.II includes
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/velocity_verlet_integrator.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining simulation general parameters. The gravitational acceleration is
  // applied along the last spatial dimension, whereas the constant force is
  // applied along the first one.
  Tensor<1, 3> g;
  g[dim - 1] = -9.81;

  // Constant force and torque applied to the particle. The torque is applied
  // around the z axis, which is the axis of rotation in 2D.
  Tensor<1, 3> applied_force;
  applied_force[0] = 1.;
  Tensor<1, 3> applied_torque;
  applied_torque[2] = 2.;

  // Time step and number of time steps
  constexpr double       dt              = 0.001;
  constexpr unsigned int number_of_steps = 5;

  // Mass and moment of inertia of the particle
  constexpr double particle_mass     = 1.;
  constexpr double moment_of_inertia = 1.;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  // Inserting one particle at the origin
  Point<dim> position1;
  int        id = 0;

  DEMSolverParameters<dim> dem_parameters;
  // Lagrangian physical properties
  Parameters::Lagrangian::LagrangianPhysicalProperties &lpp =
    dem_parameters.lagrangian_physical_properties;
  lpp.particle_type_number = 1;
  lpp.density_particle.push_back(2500);

  Particles::Particle<dim> particle1(position1, position1, id);

  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit =
    particle_handler.insert_particle(particle1, particle_cell);

  pit->get_properties()[PropertiesIndex::type] = 1;
  pit->get_properties()[PropertiesIndex::dp]   = 0.005;
  // Velocity
  pit->get_properties()[PropertiesIndex::v_x] = 0;
  pit->get_properties()[PropertiesIndex::v_y] = 0;
  pit->get_properties()[PropertiesIndex::v_z] = 0;
  // Angular velocity
  pit->get_properties()[PropertiesIndex::omega_x] = 0;
  pit->get_properties()[PropertiesIndex::omega_y] = 0;
  pit->get_properties()[PropertiesIndex::omega_z] = 0;
  // Mass
  pit->get_properties()[PropertiesIndex::mass] = particle_mass;

  std::vector<Tensor<1, 3>> torque(1);
  std::vector<Tensor<1, 3>> force(1);
  std::vector<double>       MOI(1, moment_of_inertia);

  // The integrator resets the force and the torque at the end of every step,
  // exactly like the contact force objects, which accumulate into these
  // containers, expect it to. They are therefore reapplied before every step.
  auto apply_force_and_torque = [&]() {
    force[0] += applied_force;
    torque[0] += applied_torque;
  };

  VelocityVerletIntegrator<dim, PropertiesIndex> integration_object;

  // Opening half step: the velocities are advanced by half a time step and the
  // positions by a full time step
  apply_force_and_torque();
  integration_object.integrate_start(
    particle_handler, g, dt, torque, force, MOI);

  // Regular steps
  for (unsigned int step = 1; step < number_of_steps; ++step)
    {
      apply_force_and_torque();
      integration_object.integrate(particle_handler, g, dt, torque, force, MOI);
    }

  // Closing half step: the velocities are synchronized with the positions
  apply_force_and_torque();
  integration_object.integrate_end(particle_handler, g, dt, torque, force, MOI);

  // Output
  const double final_time = number_of_steps * dt;
  deallog << "Final time: " << final_time << std::endl;

  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      const auto particle_properties = particle_iterator->get_properties();

      Tensor<1, 3> velocity;
      Tensor<1, 3> angular_velocity;
      for (int d = 0; d < 3; ++d)
        {
          velocity[d] = particle_properties[PropertiesIndex::v_x + d];
          angular_velocity[d] =
            particle_properties[PropertiesIndex::omega_x + d];
        }

      deallog << "Position: " << particle_iterator->get_location() << std::endl;
      deallog << "Velocity: " << velocity << std::endl;
      deallog << "Angular velocity: " << angular_velocity << std::endl;
    }

  // Analytical solution of a uniformly accelerated motion. The velocity Verlet
  // scheme is exact in this case, so it must be recovered to machine precision.
  // The analytical quantities are printed exactly like the ones above so that
  // both can be compared directly.
  const Tensor<1, 3> acceleration         = g + applied_force / particle_mass;
  const Tensor<1, 3> angular_acceleration = applied_torque / moment_of_inertia;

  // Only the first dim components of the acceleration contribute to the motion
  // of the particle, which lives in a dim-dimensional space
  Tensor<1, dim> analytical_position;
  for (int d = 0; d < dim; ++d)
    analytical_position[d] = 0.5 * acceleration[d] * final_time * final_time;

  deallog << "Analytical position: " << analytical_position << std::endl;
  deallog << "Analytical velocity: " << acceleration * final_time << std::endl;
  deallog << "Analytical angular velocity: "
          << angular_acceleration * final_time << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<2, DEM::DEMProperties::PropertiesIndex>();
      test<3, DEM::DEMProperties::PropertiesIndex>();
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
