// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the order of the temperature integrator by solving
 * mcdT/dt = -T for two time steps.
 */

// Tests (with common definitions)
#include <../tests/dem/test_particles_functions.h>
#include <dem/multiphysics_integrator.h>

#include <cmath>

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  const int                                 hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  const int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining simulation general parameters
  DEMSolverParameters<dim> dem_parameters;

  const double dt_1            = 0.01;
  const double dt_2            = 0.005;
  const double time_step_ratio = dt_1 / dt_2;
  const double time_final      = 1;
  const double specific_heat   = 840; // glass
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, PropertiesIndex::n_properties);

  // Inserting one particle and defining its properties
  Point<3>           position  = {0, 0, 0};
  const unsigned int id        = 0;
  const double       T_initial = 300;

  Particles::ParticleIterator<dim> pit_1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position, id);

  Tensor<1, dim>     v{{0.01, 0, 0}};
  Tensor<1, dim>     omega{{0, 0, 0}};
  const double       mass     = 1;
  const unsigned int type     = 0;
  const double       diameter = 0.005;
  set_particle_properties<dim, PropertiesIndex>(
    pit_1, type, diameter, mass, v, omega);

  pit_1->get_properties()[PropertiesIndex::T]             = T_initial;
  pit_1->get_properties()[PropertiesIndex::specific_heat] = specific_heat;

  // Initialize variables
  std::vector<double> heat_transfer_rate;
  std::vector<double> heat_source;
  heat_transfer_rate.push_back(0);
  heat_source.push_back(0);
  double T_analytical;
  double particle_temperature_error_dt_1;
  double particle_temperature_error_dt_2;
  auto   particle_properties = particle_handler.begin()->get_properties();

  // Calling temperature integrator until time_final with dt_1 time step
  double time = 0;
  while (time < time_final)
    {
      heat_transfer_rate[0] = -particle_properties[PropertiesIndex::T];
      integrate_temperature<dim, PropertiesIndex>(particle_handler,
                                                  dt_1,
                                                  heat_transfer_rate,
                                                  heat_source);
      time += dt_1;
    }

  // Calculating error with analytical solution
  T_analytical = exp(-1 / mass * 1 / specific_heat * time) * T_initial;
  particle_temperature_error_dt_1 =
    particle_properties[PropertiesIndex::T] - T_analytical;

  // Clear the particle handler and insert a new particle
  particle_handler.clear_particles();

  Particles::ParticleIterator<dim> pit_2 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position, id);

  set_particle_properties<dim, PropertiesIndex>(
    pit_2, type, diameter, mass, v, omega);
  pit_2->get_properties()[PropertiesIndex::T]             = T_initial;
  pit_2->get_properties()[PropertiesIndex::specific_heat] = specific_heat;

  // Calling temperature integrator until time_final with dt_2 time step
  time = 0;
  while (time <= time_final)
    {
      heat_transfer_rate[0] = -particle_properties[PropertiesIndex::T];
      integrate_temperature<dim, PropertiesIndex>(particle_handler,
                                                  dt_2,
                                                  heat_transfer_rate,
                                                  heat_source);
      time += dt_2;
    }

  // Calculating error with analytical solution
  T_analytical = exp(-1 / mass * 1 / specific_heat * time) * T_initial;
  particle_temperature_error_dt_2 =
    particle_properties[PropertiesIndex::T] - T_analytical;

  // Output
  deallog << "Temperature integration is a "
          << log(std::abs(particle_temperature_error_dt_1 /
                          particle_temperature_error_dt_2)) /
               log(time_step_ratio)
          << " order integration scheme." << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<3, DEM::DEMMPProperties::PropertiesIndex>();
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
