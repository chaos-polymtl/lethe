// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


/**
 * @brief This test checks the performance of the temperature integrator.
 *
 */

// Deal.ii
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>

// Tests (with common definitions)
#include <../tests/dem/multiphysics_integrator.h>

#include <../tests/tests.h>

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
  DEMSolverParameters<dim>                              dem_parameters;
  Parameters::Lagrangian::LagrangianPhysicalProperties &lagrangian_prop =
    dem_parameters.lagrangian_physical_properties;
  Parameters::Lagrangian::ModelParameters &model_param =
    dem_parameters.model_parameters;

  const double dt                                  = 0.00001;
  lagrangian_prop.particle_type_number             = 1;
  lagrangian_prop.youngs_modulus_particle[0]       = 50000000;
  lagrangian_prop.thermal_conductivity_particle[0] = 0.6;
  lagrangian_prop.heat_capacity_particle[0]        = 4;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties<PropertiesIndex>());

  // Inserting one particle and defining its properties
  Point<3> position1 = {0, 0, 0};
  int      id        = 0;

  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());

  Particles::ParticleIterator<dim> pit =
    particle_handler.insert_particle(particle1, particle_cell);

  pit->get_properties()[PropertiesIndex::type]    = 0;
  pit->get_properties()[PropertiesIndex::dp]      = 0.005;
  pit->get_properties()[PropertiesIndex::v_x]     = 0;
  pit->get_properties()[PropertiesIndex::v_y]     = 0;
  pit->get_properties()[PropertiesIndex::v_z]     = 0;
  pit->get_properties()[PropertiesIndex::omega_x] = 0;
  pit->get_properties()[PropertiesIndex::omega_y] = 0;
  pit->get_properties()[PropertiesIndex::omega_z] = 0;
  pit->get_properties()[PropertiesIndex::mass]    = 1;
  pit->get_properties()[PropertiesIndex::T]       = 300;

  std::vector<double> heat_transfer;
  double              heat_source;
  heat_transfer.push_back(0);

  // Calling temperature integrator for max_iteration iterations
  const unsigned int max_iteration    = 5000;
  double             time             = 0;
  int                output_frequency = 10;
  auto               particle         = particle_handler.begin();
  std::cout << "time temperature" << std::endl;

  for (unsigned int iteration = 0; iteration < max_iteration; ++iteration)
    {
      heat_source = 400 * cos(time);
      integrate_temperature<dim, PropertiesIndex>(
        particle_handler, dem_parameters, dt, heat_transfer, heat_source);

      if (iteration % output_frequency == 0)
        {
          std::cout << time << " "
                    << particle->get_properties()[PropertiesIndex::T]
                    << std::endl;
        }
      time += dt;
    }
  std::cout << time << " " << particle->get_properties()[PropertiesIndex::T]
            << std::endl;

  // Output
  deallog << "The temperature of the particle at time " << time
          << " seconds is: " << particle->get_properties()[PropertiesIndex::T]
          << "." << std::endl;
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
