// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the performance of linear (Hookean) particle-particle
 * contact force  is checked.
 */


#include <../tests/dem/full_contact_functions.h>



template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim> mapping(1);


  // Defining general simulation parameters
  DEMSolverParameters<dim>                              dem_parameters;
  Parameters::Lagrangian::LagrangianPhysicalProperties &lagrangian_prop =
    dem_parameters.lagrangian_physical_properties;
  Parameters::Lagrangian::ModelParameters &model_param =
    dem_parameters.model_parameters;

  Tensor<1, dim> g{{0, 0, 0}};
  double         dt                                               = 0.00001;
  double         particle_diameter                                = 0.005;
  unsigned int   output_step                                      = 10;
  lagrangian_prop.particle_type_number                            = 1;
  lagrangian_prop.youngs_modulus_particle[0]                      = 50000000;
  lagrangian_prop.poisson_ratio_particle[0]                       = 0.3;
  lagrangian_prop.restitution_coefficient_particle[0]             = 0.5;
  lagrangian_prop.friction_coefficient_particle[0]                = 0.5;
  lagrangian_prop.rolling_viscous_damping_coefficient_particle[0] = 0.5;
  lagrangian_prop.rolling_friction_coefficient_particle[0]        = 0.1;
  lagrangian_prop.surface_energy_particle[0]                      = 0.;
  lagrangian_prop.hamaker_constant_particle[0]                    = 0.;
  lagrangian_prop.density_particle[0]                             = 2500;
  model_param.rolling_resistance_method =
    RollingResistanceMethod::constant_resistance;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties<PropertiesIndex>());

  // Creating containers manager for finding cell neighbor and also broad and
  // fine particle-particle search objects
  DEMContactManager<dim, PropertiesIndex> contact_manager;

  // Finding cell neighbors list, it is required for finding the broad search
  // pairs in the contact_manager
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // Inserting two particles at distance 0
  Point<dim>                        position1 = {0.4, 0, 0};
  Point<dim>                        position2 = {0.405, 0, 0};
  Tensor<1, dim>                    v1{{0.01, 0, 0}};
  Tensor<1, dim>                    w1{{0, 0, 0}};
  Tensor<1, dim>                    v2{{0, 0, 0}};
  Tensor<1, dim>                    w2{{0, 0, 0}};
  initial_particles_properties<dim> particle_properties = {
    {position1, position2}, // initial positions
    {0, 1},                 // id
    0,                      // type
    {v1, v2},               // initial velocities
    {w1, w2},               // initial angular velocities
    1,                      // mass
    particle_diameter       // diameter
  };


  // linear

  contact_output linear_output =
    simul_full_contact<dim,
                       PropertiesIndex,
                       ParticleParticleContactForceModel::linear,
                       RollingResistanceMethod::constant_resistance>(
      triangulation,
      particle_handler,
      contact_manager,
      dem_parameters,
      particle_properties,
      g,
      dt,
      output_step,
      neighborhood_threshold,
      "linear");

  auto particle0 = particle_handler.begin();
  deallog << "The new position of particle 0 after end of contact is ("
          << particle0->get_location()[0] << " " << particle0->get_location()[1]
          << " " << particle0->get_location()[2] << ") with linear force."
          << std::endl;



  particle_handler.clear_particles();



  // hertz_mindlin_limit_overlap

  contact_output hmlo_output = simul_full_contact<
    dim,
    PropertiesIndex,
    ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
    RollingResistanceMethod::constant_resistance>(triangulation,
                                                  particle_handler,
                                                  contact_manager,
                                                  dem_parameters,
                                                  particle_properties,
                                                  g,
                                                  dt,
                                                  output_step,
                                                  neighborhood_threshold,
                                                  "hmlo");

  auto particle = particle_handler.begin();
  deallog << "The new position of particle 0 after end of contact is ("
          << particle->get_location()[0] << " " << particle->get_location()[1]
          << " " << particle->get_location()[2]
          << ") with hertz_mindlin_limit_overlap force." << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();
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
