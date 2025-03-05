// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the performance of non-linear (Hertzian)
 * particle-particle contact force  is checked.
 */


#include <../tests/dem/test_particles_functions.h>


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

  Tensor<1, dim> g{{0, 0, -9.81}};
  double         dt                                               = 0.00001;
  double         particle_diameter                                = 0.005;
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
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties<PropertiesIndex>());

  // Creating containers manager for finding cell neighbor and also broad and
  // fine particle-particle search objects
  DEMContactManager<dim, PropertiesIndex> contact_manager;

  // Finding cell neighbors
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);


  // Inserting two particles in contact
  Point<3> position1 = {0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.40499, 0, 0};
  int      id2       = 1;

  // Constructing particle iterators from particle positions (inserting
  // particles)
  Particles::ParticleIterator<dim> pit1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position1, id1);
  Particles::ParticleIterator<dim> pit2 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position2, id2);

  // Setting particle properties
  Tensor<1, dim> v1{{0.01, 0, 0}};
  Tensor<1, dim> omega1{{0, 0, 0}};
  Tensor<1, dim> v2{{0, 0, 0}};
  Tensor<1, dim> omega2{{0, 0, 0}};
  double         mass = 1;
  int            type = 0;
  set_particle_properties<dim, PropertiesIndex>(
    pit1, type, particle_diameter, mass, v1, omega1);
  set_particle_properties<dim, PropertiesIndex>(
    pit2, type, particle_diameter, mass, v2, omega2);

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  for (auto &moi_val : MOI)
    moi_val = 1;

  contact_manager.update_local_particles_in_cells(particle_handler);

  // Dummy Adaptive sparse contacts object and particle-particle broad search
  AdaptiveSparseContacts<dim, PropertiesIndex> dummy_adaptive_sparse_contacts;
  contact_manager.execute_particle_particle_broad_search(
    particle_handler, dummy_adaptive_sparse_contacts);

  // Calling fine search
  contact_manager.execute_particle_particle_fine_search(neighborhood_threshold);

  // Calling linear force
  ParticleParticleContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleParticleContactForceModel::
      hertz_mindlin_limit_overlap,
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>
    nonlinear_force_object(dem_parameters);
  nonlinear_force_object.calculate_particle_particle_contact_force(
    contact_manager.get_local_adjacent_particles(),
    contact_manager.get_ghost_adjacent_particles(),
    contact_manager.get_local_local_periodic_adjacent_particles(),
    contact_manager.get_local_ghost_periodic_adjacent_particles(),
    contact_manager.get_ghost_local_periodic_adjacent_particles(),
    dt,
    torque,
    force);

  // Output
  auto particle = particle_handler.begin();
  deallog << "The contact force vector for particle 1 is: "
          << force[particle->get_id()][0] << " " << force[particle->get_id()][1]
          << " " << force[particle->get_id()][2] << " N " << std::endl;
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
