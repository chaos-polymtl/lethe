// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the performance of non-linear (Hertzian)
 * particle-particle contact force is checked.
 */

// Deal.II

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/velocity_verlet_integrator.h>

// Tests (with common definitions)

#include <../tests/tests.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
reinitialize_contact_outcomes(
  Particles::ParticleHandler<dim>              &particle_handler,
  ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{
  const unsigned int number_of_particles =
    particle_handler.n_locally_owned_particles();
  contact_outcome.resize_interaction_containers(number_of_particles);

  for (unsigned int i = 0; i < number_of_particles; ++i)
    {
      // Reinitializing contact outcomes of particles in the system
      contact_outcome.force[i][0] = 0;
      contact_outcome.force[i][1] = 0;
      contact_outcome.force[i][2] = 0;

      contact_outcome.torque[i][0] = 0;
      contact_outcome.torque[i][1] = 0;
      contact_outcome.torque[i][2] = 0;
    }
}

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
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining general simulation parameters
  Tensor<1, 3> g{{0, 0, 0}};
  double       dt                                                    = 0.00001;
  double       particle_diameter                                     = 0.005;
  unsigned int step_end                                              = 1000;
  unsigned int output_frequency                                      = 10;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_particle[0] =
    50000000;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_particle[0] = 0.9;
  dem_parameters.lagrangian_physical_properties
    .restitution_coefficient_particle[0] = 0.9;
  dem_parameters.lagrangian_physical_properties
    .friction_coefficient_particle[0] = 0.5;
  dem_parameters.lagrangian_physical_properties
    .rolling_viscous_damping_coefficient_particle[0] = 0.5;
  dem_parameters.lagrangian_physical_properties
    .rolling_friction_coefficient_particle[0] = 0.1;
  dem_parameters.lagrangian_physical_properties.surface_energy_particle[0] = 0.;
  dem_parameters.lagrangian_physical_properties.hamaker_constant_particle[0] =
    0.;
  dem_parameters.lagrangian_physical_properties.density_particle[0] = 2500;
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, PropertiesIndex::n_properties);

  typename dem_data_structures<2>::particle_index_iterator_map
    local_particle_container;

  DEMContactManager<dim, PropertiesIndex> contact_manager;

  // Finding cell neighbors
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // Creating particle-particle force objects
  ParticleParticleContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleParticleContactForceModel::
      hertz_mindlin_limit_overlap,
    Parameters::Lagrangian::RollingResistanceMethod::constant>
    nonlinear_force_object(dem_parameters);
  VelocityVerletIntegrator<dim, PropertiesIndex> integrator_object;

  MPI_Comm communicator     = triangulation.get_communicator();
  auto     this_mpi_process = Utilities::MPI::this_mpi_process(communicator);

  // Inserting two particles in contact
  Point<2> position1 = {0, 0.003};
  int      id1       = 0;
  Point<2> position2 = {0, -0.003};
  int      id2       = 1;

  // Particle 1 is inserted in a cell owned by process1
  if (this_mpi_process == 1)
    {
      Particles::Particle<dim> particle1(position1, position1, id1);

      typename Triangulation<dim>::active_cell_iterator cell1 =
        GridTools::find_active_cell_around_point(triangulation,
                                                 particle1.get_location());
      Particles::ParticleIterator<dim> pit1 =
        particle_handler.insert_particle(particle1, cell1);
      pit1->get_properties()[PropertiesIndex::type]    = 0;
      pit1->get_properties()[PropertiesIndex::dp]      = particle_diameter;
      pit1->get_properties()[PropertiesIndex::v_x]     = 0;
      pit1->get_properties()[PropertiesIndex::v_y]     = -0.5;
      pit1->get_properties()[PropertiesIndex::v_z]     = 0;
      pit1->get_properties()[PropertiesIndex::omega_x] = 0;
      pit1->get_properties()[PropertiesIndex::omega_y] = 0;
      pit1->get_properties()[PropertiesIndex::omega_z] = 0;
      pit1->get_properties()[PropertiesIndex::mass]    = 1;
    }

  // Particle 2 is inserted in a cell owned by process0
  if (this_mpi_process == 0)
    {
      Particles::Particle<dim> particle2(position2, position2, id2);
      typename Triangulation<dim>::active_cell_iterator cell2 =
        GridTools::find_active_cell_around_point(triangulation,
                                                 particle2.get_location());
      Particles::ParticleIterator<dim> pit2 =
        particle_handler.insert_particle(particle2, cell2);
      pit2->get_properties()[PropertiesIndex::type]    = 0;
      pit2->get_properties()[PropertiesIndex::dp]      = particle_diameter;
      pit2->get_properties()[PropertiesIndex::v_x]     = 0;
      pit2->get_properties()[PropertiesIndex::v_y]     = 0.5;
      pit2->get_properties()[PropertiesIndex::v_z]     = 0;
      pit2->get_properties()[PropertiesIndex::omega_x] = 0;
      pit2->get_properties()[PropertiesIndex::omega_y] = 0;
      pit2->get_properties()[PropertiesIndex::omega_z] = 0;
      pit2->get_properties()[PropertiesIndex::mass]    = 1;
    }

  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  std::vector<double>                          MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  const unsigned int number_of_particles =
    particle_handler.get_max_local_particle_index();
  contact_outcome.resize_interaction_containers(number_of_particles);
  MOI.resize(number_of_particles);
  for (auto &moi_val : MOI)
    moi_val = 1;

  for (unsigned int iteration = 0; iteration < step_end; ++iteration)
    {
      // Reinitializing contact outcomes
      reinitialize_contact_outcomes<dim, PropertiesIndex>(particle_handler,
                                                          contact_outcome);

      particle_handler.exchange_ghost_particles();

      contact_manager.update_local_particles_in_cells(particle_handler);

      // Dummy Adaptive sparse contacts object and particle-particle broad
      // search
      AdaptiveSparseContacts<dim, PropertiesIndex>
        dummy_adaptive_sparse_contacts;
      contact_manager.execute_particle_particle_broad_search(
        particle_handler, dummy_adaptive_sparse_contacts);

      // Calling fine search
      contact_manager.execute_particle_particle_fine_search(
        neighborhood_threshold);

      // Integration
      // Calling non-linear force
      nonlinear_force_object.calculate_particle_particle_contact(
        contact_manager.get_local_adjacent_particles(),
        contact_manager.get_ghost_adjacent_particles(),
        contact_manager.get_local_local_periodic_adjacent_particles(),
        contact_manager.get_local_ghost_periodic_adjacent_particles(),
        contact_manager.get_ghost_local_periodic_adjacent_particles(),
        dt,
        contact_outcome);

      // Integration
      integrator_object.integrate(particle_handler,
                                  g,
                                  dt,
                                  contact_outcome.torque,
                                  contact_outcome.force,
                                  MOI);

      contact_manager.update_contacts();

      if (iteration % output_frequency == 0)
        {
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 1)
            {
              // Output particle 0 position
              for (auto particle = particle_handler.begin();
                   particle != particle_handler.end();
                   ++particle)
                {
                  if (particle->get_id() == 0)
                    {
                      deallog
                        << "The location of particle " << particle->get_id()
                        << " is: " << particle->get_location() << std::endl;
                    }
                }
            }
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<2, DEM::DEMProperties::PropertiesIndex>();
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
