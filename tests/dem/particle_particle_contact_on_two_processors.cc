/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019-
 */

/**
 * @brief In this test, the performance of non-linear (Hertzian)
 * particle-particle contact force is checked.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_fine_search.h>
#include <dem/particle_particle_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

void
update_contact_containers(
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    &local_adjacent_particles,
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    &ghost_adjacent_particles,
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    &cleared_local_adjacent_particles,
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    &cleared_ghost_adjacent_particles)
{
  local_adjacent_particles.clear();
  ghost_adjacent_particles.clear();

  local_adjacent_particles = cleared_local_adjacent_particles;
  ghost_adjacent_particles = cleared_ghost_adjacent_particles;
}

template <int dim>
void
update_ghost_particle_particle_contact_container_iterators(
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int,
                       particle_particle_contact_info_struct<dim>>>
    &cleared_ghost_adjacent_particles,
  const std::unordered_map<unsigned int, Particles::ParticleIterator<dim>>
    &local_particle_container)
{
  for (auto adjacent_particles_iterator =
         cleared_ghost_adjacent_particles.begin();
       adjacent_particles_iterator != cleared_ghost_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id          = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto particle_particle_map_iterator =
             pairs_in_contant_content->begin();
           particle_particle_map_iterator != pairs_in_contant_content->end();
           ++particle_particle_map_iterator)
        {
          int particle_two_id = particle_particle_map_iterator->first;
          particle_particle_map_iterator->second.particle_one =
            local_particle_container.at(particle_one_id);
          particle_particle_map_iterator->second.particle_two =
            local_particle_container.at(particle_two_id);
        }
    }
}

template <int dim>
void
update_local_particle_particle_contact_container_iterators(
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int,
                       particle_particle_contact_info_struct<dim>>>
    &cleared_local_adjacent_particles,
  const std::unordered_map<unsigned int, Particles::ParticleIterator<dim>>
    &local_particle_container)
{
  for (auto adjacent_particles_iterator =
         cleared_local_adjacent_particles.begin();
       adjacent_particles_iterator != cleared_local_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id          = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto particle_particle_map_iterator =
             pairs_in_contant_content->begin();
           particle_particle_map_iterator != pairs_in_contant_content->end();
           ++particle_particle_map_iterator)
        {
          int particle_two_id = particle_particle_map_iterator->first;

          particle_particle_map_iterator->second.particle_one =
            local_particle_container.at(particle_one_id);
          particle_particle_map_iterator->second.particle_two =
            local_particle_container.at(particle_two_id);
        }
    }
}

template <int dim>
void
update_local_particle_container(
  std::unordered_map<unsigned int, Particles::ParticleIterator<dim>>
    &                              local_particle_container,
  Particles::ParticleHandler<dim> *particle_handler)
{
  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
      local_particle_container[particle_iterator->get_id()] = particle_iterator;
    }

  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
      local_particle_container[particle_iterator->get_id()] = particle_iterator;
    }

  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
      local_particle_container[particle_iterator->get_id()] = particle_iterator;
    }
}

template <int dim>
void
locate_local_particles_in_cells(
  Particles::ParticleHandler<dim> &particle_handler,
  std::unordered_map<unsigned int, Particles::ParticleIterator<2>>
    &local_particle_container,
  std::unordered_map<unsigned int, Particles::ParticleIterator<2>>
    &ghost_particle_container,
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    &cleared_local_adjacent_particles,
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    &cleared_ghost_adjacent_particles)
{
  local_particle_container.clear();
  ghost_particle_container.clear();

  update_local_particle_container(local_particle_container, &particle_handler);

  update_local_particle_particle_contact_container_iterators(
    cleared_local_adjacent_particles, local_particle_container);

  update_ghost_particle_particle_contact_container_iterators(
    cleared_ghost_adjacent_particles, local_particle_container);
}

template <int dim>
void
reinitialize_force(Particles::ParticleHandler<dim> &particle_handler,
                   std::vector<Tensor<1, 3>> &      torque,
                   std::vector<Tensor<1, 3>> &      force)
{
  torque.resize(particle_handler.n_locally_owned_particles());
  force.resize(particle_handler.n_locally_owned_particles());

  for (unsigned int i = 0; i < torque.size(); ++i)
    {
      // Reinitializing forces and torques of particles in the system
      force[i][0] = 0;
      force[i][1] = 0;
      force[i][2] = 0;

      torque[i][0] = 0;
      torque[i][1] = 0;
      torque[i][2] = 0;
    }
}

template <int dim>
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
    .rolling_friction_coefficient_particle[0]                       = 0.1;
  dem_parameters.lagrangian_physical_properties.density_particle[0] = 2500;
  dem_parameters.model_parameters.rolling_resistance_method = Parameters::
    Lagrangian::ModelParameters::RollingResistanceMethod::constant_resistance;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties());

  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    local_adjacent_particles;
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    ghost_adjacent_particles;
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    cleared_local_adjacent_particles;
  std::unordered_map<
    unsigned int,
    std::unordered_map<unsigned int, particle_particle_contact_info_struct<2>>>
    cleared_ghost_adjacent_particles;
  std::unordered_map<unsigned int, Particles::ParticleIterator<2>>
    local_particle_container;
  std::unordered_map<unsigned int, Particles::ParticleIterator<2>>
    ghost_particle_container;

  // Finding cell neighbors
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    local_neighbor_list;
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    ghost_neighbor_list;

  FindCellNeighbors<dim> cell_neighbor_object;
  cell_neighbor_object.find_cell_neighbors(triangulation,
                                           local_neighbor_list,
                                           ghost_neighbor_list);

  // Creating broad search, fine search and particle-particle force objects
  ParticleParticleBroadSearch<dim>              broad_search_object;
  ParticleParticleFineSearch<dim>               fine_search_object;
  ParticleParticleHertzMindlinLimitOverlap<dim> nonlinear_force_object(
    dem_parameters);
  VelocityVerletIntegrator<dim> integrator_object;

  MPI_Comm communicator     = triangulation.get_communicator();
  auto     this_mpi_process = Utilities::MPI::this_mpi_process(communicator);

  // Inserting two particles in contact
  Point<2>                 position1 = {0, 0.003};
  int                      id1       = 0;
  Point<2>                 position2 = {0, -0.003};
  int                      id2       = 1;
  Particles::Particle<dim> particle1(position1, position1, id1);

  // Particle 1 is inserted in a cell owned by process1
  if (this_mpi_process == 1)
    {
      typename Triangulation<dim>::active_cell_iterator cell1 =
        GridTools::find_active_cell_around_point(triangulation,
                                                 particle1.get_location());
      Particles::ParticleIterator<dim> pit1 =
        particle_handler.insert_particle(particle1, cell1);
      pit1->get_properties()[DEM::PropertiesIndex::type]    = 0;
      pit1->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
      pit1->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
      pit1->get_properties()[DEM::PropertiesIndex::v_y]     = -0.5;
      pit1->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
      pit1->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
      pit1->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
      pit1->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
      pit1->get_properties()[DEM::PropertiesIndex::mass]    = 1;
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
      pit2->get_properties()[DEM::PropertiesIndex::type]    = 0;
      pit2->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
      pit2->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
      pit2->get_properties()[DEM::PropertiesIndex::v_y]     = 0.5;
      pit2->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
      pit2->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
      pit2->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
      pit2->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
      pit2->get_properties()[DEM::PropertiesIndex::mass]    = 1;
    }

  // Defining variables
  std::unordered_map<unsigned int, std::vector<unsigned int>>
    local_contact_pair_candidates;
  std::unordered_map<unsigned int, std::vector<unsigned int>>
    ghost_contact_pair_candidates;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  {
    unsigned int max_particle_id = 0;
    for (const auto &particle : particle_handler)
      max_particle_id = std::max(max_particle_id, particle.get_id());
    force.resize(max_particle_id + 1);
  }
#else
  force.resize(particle_handler.get_max_local_particle_index());
#endif
  torque.resize(force.size());
  MOI.resize(force.size());
  for (unsigned i = 0; i < MOI.size(); ++i)
    MOI[i] = 1;

  for (unsigned int iteration = 0; iteration < step_end; ++iteration)
    {
      // Reinitializing forces
      reinitialize_force(particle_handler, torque, force);

      particle_handler.exchange_ghost_particles();

      locate_local_particles_in_cells(particle_handler,
                                      local_particle_container,
                                      ghost_particle_container,
                                      cleared_local_adjacent_particles,
                                      cleared_ghost_adjacent_particles);

      // Calling broad search
      broad_search_object.find_particle_particle_contact_pairs(
        particle_handler,
        local_neighbor_list,
        ghost_neighbor_list,
        local_contact_pair_candidates,
        ghost_contact_pair_candidates);

      // Calling fine search
      fine_search_object.particle_particle_fine_search(
        local_contact_pair_candidates,
        ghost_contact_pair_candidates,
        cleared_local_adjacent_particles,
        cleared_ghost_adjacent_particles,
        local_particle_container,
        neighborhood_threshold);

      // Integration
      // Calling non-linear force
      nonlinear_force_object.calculate_particle_particle_contact_force(
        cleared_local_adjacent_particles,
        cleared_ghost_adjacent_particles,
        dt,
        torque,
        force);

      // Integration
      integrator_object.integrate(particle_handler, g, dt, torque, force, MOI);

      update_contact_containers(local_adjacent_particles,
                                ghost_adjacent_particles,
                                cleared_local_adjacent_particles,
                                cleared_ghost_adjacent_particles);

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
      test<2>();
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
