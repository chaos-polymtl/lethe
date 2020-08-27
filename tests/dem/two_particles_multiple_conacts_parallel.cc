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

// In this test, the performance of non-linear (Hertzian) particle-particle
// contact force  is checked

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/dem_solver_parameters.h>
#include <dem/find_cell_neighbors.h>
#include <dem/pp_broad_search.h>
#include <dem/pp_fine_search.h>
#include <dem/pp_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

void
update_contact_containers(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &ghost_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &cleared_local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &cleared_ghost_adjacent_particles)
{
  local_adjacent_particles.clear();
  ghost_adjacent_particles.clear();

  local_adjacent_particles = cleared_local_adjacent_particles;
  ghost_adjacent_particles = cleared_ghost_adjacent_particles;
}

template <int dim>
void
update_ghost_pp_contact_container_iterators(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    &cleared_ghost_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &local_particle_container)
{
  for (auto adjacent_particles_iterator =
         cleared_ghost_adjacent_particles.begin();
       adjacent_particles_iterator != cleared_ghost_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id          = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();
           ++pp_map_iterator)
        {
          int particle_two_id = pp_map_iterator->first;
          pp_map_iterator->second.particle_one =
            local_particle_container.at(particle_one_id);
          pp_map_iterator->second.particle_two =
            local_particle_container.at(particle_two_id);
        }
    }
}

template <int dim>
void
update_local_pp_contact_container_iterators(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    &cleared_local_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &local_particle_container)
{
  for (auto adjacent_particles_iterator =
         cleared_local_adjacent_particles.begin();
       adjacent_particles_iterator != cleared_local_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id          = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();
           ++pp_map_iterator)
        {
          int particle_two_id = pp_map_iterator->first;

          pp_map_iterator->second.particle_one =
            local_particle_container.at(particle_one_id);
          pp_map_iterator->second.particle_two =
            local_particle_container.at(particle_two_id);
        }
    }
}

template <int dim>
void
update_local_particle_container(
  std::unordered_map<int, Particles::ParticleIterator<dim>>
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
  std::unordered_map<int, Particles::ParticleIterator<2>>
    &local_particle_container,
  std::unordered_map<int, Particles::ParticleIterator<2>>
    &ghost_particle_container,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &cleared_local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    &cleared_ghost_adjacent_particles)
{
  local_particle_container.clear();
  ghost_particle_container.clear();

  update_local_particle_container(local_particle_container, &particle_handler);

  update_local_pp_contact_container_iterators(cleared_local_adjacent_particles,
                                              local_particle_container);

  update_ghost_pp_contact_container_iterators(cleared_ghost_adjacent_particles,
                                              local_particle_container);
}

template <int dim>
void
reinitialize_force(Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Getting properties of particle as local variable
      auto particle_properties = particle->get_properties();

      // Reinitializing forces and momentums of particles in the system
      particle_properties[DEM::PropertiesIndex::force_x] = 0;
      particle_properties[DEM::PropertiesIndex::force_y] = 0;

      particle_properties[DEM::PropertiesIndex::M_x] = 0;
      particle_properties[DEM::PropertiesIndex::M_y] = 0;

      if (dim == 3)
        {
          particle_properties[DEM::PropertiesIndex::force_z] = 0;
          particle_properties[DEM::PropertiesIndex::M_z]     = 0;
        }
    }
}

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  double                                    hyper_cube_length = 0.05;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining general simulation parameters
  const unsigned int n_properties = 21;
  Tensor<1, dim>     g{{0, 0}};
  double             dt                = 0.00001;
  double             particle_diameter = 0.005;
  int                particle_density  = 2500;
  unsigned int       step_end          = 1000;
  unsigned int       output_frequency  = 10;

  dem_parameters.physical_properties.Youngs_modulus_particle = 50000000;
  dem_parameters.physical_properties.Poisson_ratio_particle  = 0.9;
  dem_parameters.physical_properties.restitution_coefficient_particle = 0.9;
  dem_parameters.physical_properties.friction_coefficient_particle    = 0.5;
  dem_parameters.physical_properties.rolling_friction_particle        = 0.1;
  double neighborhood_threshold = 1.3 * particle_diameter;

  Particles::ParticleHandler<dim> particle_handler(triangulation,
                                                   mapping,
                                                   n_properties);

  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    local_adjacent_particles;
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    ghost_adjacent_particles;
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    cleared_local_adjacent_particles;
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    cleared_ghost_adjacent_particles;
  std::unordered_map<int, Particles::ParticleIterator<2>>
    local_particle_container;
  std::unordered_map<int, Particles::ParticleIterator<2>>
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
  PPBroadSearch<dim>            broad_search_object;
  PPFineSearch<dim>             fine_search_object;
  PPNonLinearForce<dim>         nonlinear_force_object;
  VelocityVerletIntegrator<dim> integrator_object;

  // Inserting two particles in contact
  Point<2>                 position1 = {0, 0.007};
  int                      id1       = 0;
  Point<2>                 position2 = {0, 0.001};
  int                      id2       = 1;
  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);
  pit1->get_properties()[0]  = id1;
  pit1->get_properties()[1]  = 1;
  pit1->get_properties()[2]  = particle_diameter;
  pit1->get_properties()[3]  = particle_density;
  pit1->get_properties()[4]  = 0;
  pit1->get_properties()[5]  = -0.4;
  pit1->get_properties()[6]  = 0;
  pit1->get_properties()[7]  = 0;
  pit1->get_properties()[8]  = 0;
  pit1->get_properties()[9]  = 0;
  pit1->get_properties()[10] = 0;
  pit1->get_properties()[11] = 0;
  pit1->get_properties()[12] = 0;
  pit1->get_properties()[13] = 0;
  pit1->get_properties()[14] = 0;
  pit1->get_properties()[15] = 0;
  pit1->get_properties()[16] = 1;
  pit1->get_properties()[17] = 1;

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);
  pit2->get_properties()[0]  = id2;
  pit2->get_properties()[1]  = 1;
  pit2->get_properties()[2]  = particle_diameter;
  pit2->get_properties()[3]  = particle_density;
  pit2->get_properties()[4]  = 0;
  pit2->get_properties()[5]  = 0;
  pit2->get_properties()[6]  = 0;
  pit2->get_properties()[7]  = 0;
  pit2->get_properties()[8]  = 0;
  pit2->get_properties()[9]  = 0;
  pit2->get_properties()[10] = 0;
  pit2->get_properties()[11] = 0;
  pit2->get_properties()[12] = 0;
  pit2->get_properties()[13] = 0;
  pit2->get_properties()[14] = 0;
  pit2->get_properties()[15] = 0;
  pit2->get_properties()[16] = 1;
  pit2->get_properties()[17] = 1;

  // Defining variables
  std::unordered_map<int, std::vector<int>> local_contact_pair_candidates;
  std::unordered_map<int, std::vector<int>> ghost_contact_pair_candidates;

  for (unsigned int iteration = 0; iteration < step_end; ++iteration)
    {
      // Reinitializing forces
      reinitialize_force(particle_handler);

      particle_handler.exchange_ghost_particles();

      locate_local_particles_in_cells(particle_handler,
                                      local_particle_container,
                                      ghost_particle_container,
                                      cleared_local_adjacent_particles,
                                      cleared_ghost_adjacent_particles);

      // Calling broad search
      broad_search_object.find_PP_Contact_Pairs(particle_handler,
                                                &local_neighbor_list,
                                                &ghost_neighbor_list,
                                                local_contact_pair_candidates,
                                                ghost_contact_pair_candidates);

      // Calling fine search
      fine_search_object.pp_Fine_Search(local_contact_pair_candidates,
                                        ghost_contact_pair_candidates,
                                        cleared_local_adjacent_particles,
                                        cleared_ghost_adjacent_particles,
                                        local_particle_container,
                                        neighborhood_threshold);

      // Calling non-linear force
      nonlinear_force_object.calculate_pp_contact_force(
        &cleared_local_adjacent_particles,
        &cleared_ghost_adjacent_particles,
        dem_parameters,
        dt);

      // Integration
      integrator_object.integrate(particle_handler, g, dt);

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
                        << "The exerted force on particle "
                        << particle->get_id() << " at step " << iteration
                        << " is: "
                        << particle
                             ->get_properties()[DEM::PropertiesIndex::force_y]
                        << std::endl;
                    }
                }
            }
        }
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<2>();
}
