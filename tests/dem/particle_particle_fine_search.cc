// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the performance of particle-particle fine search class
 * is evaluated.
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
#include <dem/contact_info.h>
#include <dem/dem_contact_manager.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_fine_search.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

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
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties<DEM::SolverType::cfd_dem>());

  DEMContactManager<dim> contact_manager;

  // Finding cell neighbors
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // Inserting two particles in contact
  Point<3> position1 = {0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.40499, 0, 0};
  int      id2       = 1;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());
  double       particle_diameter      = 0.005;
  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);

  pit1->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::type] =
    0;
  pit1->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::dp] =
    particle_diameter;
  pit1->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::v_x] =
    0;
  pit1->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::v_y] =
    0;
  pit1->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::v_z] =
    0;
  pit1->get_properties()
    [DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::omega_x] = 0;
  pit1->get_properties()
    [DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::omega_y] = 0;
  pit1->get_properties()
    [DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::omega_z] = 0;
  pit1->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::mass] =
    1;

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);
  pit2->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::type] =
    0;
  pit2->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::dp] =
    particle_diameter;
  pit2->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::v_x] =
    0;
  pit2->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::v_y] =
    0;
  pit2->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::v_z] =
    0;
  pit2->get_properties()
    [DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::omega_x] = 0;
  pit2->get_properties()
    [DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::omega_y] = 0;
  pit2->get_properties()
    [DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::omega_z] = 0;
  pit2->get_properties()[DEM::PropertiesIndex<DEM::SolverType::cfd_dem>::mass] =
    1;

  contact_manager.update_local_particles_in_cells(particle_handler);

  // Dummy Adaptive sparse contacts object and particle-particle broad search
  AdaptiveSparseContacts<dim> dummy_adaptive_sparse_contacts;
  contact_manager.execute_particle_particle_broad_search(
    particle_handler, dummy_adaptive_sparse_contacts);

  // Calling fine search
  contact_manager.execute_particle_particle_fine_search(neighborhood_threshold);

  // Output
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_adjacent_particles = contact_manager.get_local_adjacent_particles();

  for (auto adjacent_particles_iterator = local_adjacent_particles.begin();
       adjacent_particles_iterator != local_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      auto information_iterator =
        (&adjacent_particles_iterator->second)->begin();
      while (information_iterator !=
             (&adjacent_particles_iterator->second)->end())
        {
          deallog << "The particle pair in contact are particles: "
                  << information_iterator->second.particle_one->get_id()
                  << " and "
                  << information_iterator->second.particle_two->get_id()
                  << std::endl;
          deallog << "Tangential overlap at the beginning of contact is: "
                  << information_iterator->second.tangential_overlap[0] << " "
                  << information_iterator->second.tangential_overlap[1] << " "
                  << information_iterator->second.tangential_overlap[2]
                  << std::endl;
          ++information_iterator;
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
