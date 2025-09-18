// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Three particles are inserted manually in the x direction.
 * We check if these particles appear correctly to each other as potential
 * neighbors.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>


// Lethe
#include <dem/dem_contact_manager.h>
#include <dem/find_cell_neighbors.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Generate a cube triangulation and refine it twice globally
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);

  MappingQ1<dim> mapping;

  DEMContactManager<dim, DEM::DEMProperties::PropertiesIndex> contact_manager;
  Particles::ParticleHandler<dim> particle_handler(triangulation, mapping);

  // Finding cell neighbors list, it is required for finding the broad search
  // pairs in the contact_manager
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // inserting three particles at x = -0.4 , x = 0.4 and x = 0.8
  // which means they are inserted in three adjacent cells in x direction
  Point<3> position1 = {-0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.4, 0, 0};
  int      id2       = 1;
  Point<3> position3 = {0.8, 0, 0};
  int      id3       = 2;

  // Manually insert the three particles
  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
    pt1_info = GridTools::find_active_cell_around_point(mapping,
                                                        triangulation,
                                                        position1);
  Particles::Particle<dim>         particle1(position1, pt1_info.second, id1);
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, pt1_info.first);

  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
    pt2_info = GridTools::find_active_cell_around_point(mapping,
                                                        triangulation,
                                                        position2);
  Particles::Particle<dim>         particle2(position2, pt2_info.second, id2);
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, pt2_info.first);

  std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>
    pt3_info = GridTools::find_active_cell_around_point(mapping,
                                                        triangulation,
                                                        position3);
  Particles::Particle<dim>         particle3(position3, pt3_info.second, id3);
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, pt3_info.first);

  // Dummy Adaptive sparse contacts object for next call
  AdaptiveSparseContacts<dim, DEM::DEMProperties::PropertiesIndex>
    dummy_adaptive_sparse_contacts;

  // Calling broad search function
  contact_manager.execute_particle_particle_broad_search(
    particle_handler, dummy_adaptive_sparse_contacts);

  // Output
  typename dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_candidates =
      contact_manager.get_local_contact_pair_candidates();

  for (auto pairs_iterator = local_contact_pair_candidates.begin();
       pairs_iterator != local_contact_pair_candidates.end();
       ++pairs_iterator)
    {
      unsigned int first_particle_id = pairs_iterator->first;
      auto         candidates        = &pairs_iterator->second;
      for (auto candidate_iterator = candidates->begin();
           candidate_iterator != candidates->end();
           ++candidate_iterator)
        if (first_particle_id == 0)
          deallog << "A pair is detected: particle " << first_particle_id
                  << " and particle " << *candidate_iterator << std::endl;
    }

  for (auto pairs_iterator = local_contact_pair_candidates.begin();
       pairs_iterator != local_contact_pair_candidates.end();
       ++pairs_iterator)
    {
      unsigned int first_particle_id = pairs_iterator->first;
      auto         candidates        = &pairs_iterator->second;
      for (auto candidate_iterator = candidates->begin();
           candidate_iterator != candidates->end();
           ++candidate_iterator)
        if (first_particle_id == 1)
          deallog << "A pair is detected: particle " << first_particle_id
                  << " and particle " << *candidate_iterator << std::endl;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, dealii::numbers::invalid_unsigned_int);
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
