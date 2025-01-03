// SPDX-FileCopyrightText: Copyright (c) 2020-2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the performance of particle-wall fine search is
 * investigated.
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
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_fine_search.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
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

  // Defining general simulation parameters
  double particle_diameter = 0.005;

  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties<DEM::SolverType::cfd_dem>());

  // Inserting one particle in contact with a wall
  Point<dim>               position1 = {-0.998, 0, 0};
  int                      id1       = 0;
  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
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

  // Calling find_boundary_cells_information function to find the information of
  // boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::vector<unsigned int>     outlet_boundaries;
  boundary_cells_object.build(
    tr,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Calling particle-wall broad search
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    particle_wall_contact_list;
  find_particle_wall_contact_pairs<dim>(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);

  // Calling particle-wall fine search
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_contact_information;
  particle_wall_fine_search<dim>(particle_wall_contact_list,
                                 particle_wall_contact_information);

  // Output
  for (auto particle_wall_contact_information_iterator =
         particle_wall_contact_information.begin();
       particle_wall_contact_information_iterator !=
       particle_wall_contact_information.end();
       ++particle_wall_contact_information_iterator)
    {
      auto info_iterator =
        particle_wall_contact_information_iterator->second.begin();
      while (info_iterator !=
             particle_wall_contact_information_iterator->second.end())
        {
          deallog << "Particle " << info_iterator->second.particle->get_id()
                  << " is in contact with face " << info_iterator->first
                  << std::endl;
          deallog << "The normal vector of collision is: "
                  << info_iterator->second.normal_vector[0] << " "
                  << info_iterator->second.normal_vector[1] << " "
                  << info_iterator->second.normal_vector[2] << std::endl;
          ++info_iterator;
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
