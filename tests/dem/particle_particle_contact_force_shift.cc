// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the shift particle-particle contact model used by the
 * packed insertion method is checked for three particles. Particle 1 overlaps
 * particles 0 and 2, which do not overlap each other. With two processes,
 * particle 0 is owned by process 1 and particles 1 and 2 by process 0, so
 * particle 1 has a local-ghost and a local-local contact. The contact
 * evaluation must not move the particles, particle 1 must be displaced by the
 * sum of the displacements of its two contacts, and a single displacement
 * must separate the particles. Since the displacements do not depend on the
 * order in which the contacts are evaluated nor on the partitioning, the output
 * is the same with one and two processes.
 */

// Deal.II
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>
#include <core/parameters_lagrangian.h>

#include <../tests/dem/test_particles_functions.h>
#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/insertion_packed.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/particle_particle_contact_force.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cmath>
#include <utility>
#include <vector>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement. With two processes, the cells with y < 0
  // are owned by process 0 and the cells with y > 0 by process 1.
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  triangulation.refine_global(2);
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  DEMSolverParameters<dim> dem_parameters;
  set_default_dem_parameters(1, dem_parameters);

  const double dt                     = 0.00001;
  const double particle_diameter      = 0.005;
  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, PropertiesIndex::n_properties);

  DEMContactManager<dim, PropertiesIndex> contact_manager;

  // Finding cell neighbors
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // Creating the shift particle-particle contact model
  ParticleParticleContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleParticleContactForceModel::shift,
    Parameters::Lagrangian::RollingResistanceMethod::none>
    shift_force_object(dem_parameters);

  MPI_Comm           communicator = triangulation.get_mpi_communicator();
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(communicator);

  // Inserting three particles. Particle 1 overlaps particles 0 and 2 with a
  // normal overlap of 0.001, while particles 0 and 2 do not overlap.
  const std::vector<Point<dim>>   positions         = {Point<dim>(0, 0.002),
                                                       Point<dim>(0, -0.002),
                                                       Point<dim>(0.004, -0.002)};
  const std::vector<unsigned int> owner_of_particle = {1, 0, 0};
  const unsigned int              n_mpi_processes =
    Utilities::MPI::n_mpi_processes(communicator);

  for (unsigned int id = 0; id < positions.size(); ++id)
    {
      if (this_mpi_process != owner_of_particle[id] % n_mpi_processes)
        continue;

      Particles::Particle<dim> particle(positions[id], positions[id], id);
      typename Triangulation<dim>::active_cell_iterator cell =
        GridTools::find_active_cell_around_point(triangulation, positions[id]);
      Particles::ParticleIterator<dim> pit =
        particle_handler.insert_particle(particle, cell);
      pit->get_properties()[PropertiesIndex::type]    = 0;
      pit->get_properties()[PropertiesIndex::dp]      = particle_diameter;
      pit->get_properties()[PropertiesIndex::v_x]     = 0;
      pit->get_properties()[PropertiesIndex::v_y]     = 0;
      pit->get_properties()[PropertiesIndex::v_z]     = 0;
      pit->get_properties()[PropertiesIndex::omega_x] = 0;
      pit->get_properties()[PropertiesIndex::omega_y] = 0;
      pit->get_properties()[PropertiesIndex::omega_z] = 0;
      pit->get_properties()[PropertiesIndex::mass]    = 1;
    }

  particle_handler.sort_particles_into_subdomains_and_cells();
  const unsigned int number_of_particles =
    particle_handler.get_max_local_particle_index();

  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  contact_outcome.resize_interaction_containers(number_of_particles);
  std::vector<double> displacement(number_of_particles, 0.);

  // Contact search
  particle_handler.exchange_ghost_particles(true);
  contact_manager.update_local_particles_in_cells(particle_handler);
  AdaptiveSparseContacts<dim, PropertiesIndex> dummy_adaptive_sparse_contacts;
  contact_manager.execute_particle_particle_broad_search(
    particle_handler, dummy_adaptive_sparse_contacts);
  contact_manager.execute_particle_particle_fine_search(neighborhood_threshold);

  // Evaluate the contacts with the shift model
  auto evaluate_contacts = [&]() {
    shift_force_object.calculate_particle_particle_contact(
      contact_manager.get_local_adjacent_particles(),
      contact_manager.get_ghost_adjacent_particles(),
      contact_manager.get_local_local_periodic_adjacent_particles(),
      contact_manager.get_local_ghost_periodic_adjacent_particles(),
      contact_manager.get_ghost_local_periodic_adjacent_particles(),
      dt,
      contact_outcome);
  };

  // Return the location and accumulated displacement of a particle, gathered
  // from the process that owns it
  auto gather_particle_state = [&](const types::particle_index id) {
    Tensor<1, 3> location;
    Tensor<1, 3> accumulated_displacement;
    for (auto &particle : particle_handler)
      {
        if (particle.get_id() != id)
          continue;
        const Point<dim> particle_location = particle.get_location();
        for (unsigned int d = 0; d < dim; ++d)
          location[d] = particle_location[d];
        accumulated_displacement =
          contact_outcome.force[particle.get_local_index()];
      }
    return std::make_pair(Utilities::MPI::sum(location, communicator),
                          Utilities::MPI::sum(accumulated_displacement,
                                              communicator));
  };

  evaluate_contacts();

  for (types::particle_index id = 0; id < positions.size(); ++id)
    {
      const auto [location, accumulated_displacement] =
        gather_particle_state(id);
      deallog << "Particle " << id << " after the contact evaluation"
              << std::endl;
      deallog << "  location: " << location << std::endl;
      deallog << "  accumulated displacement: " << accumulated_displacement
              << std::endl;
    }

  // Displace the particles and update the ghost particles
  InsertionPacked<dim, PropertiesIndex>::apply_shift_displacement(
    particle_handler, particle_diameter, contact_outcome, displacement);
  particle_handler.update_ghost_particles();

  std::vector<Tensor<1, 3>> locations;
  for (types::particle_index id = 0; id < positions.size(); ++id)
    {
      locations.push_back(gather_particle_state(id).first);
      deallog << "Particle " << id
              << " location after the displacement: " << locations[id]
              << std::endl;
    }
  for (types::particle_index i = 0; i < positions.size(); ++i)
    for (types::particle_index j = i + 1; j < positions.size(); ++j)
      deallog << "Distance between particles " << i << " and " << j << ": "
              << (locations[i] - locations[j]).norm() << std::endl;

  // A second evaluation must not find any contact
  shift_force_object.reset_number_of_contacts();
  evaluate_contacts();
  deallog << "Number of contacts after the displacement: "
          << Utilities::MPI::sum(shift_force_object.get_number_of_contacts(),
                                 communicator)
          << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      mpi_initlog();
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
