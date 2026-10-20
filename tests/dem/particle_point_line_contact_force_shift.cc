// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the shift particle-point and particle-line contact
 * model used by the packed insertion method is checked in a hyper-L domain,
 * which has a reentrant vertex at the origin and a reentrant line along the z
 * axis. Particle 0 overlaps the vertex and particle 1 overlaps the line. The
 * contact evaluation must not move the particles, and a single displacement
 * must separate the particles from the boundary.
 */

// Deal.II
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/insertion_packed.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_contact_force.h>
#include <dem/particle_point_line_fine_search.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cmath>
#include <set>
#include <vector>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement. The octant [0, 1]^3 is removed from the
  // cube [-1, 1]^3. The cell [-0.5, 0]^3 only touches the boundary at the
  // origin and the cell [-0.5, 0]^2 x [0, 0.5] only touches the boundary along
  // the z axis.
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_L(triangulation, -1, 1);
  triangulation.refine_global(1);
  MappingQ<dim> mapping(1);

  const double particle_diameter      = 0.1;
  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, PropertiesIndex::n_properties);

  // Inserting two particles. Particle 0 overlaps the boundary vertex at the
  // origin and particle 1 overlaps the boundary line along the z axis.
  const std::vector<Point<dim>> positions = {Point<dim>(-0.02, -0.02, -0.02),
                                             Point<dim>(-0.03, -0.02, 0.25)};

  for (unsigned int id = 0; id < positions.size(); ++id)
    {
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

  // Construct boundary cells object and build it
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::set<types::boundary_id>  outlet_boundaries;
  boundary_cells_object.build(
    triangulation,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Particle-point and particle-line broad and fine searches
  typename DEM::dem_data_structures<dim>::particle_point_candidates
    point_candidates;
  typename DEM::dem_data_structures<dim>::particle_line_candidates
    line_candidates;
  typename DEM::dem_data_structures<dim>::particle_point_in_contact
    points_in_contact;
  typename DEM::dem_data_structures<dim>::particle_line_in_contact
    lines_in_contact;

  find_particle_point_contact_pairs<dim>(
    particle_handler,
    boundary_cells_object.get_boundary_cells_with_points(),
    point_candidates);
  find_particle_line_contact_pairs<dim>(
    particle_handler,
    boundary_cells_object.get_boundary_cells_with_lines(),
    line_candidates);

  particle_point_fine_search<dim, PropertiesIndex>(point_candidates,
                                                   neighborhood_threshold,
                                                   points_in_contact);
  particle_line_fine_search<dim, PropertiesIndex>(line_candidates,
                                                  neighborhood_threshold,
                                                  lines_in_contact);

  deallog << "Number of particle-point pairs: " << points_in_contact.size()
          << std::endl;
  deallog << "Number of particle-line pairs: " << lines_in_contact.size()
          << std::endl;

  // Evaluate the contacts with the shift model
  ParticlePointLineForce<dim, PropertiesIndex> shift_force_object;
  auto evaluate_contacts = [&]() {
    shift_force_object.calculate_particle_point_shift_displacement(
      points_in_contact, contact_outcome.force);
    shift_force_object.calculate_particle_line_shift_displacement(
      lines_in_contact, contact_outcome.force);
  };

  evaluate_contacts();

  for (const auto &particle : particle_handler)
    {
      deallog << "Particle " << particle.get_id()
              << " after the contact evaluation" << std::endl;
      deallog << "  location: " << particle.get_location() << std::endl;
      deallog << "  accumulated displacement: "
              << contact_outcome.force[particle.get_local_index()]
              << std::endl;
    }
  deallog << "Number of contacts: "
          << shift_force_object.get_number_of_contacts() << std::endl;

  // Displace the particles
  InsertionPacked<dim, PropertiesIndex>::apply_shift_displacement(
    particle_handler, particle_diameter, contact_outcome, displacement);

  for (const auto &particle : particle_handler)
    {
      const Point<dim> location = particle.get_location();
      deallog << "Particle " << particle.get_id()
              << " location after the displacement: " << location
              << std::endl;
    }

  // Particle 0 is compared with the vertex at the origin and particle 1 with
  // the line along the z axis
  for (const auto &particle : particle_handler)
    {
      const Point<dim> location = particle.get_location();
      if (particle.get_id() == 0)
        deallog << "Distance between particle 0 and the boundary vertex: "
                << location.norm() << std::endl;
      else
        deallog << "Distance between particle 1 and the boundary line: "
                << std::sqrt(location[0] * location[0] +
                             location[1] * location[1])
                << std::endl;
    }

  // A second evaluation must not find any contact
  shift_force_object.reset_number_of_contacts();
  evaluate_contacts();
  deallog << "Number of contacts after the displacement: "
          << shift_force_object.get_number_of_contacts() << std::endl;
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
