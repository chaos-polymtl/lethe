// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Inserts three batches of particles with the volume insertion and
 * checks that the random offsets differ between batches. The insertion object
 * is also checkpointed after the second batch and restored in a new object,
 * which must insert the same third batch as the original object.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

// Lethe
#include <dem/dem_solver_parameters.h>
#include <dem/insertion_volume.h>

// Tests (with common definitions)
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <../tests/tests.h>

#include <algorithm>
#include <sstream>
#include <vector>

using namespace dealii;

template <int dim, typename PropertiesIndex>
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

  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  InsertionInfo<dim>           &insert_info = dem_parameters.insertion_info;
  LagrangianPhysicalProperties &lpp =
    dem_parameters.lagrangian_physical_properties;

  // Defining simulation general parameters
  // Insertion info
  insert_info.insertion_box_point_1    = {-0.05, -0.05, -0.05};
  insert_info.insertion_box_point_2    = {0.05, 0.05, 0.05};
  insert_info.direction_sequence       = {0, 1, 2};
  insert_info.inserted_this_step       = 10;
  insert_info.distance_threshold       = 2;
  insert_info.insertion_maximum_offset = 0.75;
  insert_info.seed_for_insertion       = 19;
  insert_info.insertion_acceptance_fct =
    std::make_shared<Functions::ConstantFunction<dim>>(1.);

  // Lagrangian physical properties
  lpp.particle_type_number = 1;
  lpp.particle_average_diameter.push_back(0.005);
  lpp.distribution_type.push_back(SizeDistributionType::uniform);
  lpp.density_particle.push_back(2500);
  lpp.number.push_back(30);

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);
  // Calling uniform insertion
  std::vector<std::shared_ptr<Distribution>> distribution_object_container;
  distribution_object_container.push_back(std::make_shared<UniformDistribution>(
    dem_parameters.lagrangian_physical_properties
      .particle_average_diameter[0]));

  // Calling volume insertion
  InsertionVolume<dim, PropertiesIndex> insertion_object(
    distribution_object_container,
    tr,
    dem_parameters,
    distribution_object_container[0]->find_max_diameter());

  // Insert two batches of 10 particles
  insertion_object.insert(particle_handler, tr, dem_parameters);
  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Checkpoint the insertion object after the second batch. The archive is
  // scoped so that it is complete before the checkpoint is read back.
  std::stringstream checkpoint;
  {
    boost::archive::text_oarchive oa(checkpoint, boost::archive::no_header);
    insertion_object.serialize(oa);
  }

  // Insert the third batch with the original insertion object
  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Restore the checkpoint in a new insertion object and insert the third
  // batch with it, in a separate particle handler
  InsertionVolume<dim, PropertiesIndex> restored_insertion_object(
    distribution_object_container,
    tr,
    dem_parameters,
    distribution_object_container[0]->find_max_diameter());
  boost::archive::text_iarchive ia(checkpoint, boost::archive::no_header);
  restored_insertion_object.deserialize(ia);

  Particles::ParticleHandler<dim> restored_particle_handler(
    tr, mapping, PropertiesIndex::n_properties);
  restored_insertion_object.insert(restored_particle_handler,
                                   tr,
                                   dem_parameters);

  // Sort the particles of a particle handler by id, so that they follow the
  // insertion order
  auto sorted_particles = [](const Particles::ParticleHandler<dim> &handler) {
    std::vector<std::pair<types::particle_index, Point<dim>>> particles;
    for (const auto &particle : handler)
      particles.emplace_back(particle.get_id(), particle.get_location());
    std::sort(particles.begin(),
              particles.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    return particles;
  };

  const auto particles          = sorted_particles(particle_handler);
  const auto restored_particles = sorted_particles(restored_particle_handler);

  // Output
  for (const auto &[id, location] : particles)
    {
      deallog << "Batch " << id / insert_info.inserted_this_step + 1
              << ", particle " << id << " is inserted at: " << location[0]
              << " " << location[1] << " " << location[2] << std::endl;
    }

  // Check if the batch inserted by the restored object has the same positions
  // as the batch starting at first_index in the original particle handler. The
  // ids differ, since the restored particle handler starts from id 0.
  const double tolerance     = 1e-12;
  auto         matches_batch = [&](const unsigned int first_index) {
    if (restored_particles.size() != insert_info.inserted_this_step)
      return false;
    for (unsigned int i = 0; i < restored_particles.size(); ++i)
      if ((restored_particles[i].second - particles[first_index + i].second)
            .norm() > tolerance)
        return false;
    return true;
  };

  // The restored insertion counter must continue the sequence of the original
  // object (third batch), and must not restart it (first batch).
  deallog << "Restored third batch matches the uninterrupted third batch: "
          << (matches_batch(2 * insert_info.inserted_this_step) ? "true" :
                                                                  "false")
          << std::endl;
  deallog << "Restored third batch matches the first batch: "
          << (matches_batch(0) ? "true" : "false") << std::endl;
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
