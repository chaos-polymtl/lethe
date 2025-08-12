// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/log_collision_data.h>

#include <deal.II/base/mpi.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
log_collision_data(
  const DEMSolverParameters<dim> &parameters,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                             &particle_wall_pairs_in_contact,
  const double                current_time,
  OngoingCollisionLog<dim>   &ongoing_collision_log,
  CompletedCollisionLog<dim> &collision_event_log)
{
  // Looping over all the active particles in particle-wall pairs
  for (auto &&pairs_in_contact_content :
       particle_wall_pairs_in_contact | boost::adaptors::map_values)
    {
      // Iterating over a map which contains the required information for
      // calculation of the contact force for each particle
      for (auto &&contact_info :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          // Defining local variables which will be used within the contact
          // calculation
          auto     particle            = contact_info.particle;
          auto     particle_properties = particle->get_properties();
          Point<3> point_on_boundary   = contact_info.point_on_boundary;
          // Normal vector from the wall to the particle
          Tensor<1, 3> normal_vector = contact_info.normal_vector;
          // Getting particle 3d location
          Point<3> particle_location_3d = get_location(particle);
          // Defining a tensor which connects the point_on_boundary to the
          // center of particle
          Tensor<1, 3> point_to_particle_vector =
            particle_location_3d - point_on_boundary;
          Tensor<1, 3> projected_vector =
            find_projection(point_to_particle_vector, normal_vector);
          // Calculating the particle-wall distance using the projected vector
          double normal_overlap =
            ((particle_properties[PropertiesIndex::dp]) * 0.5) -
            (projected_vector.norm());
          types::boundary_id boundary_id = contact_info.boundary_id;

          // If we log all walls or if the boundary ID is in the list of
          // particle-wall collision boundary IDs, we log the collision
          // information.
          if (parameters.post_processing.log_collisions_with_all_walls ||
              std::find(parameters.post_processing
                          .particle_wall_collision_boundary_ids.begin(),
                        parameters.post_processing
                          .particle_wall_collision_boundary_ids.end(),
                        boundary_id) !=
                parameters.post_processing.particle_wall_collision_boundary_ids
                  .end())
            {
              unsigned int particle_id = particle->get_id();
              if (normal_overlap > 0 &&
                  !ongoing_collision_log.is_in_collision(particle_id,
                                                         boundary_id))
                {
                  collision_log<dim> start_log;
                  start_log.particle_id = particle_id;
                  start_log.dp   = particle_properties[PropertiesIndex::dp];
                  start_log.mass = particle_properties[PropertiesIndex::mass];
                  start_log.velocity[0] =
                    particle_properties[PropertiesIndex::v_x];
                  start_log.velocity[1] =
                    particle_properties[PropertiesIndex::v_y];
                  start_log.velocity[2] =
                    particle_properties[PropertiesIndex::v_z];
                  start_log.omega[0] =
                    particle_properties[PropertiesIndex::omega_x];
                  start_log.omega[1] =
                    particle_properties[PropertiesIndex::omega_y];
                  start_log.omega[2] =
                    particle_properties[PropertiesIndex::omega_z];
                  start_log.time        = current_time;
                  start_log.boundary_id = contact_info.boundary_id;
                  ongoing_collision_log.start_collision(
                    start_log); // Start logging the collision

                  // Print the start of the collision in the terminal if
                  // verbosity is set to verbose
                  if (parameters.post_processing.collision_verbosity ==
                      Parameters::Verbosity::verbose)
                    {
                      std::cout
                        << "Collision with boundary " << start_log.boundary_id
                        << " started for particle " << particle_id << std::endl;
                    }
                }

              // If the particle does not have a positive overlap anymore, the
              // collision has ended. Consequently, it now needs to be
              // removed from ongoing_collision_log since the contact has
              // reached its end.
              if (normal_overlap < 0 &&
                  ongoing_collision_log.is_in_collision(particle_id,
                                                        boundary_id))
                {
                  collision_log<dim> end_log;
                  end_log.particle_id = particle_id;
                  end_log.dp   = particle_properties[PropertiesIndex::dp];
                  end_log.mass = particle_properties[PropertiesIndex::mass];
                  end_log.velocity[0] =
                    particle_properties[PropertiesIndex::v_x];
                  end_log.velocity[1] =
                    particle_properties[PropertiesIndex::v_y];
                  end_log.velocity[2] =
                    particle_properties[PropertiesIndex::v_z];
                  end_log.omega[0] =
                    particle_properties[PropertiesIndex::omega_x];
                  end_log.omega[1] =
                    particle_properties[PropertiesIndex::omega_y];
                  end_log.omega[2] =
                    particle_properties[PropertiesIndex::omega_z];
                  end_log.time = current_time;
                  collision_log<dim> start_log;

                  // End the collision for the particle and retrieve the start
                  // log
                  ongoing_collision_log.end_collision(particle_id,
                                                      boundary_id,
                                                      start_log);
                  end_log.boundary_id = start_log.boundary_id;
                  collision_event<dim> event;
                  event.particle_id = particle_id;
                  event.start_log   = start_log;
                  event.end_log     = end_log;

                  // Add the completed collision event to the collision event
                  // log
                  collision_event_log.add_event(event);

                  // Print the end of the collision in the terminal if verbosity
                  // is set to verbose
                  if (parameters.post_processing.collision_verbosity ==
                      Parameters::Verbosity::verbose)
                    {
                      std::cout << "Collision with boundary "
                                << end_log.boundary_id << " ended for particle "
                                << particle_id << std::endl;
                    }
                }
            }
        }
    }
}

template <int dim>
void
write_collision_stats(const DEMSolverParameters<dim>   &parameters,
                      const CompletedCollisionLog<dim> &collision_event_log,
                      const MPI_Comm                   &mpi_communicator)
{
  // MPI processes information
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);

  // Separator and filename for the output file
  std::string sep;
  std::string filename = parameters.post_processing.collision_stats_file_name;

  // Check if a .csv or .dat extension is specified in the filename, if not add
  // ".csv"
  std::size_t csv_file = filename.find(".csv");
  std::size_t dat_file = filename.find(".dat");

  if ((csv_file == std::string::npos) && (dat_file == std::string::npos))
    filename += ".csv";

  // Open the file for writing or appending based on the MPI process.
  // This forces a barrier after MPI processes and will undoubtly be slow in
  // large parallel simulations. If this becomes an issue, the function should
  // be ported to use MPI I/O or a format like HDF5.
  for (unsigned int i = 0; i < n_mpi_processes; ++i)
    {
      if (this_mpi_process == i)
        {
          std::ofstream myfile;

          if (this_mpi_process == 0)
            {
              // If this is the first MPI process, we write the header
              myfile.open(filename);
              if (filename.substr(filename.find_last_of('.') + 1) == ".dat")
                {
                  myfile
                    << "particle_id diameter mass boundary_id start_time end_time start_particle_velocity_x start_particle_velocity_y start_particle_velocity_z start_particle_angular_velocity_x start_particle_angular_velocity_y start_particle_angular_velocity_z end_particle_velocity_x end_particle_velocity_y end_particle_velocity_z end_particle_angular_velocity_x end_particle_angular_velocity_y end_particle_angular_velocity_z"
                    << std::endl;
                  sep = " ";
                }
              else // .csv is default
                {
                  myfile
                    << "particle_id,diameter,mass,boundary_id,start_time,end_time,start_particle_velocity_x,start_particle_velocity_y,start_particle_velocity_z,start_particle_angular_velocity_x,start_particle_angular_velocity_y,start_particle_angular_velocity_z,end_particle_velocity_x,end_particle_velocity_y,end_particle_velocity_z,end_particle_angular_velocity_x,end_particle_angular_velocity_y,end_particle_angular_velocity_z"
                    << std::endl;
                  sep = ",";
                }
            }
          else
            {
              // If this is not the first MPI process, we open the file for
              // appending
              myfile.open(filename, std::ios::app);
              if (filename.substr(filename.find_last_of('.') + 1) == ".dat")
                sep = " ";
              else // .csv is default
                sep = ",";
            }

          // Write the collision statistics
          for (const auto &event : collision_event_log.get_events())
            {
              const auto &start = event.start_log;
              const auto &end   = event.end_log;

              // Write the collision data to the file
              myfile << start.particle_id << sep << start.dp << sep
                     << start.mass << sep << static_cast<int>(start.boundary_id)
                     << sep << start.time << sep << end.time << sep
                     << start.velocity[0] << sep << start.velocity[1] << sep
                     << start.velocity[2] << sep << start.omega[0] << sep
                     << start.omega[1] << sep << start.omega[2] << sep
                     << end.velocity[0] << sep << end.velocity[1] << sep
                     << end.velocity[2] << sep << end.omega[0] << sep
                     << end.omega[1] << sep << end.omega[2] << std::endl;
            }
          myfile.close();
        }
      // Ensure all MPI processes reach this point before continuing.
      // The barrier forces the synchronisation of the processes.
      MPI_Barrier(mpi_communicator);
    }
}

template void
log_collision_data<2, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &parameters,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<2>   &ongoing_collision_log,
  CompletedCollisionLog<2> &collision_event_log);

template void
log_collision_data<3, DEM::DEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &parameters,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<3>   &ongoing_collision_log,
  CompletedCollisionLog<3> &collision_event_log);

template void
log_collision_data<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &parameters,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<2>   &ongoing_collision_log,
  CompletedCollisionLog<2> &collision_event_log);

template void
log_collision_data<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &parameters,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<3>   &ongoing_collision_log,
  CompletedCollisionLog<3> &collision_event_log);

template void
log_collision_data<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<2> &parameters,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<2>   &ongoing_collision_log,
  CompletedCollisionLog<2> &collision_event_log);

template void
log_collision_data<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEMSolverParameters<3> &parameters,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
                           &particle_wall_pairs_in_contact,
  const double              current_time,
  OngoingCollisionLog<3>   &ongoing_collision_log,
  CompletedCollisionLog<3> &collision_event_log);

template void
write_collision_stats<2>(const DEMSolverParameters<2>   &parameters,
                         const CompletedCollisionLog<2> &collision_event_log,
                         const MPI_Comm                 &mpi_communicator);

template void
write_collision_stats<3>(const DEMSolverParameters<3>   &parameters,
                         const CompletedCollisionLog<3> &collision_event_log,
                         const MPI_Comm                 &mpi_communicator);
