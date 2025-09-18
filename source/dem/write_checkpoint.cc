// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/checkpoint_control.h>

#include <dem/write_checkpoint.h>

#include <boost/archive/text_oarchive.hpp>

#include <fstream>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
write_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim, PropertiesIndex>>        &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_objects,
  const ConditionalOStream                                &pcout,
  MPI_Comm                                                &mpi_communicator,
  const CheckpointControl &checkpoint_controller)
{
  TimerOutput::Scope timer(computing_timer, "Write checkpoint");

  pcout << "Writing restart file" << std::endl;

  std::string prefix =
    checkpoint_controller.get_filename() + "_" +
    Utilities::int_to_string(checkpoint_controller.get_next_checkpoint_id());

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      simulation_control->save(prefix);
      particles_pvdhandler.save(prefix);

      if (parameters.post_processing.lagrangian_post_processing_enabled)
        {
          grid_pvdhandler.save(prefix + "_lagrangian_postprocessing");
        }
    }

  // Prepare the particle handler for checkpointing
  particle_handler.prepare_for_serialization();

  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  oa << particle_handler;
  std::string triangulation_name = prefix + ".triangulation";
  triangulation.save(prefix + ".triangulation");

  // Write additional particle information for deserialization
  std::string   particle_filename = prefix + ".particles";
  std::ofstream output(particle_filename.c_str());
  output << oss.str() << std::endl;

  // Prepare the insertion object for checkpointing
  std::string   insertion_object_filename = prefix + ".insertion_object";
  std::ofstream oss_insertion_obj(insertion_object_filename);
  boost::archive::text_oarchive oa_insertion_obj(oss_insertion_obj,
                                                 boost::archive::no_header);
  insertion_object->serialize(oa_insertion_obj);

  // Checkpoint the serial solid objects one by one
  for (unsigned int i = 0; i < solid_objects.size(); ++i)
    {
      solid_objects[i]->write_checkpoint(prefix);
    }

  // Prepare the checkpoint controller for checkpointing
  // We don't use the same prefix, since this file needs to have the same name
  // regardless of the checkpoint id being used. This file is giving the
  // information of which checkpoint id to use when restarting.
  std::string checkpoint_controller_object_filename =
    checkpoint_controller.get_filename() + ".checkpoint_controller";
  std::ofstream oss_checkpoint_controller_obj(
    checkpoint_controller_object_filename);
  boost::archive::text_oarchive oa_checkpoint_controller_obj(
    oss_checkpoint_controller_obj, boost::archive::no_header);
  checkpoint_controller.serialize(oa_checkpoint_controller_obj);
}

template void
write_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<2>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<2> &triangulation,
  Particles::ParticleHandler<2>           &particle_handler,
  std::shared_ptr<Insertion<2, DEM::DEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_objects,
  const ConditionalOStream                        &pcout,
  MPI_Comm                                        &mpi_communicator,
  const CheckpointControl                         &checkpoint_controller);

template void
write_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<3>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<3> &triangulation,
  Particles::ParticleHandler<3>           &particle_handler,
  std::shared_ptr<Insertion<3, DEM::DEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_objects,
  const ConditionalOStream                        &pcout,
  MPI_Comm                                        &mpi_communicator,
  const CheckpointControl                         &checkpoint_controller);

template void
write_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<2>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<2> &triangulation,
  Particles::ParticleHandler<2>           &particle_handler,
  std::shared_ptr<Insertion<2, DEM::CFDDEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_objects,
  const ConditionalOStream                        &pcout,
  MPI_Comm                                        &mpi_communicator,
  const CheckpointControl                         &checkpoint_controller);

template void
write_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<3>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<3> &triangulation,
  Particles::ParticleHandler<3>           &particle_handler,
  std::shared_ptr<Insertion<3, DEM::CFDDEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_objects,
  const ConditionalOStream                        &pcout,
  MPI_Comm                                        &mpi_communicator,
  const CheckpointControl                         &checkpoint_controller);

template void
write_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<2>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<2> &triangulation,
  Particles::ParticleHandler<2>           &particle_handler,
  std::shared_ptr<Insertion<2, DEM::DEMMPProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_objects,
  const ConditionalOStream                        &pcout,
  MPI_Comm                                        &mpi_communicator,
  const CheckpointControl                         &checkpoint_controller);

template void
write_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<3>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<3> &triangulation,
  Particles::ParticleHandler<3>           &particle_handler,
  std::shared_ptr<Insertion<3, DEM::DEMMPProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_objects,
  const ConditionalOStream                        &pcout,
  MPI_Comm                                        &mpi_communicator,
  const CheckpointControl                         &checkpoint_controller);
