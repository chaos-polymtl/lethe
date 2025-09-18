// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/checkpoint_control.h>

#include <dem/dem_action_manager.h>
#include <dem/read_checkpoint.h>

#include <boost/archive/text_iarchive.hpp>

#include <fstream>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
read_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim, PropertiesIndex>>        &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_surfaces,
  CheckpointControl &checkpoint_controller)
{
  if (!DEMActionManager::get_action_manager()->check_restart_simulation())
    return;

  TimerOutput::Scope timer(computing_timer, "Read checkpoint");
  std::string        prefix = checkpoint_controller.get_filename();

  // Load checkpoint controller
  std::string checkpoint_controller_object_filename =
    prefix + ".checkpoint_controller";
  std::ifstream iss_checkpoint_controller_obj(
    checkpoint_controller_object_filename);

  assert_restart_file_exists(iss_checkpoint_controller_obj,
                             checkpoint_controller_object_filename);

  boost::archive::text_iarchive ia_checkpoint_controller_obj(
    iss_checkpoint_controller_obj, boost::archive::no_header);
  checkpoint_controller.deserialize(ia_checkpoint_controller_obj);

  // Update prefix with the checkpoint id
  prefix =
    prefix + "_" +
    Utilities::int_to_string(checkpoint_controller.get_next_checkpoint_id());

  simulation_control->read(prefix);
  particles_pvdhandler.read(prefix);

  if (parameters.post_processing.lagrangian_post_processing_enabled)
    {
      grid_pvdhandler.read(prefix + "_lagrangian_postprocessing");
    }

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());

  assert_restart_file_exists(input, particle_filename);

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());

  assert_restart_file_exists(in, filename);

  try
    {
      triangulation.load(filename.c_str());
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }


  // Unpack the information in the particle handler
  particle_handler.deserialize();


  // Load insertion object
  std::string   insertion_object_filename = prefix + ".insertion_object";
  std::ifstream iss_insertion_obj(insertion_object_filename);
  boost::archive::text_iarchive ia_insertion_obj(iss_insertion_obj,
                                                 boost::archive::no_header);
  insertion_object->deserialize(ia_insertion_obj);

  // Load solid surfaces
  for (unsigned int i = 0; i < solid_surfaces.size(); ++i)
    {
      solid_surfaces[i]->read_checkpoint(prefix);
    }
}

template void
read_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<2>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<2> &triangulation,
  Particles::ParticleHandler<2>           &particle_handler,
  std::shared_ptr<Insertion<2, DEM::DEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_surfaces,
  CheckpointControl                               &checkpoint_controller);

template void
read_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<3>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<3> &triangulation,
  Particles::ParticleHandler<3>           &particle_handler,
  std::shared_ptr<Insertion<3, DEM::DEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_surfaces,
  CheckpointControl                               &checkpoint_controller);

template void
read_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<2>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<2> &triangulation,
  Particles::ParticleHandler<2>           &particle_handler,
  std::shared_ptr<Insertion<2, DEM::CFDDEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_surfaces,
  CheckpointControl                               &checkpoint_controller);

template void
read_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<3>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<3> &triangulation,
  Particles::ParticleHandler<3>           &particle_handler,
  std::shared_ptr<Insertion<3, DEM::CFDDEMProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_surfaces,
  CheckpointControl                               &checkpoint_controller);

template void
read_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<2>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<2> &triangulation,
  Particles::ParticleHandler<2>           &particle_handler,
  std::shared_ptr<Insertion<2, DEM::DEMMPProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> &solid_surfaces,
  CheckpointControl                               &checkpoint_controller);

template void
read_checkpoint(
  TimerOutput                             &computing_timer,
  const DEMSolverParameters<3>            &parameters,
  std::shared_ptr<SimulationControl>      &simulation_control,
  PVDHandler                              &particles_pvdhandler,
  PVDHandler                              &grid_pvdhandler,
  parallel::distributed::Triangulation<3> &triangulation,
  Particles::ParticleHandler<3>           &particle_handler,
  std::shared_ptr<Insertion<3, DEM::DEMMPProperties::PropertiesIndex>>
                                                  &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> &solid_surfaces,
  CheckpointControl                               &checkpoint_controller);
