// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>
#include <dem/find_contact_detection_step.h>

using namespace dealii;

template <int dim>
void
find_particle_contact_detection_step(
  Particles::ParticleHandler<dim> &particle_handler,
  const double                     dt,
  const double                     smallest_contact_search_criterion,
  MPI_Comm                        &mpi_communicator,
  std::vector<double>             &displacement,
  const bool                       parallel_update)
{
  // Get the action manager
  auto *action_manager = DEMActionManager::get_action_manager();
  // If something else has already triggered contact search,
  // no need to do it again
  if (action_manager->check_contact_search())
    return;

  double max_displacement       = 0.;
  bool   contact_detection_step = false;

  // Updating displacement
  for (auto &particle : particle_handler)
    {
      auto  particle_properties   = particle.get_properties();
      auto &particle_displacement = displacement[particle.get_local_index()];

      // Finding displacement of each particle during last step
      particle_displacement +=
        dt * sqrt(particle_properties[DEM::PropertiesIndex::v_x] *
                    particle_properties[DEM::PropertiesIndex::v_x] +
                  particle_properties[DEM::PropertiesIndex::v_y] *
                    particle_properties[DEM::PropertiesIndex::v_y] +
                  particle_properties[DEM::PropertiesIndex::v_z] *
                    particle_properties[DEM::PropertiesIndex::v_z]);

      // Updating maximum displacement of particles
      max_displacement = std::max(max_displacement, particle_displacement);
    }

  // If the maximum displacement of particles exceeds criterion, this step
  // is a contact detection step
  contact_detection_step = max_displacement > smallest_contact_search_criterion;

  // Broadcasting contact detection step value to other processors
  if (parallel_update)
    {
      contact_detection_step =
        Utilities::MPI::logical_or(contact_detection_step, mpi_communicator);
      if (contact_detection_step)
        action_manager->contact_detection_step();
    }
}

template void
find_particle_contact_detection_step(
  Particles::ParticleHandler<2> &particle_handler,
  const double                   dt,
  const double                   smallest_contact_search_criterion,
  MPI_Comm                      &mpi_communicator,
  std::vector<double>           &displacement,
  const bool                     parallel_update);

template void
find_particle_contact_detection_step(
  Particles::ParticleHandler<3> &particle_handler,
  const double                   dt,
  const double                   smallest_contact_search_criterion,
  MPI_Comm                      &mpi_communicator,
  std::vector<double>           &displacement,
  const bool                     parallel_update);


template <int dim>
void
find_floating_mesh_mapping_step(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids)
{
  // Get the action manager
  auto *action_manager = DEMActionManager::get_action_manager();

  // If there is no solid object, no need to do anything
  if (!action_manager->check_solid_objects_enabled())
    return;

  bool floating_mesh_requires_map = false;

  for (unsigned int i_solid = 0; i_solid < solids.size(); ++i_solid)
    {
      double displacement =
        solids[i_solid]->get_max_displacement_since_mapped();
      floating_mesh_requires_map =
        floating_mesh_requires_map ||
        displacement > smallest_contact_search_criterion;
    }

  if (floating_mesh_requires_map)
    action_manager->solid_objects_search_step();
}

template void
find_floating_mesh_mapping_step<2>(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<1, 2>>> solids);

template void
find_floating_mesh_mapping_step<3>(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<2, 3>>> solids);
