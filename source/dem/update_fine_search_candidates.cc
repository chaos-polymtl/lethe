// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/data_containers.h>
#include <dem/update_fine_search_candidates.h>

using namespace dealii;

template <int dim,
          typename pairs_structure,
          typename candidates_structure,
          ContactType contact_type>
void
update_fine_search_candidates(pairs_structure      &pairs_in_contact,
                              candidates_structure &contact_candidates)
{
  for (auto pairs_in_contact_iterator = pairs_in_contact.begin();
       pairs_in_contact_iterator != pairs_in_contact.end();)
    {
      // Get the current particle id from fine search history
      auto particle_id = pairs_in_contact_iterator->first;

      // Get the reference to all the object (particle/wall/face) contact
      // candidates of the current particle from the new broad search
      auto contact_candidate_element = contact_candidates.find(particle_id);

      // Get adjacent object (particle/wall/face) content of the current
      // particle (from the fine search history)
      auto adjacent_pairs_content = &pairs_in_contact_iterator->second;

      // Loop over all the adjacent objects
      for (auto adjacent_map_iterator = adjacent_pairs_content->begin();
           adjacent_map_iterator != adjacent_pairs_content->end();)
        {
          // Get the object (2nd particle/wall/face) id in the history list
          auto object_id = adjacent_map_iterator->first;

          // The handling of the update of the fine search candidate list is
          // different when the contact is particle-particle compared to
          // particle-wall
          if constexpr (contact_type == ContactType::local_particle_particle ||
                        contact_type == ContactType::ghost_particle_particle ||
                        contact_type ==
                          ContactType::local_periodic_particle_particle ||
                        contact_type ==
                          ContactType::ghost_periodic_particle_particle ||
                        contact_type ==
                          ContactType::ghost_local_periodic_particle_particle)
            {
              if (contact_candidate_element != contact_candidates.end())
                {
                  // Check if the adjacent particle from history is a particle
                  // candidate of the current particle and get the iterator to
                  // that object if so
                  auto search_iterator =
                    std::find(contact_candidate_element->second.begin(),
                              contact_candidate_element->second.end(),
                              object_id);

                  // Particle-particle pairs are removed from the list of the
                  // new broad search candidates or from the list of the history
                  // of the fine search candidate list depending on the last
                  // checks
                  if (search_iterator !=
                      contact_candidate_element->second.end())
                    {
                      // The current 2nd is a candidate to the current particle.
                      // The 2nd particle is then removed from the new broad
                      // search particle-particle contact candidates list but
                      // kept as history in the fine search list
                      contact_candidate_element->second.erase(search_iterator);

                      ++adjacent_map_iterator;

                      // To prevent searching through the contact_candidates a
                      // second time if there is a deletion in the first search,
                      // we continue here instead of using an if, else if or
                      // else statement
                      continue;
                    }
                }

              // Get the reference to all the 2nd particle contact
              // candidates of the current particle from the new broad search
              auto object_contact_candidates =
                contact_candidates.find(object_id);

              if (object_contact_candidates != contact_candidates.end())
                {
                  // Check if the current particle is a particle candidate
                  // of the adjacent particle in history
                  auto search_object_to_particle_iterator =
                    std::find(object_contact_candidates->second.begin(),
                              object_contact_candidates->second.end(),
                              particle_id);

                  if (search_object_to_particle_iterator !=
                      object_contact_candidates->second.end())
                    {
                      if constexpr (contact_type ==
                                      ContactType::local_particle_particle ||
                                    contact_type ==
                                      ContactType::
                                        local_periodic_particle_particle)
                        {
                          // The current particle from history list is still a
                          // candidate to the 2nd particle and the contact is
                          // local The 2nd particle is then removed from the new
                          // broad search particle-particle contact candidates
                          // list but kept as history in the fine search list
                          object_contact_candidates->second.erase(
                            search_object_to_particle_iterator);
                          ++adjacent_map_iterator;
                          continue;
                        }

                      if constexpr (
                        contact_type == ContactType::ghost_particle_particle ||
                        contact_type ==
                          ContactType::ghost_periodic_particle_particle ||
                        contact_type ==
                          ContactType::ghost_local_periodic_particle_particle)
                        {
                          // The current particle from history list is still a
                          // candidate to the 2nd particle and the contact is
                          // not local The 2nd particle is then removed from the
                          // history in the fine search list particle-particle
                          // contact candidates list but kept in new broad
                          // search list
                          adjacent_map_iterator = adjacent_pairs_content->erase(
                            adjacent_map_iterator);
                          continue;
                        }
                    }
                }

              // Current particles are no longer contact candidates.
              // The particle is then removed as a candidate for fine
              // search contact candidates list
              adjacent_map_iterator =
                adjacent_pairs_content->erase(adjacent_map_iterator);
            }

          if constexpr (contact_type == ContactType::particle_wall ||
                        contact_type == ContactType::particle_floating_wall ||
                        contact_type == ContactType::particle_floating_mesh)
            {
              if (contact_candidate_element != contact_candidates.end())
                {
                  // Check if the wall/face is a candidate of the current
                  // particle and get the iterator to that object if so
                  auto search_iterator =
                    contact_candidate_element->second.find(object_id);

                  // Particle-wall/face pair is removed from the list of the new
                  // broad search candidates or from the list of the history of
                  // the fine search candidate list depending on the last check
                  if (search_iterator !=
                      contact_candidate_element->second.end())
                    {
                      // Current wall/face (from new broad search) is a
                      // candidate to the current particle. The wall/face is
                      // then removed from the new broad search
                      // particle-wall/face contact candidates list but kept as
                      // history in the fine search list
                      contact_candidate_element->second.erase(search_iterator);
                      ++adjacent_map_iterator;
                      continue;
                    }
                }

              // Current wall/face (from new broad search) is no longer a
              // candidate to the current particle. The wall/face is then
              // removed as a candidate for fine search contact candidates
              // list
              adjacent_map_iterator =
                adjacent_pairs_content->erase(adjacent_map_iterator);
            }
        }

      // If there are still particle/wall/face in the adjacent_pairs_content
      // then the pairs_in_contact_iterator remains in memory
      if (adjacent_pairs_content->size() > 0)
        ++pairs_in_contact_iterator;

      // Otherwise it is deleted. This is necessary to prevent memory inflation.
      // If this is not done, pairs_in_contact keeps on growing until it reaches
      // the total number of particles in the simulation.
      else
        pairs_in_contact_iterator =
          pairs_in_contact.erase(pairs_in_contact_iterator);
    }
}

// Local particle-particle contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::local_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &adjacent_particles,
  DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::local_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Ghost particle-particle contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::ghost_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &adjacent_particles,
  DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::ghost_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Local periodic particle-particle contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::local_periodic_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &adjacent_particles,
  DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::local_periodic_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Local-ghost particle-particle contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::ghost_periodic_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &adjacent_particles,
  DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::ghost_periodic_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Ghost-local particle-particle contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::ghost_local_periodic_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &adjacent_particles,
  DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::ghost_local_periodic_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Particle-wall contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::particle_wall_in_contact,
  DEM::dem_data_structures<2>::particle_wall_candidates,
  ContactType::particle_wall>(
  DEM::dem_data_structures<2>::particle_wall_in_contact &adjacent_particles,
  DEM::dem_data_structures<2>::particle_wall_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::particle_wall_in_contact,
  DEM::dem_data_structures<3>::particle_wall_candidates,
  ContactType::particle_wall>(
  DEM::dem_data_structures<3>::particle_wall_in_contact &adjacent_particles,
  DEM::dem_data_structures<3>::particle_wall_candidates
    &contact_pair_candidates);

// Particle-floating wall contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<2>::particle_wall_in_contact,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates,
  ContactType::particle_floating_wall>(
  DEM::dem_data_structures<2>::particle_wall_in_contact &adjacent_particles,
  DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<3>::particle_wall_in_contact,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates,
  ContactType::particle_floating_wall>(
  DEM::dem_data_structures<3>::particle_wall_in_contact &adjacent_particles,
  DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &contact_pair_candidates);

// Particle-floating mesh contacts
template void
update_fine_search_candidates<
  2,
  DEM::dem_data_structures<
    2>::particle_triangle_cell_from_mesh_potentially_in_contact,
  DEM::dem_data_structures<2>::particle_floating_wall_from_mesh_candidates,
  ContactType::particle_floating_mesh>(
  DEM::dem_data_structures<2>::
    particle_triangle_cell_from_mesh_potentially_in_contact &adjacent_particles,
  DEM::dem_data_structures<2>::particle_floating_wall_from_mesh_candidates
    &contact_pair_candidates);

template void
update_fine_search_candidates<
  3,
  DEM::dem_data_structures<
    3>::particle_triangle_cell_from_mesh_potentially_in_contact,
  DEM::dem_data_structures<3>::particle_floating_wall_from_mesh_candidates,
  ContactType::particle_floating_mesh>(
  DEM::dem_data_structures<3>::
    particle_triangle_cell_from_mesh_potentially_in_contact &adjacent_particles,
  DEM::dem_data_structures<3>::particle_floating_wall_from_mesh_candidates
    &contact_pair_candidates);
