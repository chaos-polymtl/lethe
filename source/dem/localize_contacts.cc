#include <dem/contact_type.h>
#include <dem/localize_contacts.h>

using namespace dealii;

template <int dim>
void
localize_contacts(DEMContainerManager<dim> &container_manager)
{
  // Update particle-particle contacts in local_adjacent_particles of fine
  // search step with local_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    typename DEM::dem_data_structures<dim>::particle_particle_candidates,
    ContactType::local_particle_particle>(
    container_manager.local_adjacent_particles,
    container_manager.local_contact_pair_candidates);

  // Update particle-particle contacts in global_adjacent_particles of fine
  // search step with global_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    typename DEM::dem_data_structures<dim>::particle_particle_candidates,
    ContactType::ghost_particle_particle>(
    container_manager.ghost_adjacent_particles,
    container_manager.ghost_contact_pair_candidates);

  // Update particle-wall contacts in particle_wall_pairs_in_contact of fine
  // search step with particle_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    typename DEM::dem_data_structures<dim>::particle_wall_candidates,
    ContactType::particle_wall>(container_manager.particle_wall_in_contact,
                                container_manager.particle_wall_candidates);

  // Update particle-floating wall contacts in particle_floating_wall_in_contact
  // of fine search step with particle_floating_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates,
    ContactType::particle_floating_wall>(
    container_manager.particle_floating_wall_in_contact,
    container_manager.particle_floating_wall_candidates);

  // Update particle-floating mesh contacts in particle_floating_mesh_in_contact
  // of fine search step with particle_floating_mesh_contact_candidates
  for (unsigned int solid_counter = 0;
       solid_counter <
       container_manager.particle_floating_mesh_in_contact.size();
       ++solid_counter)
    {
      update_fine_search_candidates<
        dim,
        typename DEM::dem_data_structures<
          dim>::particle_floating_wall_from_mesh_in_contact,
        typename DEM::dem_data_structures<
          dim>::particle_floating_wall_from_mesh_candidates,
        ContactType::particle_floating_mesh>(
        container_manager.particle_floating_mesh_in_contact[solid_counter],
        container_manager.particle_floating_mesh_candidates[solid_counter]);
    }
}

template <int dim,
          typename pairs_structure,
          typename candidates_structure,
          ContactType contact_type>
void
update_fine_search_candidates(pairs_structure &     pairs_in_contact,
                              candidates_structure &contact_candidates)
{
  for (auto pairs_in_contact_iterator = pairs_in_contact.begin();
       pairs_in_contact_iterator != pairs_in_contact.end();
       ++pairs_in_contact_iterator)
    {
      // Get the current particle id from fine search history
      auto particle_id = pairs_in_contact_iterator->first;

      // Get the reference to all the object (particle/wall/face) contact
      // candidates of the current particle from the new broad search
      auto contact_candidate_element = &contact_candidates[particle_id];

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
                        contact_type == ContactType::ghost_particle_particle)
            {
              // Get the reference to all the 2nd particle contact
              // candidates of the current particle from the new broad search
              auto object_contact_candidates = &contact_candidates[object_id];

              // Check if the adjacent particle from history is a particle
              // candidate of the current particle and get the iterator to that
              // object if so
              auto search_iterator =
                std::find(contact_candidate_element->begin(),
                          contact_candidate_element->end(),
                          object_id);

              // Check if the current particle is a particle candidate
              // of the adjacent particle in history
              auto search_object_to_particle_iterator =
                std::find(object_contact_candidates->begin(),
                          object_contact_candidates->end(),
                          particle_id);

              // Particle-particle pairs are removed from the list of the new
              // broad search candidates or from the list of the history of the
              // fine search candidate list depending on the last checks
              if (search_iterator != contact_candidate_element->end())
                {
                  // The current 2nd is a candidate to the current particle.
                  // The 2nd particle is then removed from the new broad search
                  // particle-particle contact candidates list but kept as
                  // history in the fine search list
                  contact_candidate_element->erase(search_iterator);
                  ++adjacent_map_iterator;
                }
              else if (search_object_to_particle_iterator !=
                       object_contact_candidates->end())
                {
                  if constexpr (contact_type ==
                                ContactType::local_particle_particle)
                    {
                      // The current particle from history list is still a
                      // candidate to the 2nd particle and the contact is local
                      // The 2nd particle is then removed from the new broad
                      // search particle-particle contact candidates list but
                      // kept as history in the fine search list
                      object_contact_candidates->erase(
                        search_object_to_particle_iterator);
                      ++adjacent_map_iterator;
                    }

                  if constexpr (contact_type ==
                                ContactType::ghost_particle_particle)
                    {
                      // The current particle from history list is still a
                      // candidate to the 2nd particle and the contact is not
                      // local The 2nd particle is then removed from the history
                      // in the fine search list particle-particle contact
                      // candidates list but kept in new broad search list
                      adjacent_pairs_content->erase(adjacent_map_iterator++);
                    }
                }
              else
                {
                  // Current particles are no longer contact candidates.
                  // The particle is then removed as a candidate for fine search
                  // contact candidates list
                  adjacent_pairs_content->erase(adjacent_map_iterator++);
                }
            }

          if constexpr (contact_type == ContactType::particle_wall ||
                        contact_type == ContactType::particle_floating_wall ||
                        contact_type == ContactType::particle_floating_mesh)
            {
              // Check if the wall/face is a candidate of the current particle
              // and get the iterator to that object if so
              auto search_iterator = contact_candidate_element->find(object_id);

              // Particle-wall/face pair is removed from the list of the new
              // broad search candidates or from the list of the history of the
              // fine search candidate list depending on the last check
              if (search_iterator != contact_candidate_element->end())
                {
                  // Current wall/face (from new broad search) is a candidate to
                  // the current particle. The wall/face is then removed from
                  // the new broad search particle-wall/face contact candidates
                  // list but kept as history in the fine search list
                  contact_candidate_element->erase(search_iterator);
                  ++adjacent_map_iterator;
                }
              else
                {
                  // Current wall/face (from new broad search) is no longer a
                  // candidate to the current particle. The wall/face is then
                  // removed as a candidate for fine search contact candidates
                  // list
                  adjacent_pairs_content->erase(adjacent_map_iterator++);
                }
            }
        }
    }
}

template void localize_contacts<2>(DEMContainerManager<2> &container_manager);

template void localize_contacts<3>(DEMContainerManager<3> &container_manager);


// Local particle-particle contacts
template void update_fine_search_candidates<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  typename DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  typename DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void update_fine_search_candidates<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  typename DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Ghost particle-particle contacts
template void update_fine_search_candidates<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  typename DEM::dem_data_structures<2>::particle_particle_candidates,
  ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  typename DEM::dem_data_structures<2>::particle_particle_candidates
    &contact_pair_candidates);

template void update_fine_search_candidates<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  typename DEM::dem_data_structures<3>::particle_particle_candidates,
  ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_particle_candidates
    &contact_candidates);

// Particle-wall contacts
template void update_fine_search_candidates<
  2,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact,
  typename DEM::dem_data_structures<2>::particle_wall_candidates,
  ContactType::particle_wall>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &adjacent_particles,
  typename DEM::dem_data_structures<2>::particle_wall_candidates
    &contact_pair_candidates);

template void update_fine_search_candidates<
  3,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact,
  typename DEM::dem_data_structures<3>::particle_wall_candidates,
  ContactType::particle_wall>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &adjacent_particles,
  typename DEM::dem_data_structures<3>::particle_wall_candidates
    &contact_pair_candidates);

// Particle-floating wall contacts
template void update_fine_search_candidates<
  2,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact,
  typename DEM::dem_data_structures<2>::particle_floating_wall_candidates,
  ContactType::particle_floating_wall>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &adjacent_particles,
  typename DEM::dem_data_structures<2>::particle_floating_wall_candidates
    &contact_pair_candidates);

template void update_fine_search_candidates<
  3,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact,
  typename DEM::dem_data_structures<3>::particle_floating_wall_candidates,
  ContactType::particle_floating_wall>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &adjacent_particles,
  typename DEM::dem_data_structures<3>::particle_floating_wall_candidates
    &contact_pair_candidates);

// Particle-floating mesh contacts
template void update_fine_search_candidates<
  2,
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_in_contact,
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_candidates,
  ContactType::particle_floating_mesh>(
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_in_contact &adjacent_particles,
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_candidates &contact_pair_candidates);

template void update_fine_search_candidates<
  3,
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_in_contact,
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_candidates,
  ContactType::particle_floating_mesh>(
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_in_contact &adjacent_particles,
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_candidates &contact_pair_candidates);