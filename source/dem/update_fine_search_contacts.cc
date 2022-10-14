#include <dem/particle_particle_contact_info_struct.h>
#include <dem/update_fine_search_contacts.h>

using namespace dealii;

template <int dim, typename pairs_structure, typename candidates_structure>
void
update_fine_search_candidates(pairs_structure &     pairs_in_contact,
                              candidates_structure &contact_candidates,
                              const ContactType     contact_type)
{
  for (auto pairs_in_contact_iterator = pairs_in_contact.begin();
       pairs_in_contact_iterator != pairs_in_contact.end();
       ++pairs_in_contact_iterator)
    {
      // Get the current particle id from fine search history
      auto particle_id = pairs_in_contact_iterator->first;

      // Get the reference to all the object (particle/wall/face) contact
      // candidates of the current particle from the new board search
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
          // Since particle-particle contacts have a specific data structure,
          // the constexpr is used to improve performance instead of ContactType
          if constexpr (std::is_same_v<
                          pairs_structure,
                          typename dem_data_containers::dem_data_structures<
                            dim>::adjacent_particle_pairs>)
            {
              // Get the reference to all the 2nd particle contact
              // candidates of the current particle from the new board search
              auto object_contact_candidates = &contact_candidates[object_id];

              // Check if the adjacent particle from history is a
              // particle candidate of the current particle and get the
              // iterator to that object if so
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
              // board search candidates or from the list of the history of the
              // fine search candidate list depending on the last checks
              if (search_iterator != contact_candidate_element->end())
                {
                  // The current 2nd is a candidate to the current particle.
                  // The 2nd particle is then removed from the new board search
                  // particle-particle contact candidates list but kept as
                  // history in the fine search list
                  contact_candidate_element->erase(search_iterator);
                  ++adjacent_map_iterator;
                }
              else if (search_object_to_particle_iterator !=
                         object_contact_candidates->end() &&
                       contact_type == ContactType::local_particle_particle)
                {
                  // The current particle from history list is still a candidate
                  // to the 2nd particle and .
                  // The 2nd particle is then removed from the new board search
                  // particle-particle contact candidates list but kept as
                  // history in the fine search list
                  object_contact_candidates->erase(
                    search_object_to_particle_iterator);
                  ++adjacent_map_iterator;
                }
              else if (search_object_to_particle_iterator !=
                         object_contact_candidates->end() &&
                       contact_type == ContactType::ghost_particle_particle)
                {
                  adjacent_pairs_content->erase(adjacent_map_iterator++);
                }
              else
                {
                  // Current particles are no longer contact candidates.
                  // The particle is then removed as a candidate for fine search
                  // contact candidates list.
                  adjacent_pairs_content->erase(adjacent_map_iterator++);
                }
            }
          else
            {
              // Check if the wall/face is a candidate of the current particle
              // and get the iterator to that object if so
              auto search_iterator = contact_candidate_element->find(object_id);

              // Particle-wall/face pair is removed from the list of the new
              // board search candidates or from the list of the history of the
              // fine search candidate list depending on the last check
              if (search_iterator != contact_candidate_element->end())
                {
                  // Current wall/face (from new board search) is a candidate to
                  // the current particle. The wall/face is then removed from
                  // the new board search particle-wall/face contact candidates
                  // list but kept as history in the fine search list
                  contact_candidate_element->erase(search_iterator);
                  ++adjacent_map_iterator;
                }
              else
                {
                  // Current wall/face (from new board search) is no longer a
                  // candidate to the current particle. The wall/face is then
                  // removed as a candidate for fine search contact candidates
                  // list.
                  adjacent_pairs_content->erase(adjacent_map_iterator++);
                }
            }
        }
    }
}

// Particle-particle contacts
template void update_fine_search_candidates<
  2,
  typename dem_data_containers::dem_data_structures<2>::adjacent_particle_pairs,
  typename dem_data_containers::dem_data_structures<
    2>::particle_particle_candidates>(
  typename dem_data_containers::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    2>::particle_particle_candidates &contact_pair_candidates,
  const ContactType                   contact_type);

template void update_fine_search_candidates<
  3,
  typename dem_data_containers::dem_data_structures<3>::adjacent_particle_pairs,
  typename dem_data_containers::dem_data_structures<
    3>::particle_particle_candidates>(
  typename dem_data_containers::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename dem_data_containers::dem_data_structures<
    3>::particle_particle_candidates &contact_candidates,
  const ContactType                   contact_type);

// Particle-wall contacts
template void update_fine_search_candidates<
  2,
  typename dem_data_containers::dem_data_structures<
    2>::particle_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    2>::particle_wall_candidates>(
  typename dem_data_containers::dem_data_structures<2>::particle_wall_in_contact
    &adjacent_particles,
  typename dem_data_containers::dem_data_structures<2>::particle_wall_candidates
    &               contact_pair_candidates,
  const ContactType contact_type);

template void update_fine_search_candidates<
  3,
  typename dem_data_containers::dem_data_structures<
    3>::particle_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    3>::particle_wall_candidates>(
  typename dem_data_containers::dem_data_structures<3>::particle_wall_in_contact
    &adjacent_particles,
  typename dem_data_containers::dem_data_structures<3>::particle_wall_candidates
    &               contact_pair_candidates,
  const ContactType contact_type);

// Particle-floating wall contacts
template void update_fine_search_candidates<
  3,
  typename dem_data_containers::dem_data_structures<
    3>::particle_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    3>::particle_floating_wall_candidates>(
  typename dem_data_containers::dem_data_structures<3>::particle_wall_in_contact
    &adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    3>::particle_floating_wall_candidates &contact_pair_candidates,
  const ContactType                        contact_type);

template void update_fine_search_candidates<
  2,
  typename dem_data_containers::dem_data_structures<
    2>::particle_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    2>::particle_floating_wall_candidates>(
  typename dem_data_containers::dem_data_structures<2>::particle_wall_in_contact
    &adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    2>::particle_floating_wall_candidates &contact_pair_candidates,
  const ContactType                        contact_type);

// Particle-floating mesh contacts
template void update_fine_search_candidates<
  2,
  std::map<typename Triangulation<1, 2>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<2>>,
           typename dem_data_containers::cut_cell_comparison<2>>,
  std::map<
    typename Triangulation<1, 2>::active_cell_iterator,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>,
    typename dem_data_containers::cut_cell_comparison<2>>>(
  std::map<typename Triangulation<1, 2>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<2>>,
           typename dem_data_containers::cut_cell_comparison<2>>
    &adjacent_particles,
  std::map<
    typename Triangulation<1, 2>::active_cell_iterator,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>,
    typename dem_data_containers::cut_cell_comparison<2>>
    &               contact_pair_candidates,
  const ContactType contact_type);

template void update_fine_search_candidates<
  3,
  std::map<typename Triangulation<2, 3>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<3>>,
           typename dem_data_containers::cut_cell_comparison<3>>,
  std::map<
    typename Triangulation<2, 3>::active_cell_iterator,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>,
    typename dem_data_containers::cut_cell_comparison<3>>>(
  std::map<typename Triangulation<2, 3>::active_cell_iterator,
           std::unordered_map<types::particle_index,
                              particle_wall_contact_info_struct<3>>,
           typename dem_data_containers::cut_cell_comparison<3>>
    &adjacent_particles,
  std::map<
    typename Triangulation<2, 3>::active_cell_iterator,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>,
    typename dem_data_containers::cut_cell_comparison<3>>
    &               contact_pair_candidates,
  const ContactType contact_type);