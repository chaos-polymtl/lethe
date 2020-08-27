#include <dem/localize_contacts.h>

using namespace dealii;

template <int dim>
void
localize_contacts(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *ghost_adjacent_particles,
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    *                                        pw_pairs_in_contact,
  std::unordered_map<int, std::vector<int>> &local_contact_pair_candidates,
  std::unordered_map<int, std::vector<int>> &ghost_contact_pair_candidates,
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    &pw_contact_candidates)

{
  for (auto adjacent_particles_iterator = local_adjacent_particles->begin();
       adjacent_particles_iterator != local_adjacent_particles->end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id = adjacent_particles_iterator->first;
      auto particle_one_contact_candidates =
        &local_contact_pair_candidates[particle_one_id];

      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();)
        {
          int  particle_two_id = pp_map_iterator->first;
          auto particle_two_contact_candidates =
            &local_contact_pair_candidates[particle_two_id];

          auto search_iterator_one =
            std::find(particle_one_contact_candidates->begin(),
                      particle_one_contact_candidates->end(),
                      particle_two_id);
          auto search_iterator_two =
            std::find(particle_two_contact_candidates->begin(),
                      particle_two_contact_candidates->end(),
                      particle_one_id);

          if (search_iterator_one != particle_one_contact_candidates->end())
            {
              particle_one_contact_candidates->erase(search_iterator_one);
              ++pp_map_iterator;
            }
          else if (search_iterator_two !=
                   particle_two_contact_candidates->end())
            {
              particle_two_contact_candidates->erase(search_iterator_two);
              ++pp_map_iterator;
            }
          else
            {
              pairs_in_contant_content->erase(pp_map_iterator++);
            }
        }
    }

  // The same for local-ghost particle containers
  for (auto adjacent_particles_iterator = ghost_adjacent_particles->begin();
       adjacent_particles_iterator != ghost_adjacent_particles->end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id = adjacent_particles_iterator->first;
      auto particle_one_contact_candidates =
        &ghost_contact_pair_candidates[particle_one_id];

      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();)
        {
          int  particle_two_id = pp_map_iterator->first;
          auto particle_two_contact_candidates =
            &ghost_contact_pair_candidates[particle_two_id];

          auto search_iterator_one =
            std::find(particle_one_contact_candidates->begin(),
                      particle_one_contact_candidates->end(),
                      particle_two_id);
          auto search_iterator_two =
            std::find(particle_two_contact_candidates->begin(),
                      particle_two_contact_candidates->end(),
                      particle_one_id);

          if (search_iterator_one != particle_one_contact_candidates->end())
            {
              particle_one_contact_candidates->erase(search_iterator_one);
              ++pp_map_iterator;
            }
          else if (search_iterator_two !=
                   particle_two_contact_candidates->end())
            {
              pairs_in_contant_content->erase(pp_map_iterator++);
            }
          else
            {
              pairs_in_contant_content->erase(pp_map_iterator++);
            }
        }
    }

  // Particle-wall contacts
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact->begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact->end();
       ++pw_pairs_in_contact_iterator)
    {
      int particle_id = pw_pairs_in_contact_iterator->first;

      auto pairs_in_contant_content = &pw_pairs_in_contact_iterator->second;

      for (auto pw_map_iterator = pairs_in_contant_content->begin();
           pw_map_iterator != pairs_in_contant_content->end();)
        {
          int  face_id = pw_map_iterator->first;
          auto pw_contact_candidate_element =
            &pw_contact_candidates[particle_id];

          auto search_iterator = pw_contact_candidate_element->find(face_id);

          if (search_iterator != pw_contact_candidate_element->end())
            {
              pw_contact_candidate_element->erase(search_iterator);
              ++pw_map_iterator;
            }
          else
            {
              pairs_in_contant_content->erase(pw_map_iterator++);
            }
        }
    }
}

template void
localize_contacts(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    *local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<2>>>
    *ghost_adjacent_particles,
  std::map<int, std::map<int, pw_contact_info_struct<2>>> *pw_pairs_in_contact,
  std::unordered_map<int, std::vector<int>> &local_contact_pair_candidates,
  std::unordered_map<int, std::vector<int>> &ghost_contact_pair_candidates,
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<2>, Tensor<1, 2>, Point<2>>>>
    &pw_contact_candidates);

template void
localize_contacts(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<3>>>
    *local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<3>>>
    *ghost_adjacent_particles,
  std::map<int, std::map<int, pw_contact_info_struct<3>>> *pw_pairs_in_contact,
  std::unordered_map<int, std::vector<int>> &local_contact_pair_candidates,
  std::unordered_map<int, std::vector<int>> &ghost_contact_pair_candidates,
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<3>, Tensor<1, 3>, Point<3>>>>
    &pw_contact_candidates);
