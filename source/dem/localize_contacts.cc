#include <dem/localize_contacts.h>

using namespace dealii;

template <int dim>
void
localize_contacts(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    *local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    *ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    *particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    *pfw_pairs_in_contact,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &local_contact_pair_candidates,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &ghost_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  unsigned int,
                                  unsigned int>>>
    &particle_wall_contact_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>>
    pfw_contact_candidates)

{
  for (auto adjacent_particles_iterator = local_adjacent_particles->begin();
       adjacent_particles_iterator != local_adjacent_particles->end();
       ++adjacent_particles_iterator)
    {
      unsigned int particle_one_id = adjacent_particles_iterator->first;
      auto         particle_one_contact_candidates =
        &local_contact_pair_candidates[particle_one_id];

      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto particle_particle_map_iterator =
             pairs_in_contant_content->begin();
           particle_particle_map_iterator != pairs_in_contant_content->end();)
        {
          unsigned int particle_two_id = particle_particle_map_iterator->first;
          auto         particle_two_contact_candidates =
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
              ++particle_particle_map_iterator;
            }
          else if (search_iterator_two !=
                   particle_two_contact_candidates->end())
            {
              particle_two_contact_candidates->erase(search_iterator_two);
              ++particle_particle_map_iterator;
            }
          else
            {
              pairs_in_contant_content->erase(particle_particle_map_iterator++);
            }
        }
    }

  // The same for local-ghost particle containers
  for (auto adjacent_particles_iterator = ghost_adjacent_particles->begin();
       adjacent_particles_iterator != ghost_adjacent_particles->end();
       ++adjacent_particles_iterator)
    {
      unsigned int particle_one_id = adjacent_particles_iterator->first;
      auto         particle_one_contact_candidates =
        &ghost_contact_pair_candidates[particle_one_id];

      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto particle_particle_map_iterator =
             pairs_in_contant_content->begin();
           particle_particle_map_iterator != pairs_in_contant_content->end();)
        {
          unsigned int particle_two_id = particle_particle_map_iterator->first;
          auto         particle_two_contact_candidates =
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
              ++particle_particle_map_iterator;
            }
          else if (search_iterator_two !=
                   particle_two_contact_candidates->end())
            {
              pairs_in_contant_content->erase(particle_particle_map_iterator++);
            }
          else
            {
              pairs_in_contant_content->erase(particle_particle_map_iterator++);
            }
        }
    }

  // Particle-wall contacts
  for (auto particle_wall_pairs_in_contact_iterator =
         particle_wall_pairs_in_contact->begin();
       particle_wall_pairs_in_contact_iterator !=
       particle_wall_pairs_in_contact->end();
       ++particle_wall_pairs_in_contact_iterator)
    {
      unsigned int particle_id = particle_wall_pairs_in_contact_iterator->first;

      auto pairs_in_contant_content =
        &particle_wall_pairs_in_contact_iterator->second;

      for (auto particle_wall_map_iterator = pairs_in_contant_content->begin();
           particle_wall_map_iterator != pairs_in_contant_content->end();)
        {
          unsigned int face_id = particle_wall_map_iterator->first;
          auto         particle_wall_contact_candidate_element =
            &particle_wall_contact_candidates[particle_id];

          auto search_iterator =
            particle_wall_contact_candidate_element->find(face_id);

          if (search_iterator != particle_wall_contact_candidate_element->end())
            {
              particle_wall_contact_candidate_element->erase(search_iterator);
              ++particle_wall_map_iterator;
            }
          else
            {
              pairs_in_contant_content->erase(particle_wall_map_iterator++);
            }
        }
    }

  // Particle-floating wall contacts
  for (auto pfw_pairs_in_contact_iterator = pfw_pairs_in_contact->begin();
       pfw_pairs_in_contact_iterator != pfw_pairs_in_contact->end();
       ++pfw_pairs_in_contact_iterator)
    {
      unsigned int particle_id = pfw_pairs_in_contact_iterator->first;

      auto pairs_in_contant_content = &pfw_pairs_in_contact_iterator->second;

      for (auto pfw_map_iterator = pairs_in_contant_content->begin();
           pfw_map_iterator != pairs_in_contant_content->end();)
        {
          unsigned int floating_wall_id = pfw_map_iterator->first;
          auto         pfw_contact_candidate_element =
            &pfw_contact_candidates[particle_id];

          auto search_iterator =
            pfw_contact_candidate_element->find(floating_wall_id);

          if (search_iterator != pfw_contact_candidate_element->end())
            {
              pfw_contact_candidate_element->erase(search_iterator);
              ++pfw_map_iterator;
            }
          else
            {
              pairs_in_contant_content->erase(pfw_map_iterator++);
            }
        }
    }
}

template void localize_contacts(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<2>>>
    *local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<2>>>
    *ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<2>>>
    *particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<2>>>
    *pfw_pairs_in_contact,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &local_contact_pair_candidates,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &ghost_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<2>,
                                  Tensor<1, 2>,
                                  Point<2>,
                                  unsigned int,
                                  unsigned int>>>
    &particle_wall_contact_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>>
    pfw_contact_candidates);

template void localize_contacts(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<3>>>
    *local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<3>>>
    *ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<3>>>
    *particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<3>>>
    *pfw_pairs_in_contact,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &local_contact_pair_candidates,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &ghost_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       std::tuple<Particles::ParticleIterator<3>,
                                  Tensor<1, 3>,
                                  Point<3>,
                                  unsigned int,
                                  unsigned int>>>
    &particle_wall_contact_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>>
    pfw_contact_candidates);
