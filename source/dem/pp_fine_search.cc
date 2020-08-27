#include <dem/pp_fine_search.h>

using namespace dealii;

template <int dim> PPFineSearch<dim>::PPFineSearch() {}

template <int dim>
void PPFineSearch<dim>::pp_Fine_Search(
    const std::unordered_map<int, std::vector<int>>
        &local_contact_pair_candidates,
    const std::unordered_map<int, std::vector<int>>
        &ghost_contact_pair_candidates,
    std::unordered_map<int,
                       std::unordered_map<int, pp_contact_info_struct<dim>>>
        &local_adjacent_particles,
    std::unordered_map<int,
                       std::unordered_map<int, pp_contact_info_struct<dim>>>
        &ghost_adjacent_particles,
    std::unordered_map<int, Particles::ParticleIterator<dim>>
        &particle_container,
    const double neighborhood_threshold) {
  // First iterating over local adjacent_particles
  for (auto adjacent_particles_iterator = local_adjacent_particles.begin();
       adjacent_particles_iterator != local_adjacent_particles.end();
       ++adjacent_particles_iterator) {
    // Each element of adjacent_particles is a map:
    auto adjacent_particles_list = &adjacent_particles_iterator->second;

    // Iterating over each map which contains the contact information
    // including particles I and II
    for (auto adjacent_particles_list_iterator =
             adjacent_particles_list->begin();
         adjacent_particles_list_iterator != adjacent_particles_list->end();) {
      // Getting contact information and particles I and II as local
      // variables
      auto adjacent_pair_information = adjacent_particles_list_iterator->second;
      auto particle_one = adjacent_pair_information.particle_one;
      auto particle_two = adjacent_pair_information.particle_two;

      // Finding the properties of the particles in contact
      Point<dim, double> particle_one_location = particle_one->get_location();
      Point<dim, double> particle_two_location = particle_two->get_location();

      // Finding distance
      const double square_distance =
          particle_one_location.distance_square(particle_two_location);
      if (square_distance > neighborhood_threshold) {
        adjacent_particles_list->erase(adjacent_particles_list_iterator++);
      } else {
        ++adjacent_particles_list_iterator;
      }
    }
  }

  // Now iterating over local_contact_pair_candidates (maps of pairs), which
  // is the output of broad search. If a pair is in vicinity (distance <
  // threshold), it is added to the local_adjacent_particles
  for (auto candidate_iterator = local_contact_pair_candidates.begin();
       candidate_iterator != local_contact_pair_candidates.end();
       ++candidate_iterator) {
    unsigned int particle_one_id = candidate_iterator->first;
    auto second_particle_container = &candidate_iterator->second;
    auto particle_one = particle_container[particle_one_id];

    for (auto second_particle_iterator = second_particle_container->begin();
         second_particle_iterator != second_particle_container->end();
         ++second_particle_iterator) {
      unsigned int particle_two_id = *second_particle_iterator;

      auto particle_two = particle_container[particle_two_id];

      // Obtaining locations of particles one and two:
      Point<dim, double> particle_one_location = particle_one->get_location();
      Point<dim, double> particle_two_location = particle_two->get_location();

      // Finding distance
      const double square_distance =
          particle_one_location.distance_square(particle_two_location);

      // If the particles distance is less than the threshold
      if (square_distance < neighborhood_threshold) {
        // Getting the particle one contact list and particle two id
        auto particle_one_contact_list =
            &local_adjacent_particles[particle_one_id];

        Tensor<1, dim> tangential_overlap;
        for (int d = 0; d < dim; ++d) {
          tangential_overlap[d] = 0;
        }

        // Initilizing the contact info and adding
        pp_contact_info_struct<dim> contact_info;
        contact_info.tangential_overlap = tangential_overlap;
        contact_info.particle_one = particle_one;
        contact_info.particle_two = particle_two;

        particle_one_contact_list->insert({particle_two_id, contact_info});
      }
    }
  }

  // Second iterating over local-ghost adjacent_particles
  for (auto adjacent_particles_iterator = ghost_adjacent_particles.begin();
       adjacent_particles_iterator != ghost_adjacent_particles.end();
       ++adjacent_particles_iterator) {
    // Each element of adjacent_particles is a map:
    auto adjacent_particles_list = &adjacent_particles_iterator->second;

    // Iterating over each map which contains the contact information
    // including particles I and II
    for (auto adjacent_particles_list_iterator =
             adjacent_particles_list->begin();
         adjacent_particles_list_iterator != adjacent_particles_list->end();) {
      // Getting contact information and particles I and II as local
      // variables
      auto adjacent_pair_information = adjacent_particles_list_iterator->second;
      auto particle_one = adjacent_pair_information.particle_one;
      auto particle_two = adjacent_pair_information.particle_two;

      // Finding the properties of the particles in contact
      Point<dim, double> particle_one_location = particle_one->get_location();
      Point<dim, double> particle_two_location = particle_two->get_location();

      // Finding distance
      const double square_distance =
          particle_one_location.distance_square(particle_two_location);
      if (square_distance > neighborhood_threshold) {
        adjacent_particles_list->erase(adjacent_particles_list_iterator++);
      } else {
        ++adjacent_particles_list_iterator;
      }
    }
  }

  // Now iterating over ghost_contact_pair_candidates (map of pairs), which
  // is the output of broad search. If a pair is in vicinity (distance <
  // threshold), it is added to the ghost_adjacent_particles
  for (auto candidate_iterator = ghost_contact_pair_candidates.begin();
       candidate_iterator != ghost_contact_pair_candidates.end();
       ++candidate_iterator) {
    unsigned int particle_one_id = candidate_iterator->first;
    auto second_particle_container = &candidate_iterator->second;

    auto particle_one = particle_container[particle_one_id];

    for (auto second_particle_iterator = second_particle_container->begin();
         second_particle_iterator != second_particle_container->end();
         ++second_particle_iterator) {
      unsigned int particle_two_id = *second_particle_iterator;

      auto particle_two = particle_container[particle_two_id];

      // Obtaining locations of particles one and two:
      Point<dim, double> particle_one_location = particle_one->get_location();
      Point<dim, double> particle_two_location = particle_two->get_location();

      // Finding distance
      const double square_distance =
          particle_one_location.distance_square(particle_two_location);

      // If the particles distance is less than the threshold
      if (square_distance < neighborhood_threshold) {
        auto particle_one_properties = particle_one->get_properties();
        auto particle_two_properties = particle_two->get_properties();

        // Getting the particle one contact list and particle two id
        auto particle_one_contact_list =
            &ghost_adjacent_particles
                [particle_one_properties[DEM::PropertiesIndex::id]];
        unsigned int particle_two_id =
            particle_two_properties[DEM::PropertiesIndex::id];

        Tensor<1, dim> tangential_overlap;
        for (int d = 0; d < dim; ++d) {
          tangential_overlap[d] = 0;
        }

        // Initilizing the contact info and adding
        pp_contact_info_struct<dim> contact_info;
        contact_info.tangential_overlap = tangential_overlap;
        contact_info.particle_one = particle_one;
        contact_info.particle_two = particle_two;

        particle_one_contact_list->insert({particle_two_id, contact_info});
      }
    }
  }
}

template class PPFineSearch<2>;
template class PPFineSearch<3>;
