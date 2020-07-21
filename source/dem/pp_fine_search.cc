#include <dem/pp_fine_search.h>

using namespace dealii;

template <int dim> PPFineSearch<dim>::PPFineSearch() {}

template <int dim>
void PPFineSearch<dim>::pp_Fine_Search(
    const std::map<std::pair<int, int>,
                   std::pair<typename Particles::ParticleIterator<dim>,
                             typename Particles::ParticleIterator<dim>>>
        &local_contact_pair_candidates,
    const std::map<std::pair<int, int>,
                   std::pair<typename Particles::ParticleIterator<dim>,
                             typename Particles::ParticleIterator<dim>>>
        &ghost_contact_pair_candidates,
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
        &local_adjacent_particles,
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
        &ghost_adjacent_particles,
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
      double distance = particle_one_location.distance(particle_two_location);
      if (distance > neighborhood_threshold) {
        adjacent_particles_list->erase(adjacent_particles_list_iterator++);
      } else {
        ++adjacent_particles_list_iterator;
      }
    }
  }

  // Now iterating over local_contact_pair_candidates (vector of pairs), which
  // is the output of broad search. If a pair is in contact (distance > 0) and
  // does not exist in the pairs_in_contact, it is added to the
  // pairs_in_contact
  for (auto candidate_iterator = local_contact_pair_candidates.begin();
       candidate_iterator != local_contact_pair_candidates.end();
       ++candidate_iterator) {
    auto particle_contact_pair = candidate_iterator->second;

    // Get particles one and two from the vector and the total array view to
    // the particle properties once to improve efficiency
    auto particle_one = particle_contact_pair.first;
    auto particle_two = particle_contact_pair.second;

    // Obtaining locations of particles one and two:
    Point<dim, double> particle_one_location = particle_one->get_location();
    Point<dim, double> particle_two_location = particle_two->get_location();

    // Finding distance
    double distance = particle_one_location.distance(particle_two_location);

    // If the particles distance is less than the threshold
    if (distance < neighborhood_threshold) {
      auto particle_one_properties = particle_one->get_properties();
      auto particle_two_properties = particle_two->get_properties();

      // Getting the particle one contact list and particle two id
      auto particle_one_contact_list =
          &local_adjacent_particles
              [particle_one_properties[DEM::PropertiesIndex::id]];
      unsigned int particle_two_id =
          particle_two_properties[DEM::PropertiesIndex::id];

      // CORRECT**
      //     auto particle_two_contact_list =
      //         &local_adjacent_particles
      //            [particle_two_properties[DEM::PropertiesIndex::id]];
      //     unsigned int particle_one_id =
      //         particle_one_properties[DEM::PropertiesIndex::id];

      // If the pair does not exist in the local_adjacent_particles
      //   if (particle_one_contact_list->count(particle_two_id) <= 0 &&
      //     particle_two_contact_list->count(particle_one_id) <= 0) {
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
      //    }
    }
  }

  // Second iterating over ghost adjacent_particles
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
      double distance = particle_one_location.distance(particle_two_location);
      if (distance > neighborhood_threshold) {
        adjacent_particles_list->erase(adjacent_particles_list_iterator++);
      } else {
        ++adjacent_particles_list_iterator;
      }
    }
  }

  // Now iterating over ghost_contact_pair_candidates (vector of pairs), which
  // is the output of broad search. If a pair is in contact (distance > 0) and
  // does not exist in the pairs_in_contact, it is added to the
  // pairs_in_contact
  for (auto candidate_iterator = ghost_contact_pair_candidates.begin();
       candidate_iterator != ghost_contact_pair_candidates.end();
       ++candidate_iterator) {
    auto particle_contact_pair = candidate_iterator->second;

    // Get particles one and two from the vector and the total array view to
    // the particle properties once to improve efficiency
    auto particle_one = particle_contact_pair.first;
    auto particle_two = particle_contact_pair.second;

    // Obtaining locations of particles one and two:
    Point<dim, double> particle_one_location = particle_one->get_location();
    Point<dim, double> particle_two_location = particle_two->get_location();

    // Finding distance
    double distance = particle_one_location.distance(particle_two_location);

    // If the particles distance is less than the threshold
    if (distance < neighborhood_threshold) {
      auto particle_one_properties = particle_one->get_properties();
      auto particle_two_properties = particle_two->get_properties();
      // Getting the particle one contact list and particle two id
      auto particle_one_contact_list =
          &ghost_adjacent_particles
              [particle_one_properties[DEM::PropertiesIndex::id]];
      unsigned int particle_two_id =
          particle_two_properties[DEM::PropertiesIndex::id];

      // CORRECT**
      //    auto particle_two_contact_list =
      //        &ghost_adjacent_particles
      //            [particle_two_properties[DEM::PropertiesIndex::id]];
      //    unsigned int particle_one_id =
      //        particle_one_properties[DEM::PropertiesIndex::id];

      // If the pair does not exist in the ghost_adjacent_particles
      //   if (particle_one_contact_list->count(particle_two_id) <= 0 &&
      //      particle_two_contact_list->count(particle_one_id) <= 0) {
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
      //     }
    }
  }
}

template class PPFineSearch<2>;
template class PPFineSearch<3>;
