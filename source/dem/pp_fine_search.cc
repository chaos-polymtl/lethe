#include <dem/pp_fine_search.h>

using namespace dealii;

template <int dim>
PPFineSearch<dim>::PPFineSearch()
{}

template <int dim>
void
PPFineSearch<dim>::particle_particle_fine_search(
  const std::unordered_map<types::particle_index,
                           std::vector<types::particle_index>>
    &local_contact_pair_candidates,
  const std::unordered_map<types::particle_index,
                           std::vector<types::particle_index>>
    &ghost_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
    &ghost_adjacent_particles,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &          particle_container,
  const double neighborhood_threshold)
{
  // First iterating over local adjacent_particles
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)

    {
      // Iterating over each map which contains the contact information
      // including particles I and II
      for (auto adjacent_particles_list_iterator =
             adjacent_particles_list.begin();
           adjacent_particles_list_iterator != adjacent_particles_list.end();)
        {
          // Getting contact information and particles I and II as local
          // variables
          auto adjacent_pair_information =
            adjacent_particles_list_iterator->second;
          auto particle_one = adjacent_pair_information.particle_one;
          auto particle_two = adjacent_pair_information.particle_two;

          // Finding the properties of the particles in contact
          Point<dim, double> particle_one_location =
            particle_one->get_location();
          Point<dim, double> particle_two_location =
            particle_two->get_location();

          // Finding distance
          const double square_distance =
            particle_one_location.distance_square(particle_two_location);
          if (square_distance > neighborhood_threshold)
            {
              adjacent_particles_list.erase(adjacent_particles_list_iterator++);
            }
          else
            {
              ++adjacent_particles_list_iterator;
            }
        }
    }

  // Now iterating over local_contact_pair_candidates (maps of pairs), which
  // is the output of broad search. If a pair is in vicinity (distance <
  // threshold), it is added to the local_adjacent_particles
  for (auto const &[particle_one_id, second_particle_container] :
       local_contact_pair_candidates)
    {
      auto particle_one = particle_container[particle_one_id];

      for (const unsigned int &particle_two_id : second_particle_container)
        {
          auto particle_two = particle_container[particle_two_id];

          // Obtaining locations of particles one and two:
          Point<dim, double> particle_one_location =
            particle_one->get_location();
          Point<dim, double> particle_two_location =
            particle_two->get_location();

          // Finding distance
          const double square_distance =
            particle_one_location.distance_square(particle_two_location);

          // If the particles distance is less than the threshold
          if (square_distance < neighborhood_threshold)
            {
              // Getting the particle one contact list and particle two id
              auto particle_one_contact_list =
                &local_adjacent_particles[particle_one_id];

              Tensor<1, dim> tangential_overlap;
              for (int d = 0; d < dim; ++d)
                {
                  tangential_overlap[d] = 0;
                }

              // Initilizing the contact info and adding
              pp_contact_info_struct<dim> contact_info;
              contact_info.tangential_overlap = tangential_overlap;
              contact_info.particle_one       = particle_one;
              contact_info.particle_two       = particle_two;

              particle_one_contact_list->insert(
                {particle_two_id, contact_info});
            }
        }
    }

  // Second iterating over local-ghost adjacent_particles
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      // Iterating over each map which contains the contact information
      // including particles I and II
      for (auto adjacent_particles_list_iterator =
             adjacent_particles_list.begin();
           adjacent_particles_list_iterator != adjacent_particles_list.end();)
        {
          // Getting contact information and particles I and II as local
          // variables
          auto adjacent_pair_information =
            adjacent_particles_list_iterator->second;
          auto particle_one = adjacent_pair_information.particle_one;
          auto particle_two = adjacent_pair_information.particle_two;

          // Finding the properties of the particles in contact
          Point<dim, double> particle_one_location =
            particle_one->get_location();
          Point<dim, double> particle_two_location =
            particle_two->get_location();

          // Finding distance
          const double square_distance =
            particle_one_location.distance_square(particle_two_location);
          if (square_distance > neighborhood_threshold)
            {
              adjacent_particles_list.erase(adjacent_particles_list_iterator++);
            }
          else
            {
              ++adjacent_particles_list_iterator;
            }
        }
    }

  // Now iterating over ghost_contact_pair_candidates (map of pairs), which
  // is the output of broad search. If a pair is in vicinity (distance <
  // threshold), it is added to the ghost_adjacent_particles
  for (auto const &[particle_one_id, second_particle_container] :
       ghost_contact_pair_candidates)
    {
      auto particle_one = particle_container[particle_one_id];

      for (const unsigned int &particle_two_id : second_particle_container)
        {
          auto particle_two = particle_container[particle_two_id];

          // Obtaining locations of particles one and two:
          Point<dim, double> particle_one_location =
            particle_one->get_location();
          Point<dim, double> particle_two_location =
            particle_two->get_location();

          // Finding distance
          const double square_distance =
            particle_one_location.distance_square(particle_two_location);

          // If the particles distance is less than the threshold
          if (square_distance < neighborhood_threshold)
            {
              // Getting the particle one contact list and particle two id
              auto particle_one_contact_list =
                &ghost_adjacent_particles[particle_one->get_id()];
              unsigned int particle_two_id = particle_two->get_id();

              Tensor<1, dim> tangential_overlap;
              for (int d = 0; d < dim; ++d)
                {
                  tangential_overlap[d] = 0;
                }

              // Initilizing the contact info and adding
              pp_contact_info_struct<dim> contact_info;
              contact_info.tangential_overlap = tangential_overlap;
              contact_info.particle_one       = particle_one;
              contact_info.particle_two       = particle_two;

              particle_one_contact_list->insert(
                {particle_two_id, contact_info});
            }
        }
    }
}

template class PPFineSearch<2>;
template class PPFineSearch<3>;
