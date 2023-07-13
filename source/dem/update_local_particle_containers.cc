#include <dem/update_local_particle_containers.h>

using namespace dealii;

template <int dim>
void
update_particle_container(
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &                                    particle_container,
  const Particles::ParticleHandler<dim> *particle_handler)
{
  // Clear the current particle container
  particle_container.clear();

  // Update the local particle container with local particles
  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }

  // Update the local particle container with ghost particles
  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }
}

template <int dim, typename pairs_structure, ContactType contact_type>
void
update_contact_container_iterators(
  pairs_structure &pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures)
{
  // Loop over particle-object (particle/wall/line/point) pairs in contact
  for (auto pairs_in_contact_iterator = pairs_in_contact.begin();
       pairs_in_contact_iterator != pairs_in_contact.end();
       ++pairs_in_contact_iterator)
    {
      // Get the adjacent objects content
      auto adjacent_pairs_content = &pairs_in_contact_iterator->second;

      if constexpr (contact_type == ContactType::local_particle_particle ||
                    contact_type == ContactType::ghost_particle_particle ||
                    contact_type ==
                      ContactType::local_periodic_particle_particle ||
                    contact_type ==
                      ContactType::ghost_periodic_particle_particle ||
                    contact_type ==
                      ContactType::ghost_local_periodic_particle_particle ||
                    contact_type == ContactType::particle_wall ||
                    contact_type == ContactType::particle_floating_wall)
        {
          // Get current particle id
          unsigned int particle_id = pairs_in_contact_iterator->first;

          // Loop over all the other objects of contact_type in contact
          for (auto adjacent_map_iterator = adjacent_pairs_content->begin();
               adjacent_map_iterator != adjacent_pairs_content->end();
               ++adjacent_map_iterator)
            {
              if constexpr (contact_type ==
                              ContactType::local_particle_particle ||
                            contact_type ==
                              ContactType::ghost_particle_particle ||
                            contact_type ==
                              ContactType::local_periodic_particle_particle ||
                            contact_type ==
                              ContactType::ghost_periodic_particle_particle ||
                            contact_type ==
                              ContactType::
                                ghost_local_periodic_particle_particle)
                {
                  // For particle-particle contacts, iterators of both particles
                  // must be updated
                  adjacent_map_iterator->second.particle_one =
                    particle_container[particle_id];

                  unsigned int particle_two_id = adjacent_map_iterator->first;
                  adjacent_map_iterator->second.particle_two =
                    particle_container[particle_two_id];

                  if (clear_contact_structures)
                    {
                      adjacent_map_iterator->second.tangential_overlap.clear();
                      adjacent_map_iterator->second.tangential_relative_velocity
                        .clear();
                    }
                }

              if constexpr (contact_type == ContactType::particle_wall ||
                            contact_type == ContactType::particle_floating_wall)
                {
                  // Particle iterator is updated
                  adjacent_map_iterator->second.particle =
                    particle_container[particle_id];
                }
            }
        }

      if constexpr (contact_type == ContactType::particle_floating_mesh)
        {
          for (auto adjacent_map_iterator = adjacent_pairs_content->begin();
               adjacent_map_iterator != adjacent_pairs_content->end();
               ++adjacent_map_iterator)
            {
              // Get the particle in contact with the solid of the floating mesh
              unsigned int particle_id = adjacent_map_iterator->first;

              // Update the particle of the contact
              adjacent_map_iterator->second.particle =
                particle_container[particle_id];

              // Note : particle_floating_wall_from_mesh_in_contact is stored as
              // <cell iterator, <particle id, particle-wall info>> and
              // particle_wall_in_contact as is stored as
              // <particle id, <face id, particle-wall info>>,
              // which restrains to use the same iterator loop logic
            }
        }

      if constexpr (contact_type == ContactType::particle_line ||
                    contact_type == ContactType::particle_point)
        {
          // Get current particle id and particle iterator is updated
          unsigned int particle_id         = pairs_in_contact_iterator->first;
          adjacent_pairs_content->particle = particle_container[particle_id];
        }
    }
}

// Particle container
template void
update_particle_container(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &                                  particle_container,
  const Particles::ParticleHandler<2> *particle_handler);

template void
update_particle_container(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &                                  particle_container,
  const Particles::ParticleHandler<3> *particle_handler);

// Local particle-particle contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Ghost particle-particle contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Local-local particle-particle periodic contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::local_periodic_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::local_periodic_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Local-ghost particle-particle periodic contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_periodic_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_periodic_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Ghost-local particle-particle periodic contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_local_periodic_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_local_periodic_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Particle-wall contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact,
  ContactType::particle_wall>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact,
  ContactType::particle_wall>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Particle-floating wall contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact,
  ContactType::particle_floating_wall>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact,
  ContactType::particle_floating_wall>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Particle-floating mesh contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_in_contact,
  ContactType::particle_floating_mesh>(
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_in_contact &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_in_contact,
  ContactType::particle_floating_mesh>(
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_in_contact &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Particle-point contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info,
  ContactType::particle_point>(
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info,
  ContactType::particle_point>(
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

// Particle-line contact container
template void
update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info,
  ContactType::particle_line>(
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);

template void
update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info,
  ContactType::particle_line>(
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures);
