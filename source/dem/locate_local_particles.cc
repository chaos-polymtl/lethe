#include <dem/data_containers.h>
#include <dem/locate_local_particles.h>

using namespace dealii;

template <int dim>
void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
    &ghost_adjacent_particles,
  typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_floating_wall_pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
    &particle_floating_mesh_pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
    &particle_points_in_contact,
  typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
    &particle_lines_in_contact)
{
  // Update the iterators to local particles in a map of particles
  update_particle_container<dim>(particle_container, &particle_handler);

  // Update contact containers for local particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::local_particle_particle>(local_adjacent_particles,
                                          particle_container);

  // Update contact containers for local particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::ghost_particle_particle>(ghost_adjacent_particles,
                                          particle_container);

  // Update contact containers for particle-wall pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    ContactType::particle_wall>(particle_wall_pairs_in_contact,
                                particle_container);

  // Update contact containers for particle-floating wall pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    ContactType::particle_floating_wall>(
    particle_floating_wall_pairs_in_contact, particle_container);

  // Update contact containers for particle-floating mesh pairs in contact for
  // every solid objects of mesh
  for (unsigned int solid_counter = 0;
       solid_counter < particle_floating_mesh_pairs_in_contact.size();
       ++solid_counter)
    {
      update_contact_container_iterators<
        dim,
        typename DEM::dem_data_structures<
          dim>::particle_floating_wall_from_mesh_in_contact,
        ContactType::particle_floating_mesh>(
        particle_floating_mesh_pairs_in_contact[solid_counter],
        particle_container);
    }

  // Update contact containers for particle-line pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_point_line_contact_info,
    ContactType::particle_point>(particle_lines_in_contact, particle_container);

  // Update contact containers for particle-point pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_point_line_contact_info,
    ContactType::particle_point>(particle_points_in_contact,
                                 particle_container);
}

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
    &particle_container)
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
                              ContactType::ghost_particle_particle)
                {
                  // For particle-particle contacts, iterators of both particles
                  // must be updated
                  adjacent_map_iterator->second.particle_one =
                    particle_container[particle_id];

                  unsigned int particle_two_id = adjacent_map_iterator->first;
                  adjacent_map_iterator->second.particle_two =
                    particle_container[particle_two_id];
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

template void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<2> &particle_handler,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &ghost_adjacent_particles,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &particle_floating_wall_pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_floating_mesh_in_contact
    &particle_floating_mesh_pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &particle_points_in_contact,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &particle_lines_in_contact);

template void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<3> &particle_handler,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &ghost_adjacent_particles,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &particle_floating_wall_pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_floating_mesh_in_contact
    &particle_floating_mesh_pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &particle_points_in_contact,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &particle_lines_in_contact);


// Particle container
template void update_particle_container(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &                                  particle_container,
  const Particles::ParticleHandler<2> *particle_handler);

template void update_particle_container(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &                                  particle_container,
  const Particles::ParticleHandler<3> *particle_handler);

// Local particle-particle contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Ghost particle-particle contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-wall contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact,
  ContactType::particle_wall>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact,
  ContactType::particle_wall>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-floating wall contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_wall_in_contact,
  ContactType::particle_floating_wall>(
  typename DEM::dem_data_structures<2>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_wall_in_contact,
  ContactType::particle_floating_wall>(
  typename DEM::dem_data_structures<3>::particle_wall_in_contact
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-floating mesh contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_in_contact,
  ContactType::particle_floating_mesh>(
  typename DEM::dem_data_structures<
    2>::particle_floating_wall_from_mesh_in_contact &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_in_contact,
  ContactType::particle_floating_mesh>(
  typename DEM::dem_data_structures<
    3>::particle_floating_wall_from_mesh_in_contact &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-point contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info,
  ContactType::particle_point>(
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info,
  ContactType::particle_point>(
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-line contact container
template void update_contact_container_iterators<
  2,
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info,
  ContactType::particle_line>(
  typename DEM::dem_data_structures<2>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void update_contact_container_iterators<
  3,
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info,
  ContactType::particle_line>(
  typename DEM::dem_data_structures<3>::particle_point_line_contact_info
    &pairs_in_contact,
  typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);
