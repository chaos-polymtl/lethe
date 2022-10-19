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
  update_particle_container<dim>(particle_container, &particle_handler);

  update_local_particle_particle_contact_container_iterators<dim>(
    local_adjacent_particles, particle_container);

  update_ghost_particle_particle_contact_container_iterators<dim>(
    ghost_adjacent_particles, particle_container);

  update_particle_wall_contact_container_iterators<dim>(
    particle_wall_pairs_in_contact, particle_container);

  // Calling the same function for floating walls
  update_particle_wall_contact_container_iterators<dim>(
    particle_floating_wall_pairs_in_contact, particle_container);

  // Calling the same function for floating mesh
  update_particle_floating_mesh_contact_container_iterators<dim>(
    particle_floating_mesh_pairs_in_contact, particle_container);

  update_particle_point_line_contact_container_iterators<dim>(
    particle_points_in_contact, particle_lines_in_contact, particle_container);
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
