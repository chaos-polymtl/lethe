#include <dem/locate_local_particles.h>

using namespace dealii;

template <int dim>
void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &particle_container,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    &pfw_pairs_in_contact,
  std::vector<
    std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
             std::unordered_map<types::particle_index,
                                particle_wall_contact_info_struct<dim>>,
             dem_data_containers::cut_cell_comparison<dim>>>
    &particle_floating_mesh_in_contact,
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<dim>>
    &particle_points_in_contact,
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<dim>>
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
  update_particle_wall_contact_container_iterators<dim>(pfw_pairs_in_contact,
                                                        particle_container);

  // Calling the same function for floating mesh
  update_particle_floating_wall_contact_container_iterators<dim>(
    particle_floating_mesh_in_contact, particle_container);

  update_particle_point_line_contact_container_iterators<dim>(
    particle_points_in_contact, particle_lines_in_contact, particle_container);
}

template void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<2> &particle_handler,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<2>>
    &particle_container,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<2>>>
    &ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<2>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<2>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<2>>>
    &pfw_pairs_in_contact,
  std::vector<std::map<typename Triangulation<1, 2>::active_cell_iterator,
                       std::unordered_map<types::particle_index,
                                          particle_wall_contact_info_struct<2>>,
                       dem_data_containers::cut_cell_comparison<2>>>
    &particle_floating_mesh_in_contact,
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<2>>
    &particle_points_in_contact,
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<2>>
    &particle_lines_in_contact);

template void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<3> &particle_handler,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<3>>
    &particle_container,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<3>>>
    &ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<3>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<3>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<3>>>
    &pfw_pairs_in_contact,
  std::vector<std::map<typename Triangulation<2, 3>::active_cell_iterator,
                       std::unordered_map<types::particle_index,
                                          particle_wall_contact_info_struct<3>>,
                       dem_data_containers::cut_cell_comparison<3>>>
    &particle_floating_mesh_in_contact,
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<3>>
    &particle_points_in_contact,
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<3>>
    &particle_lines_in_contact);
