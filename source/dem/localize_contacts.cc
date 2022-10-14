#include <dem/contact_type.h>
#include <dem/localize_contacts.h>

using namespace dealii;

template <int dim>
void
localize_contacts(
  typename dem_data_containers::dem_data_structures<
    dim>::adjacent_particle_pairs &local_adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    dim>::adjacent_particle_pairs &ghost_adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_wall_in_contact &particle_wall_pairs_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_wall_in_contact &particle_floating_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_mesh_in_contact &particle_floating_mesh_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_particle_candidates &local_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_particle_candidates &ghost_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_wall_candidates &particle_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_wall_candidates
    &particle_floating_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates)
{
  // Update particle-particle contacts in local_adjacent_particles of fine
  // search step with local_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_containers::dem_data_structures<
      dim>::adjacent_particle_pairs,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_particle_candidates>(local_adjacent_particles,
                                          local_contact_pair_candidates,
                                          ContactType::local_particle_particle);

  // Update particle-particle contacts in global_adjacent_particles of fine
  // search step with global_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_containers::dem_data_structures<
      dim>::adjacent_particle_pairs,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_particle_candidates>(ghost_adjacent_particles,
                                          ghost_contact_pair_candidates,
                                          ContactType::ghost_particle_particle);

  // Update particle-wall contacts in particle_wall_pairs_in_contact of fine
  // search step with particle_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_in_contact,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_candidates>(particle_wall_pairs_in_contact,
                                      particle_wall_contact_candidates,
                                      ContactType::particle_wall);

  // Update particle-floating wall contacts in particle_floating_wall_in_contact
  // of fine search step with particle_floating_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_in_contact,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_floating_wall_candidates>(
    particle_floating_wall_in_contact,
    particle_floating_wall_contact_candidates,
    ContactType::particle_floating_wall);

  // Update particle-floating mesh contacts in particle_floating_mesh_in_contact
  // of fine search step with particle_floating_mesh_contact_candidates
  for (unsigned int solid_counter = 0;
       solid_counter < particle_floating_mesh_in_contact.size();
       ++solid_counter)
    {
      update_fine_search_candidates<
        dim,
        std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
                 std::unordered_map<types::particle_index,
                                    particle_wall_contact_info_struct<dim>>,
                 typename dem_data_containers::cut_cell_comparison<dim>>,
        std::map<typename Triangulation<dim - 1, dim>::active_cell_iterator,
                 std::unordered_map<types::particle_index,
                                    Particles::ParticleIterator<dim>>,
                 typename dem_data_containers::cut_cell_comparison<dim>>>(
        particle_floating_mesh_in_contact[solid_counter],
        particle_floating_mesh_contact_candidates[solid_counter],
        ContactType::particle_floating_mesh);
    }
}

template void localize_contacts<2>(
  typename dem_data_containers::dem_data_structures<2>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename dem_data_containers::dem_data_structures<2>::adjacent_particle_pairs
    &ghost_adjacent_particles,
  typename dem_data_containers::dem_data_structures<2>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename dem_data_containers::dem_data_structures<2>::particle_wall_in_contact
    &particle_floating_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    2>::particle_floating_mesh_in_contact &particle_floating_mesh_in_contact,
  typename dem_data_containers::dem_data_structures<
    2>::particle_particle_candidates &local_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<
    2>::particle_particle_candidates &ghost_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<2>::particle_wall_candidates
    &particle_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    2>::particle_floating_wall_candidates
    &particle_floating_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    2>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates);

template void localize_contacts<3>(
  typename dem_data_containers::dem_data_structures<3>::adjacent_particle_pairs
    &local_adjacent_particles,
  typename dem_data_containers::dem_data_structures<3>::adjacent_particle_pairs
    &ghost_adjacent_particles,
  typename dem_data_containers::dem_data_structures<3>::particle_wall_in_contact
    &particle_wall_pairs_in_contact,
  typename dem_data_containers::dem_data_structures<3>::particle_wall_in_contact
    &particle_floating_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    3>::particle_floating_mesh_in_contact &particle_floating_mesh_in_contact,
  typename dem_data_containers::dem_data_structures<
    3>::particle_particle_candidates &local_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<
    3>::particle_particle_candidates &ghost_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<3>::particle_wall_candidates
    &particle_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    3>::particle_floating_wall_candidates
    &particle_floating_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    3>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates);