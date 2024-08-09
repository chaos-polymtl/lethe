#include <dem/particle_point_line_broad_search.h>

using namespace dealii;
using namespace DEM;

// This function finds all the particle-point contact candidates
template <int dim>
void
find_particle_point_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<dim>::particle_point_candidates
    &particle_point_contact_candidates)
{
  // Clear the candidate map
  particle_point_contact_candidates.clear();

  // Defining and resetting a local particle-point candidate counter. This is
  // used as a key to the output map
  int contact_candidate_counter = 0;

  // Iterating over the boundary_cells_with_points which is the output of
  // the find_boundary_cells_information class.
  // This vector contains all the required information of the boundary
  // cells with points. In this loop we find the particles located in each of
  // these boundary cells with points
  for (auto map_iterator = boundary_cells_with_points.begin();
       map_iterator != boundary_cells_with_points.end();
       ++map_iterator)
    {
      auto cells_with_boundary_points_information = &map_iterator->second;

      // Getting the cell and boundary vertex as local variables
      auto cell_with_boundary_point =
        cells_with_boundary_points_information->first;

      Point<dim> vertex_location =
        cells_with_boundary_points_information->second;

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell =
          particle_handler.particles_in_cell(cell_with_boundary_point);

      for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Making the pair and adding it to the
          // particle_point_contact_candidates map. This map is the output of
          // this function
          particle_point_contact_candidates.insert(
            {contact_candidate_counter,
             std::make_pair(particles_in_cell_iterator, vertex_location)});
          ++contact_candidate_counter;
        }
    }
}

template <int dim>
void
find_particle_point_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<dim>::particle_point_candidates
                                    &particle_point_contact_candidates,
  const AdaptiveSparseContacts<dim> &sparse_contacts_object)
{
  // Clear the candidate map
  particle_point_contact_candidates.clear();

  // Defining and resetting a local particle-point candidate counter. This is
  // used as a key to the output map
  int contact_candidate_counter = 0;

  // Iterating over the boundary_cells_with_points which is the output of
  // the find_boundary_cells_information class.
  // This vector contains all the required information of the boundary
  // cells with points. In this loop we find the particles located in each of
  // these boundary cells with points
  for (auto map_iterator = boundary_cells_with_points.begin();
       map_iterator != boundary_cells_with_points.end();
       ++map_iterator)
    {
      auto cells_with_boundary_points_information = &map_iterator->second;

      // Getting the cell and boundary vertex as local variables
      auto cell_with_boundary_point =
        cells_with_boundary_points_information->first;

      // If main cell has status other than mobile, skip to next cell
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(cell_with_boundary_point);
      if (main_cell_mobility_status != AdaptiveSparseContacts<dim>::mobile)
        continue;

      Point<dim> vertex_location =
        cells_with_boundary_points_information->second;

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell =
          particle_handler.particles_in_cell(cell_with_boundary_point);

      for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Making the pair and adding it to the
          // particle_point_contact_candidates map. This map is the output of
          // this function
          particle_point_contact_candidates.insert(
            {contact_candidate_counter,
             std::make_pair(particles_in_cell_iterator, vertex_location)});
          ++contact_candidate_counter;
        }
    }
}

// This function finds all the particle-line contact candidates
template <int dim>
void
find_particle_line_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<dim>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<dim>::particle_line_candidates
    &particle_line_contact_candidates)
{
  // Clear the candidates map
  particle_line_contact_candidates.clear();

  // Defining and resetting a local particle-line candidate counter. This is
  // used as a key to the output map
  unsigned int contact_candidate_counter = 0;

  // Iterating over the  boundary_cells_with_lines which is the output of
  // the find_boundary_cells_information class.
  // This vector contains all the required information of the boundary
  // cells with lines. In this loop we find the particles located in each of
  // these boundary cells with lines
  for (auto map_iterator = boundary_cells_with_lines.begin();
       map_iterator != boundary_cells_with_lines.end();
       ++map_iterator)
    {
      const cell_line_info<dim> &cells_with_boundary_lines_info =
        map_iterator->second;

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(
          cells_with_boundary_lines_info.cell);

      for (auto particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Adding the particle line contact info to the
          // particle_line_contact_candidates map.
          particle_line_contact_candidates.emplace(
            contact_candidate_counter,
            particle_line_contact_info<dim>{
              particles_in_cell_iterator,
              cells_with_boundary_lines_info.point_one,
              cells_with_boundary_lines_info.point_two});
          ++contact_candidate_counter;
        }
    }
}

template <int dim>
void
find_particle_line_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<dim>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<dim>::particle_line_candidates
                                    &particle_line_contact_candidates,
  const AdaptiveSparseContacts<dim> &sparse_contacts_object)
{
  // Clear the candidates map
  particle_line_contact_candidates.clear();

  // Defining and resetting a local particle-line candidate counter. This is
  // used as a key to the output map
  unsigned int contact_candidate_counter = 0;

  // Iterating over the  boundary_cells_with_lines which is the output of
  // the find_boundary_cells_information class.
  // This vector contains all the required information of the boundary
  // cells with lines. In this loop we find the particles located in each of
  // these boundary cells with lines
  for (auto map_iterator = boundary_cells_with_lines.begin();
       map_iterator != boundary_cells_with_lines.end();
       ++map_iterator)
    {
      const cell_line_info<dim> &cells_with_boundary_lines_info =
        map_iterator->second;

      // If main cell has status other than mobile, skip to next cell
      unsigned int main_cell_mobility_status =
        sparse_contacts_object.check_cell_mobility(
          cells_with_boundary_lines_info.cell);
      if (main_cell_mobility_status != AdaptiveSparseContacts<dim>::mobile)
        continue;

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(
          cells_with_boundary_lines_info.cell);

      for (auto particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Adding the particle line contact info to the
          // particle_line_contact_candidates map.
          particle_line_contact_candidates.emplace(
            contact_candidate_counter,
            particle_line_contact_info<dim>{
              particles_in_cell_iterator,
              cells_with_boundary_lines_info.point_one,
              cells_with_boundary_lines_info.point_two});
          ++contact_candidate_counter;
        }
    }
}

template void
find_particle_line_contact_pairs<2>(
  const Particles::ParticleHandler<2> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<2>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<2>::particle_line_candidates
    &particle_line_contact_candidates);

template void
find_particle_line_contact_pairs<3>(
  const Particles::ParticleHandler<3> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<3>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<3>::particle_line_candidates
    &particle_line_contact_candidates);

template void
find_particle_line_contact_pairs<2>(
  const Particles::ParticleHandler<2> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<2>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<2>::particle_line_candidates
                                  &particle_line_contact_candidates,
  const AdaptiveSparseContacts<2> &sparse_contacts_object);

template void
find_particle_line_contact_pairs<3>(
  const Particles::ParticleHandler<3> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<3>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<3>::particle_line_candidates
                                  &particle_line_contact_candidates,
  const AdaptiveSparseContacts<3> &sparse_contacts_object);

template void
find_particle_point_contact_pairs<2>(
  const Particles::ParticleHandler<2> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<2>::active_cell_iterator, Point<2>>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<2>::particle_point_candidates
    &particle_point_contact_candidates);

template void
find_particle_point_contact_pairs<3>(
  const Particles::ParticleHandler<3> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<3>::active_cell_iterator, Point<3>>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<3>::particle_point_candidates
    &particle_point_contact_candidates);

template void
find_particle_point_contact_pairs<2>(
  const Particles::ParticleHandler<2> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<2>::active_cell_iterator, Point<2>>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<2>::particle_point_candidates
                                  &particle_point_contact_candidates,
  const AdaptiveSparseContacts<2> &sparse_contacts_object);

template void
find_particle_point_contact_pairs<3>(
  const Particles::ParticleHandler<3> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<3>::active_cell_iterator, Point<3>>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<3>::particle_point_candidates
                                  &particle_point_contact_candidates,
  const AdaptiveSparseContacts<3> &sparse_contacts_object);
