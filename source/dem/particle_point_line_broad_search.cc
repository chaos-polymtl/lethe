#include <dem/particle_point_line_broad_search.h>

using namespace dealii;

// Constructor
template <int dim>
ParticlePointLineBroadSearch<dim>::ParticlePointLineBroadSearch()
{}

// This function finds all the particle-point contact candidates
template <int dim>
typename DEM::dem_data_structures<dim>::particle_point_candidates
ParticlePointLineBroadSearch<dim>::find_particle_point_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
    &boundary_cells_with_points)
{
  std::unordered_map<types::particle_index,
                     std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
    particle_point_contact_candidates;

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
  return particle_point_contact_candidates;
}

template <int dim>
typename DEM::dem_data_structures<dim>::particle_point_candidates
ParticlePointLineBroadSearch<dim>::find_particle_point_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<
    std::string,
    std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
    &                         boundary_cells_with_points,
  const DisableContacts<dim> &disable_contacts_object)
{
  std::unordered_map<types::particle_index,
                     std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
    particle_point_contact_candidates;

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
        disable_contacts_object.check_cell_mobility(cell_with_boundary_point);
      if (main_cell_mobility_status != DisableContacts<dim>::mobile)
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
  return particle_point_contact_candidates;
}

// This function finds all the particle-line contact candidates
template <int dim>
typename DEM::dem_data_structures<dim>::particle_line_candidates
ParticlePointLineBroadSearch<dim>::find_particle_line_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<
    std::string,
    std::tuple<typename Triangulation<dim>::active_cell_iterator,
               Point<dim>,
               Point<dim>>> &boundary_cells_with_lines)
{
  std::unordered_map<
    types::particle_index,
    std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
    particle_line_contact_candidates;

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
      auto cells_with_boundary_lines_information = &map_iterator->second;

      // Getting the cell and the locations of boundary vertices (beginning and
      // ending points of the boundary line) as local variables
      auto cell_with_boundary_line =
        std::get<0>(*cells_with_boundary_lines_information);

      Point<dim> first_vertex_location =
        std::get<1>(*cells_with_boundary_lines_information);
      Point<dim> second_vertex_location =
        std::get<2>(*cells_with_boundary_lines_information);

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell =
          particle_handler.particles_in_cell(cell_with_boundary_line);

      for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Making the tuple (particle, beginning and ending points of the
          // boundary line) and adding it to the
          // particle_line_contact_candidates map. This map is the output of
          // this function
          particle_line_contact_candidates.insert(
            {contact_candidate_counter,
             std::make_tuple(particles_in_cell_iterator,
                             first_vertex_location,
                             second_vertex_location)});
          ++contact_candidate_counter;
        }
    }
  return particle_line_contact_candidates;
}

template <int dim>
typename DEM::dem_data_structures<dim>::particle_line_candidates
ParticlePointLineBroadSearch<dim>::find_particle_line_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<
    std::string,
    std::tuple<typename Triangulation<dim>::active_cell_iterator,
               Point<dim>,
               Point<dim>>> & boundary_cells_with_lines,
  const DisableContacts<dim> &disable_contacts_object)
{
  std::unordered_map<
    types::particle_index,
    std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
    particle_line_contact_candidates;

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
      auto cells_with_boundary_lines_information = &map_iterator->second;

      // Getting the cell and the locations of boundary vertices (beginning and
      // ending points of the boundary line) as local variables
      auto cell_with_boundary_line =
        std::get<0>(*cells_with_boundary_lines_information);

      // If main cell has status other than mobile, skip to next cell
      unsigned int main_cell_mobility_status =
        disable_contacts_object.check_cell_mobility(cell_with_boundary_line);
      if (main_cell_mobility_status != DisableContacts<dim>::mobile)
        continue;

      Point<dim> first_vertex_location =
        std::get<1>(*cells_with_boundary_lines_information);
      Point<dim> second_vertex_location =
        std::get<2>(*cells_with_boundary_lines_information);

      // Finding particles located in the corresponding cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell =
          particle_handler.particles_in_cell(cell_with_boundary_line);

      for (auto particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Making the tuple (particle, beginning and ending points of the
          // boundary line) and adding it to the
          // particle_line_contact_candidates map. This map is the output of
          // this function
          particle_line_contact_candidates.insert(
            {contact_candidate_counter,
             std::make_tuple(particles_in_cell_iterator,
                             first_vertex_location,
                             second_vertex_location)});
          ++contact_candidate_counter;
        }
    }
  return particle_line_contact_candidates;
}

template class ParticlePointLineBroadSearch<2>;
template class ParticlePointLineBroadSearch<3>;
