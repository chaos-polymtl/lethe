#include <dem/particle_point_line_fine_search.h>

using namespace dealii;

template <int dim>
ParticlePointLineFineSearch<dim>::ParticlePointLineFineSearch()
{}

// In this function, the output of particle-point broad search is investigated
// to see if the pairs are in contact or not. If they are in contact, the normal
// overlap, normal vector and normal relative velocity of the contact are
// calculated. The output of this function is used for calculation of the
// contact force
template <int dim>
std::unordered_map<int, particle_point_line_contact_info_struct<dim>>
ParticlePointLineFineSearch<dim>::particle_point_fine_search(
  const std::
    unordered_map<int, std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
      &         particle_point_contact_candidates,
  const double &neighborhood_threshold)
{
  std::unordered_map<int, particle_point_line_contact_info_struct<dim>>
    particle_point_pairs_in_contact;

  // Iterating over contact candidates from broad search. If a particle-point
  // pair is in contact (distance > 0) it is inserted into the output of this
  // function (particle_point_pairs_in_contact)
  for (auto contact_pair_candidates_iterator =
         particle_point_contact_candidates.begin();
       contact_pair_candidates_iterator !=
       particle_point_contact_candidates.end();
       ++contact_pair_candidates_iterator)
    {
      // Get the value of the map (pair candidate) from the
      // contact_pair_candidates_iterator
      auto pair_candidates = &contact_pair_candidates_iterator->second;

      // Get the particle, particle properties and boundary vertex location once
      // to improve efficiency
      auto             particle            = std::get<0>(*pair_candidates);
      auto             particle_properties = particle->get_properties();
      const Point<dim> vertex_location     = std::get<1>(*pair_candidates);

      // Calculation of the square_distance between the particle and boundary
      // vertex
      const double square_distance =
        ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
        vertex_location.distance_square(particle->get_location());

      // If the distance is larger than neighberhood threshold, then the
      // particle-point pair are in contact
      if (square_distance > neighborhood_threshold)
        {
          // Creating a sample from the particle_point_line_contact_info_struct
          // and adding contact info to the sample
          particle_point_line_contact_info_struct<dim> contact_info;
          contact_info.particle  = particle;
          contact_info.point_one = vertex_location;

          particle_point_pairs_in_contact.insert(
            {particle->get_id(), contact_info});
        }
    }
  return particle_point_pairs_in_contact;
}

// In this function, the output of particle-line broad search is investigated
// to see if the pairs are in contact or not. If they are in contact, the normal
// overlap, normal vector and normal relative velocity of the contact are
// calculated. The output of this function is used for calculation of the
// contact force
template <int dim>
std::unordered_map<int, particle_point_line_contact_info_struct<dim>>
ParticlePointLineFineSearch<dim>::particle_line_fine_search(
  const std::unordered_map<
    int,
    std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
    &           particle_line_contact_candidates,
  const double &neighborhood_threshold)
{
  std::unordered_map<int, particle_point_line_contact_info_struct<dim>>
    particle_line_pairs_in_contact;

  // Iterating over contact candidates from broad search. If a particle-line
  // pair is in contact (distance > 0) it is inserted into the output of this
  // function (particle_line_pairs_in_contact)
  for (auto contact_pair_candidates_iterator =
         particle_line_contact_candidates.begin();
       contact_pair_candidates_iterator !=
       particle_line_contact_candidates.end();
       ++contact_pair_candidates_iterator)
    {
      // Get the value of the map (pair candidate) from the
      // contact_pair_candidates_iterator
      auto pair_candidates = &contact_pair_candidates_iterator->second;

      // Get the particle, particle properties and the locations of beginning
      // and ending vertices of the boundary line once to improve efficiency
      auto             particle            = std::get<0>(*pair_candidates);
      auto             particle_properties = particle->get_properties();
      const Point<dim> vertex_one_location = std::get<1>(*pair_candidates);
      const Point<dim> vertex_two_location = std::get<2>(*pair_candidates);

      // For finding the particle-line distance, the projection of the particle
      // on the line should be obtained
      Point<dim> projection = find_projection_point(particle->get_location(),
                                                    vertex_one_location,
                                                    vertex_two_location);

      // Calculation of the distance between the particle and boundary line
      const double square_distance =
        ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
        projection.distance_square(particle->get_location());

      // If the distance is positive, then the particle-line pair are in
      // contact
      if (square_distance > neighborhood_threshold)
        {
          // Creating a sample from the particle_point_line_contact_info_struct
          // and adding contact info to the sample
          particle_point_line_contact_info_struct<dim> contact_info;
          contact_info.particle  = particle;
          contact_info.point_one = vertex_one_location;
          contact_info.point_two = vertex_two_location;

          particle_line_pairs_in_contact.insert(
            {particle->get_id(), contact_info});
        }
    }
  return particle_line_pairs_in_contact;
}

template <int dim>
Point<dim>
ParticlePointLineFineSearch<dim>::find_projection_point(
  const Point<dim> &point_p,
  const Point<dim> &point_a,
  const Point<dim> &point_b)
{
  Tensor<1, dim> vector_ab = point_b - point_a;
  Tensor<1, dim> vector_ap = point_p - point_a;

  Point<dim> projection =
    point_a + ((vector_ap * vector_ab) / (vector_ab * vector_ab)) * vector_ab;

  return projection;
}

template class ParticlePointLineFineSearch<2>;
template class ParticlePointLineFineSearch<3>;
