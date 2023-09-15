#include <core/tensors_and_points_dimension_manipulation.h>

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
typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
ParticlePointLineFineSearch<dim>::particle_point_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold)
{
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<dim>>
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
      auto     particle            = std::get<0>(*pair_candidates);
      auto     particle_properties = particle->get_properties();
      Point<3> vertex_location_3d;

      if constexpr (dim == 3)
        vertex_location_3d = std::get<1>(*pair_candidates);

      if constexpr (dim == 2)
        vertex_location_3d = point_nd_to_3d(std::get<1>(*pair_candidates));

      Point<3> particle_location_3d;
      if constexpr (dim == 3)
        particle_location_3d = particle->get_location();

      if constexpr (dim == 2)
        particle_location_3d = point_nd_to_3d(particle->get_location());

      // Calculation of the square_distance between the particle and boundary
      // vertex
      const double square_distance =
        ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
        vertex_location_3d.distance_square(particle_location_3d);

      // If the distance is larger than neighberhood threshold, then the
      // particle-point pair are in contact
      if (square_distance > neighborhood_threshold)
        {
          // Creating a sample from the particle_point_line_contact_info_struct
          // and adding contact info to the sample
          particle_point_line_contact_info_struct<dim> contact_info;
          contact_info.particle  = particle;
          contact_info.point_one = vertex_location_3d;

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
typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
ParticlePointLineFineSearch<dim>::particle_line_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold)
{
  std::unordered_map<types::particle_index,
                     particle_point_line_contact_info_struct<dim>>
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
      auto     particle            = std::get<0>(*pair_candidates);
      auto     particle_properties = particle->get_properties();
      Point<3> vertex_one_location_3d;
      Point<3> vertex_two_location_3d;


      if constexpr (dim == 3)
        {
          vertex_one_location_3d = std::get<1>(*pair_candidates);
          vertex_two_location_3d = std::get<2>(*pair_candidates);
        }

      if constexpr (dim == 2)
        {
          vertex_one_location_3d =
            point_nd_to_3d(std::get<1>(*pair_candidates));
          vertex_two_location_3d =
            point_nd_to_3d(std::get<2>(*pair_candidates));
        }

      Point<3> particle_location_3d;
      if constexpr (dim == 3)
        particle_location_3d = particle->get_location();

      if constexpr (dim == 2)
        particle_location_3d = point_nd_to_3d(particle->get_location());

      // For finding the particle-line distance, the projection of the particle
      // on the line should be obtained
      Point<3> projection = find_projection_point(particle_location_3d,
                                                  vertex_one_location_3d,
                                                  vertex_two_location_3d);

      // Calculation of the distance between the particle and boundary line
      const double square_distance =
        ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
        projection.distance_square(particle_location_3d);

      // If the distance is positive, then the particle-line pair are in
      // contact
      if (square_distance > neighborhood_threshold)
        {
          // Creating a sample from the particle_point_line_contact_info_struct
          // and adding contact info to the sample
          particle_point_line_contact_info_struct<dim> contact_info;
          contact_info.particle  = particle;
          contact_info.point_one = vertex_one_location_3d;
          contact_info.point_two = vertex_two_location_3d;

          particle_line_pairs_in_contact.insert(
            {particle->get_id(), contact_info});
        }
    }
  return particle_line_pairs_in_contact;
}

template <int dim>
Point<3>
ParticlePointLineFineSearch<dim>::find_projection_point(const Point<3> &point_p,
                                                        const Point<3> &point_a,
                                                        const Point<3> &point_b)
{
  Tensor<1, 3> vector_ab = point_b - point_a;
  Tensor<1, 3> vector_ap = point_p - point_a;

  Point<3> projection =
    point_a + ((vector_ap * vector_ab) / (vector_ab * vector_ab)) * vector_ab;

  return projection;
}

template class ParticlePointLineFineSearch<2>;
template class ParticlePointLineFineSearch<3>;
