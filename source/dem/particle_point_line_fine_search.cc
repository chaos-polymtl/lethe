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
std::map<int, particle_point_line_contact_info_struct<dim>>
ParticlePointLineFineSearch<dim>::Particle_Point_Fine_Search(
  const std::map<int, std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
    &particle_point_contact_candidates)
{
  std::map<int, particle_point_line_contact_info_struct<dim>>
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
      auto       particle            = std::get<0>(*pair_candidates);
      auto       particle_properties = particle->get_properties();
      Point<dim> vertex_location     = std::get<1>(*pair_candidates);

      // Calculation of the distance between the particle and boundary vertex
      double distance = ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
                        vertex_location.distance(particle->get_location());

      // If the distance is positive, then the particle-point pair are in
      // contact
      if (distance > 0)
        {
          // Calculation of normal vector, particle velocity and contact
          // relative velocity
          Tensor<1, dim> point_to_particle_vector =
            particle->get_location() - vertex_location;
          Tensor<1, dim> normal_vector =
            point_to_particle_vector / point_to_particle_vector.norm();

          Tensor<1, dim> particle_velocity;
          particle_velocity[0] = particle_properties[DEM::PropertiesIndex::v_x];
          particle_velocity[1] = particle_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3)
            {
              particle_velocity[2] =
                particle_properties[DEM::PropertiesIndex::v_z];
            }

          Tensor<1, dim> particle_omega;
          particle_omega[0] =
            particle_properties[DEM::PropertiesIndex::omega_x];
          particle_omega[1] =
            particle_properties[DEM::PropertiesIndex::omega_y];
          if (dim == 3)
            {
              particle_omega[2] =
                particle_properties[DEM::PropertiesIndex::omega_z];
            }

          // Defining relative contact velocity
          Tensor<1, dim> contact_relative_velocity;
          if (dim == 3)
            {
              contact_relative_velocity =
                particle_velocity +
                cross_product_3d(
                  (((particle_properties[DEM::PropertiesIndex::dp]) / 2) *
                   particle_omega),
                  normal_vector);
            }

          if (dim == 2)
            {
              contact_relative_velocity = particle_velocity;
            }

          // Calculation of normal relative velocity
          double normal_relative_velocity_value =
            contact_relative_velocity * normal_vector;

          // Creating a sample from the particle_point_line_contact_info_struct
          // and adding contact info to the sample
          particle_point_line_contact_info_struct<dim> contact_info;
          contact_info.particle       = particle;
          contact_info.normal_overlap = distance;
          contact_info.normal_vector  = normal_vector;
          contact_info.normal_relative_velocity =
            normal_relative_velocity_value;

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
std::map<int, particle_point_line_contact_info_struct<dim>>
ParticlePointLineFineSearch<dim>::Particle_Line_Fine_Search(
  const std::map<
    int,
    std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
    &particle_line_contact_candidates)
{
  std::map<int, particle_point_line_contact_info_struct<dim>>
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
      auto       particle            = std::get<0>(*pair_candidates);
      auto       particle_properties = particle->get_properties();
      Point<dim> vertex_one_location = std::get<1>(*pair_candidates);
      Point<dim> vertex_two_location = std::get<2>(*pair_candidates);

      // For finding the particle-line distance, the projection of the particle
      // on the line should be obtained
      Point<dim> projection = find_projection_point(particle->get_location(),
                                                    vertex_one_location,
                                                    vertex_two_location);

      // Calculation of the distance between the particle and boundary line
      double distance = ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
                        projection.distance(particle->get_location());

      // If the distance is positive, then the particle-line pair are in
      // contact
      if (distance > 0)
        {
          // Calculation of normal vector, particle velocity and contact
          // relative velocity
          Tensor<1, dim> point_to_particle_vector =
            particle->get_location() - projection;
          Tensor<1, dim> normal_vector =
            point_to_particle_vector / point_to_particle_vector.norm();

          Tensor<1, dim> particle_velocity;
          particle_velocity[0] = particle_properties[DEM::PropertiesIndex::v_x];
          particle_velocity[1] = particle_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3)
            {
              particle_velocity[2] =
                particle_properties[DEM::PropertiesIndex::v_z];
            }

          Tensor<1, dim> particle_omega;
          particle_omega[0] =
            particle_properties[DEM::PropertiesIndex::omega_x];
          particle_omega[1] =
            particle_properties[DEM::PropertiesIndex::omega_y];
          if (dim == 3)
            {
              particle_omega[2] =
                particle_properties[DEM::PropertiesIndex::omega_z];
            }

          // Defining relative contact velocity
          Tensor<1, dim> contact_relative_velocity;
          if (dim == 3)
            {
              contact_relative_velocity =
                particle_velocity +
                cross_product_3d(
                  (((particle_properties[DEM::PropertiesIndex::dp]) / 2) *
                   particle_omega),
                  normal_vector);
            }

          if (dim == 2)
            {
              contact_relative_velocity = particle_velocity;
            }

          // Calculation of normal relative velocity
          double normal_relative_velocity_value =
            contact_relative_velocity * normal_vector;

          // Creating a sample from the particle_point_line_contact_info_struct
          // and adding contact info to the sample
          particle_point_line_contact_info_struct<dim> contact_info;
          contact_info.particle       = particle;
          contact_info.normal_overlap = distance;
          contact_info.normal_vector  = normal_vector;
          contact_info.normal_relative_velocity =
            normal_relative_velocity_value;

          particle_line_pairs_in_contact.insert(
            {particle->get_id(), contact_info});
        }
    }
  return particle_line_pairs_in_contact;
}

template <int dim>
Point<dim>
ParticlePointLineFineSearch<dim>::find_projection_point(Point<dim> point_p,
                                                        Point<dim> point_a,
                                                        Point<dim> point_b)
{
  Tensor<1, dim> vector_ab = point_b - point_a;
  Tensor<1, dim> vector_ap = point_p - point_a;

  Point<dim> projection =
    point_a + ((vector_ap * vector_ab) / (vector_ab * vector_ab)) * vector_ab;

  return projection;
}

template class ParticlePointLineFineSearch<2>;
template class ParticlePointLineFineSearch<3>;
