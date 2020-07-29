#include <dem/pw_fine_search.h>

using namespace dealii;

template <int dim>
PWFineSearch<dim>::PWFineSearch()
{}

template <int dim>
void
PWFineSearch<dim>::pw_Fine_Search(
  std::map<
    std::pair<int, int>,
    std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>
    &pw_contact_pair_candidates,
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    &    pw_pairs_in_contact,
  double dt)
{
  // Iterating over all element of pw_pairs_in_contact vector, i.e. lterating
  // over all the particles. The size of this vector (pw_pairs_in_contact) is
  // equal to the number of particles and element i of this vector corresponds
  // to particle i. Each element of this vector (i.e. each particle) consists of
  // a map container. The key of this map is the id of the boundary which is in
  // contact with the corresponding particle. The contact information are stored
  // in this map
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact.begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact.end();
       ++pw_pairs_in_contact_iterator)
    {
      // Defining element (particle->get_id())th of the pw_pairs_in_contact
      // vector as a local map (pw_contact_map)
      auto pw_contact_map = &pw_pairs_in_contact_iterator->second;

      // Iterating (with defining the iterator contact_pairs_iterator) over
      // elemens of pw_pairs_in_contact_iterator which is a pointer to a map
      for (auto contact_pairs_iterator = pw_contact_map->begin();
           contact_pairs_iterator != pw_contact_map->end();)
        {
          // For each contact, the particle is taken
          // from the iterator and defined as local parameters. Similarly the
          // information tuple is also defined as a local variable
          auto information_tuple   = &contact_pairs_iterator->second;
          auto particle            = information_tuple->particle;
          auto particle_properties = particle->get_properties();

          // Normal vector of the boundary and a point on the boudary are
          // defined as local parameters
          Tensor<1, dim> normal_vector = information_tuple->normal_vector;
          Point<dim> point_on_boundary = information_tuple->point_on_boundary;

          // A vector (point_to_particle_vector) is defined which connects the
          // center of particle to the point_on_boundary. This vector will then
          // be projected on the normal vector of the boundary to obtain the
          // particle-wall distance
          Tensor<1, dim> point_to_particle_vector =
            particle->get_location() - point_on_boundary;

          // Finding the projected vector on the normal vector of the boundary.
          // Here we have used the private function find_projection. Using this
          // projected vector, the particle-wall distance is calculated
          Tensor<1, dim> projected_vector =
            find_projection(point_to_particle_vector, normal_vector);
          double distance =
            ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
            (projected_vector.norm());

          // Check to see if particle-wall pair is still in contact
          if (distance > 0.0)
            {
              // If they are still in contact

              // Using velocity and angular velocity of particle as
              // local vectors
              Tensor<1, dim> particle_velocity;
              particle_velocity[0] =
                particle_properties[DEM::PropertiesIndex::v_x];
              particle_velocity[1] =
                particle_properties[DEM::PropertiesIndex::v_y];
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
              Tensor<1, dim> normal_relative_velocity =
                normal_relative_velocity_value * normal_vector;

              // Calculation of tangential relative velocity
              Tensor<1, dim> tangential_relative_velocity =
                contact_relative_velocity - normal_relative_velocity;

              // Calculation of new tangential_overlap, since this value is
              // history-dependent, it needs the value at previous time-step
              Tensor<1, dim> tangential_overlap =
                information_tuple->tangential_overlap -
                (information_tuple->tangential_overlap * normal_vector) *
                  normal_vector;

              Tensor<1, dim> modified_tangential_overlap;
              if (tangential_overlap.norm() != 0.0)
                {
                  modified_tangential_overlap =
                    (information_tuple->tangential_overlap.norm() /
                     tangential_overlap.norm()) *
                      tangential_overlap +
                    information_tuple->tangential_relative_velocity * dt;
                }
              else
                {
                  modified_tangential_overlap =
                    tangential_overlap +
                    information_tuple->tangential_relative_velocity * dt;
                }

              // Creating a sample from the pw_contact_info_struct and adding
              // contact info to the sample
              pw_contact_info_struct<dim> contact_info;
              contact_info.particle          = particle;
              contact_info.normal_vector     = normal_vector;
              contact_info.point_on_boundary = point_on_boundary;
              contact_info.normal_overlap    = distance;
              contact_info.normal_relative_velocity =
                normal_relative_velocity_value;
              contact_info.tangential_overlap = modified_tangential_overlap;
              contact_info.tangential_relative_velocity =
                tangential_relative_velocity;

              // pw_contact_map[face_id] = contact_info;
              *information_tuple = contact_info;
              ++contact_pairs_iterator;
            }

          // If the particle-wall pair is not in contact anymore (i.e. the
          // contact is finished and distance <= 0), this element should be
          // erased from pw_pairs_in_contact
          else
            {
              //(pw_pairs_in_contact[particle->get_id()]).erase(contact_pairs_iterator->first);
              pw_contact_map->erase(contact_pairs_iterator++);
            }
        }
    }

  // Now iterating over contact candidates from broad search. If a particle-wall
  // pair is in contact (distance > 0) and does not exist in the
  // pw_pairs_in_contact, it is added to the pw_pairs_in_contact
  for (auto map_iterator = pw_contact_pair_candidates.begin();
       map_iterator != pw_contact_pair_candidates.end();
       ++map_iterator)
    {
      auto contact_ids              = map_iterator->first;
      auto particle_pair_candidates = &map_iterator->second;

      // Get the particle and face id from the vector and the total array
      // view to the particle properties once to improve efficiency
      auto particle            = std::get<0>(*particle_pair_candidates);
      auto particle_properties = particle->get_properties();
      int  particle_id         = contact_ids.first;
      int  face_id             = contact_ids.second;

      // Normal vector of the boundary and a point on the boudary are defined as
      // local parameters
      Tensor<1, dim> normal_vector     = std::get<1>(*particle_pair_candidates);
      Point<dim>     point_on_boundary = std::get<2>(*particle_pair_candidates);

      // A vector (point_to_particle_vector) is defined which connects the
      // center of particle to the point_on_boundary. This vector will then be
      // projected on the normal vector of the boundary to obtain the
      // particle-wall distance
      Tensor<1, dim> point_to_particle_vector =
        particle->get_location() - point_on_boundary;

      // Finding the projected vector on the normal vector of the boundary. Here
      // we have used the private function find_projection. Using this projected
      // vector, the particle-wall distance is calculated
      Tensor<1, dim> projected_vector =
        find_projection(point_to_particle_vector, normal_vector);
      double distance = ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
                        (projected_vector.norm());

      // Check to see if the particle-wall pair is in contact
      if (distance > 0)
        {
          // Check to see if in the (particle->get_id())th element of the
          // pw_pairs_in_contact vector, an element with the same key as the
          // boundary exists or not. If there exist an element with this key
          // value, it shows that this particle-wall pair has already been found
          // and there is no need to calculate its properties again
          //  if (pw_pairs_in_contact[particle->get_id()].count(face_id) <= 0) {
          // If the pair is in contact (distance>0) and the pair does not
          // exist in the pw_pairs_in_contact vector, the contact properties
          // should be obtained and added to the pw_pairs_in_contact vector
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
          Tensor<1, dim> normal_relative_velocity =
            normal_relative_velocity_value * normal_vector;

          // Calculation of tangential relative velocity
          Tensor<1, dim> tangential_relative_velocity =
            contact_relative_velocity - normal_relative_velocity;

          // Setting tangential overlap of the new particle-wall contact
          // pair equal to zero
          Tensor<1, dim> tangential_overlap;
          tangential_overlap[0] = 0.0;
          tangential_overlap[1] = 0.0;
          if (dim == 3)
            {
              tangential_overlap[2] = 0.0;
            }

          // Creating a sample from the pw_contact_info_struct and adding
          // contact info to the sample
          pw_contact_info_struct<dim> contact_info;
          contact_info.particle          = particle;
          contact_info.normal_vector     = normal_vector;
          contact_info.point_on_boundary = point_on_boundary;
          contact_info.normal_overlap    = distance;
          contact_info.normal_relative_velocity =
            normal_relative_velocity_value;
          contact_info.tangential_overlap = tangential_overlap;
          contact_info.tangential_relative_velocity =
            tangential_relative_velocity;

          pw_pairs_in_contact[particle_id].insert({face_id, contact_info});
          //  }
        }
    }
}

template <int dim>
Tensor<1, dim> PWFineSearch<dim>::find_projection(Tensor<1, dim> vector_a,
                                                  Tensor<1, dim> vector_b)
{
  Tensor<1, dim> vector_c;
  vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

  return vector_c;
}

template class PWFineSearch<2>;
template class PWFineSearch<3>;
