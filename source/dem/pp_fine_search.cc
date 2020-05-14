#include <dem/pp_fine_search.h>

using namespace dealii;

template <int dim> PPFineSearch<dim>::PPFineSearch() {}

template <int dim>
void PPFineSearch<dim>::pp_Fine_Search(
    std::vector<std::pair<typename Particles::ParticleIterator<dim>,
                          typename Particles::ParticleIterator<dim>>>
        &contact_pair_candidates,
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
        &pairs_in_contact_info,
    double dt) {
  // Iterating over pairs_in_contact_info, which is equivalent to iteration over
  // all the particles
  {
    for (auto pairs_in_contact_info_iterator = pairs_in_contact_info.begin();
         pairs_in_contact_info_iterator != pairs_in_contact_info.end();
         ++pairs_in_contact_info_iterator) {
      // Each element of particle_in_contact_info is a map:
      auto contact_partners_map = &pairs_in_contact_info_iterator->second;

      // Iterating over each map which contains the contact information
      // including particles I and II

      for (auto map_iterator = contact_partners_map->begin();
           map_iterator != contact_partners_map->end();) {
        // Getting contact information and particles I and II as local
        // variables
        auto contact_information = map_iterator->second;
        auto particle_one = contact_information.particle_one;
        auto particle_two = contact_information.particle_two;

        // Finding the locations of the particles in contact
        Point<dim, double> particle_one_location = particle_one->get_location();
        Point<dim, double> particle_two_location = particle_two->get_location();

        // Get the total array view to the particle properties once to
        // improve efficiency
        auto particle_one_properties = particle_one->get_properties();
        auto particle_two_properties = particle_two->get_properties();

        // Finding the distance of particles based on their new position, if
        // the particles are still in contact, the distance will be equal to
        // the normal overlap
        double overlap_distance =
            0.5 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                   particle_two_properties[DEM::PropertiesIndex::dp]) -
            particle_one_location.distance(particle_two_location);

        // If the pair is still in contact
        if (overlap_distance > 0) {
          // contact_vector shows a vector from location of particle_one
          // to location of particle_two
          Tensor<1, dim> contact_vector =
              (particle_two_location - particle_one_location);

          // Using contact_vector, the contact normal vector is obtained
          Tensor<1, dim> normal_unit_vector =
              contact_vector / contact_vector.norm();

          // Using velocities and angular velocities of particles one and
          // two as vectors
          Tensor<1, dim> particle_one_velocity;
          particle_one_velocity[0] =
              particle_one_properties[DEM::PropertiesIndex::v_x];
          particle_one_velocity[1] =
              particle_one_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3) {
            particle_one_velocity[2] =
                particle_one_properties[DEM::PropertiesIndex::v_z];
          }

          Tensor<1, dim> particle_two_velocity;
          particle_two_velocity[0] =
              particle_two_properties[DEM::PropertiesIndex::v_x];
          particle_two_velocity[1] =
              particle_two_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3) {
            particle_two_velocity[2] =
                particle_two_properties[DEM::PropertiesIndex::v_z];
          }

          Tensor<1, dim> particle_one_omega;
          particle_one_omega[0] =
              particle_one_properties[DEM::PropertiesIndex::omega_x];
          particle_one_omega[1] =
              particle_one_properties[DEM::PropertiesIndex::omega_y];
          if (dim == 3) {
            particle_one_omega[2] =
                particle_one_properties[DEM::PropertiesIndex::omega_z];
          }

          Tensor<1, dim> particle_two_omega;
          particle_two_omega[0] =
              particle_two_properties[DEM::PropertiesIndex::omega_x];
          particle_two_omega[1] =
              particle_two_properties[DEM::PropertiesIndex::omega_y];
          if (dim == 3) {
            particle_two_omega[2] =
                particle_two_properties[DEM::PropertiesIndex::omega_z];
          }

          // Defining relative contact velocity
          Tensor<1, dim> contact_relative_velocity;
          if (dim == 3) {
            contact_relative_velocity =
                (particle_one_velocity - particle_two_velocity) +
                (cross_product_3d(
                    0.5 * (particle_one_properties[DEM::PropertiesIndex::dp] *
                               particle_one_omega +
                           particle_two_properties[DEM::PropertiesIndex::dp] *
                               particle_two_omega),
                    normal_unit_vector));
          } else {
            contact_relative_velocity =
                particle_one_velocity - particle_two_velocity;
          }

          // Calculation of normal relative velocity. Note that in the
          // following line the product acts as inner product since both
          // sides are vectors, while in the second line the product is
          // scalar and vector product
          double normal_relative_velocity_value =
              contact_relative_velocity * normal_unit_vector;
          Tensor<1, dim> normal_relative_velocity =
              normal_relative_velocity_value * normal_unit_vector;

          // Calculation of tangential relative velocity
          Tensor<1, dim> tangential_relative_velocity =
              contact_relative_velocity - normal_relative_velocity;

          // Calculation of new tangential_overlap, since this value is
          // history-dependent it needs the value at previous time-step
          // This variable is the main reason that we have iteration over
          // two different vectors (pairs_in_contact and
          // contact_pair_candidates): tangential_overlap of the particles
          // which were already in contact (pairs_in_contact) needs to
          // modified using its history, while the tangential_overlaps of
          // new particles are equal to zero
          Tensor<1, dim> last_step_tangential_overlap =
              contact_information.tangential_overlap;
          Tensor<1, dim> tangential_overlap =
              last_step_tangential_overlap -
              (last_step_tangential_overlap * normal_unit_vector) *
                  normal_unit_vector;

          // Adding a small value to the tangential_overlap_norm to avoid
          // 0/0 occurance
          double tangential_overlap_norm = tangential_overlap.norm() + DBL_MIN;

          Tensor<1, dim> modified_tangential_overlap =
              (last_step_tangential_overlap.norm() / tangential_overlap_norm) *
                  tangential_overlap +
              contact_information.tangential_relative_velocity * dt;

          // Creating a sample from the pp_contact_info_struct and adding
          // contact info to the sample
          pp_contact_info_struct<dim> contact_info;

          contact_info.normal_overlap = overlap_distance;
          contact_info.normal_vector = normal_unit_vector;
          contact_info.normal_relative_velocity =
              normal_relative_velocity_value;
          contact_info.tangential_relative_velocity =
              tangential_relative_velocity;
          contact_info.tangential_overlap = modified_tangential_overlap;
          contact_info.particle_one = particle_one;
          contact_info.particle_two = particle_two;

          contact_partners_map->insert_or_assign(particle_two->get_id(),
                                                 contact_info);
          ++map_iterator;
        } else {
          // If the particles are not in contact anymore (i.e. the contact
          // is finished and distance <= 0), this element should be erased
          // from pairs_in_contact_info:
          contact_partners_map->erase(map_iterator++);
        }
      }
    }
  }

  {
    // Now iterating over contact_pair_candidates (vector of pairs), which is
    // the output of broad search. If a pair is in contact (distance > 0) and
    // does not exist in the pairs_in_contact, it is added to the
    // pairs_in_contact
    for (auto candidate_iterator = contact_pair_candidates.begin();
         candidate_iterator != contact_pair_candidates.end();
         ++candidate_iterator) {
      // Get particles one and two from the vector and the total array view to
      // the particle properties once to improve efficiency
      auto particle_one = candidate_iterator->first;
      auto particle_two = candidate_iterator->second;
      auto particle_one_properties = particle_one->get_properties();
      auto particle_two_properties = particle_two->get_properties();

      // Obtaining locations of particles one and two:
      Point<dim, double> particle_one_location = particle_one->get_location();
      Point<dim, double> particle_two_location = particle_two->get_location();

      // Calculation of the distance between particles one and two:
      double size_distance =
          0.5 * ((particle_one_properties[DEM::PropertiesIndex::dp] +
                  particle_two_properties[DEM::PropertiesIndex::dp]));
      auto contact_vector = particle_two_location - particle_one_location;
      double distance = size_distance * size_distance;
      for (int d = 0; d < dim; ++d) {
        distance -= contact_vector[d] * contact_vector[d];
      }

      // Check to see if particle pair is in contact:
      if (distance > 0) {
        // Check to see if the pair already exists in pairs_in_contact
        // vector or not. Note that the pair shoule be searched in the
        // (particle_one_properties[DEM::PropertiesIndex::id])th element of
        // the vector as well as
        // (particle_two_properties[DEM::PropertiesIndex::id])th element of
        // the vector. If the pair does not exist in these two locations, it
        // should be added to pairs_in_contact

        // Since insert already checks if the element with the same key
        // exists, only the element with key = particle one id, is checked
        if (pairs_in_contact_info
                [particle_one_properties[DEM::PropertiesIndex::id]]
                    .count(particle_two_properties[DEM::PropertiesIndex::id]) <=
            0) {
          // Finding normal overlap:
          double normal_overlap =
              0.5 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                     particle_two_properties[DEM::PropertiesIndex::dp]) -
              particle_one_location.distance(particle_two_location);

          // Using contact_vector, the contact normal vector is obtained
          Tensor<1, dim> normal_vector = contact_vector / contact_vector.norm();

          // Using velocities and angular velocities of particles one and
          // two as vectors
          Tensor<1, dim> particle_one_velocity;
          particle_one_velocity[0] =
              particle_one_properties[DEM::PropertiesIndex::v_x];
          particle_one_velocity[1] =
              particle_one_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3) {
            particle_one_velocity[2] =
                particle_one_properties[DEM::PropertiesIndex::v_z];
          }

          Tensor<1, dim> particle_two_velocity;
          particle_two_velocity[0] =
              particle_two_properties[DEM::PropertiesIndex::v_x];
          particle_two_velocity[1] =
              particle_two_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3) {
            particle_two_velocity[2] =
                particle_two_properties[DEM::PropertiesIndex::v_z];
          }

          Tensor<1, dim> particle_one_omega;
          particle_one_omega[0] =
              particle_one_properties[DEM::PropertiesIndex::omega_x];
          particle_one_omega[1] =
              particle_one_properties[DEM::PropertiesIndex::omega_y];
          if (dim == 3) {
            particle_one_omega[2] =
                particle_one_properties[DEM::PropertiesIndex::omega_z];
          }

          Tensor<1, dim> particle_two_omega;
          particle_two_omega[0] =
              particle_two_properties[DEM::PropertiesIndex::omega_x];
          particle_two_omega[1] =
              particle_two_properties[DEM::PropertiesIndex::omega_y];
          if (dim == 3) {
            particle_two_omega[2] =
                particle_two_properties[DEM::PropertiesIndex::omega_z];
          }

          // Defining relative contact velocity
          Tensor<1, dim> contact_relative_velocity;
          if (dim == 3) {
            contact_relative_velocity =
                (particle_one_velocity - particle_two_velocity) +
                (cross_product_3d(
                    0.5 * (particle_one_properties[DEM::PropertiesIndex::dp] *
                               particle_one_omega +
                           particle_two_properties[DEM::PropertiesIndex::dp] *
                               particle_two_omega),
                    normal_vector));
          } else {
            contact_relative_velocity =
                particle_one_velocity - particle_two_velocity;
          }

          // Calculation of normal relative velocity. Note that in the
          // following line the product acts as inner product since both
          // sides are vectors, while in the second line the product is
          // scalar and vector product
          double normal_relative_velocity_value =
              contact_relative_velocity * normal_vector;
          Tensor<1, dim> normal_relative_velocity =
              normal_relative_velocity_value * normal_vector;

          // Calculation of tangential relative velocity
          Tensor<1, dim> tangential_relative_velocity =
              contact_relative_velocity - normal_relative_velocity;

          // For new pairs added to pairs_in_contact, the tangential
          // overlap is equal to zero
          Tensor<1, dim> tangential_overlap;
          tangential_overlap[0] = 0.0;
          tangential_overlap[1] = 0.0;
          if (dim == 3) {
            tangential_overlap[2] = 0.0;
          }

          // Creating a sample from the contact_info_struct and adding
          // contact info to the sample
          {
            pp_contact_info_struct<dim> contact_info;

            contact_info.normal_overlap = normal_overlap;
            contact_info.normal_vector = normal_vector;
            contact_info.normal_relative_velocity =
                normal_relative_velocity_value;
            contact_info.tangential_relative_velocity =
                tangential_relative_velocity;
            contact_info.tangential_overlap = tangential_overlap;
            contact_info.particle_one = particle_one;
            contact_info.particle_two = particle_two;

            {
              pairs_in_contact_info
                  [particle_one_properties[DEM::PropertiesIndex::id]]
                      .insert(
                          {particle_two_properties[DEM::PropertiesIndex::id],
                           contact_info});
            }
          }
        }
      }
    }
  }
}

template class PPFineSearch<2>;
template class PPFineSearch<3>;
