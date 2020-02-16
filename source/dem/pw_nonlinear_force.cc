#include <dem/pw_nonlinear_force.h>

using namespace dealii;

template <int dim, int spacedim>
void PWNonLinearForce<dim, spacedim>::calculate_pw_contact_force(
    std::vector<std::map<int, pw_contact_info_struct<dim, spacedim>>>
        &pw_pairs_in_contact,
    const physical_info_struct<dim> &physical_properties) {

  // Looping over pw_pairs_in_contact, which means looping over all the active
  // particles with iterator pw_pairs_in_contact_iterator
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact.begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact.end();
       ++pw_pairs_in_contact_iterator) {

    // Now an iterator (pw_contact_information_iterator) on each element of the
    // pw_pairs_in_contact vector is defined. This iterator iterates over a
    // map which contains the required information for calculation of the
    // contact force for each particle
    auto pw_contact_information_iterator =
        pw_pairs_in_contact_iterator->begin();
    while (pw_contact_information_iterator !=
           pw_pairs_in_contact_iterator->end()) {
      // Defining the iterator's second value (map value) as a local
      // parameter
      auto contact_information = pw_contact_information_iterator->second;

      // Defining total force of the contact, properties of particle as local
      // parameters
      Tensor<1, dim> total_force;
      auto particle = contact_information.particle;
      auto particle_properties = particle->get_properties();

      // Calculation of effective Young's modulus and shear
      // modulus of the contact
      double effective_youngs_modulus =
          physical_properties.Young_modulus_wall /
          (2.0 * (1.0 - pow(physical_properties.Poisson_ratio_wall, 2.0)));
      double effective_shear_modulus =
          physical_properties.Young_modulus_wall /
          (4.0 * (2.0 - physical_properties.Poisson_ratio_wall) *
           (1.0 + physical_properties.Poisson_ratio_wall));

      // Calculation of model parameters (betha, sn and st). These values
      // are used to consider non-linear relation of the contact force to
      // the normal overlap
      double model_parameter_betha =
          log(physical_properties.Poisson_ratio_wall) /
          sqrt(pow(log(physical_properties.Poisson_ratio_wall), 2.0) + 9.8696);
      double model_parameter_sn =
          2.0 * effective_youngs_modulus *
          sqrt(particle_properties[DEM::PropertiesIndex::dp] *
               contact_information.normal_overlap);
      double model_parameter_st =
          8.0 * effective_shear_modulus *
          sqrt(particle_properties[DEM::PropertiesIndex::dp] *
               contact_information.normal_overlap);

      // Calculation of normal and tangential spring and dashpot constants
      // using particle and wall properties
      double normal_spring_constant =
          1.3333 * effective_youngs_modulus *
          sqrt(particle_properties[DEM::PropertiesIndex::dp] *
               contact_information.normal_overlap);
      double normal_damping_constant =
          -1.8257 * model_parameter_betha *
          sqrt(model_parameter_sn *
               particle_properties[DEM::PropertiesIndex::mass]);
      double tangential_spring_constant =
          8.0 * effective_shear_modulus *
          sqrt(particle_properties[DEM::PropertiesIndex::dp] *
               contact_information.normal_overlap);
      double tangential_damping_constant =
          -1.8257 * model_parameter_betha *
          sqrt(model_parameter_st *
               particle_properties[DEM::PropertiesIndex::mass]);

      // Calculation of normal force using spring and dashpot normal forces
      Tensor<1, dim> spring_normal_force =
          (normal_spring_constant * contact_information.normal_overlap) *
          contact_information.normal_vector;

      Tensor<1, dim> dashpot_normal_force =
          (normal_damping_constant *
           contact_information.normal_relative_velocity) *
          contact_information.normal_vector;

      Tensor<1, dim> normal_force = spring_normal_force - dashpot_normal_force;

      // Calculation of tangential force using spring and dashpot tangential
      // forces
      Tensor<1, dim> spring_tangential_force =
          (tangential_spring_constant *
           contact_information.tangential_overlap) *
          contact_information.tangential_vector;
      Tensor<1, dim> dashpot_tangential_force =
          (tangential_damping_constant *
           contact_information.tangential_relative_velocity) *
          contact_information.tangential_vector;
      Tensor<1, dim> tangential_force =
          spring_tangential_force - dashpot_tangential_force;

      double maximum_tangential_overlap =
          contact_information.normal_overlap *
          physical_properties.friction_coefficient_wall *
          ((2 - physical_properties.Poisson_ratio_wall) /
           (2 * (1 - physical_properties.Poisson_ratio_wall)));

      // Check for gross sliding
      if (contact_information.tangential_overlap <=
          maximum_tangential_overlap) {
        // No gross sliding here
        total_force = normal_force + tangential_force;
        // Tangential overlap is not changed
      } else {
        // Gross sliding occurs and the tangential overlap and tangnetial force
        // are limited to Coulumb's criterion
        contact_information.tangential_overlap =
            maximum_tangential_overlap *
            boost::math::sign(contact_information.tangential_overlap);

        Tensor<1, dim> coulumb_tangential_force =
            (-1.0 * physical_properties.friction_coefficient_wall *
             normal_force.norm() *
             boost::math::sign(contact_information.tangential_overlap)) *
            contact_information.tangential_vector;
        total_force = normal_force + coulumb_tangential_force;
      }

      // Updating the body force of particles in the particle handler
      for (int d = 0; d < dim; ++d) {
        particle_properties[DEM::PropertiesIndex::force_x + d] =
            particle_properties[DEM::PropertiesIndex::force_x + d] +
            total_force[d];
      }

      // Calculation of torque
      /*
       Point<dim> torqueTi;
       torqueTi =
  ((particle_properties[2])/2.0) *
  cross_product_3d( contact_information.normal_vector , total_force);
  Point<dim> omegai = {particle_properties[16] ,
  particle_properties[17] ,
  particle_properties[18]};

      Point<dim> omegaiw = {0.0, 0.0, 0.0};
      double omegaNorm = omegai.norm();
      if(omegaNorm != 0)
      {omegaiw = omegai / omegaNorm ;}
      Point<dim> torquer;
     torquer = -1.0 * physical_properties.rolling_friction_coefficient_wall *
  ((particle_properties[2])/2.0) *
  normal_force.norm() * omegaiw;

     particle_properties[21] =
  particle_properties[21] + torqueTi[0] +
  torquer[0]; particle_properties[22] =
  particle_properties[22] + torqueTi[1] +
  torquer[1]; particle_properties[23] =
  particle_properties[23] + torqueTi[2] +
  torquer[2];
  */

      ++pw_contact_information_iterator;
    }
  }
}

template class PWNonLinearForce<2, 2>;
template class PWNonLinearForce<3, 3>;
