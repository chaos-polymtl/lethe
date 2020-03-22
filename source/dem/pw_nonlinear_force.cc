#include <dem/pw_nonlinear_force.h>

using namespace dealii;

template <int dim>
void PWNonLinearForce<dim>::calculate_pw_contact_force(
    const std::vector<std::map<int, pw_contact_info_struct<dim>>>
        *pw_pairs_in_contact,
    const DEMSolverParameters<dim> &dem_parameters) {
  // Defining physical properties as local variable
  const auto physical_properties = dem_parameters.physicalProperties;

  // Looping over pw_pairs_in_contact, which means looping over all the active
  // particles with iterator pw_pairs_in_contact_iterator
  for (auto pw_pairs_in_contact_iterator = pw_pairs_in_contact->begin();
       pw_pairs_in_contact_iterator != pw_pairs_in_contact->end();
       ++pw_pairs_in_contact_iterator) {
    // Now an iterator (pw_contact_information_iterator) on each element of
    // the pw_pairs_in_contact vector is defined. This iterator iterates over
    // a map which contains the required information for calculation of the
    // contact force for each particle

    // auto pw_contact_information_iterator =
    //    pw_pairs_in_contact_iterator->begin();
    for (auto pw_contact_information_iterator =
             pw_pairs_in_contact_iterator->begin();
         pw_contact_information_iterator != pw_pairs_in_contact_iterator->end();
         ++pw_contact_information_iterator) {
      // Defining the iterator's second value (map value) as a local
      // parameter
      auto contact_information = pw_contact_information_iterator->second;

      // Defining total force of the contact, properties of particle as
      // local parameters
      Tensor<1, dim> total_force;
      auto particle = contact_information.particle;
      auto particle_properties = particle->get_properties();

      // Calculation of effective Young's modulus and shear
      // modulus of the contact
      double effective_youngs_modulus =
          physical_properties.Youngs_modulus_wall /
          (2.0 * (1.0 - pow(physical_properties.Poisson_ratio_wall, 2.0)));
      double effective_shear_modulus =
          physical_properties.Youngs_modulus_wall /
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

      double coulumb_force_value =
          physical_properties.friction_coefficient_wall * normal_force.norm();

      // Check for gross sliding
      if (tangential_force.norm() <= coulumb_force_value) {
        // No gross sliding here
        total_force = normal_force + tangential_force;
        // Tangential overlap is not changed
      } else {
        // Gross sliding occurs and the tangential overlap and tangnetial
        // force are limited to Coulumb's criterion
        contact_information.tangential_overlap =
            tangential_force.norm() / tangential_spring_constant;

        Tensor<1, dim> coulumb_tangential_force =
            (-1.0 * coulumb_force_value) *
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
      // First calculation of torque due to tangential force acting on
      // particle
      Tensor<1, dim> tangential_toruqe =
          ((particle_properties[DEM::PropertiesIndex::dp]) / 2.0) *
          cross_product_3d(contact_information.normal_vector, total_force);

      // Getting the angular velocity of particle in the vector format
      Tensor<1, dim> angular_velocity;
      for (int d = 0; d < dim; ++d) {
        angular_velocity[d] =
            particle_properties[DEM::PropertiesIndex::omega_x + d];
      }

      // Calculation of particle-wall angular velocity (norm of the

      // particle angular velocity)
      Tensor<1, dim> pw_angular_velocity;
      for (int d = 0; d < dim; ++d) {
        pw_angular_velocity[d] = 0.0;
      }
      double omegaNorm = angular_velocity.norm();
      if (omegaNorm != 0) {
        pw_angular_velocity = angular_velocity / omegaNorm;
      }

      // Calcualation of rolling resistance torque
      Tensor<1, dim> rolling_resistance_torque =
          -1.0 * physical_properties.rolling_friction_wall *
          ((particle_properties[DEM::PropertiesIndex::dp]) / 2.0) *
          normal_force.norm() * pw_angular_velocity;

      // Updating the acting toruqe on the particle
      for (int d = 0; d < dim; ++d) {
        particle_properties[DEM::PropertiesIndex::M_x + d] =
            particle_properties[DEM::PropertiesIndex::M_x + d] +
            tangential_toruqe[d] + rolling_resistance_torque[d];
      }

      //++pw_contact_information_iterator;
    }
  }
}

template class PWNonLinearForce<2>;
template class PWNonLinearForce<3>;
