#include <dem/pp_nonlinear_force.h>

using namespace dealii;

template <int dim>
void
PPNonLinearForce<dim>::calculate_pp_contact_force(
  const std::vector<std::map<int, pp_contact_info_struct<dim>>>
    *                             pairs_in_contact_info,
  const DEMSolverParameters<dim> &dem_parameters)
{
  // Defining physical properties as local variable
  const auto physical_properties = dem_parameters.physicalProperties;

  // Looping over pairs_in_contact_info, which means looping over all the active
  // particles with iterator pairs_in_contact_info_iterator
  for (auto pairs_in_contact_info_iterator = pairs_in_contact_info->begin();
       pairs_in_contact_info_iterator != pairs_in_contact_info->end();
       ++pairs_in_contact_info_iterator)
    {
      // Now an iterator (contact_information_iterator) on each element of the
      // pairs_in_contact_info vector is defined. This iterator iterates over a
      // map which contains the required information for calculation of the
      // contact  Pointforce for each particle (i.e. each element of the
      // pairs_in_contact_info vector)

      // auto contact_information_iterator =
      // pairs_in_contact_info_iterator->begin();
      for (auto contact_information_iterator =
             pairs_in_contact_info_iterator->begin();
           contact_information_iterator !=
           pairs_in_contact_info_iterator->end();
           ++contact_information_iterator)
        {
          // Defining the iterator's second value (map value) as a local
          // parameter
          auto contact_information_iterator_second =
            contact_information_iterator->second;

          // Defining total force of the contact, properties of particles one
          // and two as local parameters
          Tensor<1, dim> total_force;
          auto           particle_one_properties =
            contact_information_iterator_second.particle_one->get_properties();
          auto particle_two_properties =
            contact_information_iterator_second.particle_two->get_properties();

          // Calculation of effective mass, radius, Young's modulus and shear
          // modulus of the contact
          double effective_mass =
            (particle_one_properties[DEM::PropertiesIndex::mass] *
             particle_two_properties[DEM::PropertiesIndex::mass]) /
            (particle_one_properties[DEM::PropertiesIndex::mass] +
             particle_two_properties[DEM::PropertiesIndex::mass]);
          double effective_radius =
            (particle_one_properties[DEM::PropertiesIndex::dp] *
             particle_two_properties[DEM::PropertiesIndex::dp]) /
            (2.0 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                    particle_two_properties[DEM::PropertiesIndex::dp]));
          double effective_youngs_modulus =
            physical_properties.Youngs_modulus_particle /
            (2.0 *
             (1.0 - pow(physical_properties.Poisson_ratio_particle, 2.0)));
          double effective_shear_modulus =
            physical_properties.Youngs_modulus_particle /
            (4.0 * (2.0 - physical_properties.Poisson_ratio_particle) *
             (1.0 + physical_properties.Poisson_ratio_particle));

          // Calculation of model parameters (betha, sn and st). These values
          // are used to consider non-linear relation of the contact force to
          // the normal overlap
          double model_parameter_betha =
            log(physical_properties.Poisson_ratio_particle) /
            sqrt(pow(log(physical_properties.Poisson_ratio_particle), 2.0) +
                 9.8696);
          double model_parameter_sn =
            2.0 * effective_youngs_modulus *
            sqrt(effective_radius *
                 contact_information_iterator_second.normal_overlap);
          double model_parameter_st =
            8.0 * effective_shear_modulus *
            sqrt(effective_radius *
                 contact_information_iterator_second.normal_overlap);

          // Calculation of normal and tangential spring and dashpot constants
          // using particle properties
          double normal_spring_constant =
            1.3333 * effective_youngs_modulus *
            sqrt(effective_radius *
                 contact_information_iterator_second.normal_overlap);
          double normal_damping_constant =
            -1.8257 * model_parameter_betha *
            sqrt(model_parameter_sn * effective_mass);
          double tangential_spring_constant =
            8.0 * effective_shear_modulus *
            sqrt(effective_radius *
                 contact_information_iterator_second.normal_overlap);
          double tangential_damping_constant =
            -1.8257 * model_parameter_betha *
            sqrt(model_parameter_st * effective_mass);

          // Calculation of normal force using spring and dashpot normal forces
          Tensor<1, dim> spring_normal_force =
            (normal_spring_constant *
             contact_information_iterator_second.normal_overlap) *
            contact_information_iterator_second.normal_vector;
          Tensor<1, dim> dashpot_normal_force =
            (normal_damping_constant *
             contact_information_iterator_second.normal_relative_velocity) *
            contact_information_iterator_second.normal_vector;
          Tensor<1, dim> normal_force;
          normal_force = spring_normal_force + dashpot_normal_force;

          double maximum_tangential_overlap =
            physical_properties.friction_coefficient_particle *
            normal_force.norm() / tangential_spring_constant;

          // Check for gross sliding
          if (contact_information_iterator_second.tangential_overlap.norm() >
              maximum_tangential_overlap)
            {
              // Gross sliding occurs and the tangential overlap and tangnetial
              // force are limited to Coulumb's criterion
              contact_information_iterator_second.tangential_overlap =
                maximum_tangential_overlap *
                (contact_information_iterator_second.tangential_overlap /
                 contact_information_iterator_second.tangential_overlap.norm());
            }
          // Calculation of tangential force using spring and dashpot tangential
          // forces
          Tensor<1, dim> spring_tangential_force =
            tangential_spring_constant *
            contact_information_iterator_second.tangential_overlap;
          Tensor<1, dim> dashpot_tangential_force =
            tangential_damping_constant *
            contact_information_iterator_second.tangential_relative_velocity;
          Tensor<1, dim> tangential_force;
          tangential_force = spring_tangential_force + dashpot_tangential_force;

          // Calculation of total force
          total_force = normal_force + tangential_force;

          // Updating the force of particles in the particle handler
          for (int d = 0; d < dim; ++d)
            {
              particle_one_properties[DEM::PropertiesIndex::force_x + d] =
                particle_one_properties[DEM::PropertiesIndex::force_x + d] -
                total_force[d];
              particle_two_properties[DEM::PropertiesIndex::force_x + d] =
                particle_two_properties[DEM::PropertiesIndex::force_x + d] +
                total_force[d];
            }
          /*
                // Calculation of torque
                // Torque caused by tangential force (tangential_torque)
                Tensor<1, dim> tangential_torque_particle_one,
                    tangential_torque_particle_two;

                tangential_torque_particle_one =
                    (particle_one_properties[DEM::PropertiesIndex::dp] / 2.0) *
                    cross_product_3d(contact_information_iterator_second.normal_vector,
                                     (-1.0 * total_force));
                tangential_torque_particle_two =
                    (particle_one_properties[DEM::PropertiesIndex::dp] / 2.0) *
                    cross_product_3d(contact_information_iterator_second.normal_vector,
                                     total_force);
                tangential_torque_particle_one = -1.0 *
             tangential_torque_particle_two;

                // Rolling resistance torque
                // For calculation of rolling resistance torque, we need to
        obtain
                // omega_ij using rotational velocities of particles one and two
                Tensor<1, dim> particle_one_angular_velocity,
                    particle_two_angular_velocity, omega_ij;
                for (int d = 0; d < dim; ++d) {
                  particle_one_angular_velocity[d] =
                      particle_one_properties[DEM::PropertiesIndex::omega_x +
        d]; particle_two_angular_velocity[d] =
                      particle_two_properties[DEM::PropertiesIndex::omega_x +
        d]; omega_ij[d] = 0.0;
                }

                double omega_value =
                    (particle_one_angular_velocity -
        particle_two_angular_velocity) .norm(); if (omega_value != 0) { omega_ij
        = (particle_one_angular_velocity - particle_two_angular_velocity) /
        omega_value; omega_ij = (particle_two_angular_velocity -
             particle_one_angular_velocity) / omega_value;
                }

                // Calculation of rolling resistance torque

                Tensor<1, dim> rolling_resistance_torque =
                    -1.0 * physical_properties.rolling_friction_particle *
                    effective_radius * normal_force.norm() * omega_ij;

                // Updating the torque acting on particles
                for (int d = 0; d < dim; ++d) {
                  particle_one_properties[DEM::PropertiesIndex::M_x + d] =
                      particle_one_properties[DEM::PropertiesIndex::M_x + d] +
                      tangential_torque_particle_one[d] +
             rolling_resistance_torque[d];
                  particle_two_properties[DEM::PropertiesIndex::M_x + d] =
                      particle_two_properties[DEM::PropertiesIndex::M_x + d] +
                      tangential_torque_particle_two[d] +
             rolling_resistance_torque[d];

        }
          */
        }
    }
}

template class PPNonLinearForce<2>;
template class PPNonLinearForce<3>;
