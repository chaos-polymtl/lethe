#include <dem/pp_nonlinear_force.h>

using namespace DEM;

template <int dim>
void
PPNonLinearForce<dim>::calculate_pp_contact_force(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *                                               ghost_adjacent_particles,
  const Parameters::Lagrangian::PhysicalProperties &physical_properties,
  const double &                                    dt)
{
  // Updating contact force of particles in local-local and local-ghost contact
  // pairs are differnet. Consequently, contact forces of local-local and
  // local-ghost particle pairs are performed in separate loops

  // Looping over local_adjacent_particles with iterator
  // adjacent_particles_iterator
  for (auto adjacent_particles_iterator = local_adjacent_particles->begin();
       adjacent_particles_iterator != local_adjacent_particles->end();
       ++adjacent_particles_iterator)
    {
      // Now an iterator (adjacent_particles_list_iterator) on each element of
      // the adjacent_particles_iterator map is defined. This iterator iterates
      // over another map which contains the required information for
      // calculation of the contact force
      auto adjacent_particles_list = &adjacent_particles_iterator->second;

      if (!adjacent_particles_list->empty())
        {
          for (auto adjacent_particles_list_iterator =
                 adjacent_particles_list->begin();
               adjacent_particles_list_iterator !=
               adjacent_particles_list->end();
               ++adjacent_particles_list_iterator)
            {
              // Defining the iterator's second value (map value) as a local
              // parameter
              auto contact_info = &adjacent_particles_list_iterator->second;

              // Getting information (location and propertis) of particle one
              // and two in contact
              auto             particle_one = contact_info->particle_one;
              auto             particle_two = contact_info->particle_two;
              const Point<dim> particle_one_location =
                particle_one->get_location();
              const Point<dim> particle_two_location =
                particle_two->get_location();
              auto particle_one_properties = particle_one->get_properties();
              auto particle_two_properties = particle_two->get_properties();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0)
                // This means that the adjacent particles are in contact
                {
                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    *contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  this->calculate_nonlinear_contact_force_and_torque(
                    physical_properties,
                    *contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    normal_overlap,
                    particle_one_properties,
                    particle_two_properties,
                    normal_force,
                    tangential_force,
                    tangential_torque,
                    rolling_resistance_torque);

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_real(particle_one_properties,
                                                    particle_two_properties,
                                                    normal_force,
                                                    tangential_force,
                                                    tangential_torque,
                                                    rolling_resistance_torque);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info->tangential_overlap[d] = 0;
                    }
                }
            }
        }
    }

  // Doing the same calculations for local-ghost particle pairs

  // Looping over ghost_adjacent_particles with iterator
  // adjacent_particles_iterator
  for (auto adjacent_particles_iterator = ghost_adjacent_particles->begin();
       adjacent_particles_iterator != ghost_adjacent_particles->end();
       ++adjacent_particles_iterator)
    {
      // Now an iterator (adjacent_particles_list_iterator) on each element of
      // the adjacent_particles_iterator map is defined. This iterator iterates
      // over another map which contains the required information for
      // calculation of the contact force
      auto adjacent_particles_list = &adjacent_particles_iterator->second;

      if (!adjacent_particles_list->empty())
        {
          for (auto adjacent_particles_list_iterator =
                 adjacent_particles_list->begin();
               adjacent_particles_list_iterator !=
               adjacent_particles_list->end();
               ++adjacent_particles_list_iterator)
            {
              // Defining the iterator's second value (map value) as a local
              // parameter
              auto contact_info = &adjacent_particles_list_iterator->second;

              // Getting information (location and propertis) of particle one
              // and two in contact
              auto             particle_one = contact_info->particle_one;
              auto             particle_two = contact_info->particle_two;
              const Point<dim> particle_one_location =
                particle_one->get_location();
              const Point<dim> particle_two_location =
                particle_two->get_location();
              auto particle_one_properties = particle_one->get_properties();
              auto particle_two_properties = particle_two->get_properties();

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    *contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  this->calculate_nonlinear_contact_force_and_torque(
                    physical_properties,
                    *contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    normal_overlap,
                    particle_one_properties,
                    particle_two_properties,
                    normal_force,
                    tangential_force,
                    tangential_torque,
                    rolling_resistance_torque);

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_ghost(particle_one_properties,
                                                     normal_force,
                                                     tangential_force,
                                                     tangential_torque,
                                                     rolling_resistance_torque);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info->tangential_overlap[d] = 0;
                    }
                }
            }
        }
    }
}

// Calculates nonlinear contact force and torques
template <int dim>
void
PPNonLinearForce<dim>::calculate_nonlinear_contact_force_and_torque(
  const Parameters::Lagrangian::PhysicalProperties &physical_properties,
  pp_contact_info_struct<dim> &                     contact_info,
  const double &                 normal_relative_velocity_value,
  const Tensor<1, dim> &         normal_unit_vector,
  const double &                 normal_overlap,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  Tensor<1, dim> &               normal_force,
  Tensor<1, dim> &               tangential_force,
  Tensor<1, dim> &               tangential_torque,
  Tensor<1, dim> &               rolling_resistance_torque)
{
  // Calculation of effective mass, radius, Young's modulus and shear
  // modulus of the contact
  const double effective_mass =
    (particle_one_properties[PropertiesIndex::mass] *
     particle_two_properties[PropertiesIndex::mass]) /
    (particle_one_properties[PropertiesIndex::mass] +
     particle_two_properties[PropertiesIndex::mass]);
  const double effective_radius =
    (particle_one_properties[PropertiesIndex::dp] *
     particle_two_properties[PropertiesIndex::dp]) /
    (2 * (particle_one_properties[PropertiesIndex::dp] +
          particle_two_properties[PropertiesIndex::dp]));
  const double effective_youngs_modulus =
    physical_properties.youngs_modulus_particle /
    (2 * (1 - (physical_properties.poisson_ratio_particle *
               physical_properties.poisson_ratio_particle)));
  const double effective_shear_modulus =
    physical_properties.youngs_modulus_particle /
    (4 * (2 - physical_properties.poisson_ratio_particle) *
     (1 + physical_properties.poisson_ratio_particle));

  // Calculation of model parameters (beta, sn and st). These values
  // are used to consider non-linear relation of the contact force to
  // the normal overlap
  const double restitution_coefficient_particle_log =
    std::log(physical_properties.restitution_coefficient_particle);
  const double model_parameter_beta =
    restitution_coefficient_particle_log /
    sqrt(restitution_coefficient_particle_log *
           restitution_coefficient_particle_log +
         9.8696);
  const double radius_times_overlap_sqrt =
    sqrt(effective_radius * normal_overlap);
  const double model_parameter_sn =
    2 * effective_youngs_modulus * radius_times_overlap_sqrt;
  double model_parameter_st =
    8 * effective_shear_modulus * radius_times_overlap_sqrt;

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties
  double normal_spring_constant =
    1.3333 * effective_youngs_modulus * radius_times_overlap_sqrt;
  double normal_damping_constant =
    -1.8257 * model_parameter_beta * sqrt(model_parameter_sn * effective_mass);
  double tangential_spring_constant =
    8 * effective_shear_modulus * radius_times_overlap_sqrt + DBL_MIN;
  double tangential_damping_constant =
    -1.8257 * model_parameter_beta * sqrt(model_parameter_st * effective_mass);

  // Calculation of normal force using spring and dashpot normal forces
  normal_force =
    ((normal_spring_constant * normal_overlap) * normal_unit_vector) +
    ((normal_damping_constant * normal_relative_velocity_value) *
     normal_unit_vector);

  // Calculation of tangential force using spring and dashpot tangential
  // forces. Since we need dashpot tangential force in the gross sliding again,
  // we define it as a separate variable
  Tensor<1, dim> dashpot_tangential_force =
    (tangential_damping_constant * contact_info.tangential_relative_velocity);
  tangential_force =
    (tangential_spring_constant * contact_info.tangential_overlap) +
    dashpot_tangential_force;

  double coulomb_threshold =
    physical_properties.friction_coefficient_particle * normal_force.norm();
  // Check for gross sliding
  if (tangential_force.norm() > coulomb_threshold)
    {
      // Gross sliding occurs and the tangential overlap and tangnetial
      // force are limited to Coulumb's criterion
      tangential_force =
        coulomb_threshold * (tangential_force / tangential_force.norm());

      contact_info.tangential_overlap =
        (tangential_force - dashpot_tangential_force) /
        (tangential_spring_constant + DBL_MIN);
    }

  // Calculation of torque
  // Torque caused by tangential force (tangential_torque)
  if (dim == 3)
    {
      tangential_torque =
        cross_product_3d((effective_radius * normal_unit_vector),
                         tangential_force);
    }

  // Rolling resistance torque
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, dim> particle_one_angular_velocity, particle_two_angular_velocity,
    omega_ij, omega_ij_direction;
  for (int d = 0; d < dim; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[PropertiesIndex::omega_x + d];
    }

  omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
  omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  // Calculation of rolling resistance torque
  rolling_resistance_torque = -physical_properties.rolling_friction_particle *
                              effective_radius * normal_force.norm() *
                              omega_ij_direction;
}

template class PPNonLinearForce<2>;
template class PPNonLinearForce<3>;
