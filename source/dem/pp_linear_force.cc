#include <dem/pp_linear_force.h>

using namespace dealii;

template <int dim>
void
PPLinearForce<dim>::calculate_pp_contact_force(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *                             ghost_adjacent_particles,
  const DEMSolverParameters<dim> &dem_parameters,
  const double &                  dt)
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
      for (auto adjacent_particles_list_iterator =
             adjacent_particles_list->begin();
           adjacent_particles_list_iterator != adjacent_particles_list->end();
           ++adjacent_particles_list_iterator)
        {
          // Defining the iterator's second value (map value) as a local
          // parameter
          auto contact_info = &adjacent_particles_list_iterator->second;

          // Getting information (location and propertis) of particle one and
          // two in contact
          auto       particle_one            = contact_info->particle_one;
          auto       particle_two            = contact_info->particle_two;
          Point<dim> particle_one_location   = particle_one->get_location();
          Point<dim> particle_two_location   = particle_two->get_location();
          auto       particle_one_properties = particle_one->get_properties();
          auto       particle_two_properties = particle_two->get_properties();

          // Calculation of normal overlap
          double normal_overlap =
            0.5 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                   particle_two_properties[DEM::PropertiesIndex::dp]) -
            particle_one_location.distance(particle_two_location);

          if (normal_overlap > 0)
            {
              // This means that the adjacent particles are in contact

              // Defining physical properties as local variable
              const auto physical_properties =
                dem_parameters.physical_properties;

              // Since the normal overlap is already calculated we update this
              // element of the container here. The rest of information are
              // updated using the following function
              contact_info->normal_overlap = normal_overlap;
              this->update_contact_information(*contact_info,
                                               particle_one_properties,
                                               particle_two_properties,
                                               particle_one_location,
                                               particle_two_location,
                                               dt);

              // This tuple (forces and torques) contains four elements which
              // are: 1, normal force, 2, tangential force, 3, tangential torque
              // and 4, rolling resistance torque, respectively
              std::tuple<Tensor<1, dim>,
                         Tensor<1, dim>,
                         Tensor<1, dim>,
                         Tensor<1, dim>>
                forces_and_torques =
                  this->calculate_linear_contact_force_and_torque(
                    physical_properties,
                    *contact_info,
                    particle_one_properties,
                    particle_two_properties);

              // Apply the calculated forces and torques on the particle pair
              this->apply_force_and_torque_real(particle_one_properties,
                                                particle_two_properties,
                                                forces_and_torques);
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
      for (auto adjacent_particles_list_iterator =
             adjacent_particles_list->begin();
           adjacent_particles_list_iterator != adjacent_particles_list->end();
           ++adjacent_particles_list_iterator)
        {
          // Defining the iterator's second value (map value) as a local
          // parameter
          auto contact_info = &adjacent_particles_list_iterator->second;

          // Getting information (location and propertis) of particle one and
          // two in contact
          auto       particle_one            = contact_info->particle_one;
          auto       particle_two            = contact_info->particle_two;
          Point<dim> particle_one_location   = particle_one->get_location();
          Point<dim> particle_two_location   = particle_two->get_location();
          auto       particle_one_properties = particle_one->get_properties();
          auto       particle_two_properties = particle_two->get_properties();

          // Calculation of normal overlap
          double normal_overlap =
            0.5 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                   particle_two_properties[DEM::PropertiesIndex::dp]) -
            particle_one_location.distance(particle_two_location);

          if (normal_overlap > 0)
            {
              // This means that the adjacent particles are in contact

              // Defining physical properties as local variable
              const auto physical_properties =
                dem_parameters.physical_properties;

              // Since the normal overlap is already calculated we update this
              // element of the container here. The rest of information are
              // updated using the following function
              contact_info->normal_overlap = normal_overlap;
              this->update_contact_information(*contact_info,
                                               particle_one_properties,
                                               particle_two_properties,
                                               particle_one_location,
                                               particle_two_location,
                                               dt);

              // This tuple (forces and torques) contains four elements which
              // are: 1, normal force, 2, tangential force, 3, tangential torque
              // and 4, rolling resistance torque, respectively
              std::tuple<Tensor<1, dim>,
                         Tensor<1, dim>,
                         Tensor<1, dim>,
                         Tensor<1, dim>>
                forces_and_torques =
                  this->calculate_linear_contact_force_and_torque(
                    physical_properties,
                    *contact_info,
                    particle_one_properties,
                    particle_two_properties);

              // Apply the calculated forces and torques on the particle pair
              this->apply_force_and_torque_ghost(particle_one_properties,
                                                 forces_and_torques);
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

// Calculates nonlinear contact force and torques
template <int dim>
std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
PPLinearForce<dim>::calculate_linear_contact_force_and_torque(
  const Parameters::Lagrangian::PhysicalProperties &physical_properties,
  pp_contact_info_struct<dim> &                     contact_info,
  const ArrayView<const double> &                   particle_one_properties,
  const ArrayView<const double> &                   particle_two_properties)
{
  // Calculation of effective mass, radius and Young's modulus of the
  // contact
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
    (2.0 * (1.0 - pow(physical_properties.Poisson_ratio_particle, 2.0)));

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties
  double normal_spring_constant =
    1.0667 * sqrt(effective_radius) * effective_youngs_modulus *
    pow((1.0667 * effective_mass * contact_info.normal_relative_velocity *
         contact_info.normal_relative_velocity /
         (sqrt(effective_radius) * effective_youngs_modulus)),
        0.2);
  double tangential_spring_constant =
    1.0667 * sqrt(effective_radius) * effective_youngs_modulus *
      pow((1.0667 * effective_mass * contact_info.tangential_relative_velocity *
           contact_info.tangential_relative_velocity /
           (sqrt(effective_radius) * effective_youngs_modulus)),
          0.2) +
    DBL_MIN;
  double normal_damping_constant = sqrt(
    (4 * effective_mass * normal_spring_constant) /
    (1 +
     pow((3.1415 / (log(physical_properties.restitution_coefficient_particle) +
                    DBL_MIN)),
         2)));
  double tangential_damping_constant = sqrt(
    (4 * effective_mass * tangential_spring_constant) /
    (1 +
     pow((3.1415 / (log(physical_properties.restitution_coefficient_particle) +
                    DBL_MIN)),
         2)));

  // Calculation of normal force using spring and dashpot normal forces
  Tensor<1, dim> spring_normal_force =
    (normal_spring_constant * contact_info.normal_overlap) *
    contact_info.normal_unit_vector;
  Tensor<1, dim> dashpot_normal_force =
    (normal_damping_constant * contact_info.normal_relative_velocity) *
    contact_info.normal_unit_vector;
  Tensor<1, dim> normal_force = spring_normal_force + dashpot_normal_force;

  double maximum_tangential_overlap =
    physical_properties.friction_coefficient_particle * normal_force.norm() /
    tangential_spring_constant;

  // Check for gross sliding
  if (contact_info.tangential_overlap.norm() > maximum_tangential_overlap)
    {
      // Gross sliding occurs and the tangential overlap and tangnetial
      // force are limited to Coulumb's criterion
      contact_info.tangential_overlap =
        maximum_tangential_overlap * (contact_info.tangential_overlap /
                                      contact_info.tangential_overlap.norm());
    }
  // Calculation of tangential force using spring and dashpot tangential
  // forces
  Tensor<1, dim> spring_tangential_force =
    tangential_spring_constant * contact_info.tangential_overlap;
  Tensor<1, dim> dashpot_tangential_force =
    tangential_damping_constant * contact_info.tangential_relative_velocity;
  Tensor<1, dim> tangential_force =
    -1.0 * spring_tangential_force + dashpot_tangential_force;

  // Calculation of torque
  // Torque caused by tangential force (tangential_torque)
  Tensor<1, dim> tangential_torque;

  if (dim == 3)
    {
      tangential_torque =
        cross_product_3d((0.5 *
                          particle_one_properties[DEM::PropertiesIndex::dp] *
                          contact_info.normal_unit_vector),
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
        particle_one_properties[DEM::PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex::omega_x + d];
    }

  omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
  double omega_ij_value = omega_ij.norm() + DBL_MIN;
  omega_ij_direction    = omega_ij / omega_ij_value;

  // Calculation of rolling resistance torque
  Tensor<1, dim> rolling_resistance_torque =
    -1.0 * physical_properties.rolling_friction_particle * effective_radius *
    normal_force.norm() * omega_ij_direction;

  return std::make_tuple(normal_force,
                         tangential_force,
                         tangential_torque,
                         rolling_resistance_torque);
}

template class PPLinearForce<2>;
template class PPLinearForce<3>;
