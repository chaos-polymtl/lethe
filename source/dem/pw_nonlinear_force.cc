#include <dem/pw_nonlinear_force.h>

using namespace dealii;

template <int dim>
PWNonLinearForce<dim>::PWNonLinearForce(
  const std::unordered_map<int, Tensor<1, dim>> boundary_translational_velocity,
  const std::unordered_map<int, double>         boundary_rotational_speed,
  const std::unordered_map<int, Tensor<1, dim>> boundary_rotational_vector,
  const double                                  triangulation_radius,
  const DEMSolverParameters<dim> &              dem_parameters)
{
  this->boundary_translational_velocity_map = boundary_translational_velocity;
  this->boundary_rotational_speed_map       = boundary_rotational_speed;
  this->boundary_rotational_vector          = boundary_rotational_vector;
  this->triangulation_radius                = triangulation_radius;

  const double wall_youngs_modulus =
    dem_parameters.physical_properties.youngs_modulus_wall;
  const double wall_poisson_ratio =
    dem_parameters.physical_properties.poisson_ratio_wall;
  const double wall_restitution_coefficient =
    dem_parameters.physical_properties.restitution_coefficient_wall;
  const double wall_friction_coefficient =
    dem_parameters.physical_properties.friction_coefficient_wall;
  const double wall_rolling_friction_coefficient =
    dem_parameters.physical_properties.rolling_friction_wall;

  for (unsigned int i = 0;
       i < dem_parameters.physical_properties.particle_type_number;
       ++i)
    {
      const double particle_youngs_modulus =
        dem_parameters.physical_properties.youngs_modulus_particle.at(i);
      const double particle_poisson_ratio =
        dem_parameters.physical_properties.poisson_ratio_particle.at(i);
      const double particle_restitution_coefficient =
        dem_parameters.physical_properties.restitution_coefficient_particle.at(
          i);
      const double particle_friction_coefficient =
        dem_parameters.physical_properties.friction_coefficient_particle.at(i);
      const double particle_rolling_friction_coefficient =
        dem_parameters.physical_properties.rolling_friction_coefficient_particle
          .at(i);

      this->effective_youngs_modulus.insert(
        {i,
         (particle_youngs_modulus * wall_youngs_modulus) /
           (wall_youngs_modulus *
              (1 - particle_poisson_ratio * particle_poisson_ratio) +
            particle_youngs_modulus *
              (1 - wall_poisson_ratio * wall_poisson_ratio))});

      this->effective_shear_modulus.insert(
        {i,
         (particle_youngs_modulus * wall_youngs_modulus) /
           ((2 * wall_youngs_modulus * (2 - particle_poisson_ratio) *
             (1 + particle_poisson_ratio)) +
            (2 * particle_youngs_modulus * (2 - wall_poisson_ratio) *
             (1 + wall_poisson_ratio)))});

      this->effective_coefficient_of_restitution.insert(
        {i,
         2 * particle_restitution_coefficient * wall_restitution_coefficient /
           (particle_restitution_coefficient + wall_restitution_coefficient)});

      this->effective_coefficient_of_friction.insert(
        {i,
         2 * particle_friction_coefficient * wall_friction_coefficient /
           (particle_friction_coefficient + wall_friction_coefficient)});

      this->effective_coefficient_of_rolling_friction.insert(
        {i,
         2 * particle_rolling_friction_coefficient *
           wall_rolling_friction_coefficient /
           (particle_rolling_friction_coefficient +
            wall_rolling_friction_coefficient)});
    }
}

template <int dim>
void
PWNonLinearForce<dim>::calculate_pw_contact_force(
  std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
    &           pw_pairs_in_contact,
  const double &dt)
{
  // Looping over pw_pairs_in_contact, which means looping over all the active
  // particles with iterator pw_pairs_in_contact_iterator
  for (auto &&pairs_in_contact_content :
       pw_pairs_in_contact | boost::adaptors::map_values)
    {
      // Now an iterator (pw_contact_information_iterator) on each element of
      // the pw_pairs_in_contact vector is defined. This iterator iterates over
      // a map which contains the required information for calculation of the
      // contact force for each particle
      for (auto &&contact_information :
           pairs_in_contact_content | boost::adaptors::map_values)
        {
          // Defining the total force of contact, properties of particle as
          // local parameters
          auto particle            = contact_information.particle;
          auto particle_properties = particle->get_properties();

          auto normal_vector     = contact_information.normal_vector;
          auto point_on_boundary = contact_information.point_on_boundary;

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
            this->find_projection(point_to_particle_vector, normal_vector);
          double normal_overlap =
            ((particle_properties[DEM::PropertiesIndex::dp]) / 2) -
            (projected_vector.norm());

          if (normal_overlap > 0)
            {
              contact_information.normal_overlap = normal_overlap;

              this->update_contact_information(contact_information,
                                               particle_properties,
                                               dt);

              // This tuple (forces and torques) contains four elements which
              // are: 1, normal force, 2, tangential force, 3, tangential torque
              // and 4, rolling resistance torque, respectively
              std::tuple<Tensor<1, dim>,
                         Tensor<1, dim>,
                         Tensor<1, dim>,
                         Tensor<1, dim>>
                forces_and_torques =
                  this->calculate_nonlinear_contact_force_and_torque(
                    contact_information, particle_properties);

              // Apply the calculated forces and torques on the particle pair
              this->apply_force_and_torque(particle_properties,
                                           forces_and_torques);
            }
          else
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_information.normal_overlap        = 0;
                  contact_information.tangential_overlap[d] = 0;
                }
            }
        }
    }
}

// Calculates nonlinear contact force and torques
template <int dim>
std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
PWNonLinearForce<dim>::calculate_nonlinear_contact_force_and_torque(
  pw_contact_info_struct<dim> &  contact_info,
  const ArrayView<const double> &particle_properties)
{
  const unsigned int particle_type =
    particle_properties[DEM::PropertiesIndex::type];

  // Calculation of model parameters (beta, sn and st). These values
  // are used to consider non-linear relation of the contact force to
  // the normal overlap
  double model_parameter_beta =
    log(this->effective_coefficient_of_restitution[particle_type]) /
    sqrt(
      pow(log(this->effective_coefficient_of_restitution[particle_type]), 2) +
      9.8696);
  double model_parameter_sn =
    2 * this->effective_youngs_modulus[particle_type] *
    sqrt(particle_properties[DEM::PropertiesIndex::dp] *
         contact_info.normal_overlap);

  // Calculation of normal and tangential spring and dashpot constants
  // using particle and wall properties
  double normal_spring_constant =
    1.3333 * this->effective_youngs_modulus[particle_type] *
    sqrt(particle_properties[DEM::PropertiesIndex::dp] / 2 *
         contact_info.normal_overlap);
  double normal_damping_constant =
    -1.8257 * model_parameter_beta *
    sqrt(model_parameter_sn * particle_properties[DEM::PropertiesIndex::mass]);
  double tangential_spring_constant =
    8 * this->effective_shear_modulus[particle_type] *
      sqrt(particle_properties[DEM::PropertiesIndex::dp] / 2 *
           contact_info.normal_overlap) +
    DBL_MIN;

  // Calculation of normal force using spring and dashpot normal forces
  Tensor<1, dim> spring_normal_force =
    (normal_spring_constant * contact_info.normal_overlap) *
    contact_info.normal_vector;
  Tensor<1, dim> dashpot_normal_force =
    (normal_damping_constant * contact_info.normal_relative_velocity) *
    contact_info.normal_vector;
  Tensor<1, dim> normal_force = spring_normal_force - dashpot_normal_force;

  // Calculation of tangential force using spring and dashpot tangential
  // forces
  Tensor<1, dim> spring_tangential_force =
    tangential_spring_constant * contact_info.tangential_overlap;
  Tensor<1, dim> tangential_force = -spring_tangential_force;

  double coulomb_threshold =
    this->effective_coefficient_of_friction[particle_type] *
    normal_force.norm();

  // Check for gross sliding
  if (tangential_force.norm() > coulomb_threshold)
    {
      // Gross sliding occurs and the tangential overlap and tangnetial
      // force are limited to Coulumb's criterion
      tangential_force =
        coulomb_threshold * (tangential_force / tangential_force.norm());

      contact_info.tangential_overlap =
        -tangential_force / (tangential_spring_constant + DBL_MIN);
    }

  // Calculation of torque
  // Torque caused by tangential force (tangential_torque)
  Tensor<1, dim> tangential_torque;

  if (dim == 3)
    {
      tangential_torque =
        cross_product_3d((0.5 * particle_properties[DEM::PropertiesIndex::dp] *
                          contact_info.normal_vector),
                         tangential_force);
    }

  // Rolling resistance torque
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, dim> angular_velocity, pw_angular_velocity;

  for (int d = 0; d < dim; ++d)
    {
      angular_velocity[d] =
        particle_properties[DEM::PropertiesIndex::omega_x + d];
    }

  double angular_velocity_value = angular_velocity.norm() + DBL_MIN;
  pw_angular_velocity           = angular_velocity / angular_velocity_value;

  // Calcualation of rolling resistance torque
  Tensor<1, dim> rolling_resistance_torque =
    -this->effective_coefficient_of_rolling_friction[particle_type] *
    ((particle_properties[DEM::PropertiesIndex::dp]) / 2) *
    normal_force.norm() * pw_angular_velocity;

  return std::make_tuple(normal_force,
                         tangential_force,
                         tangential_torque,
                         rolling_resistance_torque);
}

template class PWNonLinearForce<2>;
template class PWNonLinearForce<3>;
