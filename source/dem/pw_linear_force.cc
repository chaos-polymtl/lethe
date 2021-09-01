#include <dem/pw_linear_force.h>

using namespace dealii;

template <int dim>
PWLinearForce<dim>::PWLinearForce(
  const std::unordered_map<unsigned int, Tensor<1, dim>>
                                                 boundary_translational_velocity,
  const std::unordered_map<unsigned int, double> boundary_rotational_speed,
  const std::unordered_map<unsigned int, Tensor<1, dim>>
                                        boundary_rotational_vector,
  const double                          triangulation_radius,
  const DEMSolverParameters<dim> &      dem_parameters,
  const std::vector<types::boundary_id> boundary_index)
{
  this->boundary_translational_velocity_map = boundary_translational_velocity;
  this->boundary_rotational_speed_map       = boundary_rotational_speed;
  this->boundary_rotational_vector          = boundary_rotational_vector;
  this->triangulation_radius                = triangulation_radius;

  const double wall_youngs_modulus =
    dem_parameters.lagrangian_physical_properties.youngs_modulus_wall;
  const double wall_poisson_ratio =
    dem_parameters.lagrangian_physical_properties.poisson_ratio_wall;
  const double wall_restitution_coefficient =
    dem_parameters.lagrangian_physical_properties.restitution_coefficient_wall;
  const double wall_friction_coefficient =
    dem_parameters.lagrangian_physical_properties.friction_coefficient_wall;
  const double wall_rolling_friction_coefficient =
    dem_parameters.lagrangian_physical_properties.rolling_friction_wall;

  for (unsigned int i = 0;
       i < dem_parameters.lagrangian_physical_properties.particle_type_number;
       ++i)
    {
      const double particle_youngs_modulus =
        dem_parameters.lagrangian_physical_properties.youngs_modulus_particle
          .at(i);
      const double particle_poisson_ratio =
        dem_parameters.lagrangian_physical_properties.poisson_ratio_particle.at(
          i);
      const double particle_restitution_coefficient =
        dem_parameters.lagrangian_physical_properties
          .restitution_coefficient_particle.at(i);
      const double particle_friction_coefficient =
        dem_parameters.lagrangian_physical_properties
          .friction_coefficient_particle.at(i);
      const double particle_rolling_friction_coefficient =
        dem_parameters.lagrangian_physical_properties
          .rolling_friction_coefficient_particle.at(i);


      this->effective_youngs_modulus[i] =
        (particle_youngs_modulus * wall_youngs_modulus) /
        (wall_youngs_modulus *
           (1 - particle_poisson_ratio * particle_poisson_ratio) +
         particle_youngs_modulus *
           (1 - wall_poisson_ratio * wall_poisson_ratio) +
         DBL_MIN);

      this->effective_coefficient_of_restitution.insert(
        {i,
         2 * particle_restitution_coefficient * wall_restitution_coefficient /
           (particle_restitution_coefficient + wall_restitution_coefficient +
            DBL_MIN)});

      this->effective_coefficient_of_friction.insert(
        {i,
         2 * particle_friction_coefficient * wall_friction_coefficient /
           (particle_friction_coefficient + wall_friction_coefficient +
            DBL_MIN)});

      this->effective_coefficient_of_rolling_friction.insert(
        {i,
         2 * particle_rolling_friction_coefficient *
           wall_rolling_friction_coefficient /
           (particle_rolling_friction_coefficient +
            wall_rolling_friction_coefficient + DBL_MIN)});
    }

  if (dem_parameters.model_parameters.rolling_resistance_method ==
      Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
        no_resistance)
    {
      calculate_rolling_resistance_torque = &PWLinearForce<dim>::no_resistance;
    }
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
             constant_resistance)
    {
      calculate_rolling_resistance_torque =
        &PWLinearForce<dim>::constant_resistance;
    }
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
             viscous_resistance)
    {
      calculate_rolling_resistance_torque =
        &PWLinearForce<dim>::viscous_resistance;
    }

  this->calculate_force_torque_on_boundary =
    dem_parameters.forces_torques.calculate_force_torque;
  this->center_mass_container =
    dem_parameters.forces_torques.triangulation_center_mass;
  this->triangulation_mass   = dem_parameters.forces_torques.triangulation_mass;
  this->inclined_plane_angle = dem_parameters.grid_motion.inclined_plane_angle;
  this->boundary_index       = boundary_index;
  for (unsigned int d = 0; d < dim; ++d)
    this->gravity[d] = dem_parameters.physical_properties.g[d];
  this->force_on_walls  = this->initialize();
  this->torque_on_walls = this->initialize();
}

template <int dim>
void
PWLinearForce<dim>::calculate_pw_contact_force(
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, pw_contact_info_struct<dim>>>
    &                          pw_pairs_in_contact,
  const double &               dt,
  std::vector<Tensor<1, dim>> &momentum,
  std::vector<Tensor<1, dim>> &force)
{
  PWContactForce<dim>::force_on_walls  = PWContactForce<dim>::initialize();
  PWContactForce<dim>::torque_on_walls = PWContactForce<dim>::initialize();
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
          // Defining total force of the contact, properties of particle as
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
                  this->calculate_linear_contact_force_and_torque(
                    contact_information, particle_properties);

              // Getting particle's momentum and force
              unsigned int    particle_id       = particle->get_id();
              Tensor<1, dim> &particle_momentum = momentum[particle_id];
              Tensor<1, dim> &particle_force    = force[particle_id];

              // Apply the calculated forces and torques on the particle pair
              this->apply_force_and_torque(forces_and_torques,
                                           particle_momentum,
                                           particle_force,
                                           point_on_boundary,
                                           contact_information.boundary_id);
            }
          else
            {
              for (int d = 0; d < dim; ++d)
                {
                  contact_information.tangential_overlap[d] = 0;
                }
            }
        }
    }
  this->mpi_correction_over_calculation_of_forces_and_torques();
}

// Calculates linear contact force and torques
template <int dim>
std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
PWLinearForce<dim>::calculate_linear_contact_force_and_torque(
  pw_contact_info_struct<dim> &  contact_info,
  const ArrayView<const double> &particle_properties)
{
  const unsigned int particle_type =
    particle_properties[DEM::PropertiesIndex::type];

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties
  double rp_sqrt = sqrt(particle_properties[DEM::PropertiesIndex::dp] * 0.5);

  double normal_spring_constant =
    1.0667 * rp_sqrt * this->effective_youngs_modulus[particle_type] *
    pow((0.9375 * particle_properties[DEM::PropertiesIndex::mass] *
         contact_info.normal_relative_velocity *
         contact_info.normal_relative_velocity /
         (rp_sqrt * this->effective_youngs_modulus[particle_type])),
        0.2);
  double tangential_spring_constant =
    -1.0667 * rp_sqrt * this->effective_youngs_modulus[particle_type] *
      pow((0.9375 * particle_properties[DEM::PropertiesIndex::mass] *
           contact_info.tangential_relative_velocity *
           contact_info.tangential_relative_velocity /
           (rp_sqrt * this->effective_youngs_modulus[particle_type])),
          0.2) +
    DBL_MIN;
  double normal_damping_constant = sqrt(
    (4 * particle_properties[DEM::PropertiesIndex::mass] *
     normal_spring_constant) /
    (1 + pow((M_PI /
              (log(this->effective_coefficient_of_restitution[particle_type]) +
               DBL_MIN)),
             2)));

  // Calculation of normal force using spring and dashpot normal forces
  Tensor<1, dim> normal_force =
    (normal_spring_constant * contact_info.normal_overlap -
     normal_damping_constant * contact_info.normal_relative_velocity) *
    contact_info.normal_vector;
  ;

  // Calculation of tangential force
  Tensor<1, dim> tangential_force =
    tangential_spring_constant * contact_info.tangential_overlap;

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
        tangential_force / (tangential_spring_constant + DBL_MIN);
    }

  // Calculation of torque
  // First calculation of torque due to tangential force acting on
  // particle
  Tensor<1, dim> tangential_torque;
  if (dim == 3)
    {
      tangential_torque =
        cross_product_3d((0.5 * particle_properties[DEM::PropertiesIndex::dp] *
                          contact_info.normal_vector),
                         tangential_force);
    }

  // Rolling resistance torque
  Tensor<1, dim> rolling_resistance_torque =
    (this->*calculate_rolling_resistance_torque)(
      particle_properties,
      this->effective_coefficient_of_rolling_friction[particle_type],
      normal_force.norm(),
      contact_info.normal_vector);

  return std::make_tuple(normal_force,
                         tangential_force,
                         tangential_torque,
                         rolling_resistance_torque);
}

template class PWLinearForce<2>;
template class PWLinearForce<3>;
