#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/particle_particle_linear_force.h>

using namespace DEM;

template <int dim>
ParticleParticleLinearForce<dim>::ParticleParticleLinearForce(
  const DEMSolverParameters<dim> &dem_parameters)
{
  for (unsigned int i = 0;
       i < dem_parameters.lagrangian_physical_properties.particle_type_number;
       ++i)
    {
      const double youngs_modulus_i =
        dem_parameters.lagrangian_physical_properties.youngs_modulus_particle
          .at(i);
      const double poisson_ratio_i =
        dem_parameters.lagrangian_physical_properties.poisson_ratio_particle.at(
          i);
      const double restitution_coefficient_i =
        dem_parameters.lagrangian_physical_properties
          .restitution_coefficient_particle.at(i);
      const double friction_coefficient_i =
        dem_parameters.lagrangian_physical_properties
          .friction_coefficient_particle.at(i);
      const double rolling_friction_coefficient_i =
        dem_parameters.lagrangian_physical_properties
          .rolling_friction_coefficient_particle.at(i);

      for (unsigned int j = 0;
           j <
           dem_parameters.lagrangian_physical_properties.particle_type_number;
           ++j)
        {
          const double youngs_modulus_j =
            dem_parameters.lagrangian_physical_properties
              .youngs_modulus_particle.at(j);
          const double poisson_ratio_j =
            dem_parameters.lagrangian_physical_properties.poisson_ratio_particle
              .at(j);
          const double restitution_coefficient_j =
            dem_parameters.lagrangian_physical_properties
              .restitution_coefficient_particle.at(j);
          const double friction_coefficient_j =
            dem_parameters.lagrangian_physical_properties
              .friction_coefficient_particle.at(j);
          const double rolling_friction_coefficient_j =
            dem_parameters.lagrangian_physical_properties
              .rolling_friction_coefficient_particle.at(j);

          this->effective_youngs_modulus[i][j] =
            (youngs_modulus_i * youngs_modulus_j) /
            ((youngs_modulus_j * (1.0 - poisson_ratio_i * poisson_ratio_i)) +
             (youngs_modulus_i * (1.0 - poisson_ratio_j * poisson_ratio_j)) +
             DBL_MIN);

          this->effective_shear_modulus[i][j] =
            (youngs_modulus_i * youngs_modulus_j) /
            (2.0 * ((youngs_modulus_j * (2.0 - poisson_ratio_i) *
                     (1.0 + poisson_ratio_i)) +
                    (youngs_modulus_i * (2.0 - poisson_ratio_j) *
                     (1.0 + poisson_ratio_j))) +
             DBL_MIN);

          this->effective_coefficient_of_restitution[i][j] =
            harmonic_mean(restitution_coefficient_i, restitution_coefficient_j);

          this->effective_coefficient_of_friction[i][j] =
            harmonic_mean(friction_coefficient_i, friction_coefficient_j);

          this->effective_coefficient_of_rolling_friction[i][j] =
            harmonic_mean(rolling_friction_coefficient_i,
                          rolling_friction_coefficient_j);
        }
    }
  if (dem_parameters.model_parameters.rolling_resistance_method ==
      Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
        no_resistance)
    rolling_reistance_model =
      RollingResistanceTorqueModel::no_rolling_resistance;
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
             constant_resistance)
    rolling_reistance_model =
      RollingResistanceTorqueModel::constant_rolling_resistance;
  else if (dem_parameters.model_parameters.rolling_resistance_method ==
           Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
             viscous_resistance)
    rolling_reistance_model =
      RollingResistanceTorqueModel::viscous_rolling_resistance;
}

template <int dim>
void
ParticleParticleLinearForce<dim>::calculate_particle_particle_contact_force(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &                        ghost_adjacent_particles,
  const double &             dt,
  std::vector<Tensor<1, 3>> &torque,
  std::vector<Tensor<1, 3>> &force)
{
  // Contact forces calculations of local-local and local-ghost particle
  // pairs are performed in separate loops

  // Looping over local_adjacent_particles values with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          for (auto &&contact_info :
               adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle one
              // and two in contact
              auto     particle_one = contact_info.particle_one;
              auto     particle_two = contact_info.particle_two;
              Point<3> particle_one_location;
              Point<3> particle_two_location;
              auto     particle_one_properties = particle_one->get_properties();
              auto     particle_two_properties = particle_two->get_properties();

              if constexpr (dim == 3)
                {
                  particle_one_location = particle_one->get_location();
                  particle_two_location = particle_two->get_location();
                }

              if constexpr (dim == 2)
                {
                  particle_one_location =
                    point_nd_to_3d(particle_one->get_location());
                  particle_two_location =
                    point_nd_to_3d(particle_two->get_location());
                }


              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  this->calculate_linear_contact_force_and_torque(
                    contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    normal_overlap,
                    particle_one_properties,
                    particle_two_properties,
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    particle_two_tangential_torque,
                    rolling_resistance_torque);

                  // Getting particles' torque and force
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
                  types::particle_index particle_one_id =
                    particle_one->get_id();
                  types::particle_index particle_two_id =
                    particle_two->get_id();
#else
                  types::particle_index particle_one_id =
                    particle_one->get_local_index();
                  types::particle_index particle_two_id =
                    particle_two->get_local_index();
#endif
                  Tensor<1, 3> &particle_one_torque = torque[particle_one_id];
                  Tensor<1, 3> &particle_two_torque = torque[particle_two_id];
                  Tensor<1, 3> &particle_one_force  = force[particle_one_id];
                  Tensor<1, 3> &particle_two_force  = force[particle_two_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_local_particles(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    particle_two_tangential_torque,
                    rolling_resistance_torque,
                    particle_one_torque,
                    particle_two_torque,
                    particle_one_force,
                    particle_two_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d] = 0.0;
                    }
                }
            }
        }
    }

  // Doing the same calculations for local-ghost particle pairs

  // Looping over ghost_adjacent_particles with iterator
  // adjacent_particles_iterator
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          for (auto &&contact_info :
               adjacent_particles_list | boost::adaptors::map_values)
            {
              // Getting information (location and properties) of particle one
              // and two in contact
              auto     particle_one = contact_info.particle_one;
              auto     particle_two = contact_info.particle_two;
              Point<3> particle_one_location;
              Point<3> particle_two_location;
              auto     particle_one_properties = particle_one->get_properties();
              auto     particle_two_properties = particle_two->get_properties();

              if constexpr (dim == 3)
                {
                  particle_one_location = particle_one->get_location();
                  particle_two_location = particle_two->get_location();
                }

              if constexpr (dim == 2)
                {
                  particle_one_location =
                    point_nd_to_3d(particle_one->get_location());
                  particle_two_location =
                    point_nd_to_3d(particle_two->get_location());
                }

              // Calculation of normal overlap
              double normal_overlap =
                0.5 * (particle_one_properties[PropertiesIndex::dp] +
                       particle_two_properties[PropertiesIndex::dp]) -
                particle_one_location.distance(particle_two_location);

              if (normal_overlap > 0.0)
                {
                  // This means that the adjacent particles are in contact

                  // Since the normal overlap is already calculated we update
                  // this element of the container here. The rest of information
                  // are updated using the following function
                  this->update_contact_information(
                    contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  this->calculate_linear_contact_force_and_torque(
                    contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    normal_overlap,
                    particle_one_properties,
                    particle_two_properties,
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    particle_two_tangential_torque,
                    rolling_resistance_torque);

                  // Getting torque and force of particle one
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
                  types::particle_index particle_one_id =
                    particle_one->get_id();
#else
                  types::particle_index particle_one_id =
                    particle_one->get_local_index();
#endif
                  Tensor<1, 3> &particle_one_torque = torque[particle_one_id];
                  Tensor<1, 3> &particle_one_force  = force[particle_one_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_on_ghost_particles(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    rolling_resistance_torque,
                    particle_one_torque,
                    particle_one_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d] = 0.0;
                    }
                }
            }
        }
    }
}

template <int dim>
void
ParticleParticleLinearForce<dim>::calculate_IB_particle_particle_contact_force(
  const double &                              normal_overlap,
  particle_particle_contact_info_struct<dim> &contact_info,
  Tensor<1, 3> &                              normal_force,
  Tensor<1, 3> &                              tangential_force,
  Tensor<1, 3> &                              particle_one_tangential_torque,
  Tensor<1, 3> &                              particle_two_tangential_torque,
  Tensor<1, 3> &                              rolling_resistance_torque,
  IBParticle<dim> &                           particle_one,
  IBParticle<dim> &                           particle_two,
  const Point<dim> &                          particle_one_location,
  const Point<dim> &                          particle_two_location,
  const double &                              dt,
  const double &                              particle_one_radius,
  const double &                              particle_two_radius,
  const double &                              particle_one_mass,
  const double &                              particle_two_mass)
{
  Point<3> particle_one_location_3d;
  Point<3> particle_two_location_3d;

  if constexpr (dim == 3)
    {
      particle_one_location_3d = particle_one_location;
      particle_two_location_3d = particle_two_location;
    }

  if constexpr (dim == 2)
    {
      particle_one_location_3d = point_nd_to_3d(particle_one_location);
      particle_two_location_3d = point_nd_to_3d(particle_two_location);
    }
  auto particle_one_properties = particle_one.get_properties();
  particle_one_properties[PropertiesIndex::mass] = particle_one_mass;
  particle_one_properties[PropertiesIndex::type] = 0;
  particle_one_properties[PropertiesIndex::dp]   = 2.0 * particle_one_radius;

  auto particle_two_properties = particle_one.get_properties();
  particle_two_properties[PropertiesIndex::mass] = particle_two_mass;
  particle_two_properties[PropertiesIndex::type] = 0;
  particle_two_properties[PropertiesIndex::dp]   = 2.0 * particle_two_radius;
  // DEM::PropertiesIndex::type is the first (0) property of particles in the
  // DEM solver. For the IB particles, the first property is ID. For force and
  // torque calculations, we need pair-wise properties (such as effective
  // Young's modulus, effective coefficient of restitution, etc.) We rewrite
  // these pair-wise properties by using the ID of IB particles (using
  // DEM::PropertiesIndex::type) and use them in force calculations.
  const unsigned int particle_one_type =
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];



  this->effective_youngs_modulus[particle_one_type][particle_two_type] =
    (particle_one.youngs_modulus * particle_two.youngs_modulus) /
    ((particle_two.youngs_modulus *
      (1.0 - particle_one.poisson_ratio * particle_one.poisson_ratio)) +
     (particle_one.youngs_modulus *
      (1.0 - particle_two.poisson_ratio * particle_two.poisson_ratio)) +
     DBL_MIN);

  this->effective_shear_modulus[particle_one_type][particle_two_type] =
    (particle_one.youngs_modulus * particle_two.youngs_modulus) /
    (2.0 * ((particle_two.youngs_modulus * (2.0 - particle_one.poisson_ratio) *
             (1.0 + particle_one.poisson_ratio)) +
            (particle_one.youngs_modulus * (2.0 - particle_two.poisson_ratio) *
             (1.0 + particle_two.poisson_ratio))) +
     DBL_MIN);

  this->effective_coefficient_of_restitution[particle_one_type]
                                            [particle_two_type] =
    harmonic_mean(particle_one.restitution_coefficient,
                  particle_two.restitution_coefficient);
  this
    ->effective_coefficient_of_friction[particle_one_type][particle_two_type] =
    harmonic_mean(particle_one.friction_coefficient,
                  particle_two.friction_coefficient);
  this->effective_coefficient_of_rolling_friction[particle_one_type]
                                                 [particle_two_type] =
    harmonic_mean(particle_one.rolling_friction_coefficient,
                  particle_two.rolling_friction_coefficient);

  // Since the normal overlap is already calculated we update
  // this element of the container here. The rest of information
  // are updated using the following function
  this->update_contact_information(contact_info,
                                   normal_relative_velocity_value,
                                   normal_unit_vector,
                                   particle_one_properties,
                                   particle_two_properties,
                                   particle_one_location_3d,
                                   particle_two_location_3d,
                                   dt);

  calculate_linear_contact_force_and_torque(contact_info,
                                            normal_relative_velocity_value,
                                            normal_unit_vector,
                                            normal_overlap,
                                            particle_one_properties,
                                            particle_two_properties,
                                            normal_force,
                                            tangential_force,
                                            particle_one_tangential_torque,
                                            particle_two_tangential_torque,
                                            rolling_resistance_torque);
}

// Calculates linear contact force and torques
template <int dim>
void
ParticleParticleLinearForce<dim>::calculate_linear_contact_force_and_torque(
  particle_particle_contact_info_struct<dim> &contact_info,
  const double &                              normal_relative_velocity_value,
  const Tensor<1, 3> &                        normal_unit_vector,
  const double &                              normal_overlap,
  const ArrayView<const double> &             particle_one_properties,
  const ArrayView<const double> &             particle_two_properties,
  Tensor<1, 3> &                              normal_force,
  Tensor<1, 3> &                              tangential_force,
  Tensor<1, 3> &                              particle_one_tangential_torque,
  Tensor<1, 3> &                              particle_two_tangential_torque,
  Tensor<1, 3> &                              rolling_resistance_torque)
{
  // Calculation of effective radius and mass
  this->find_effective_radius_and_mass(particle_one_properties,
                                       particle_two_properties);
  const unsigned int particle_one_type =
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties
  double normal_spring_constant =
    1.0667 * sqrt(this->effective_radius) *
    this->effective_youngs_modulus
      [particle_one_type][particle_two_properties[PropertiesIndex::type]] *
    pow(
      (0.9375 * this->effective_mass * normal_relative_velocity_value *
       normal_relative_velocity_value /
       (sqrt(this->effective_radius) *
        this->effective_youngs_modulus[particle_one_type][particle_two_type])),
      0.2);
  double tangential_spring_constant =
    1.0667 * sqrt(this->effective_radius) *
      this->effective_youngs_modulus[particle_one_type][particle_two_type] *
      pow((0.9375 * this->effective_mass *
           contact_info.tangential_relative_velocity *
           contact_info.tangential_relative_velocity /
           (sqrt(this->effective_radius) *
            this->effective_youngs_modulus[particle_one_type]
                                          [particle_two_type])),
          0.2) +
    DBL_MIN;
  double normal_damping_constant = sqrt(
    (4.0 * this->effective_mass * normal_spring_constant) /
    (1.0 +
     (M_PI /
      (log(this->effective_coefficient_of_restitution[particle_one_type]
                                                     [particle_two_type]) +
       DBL_MIN)) *
       (M_PI /
        (log(this->effective_coefficient_of_restitution[particle_one_type]
                                                       [particle_two_type]) +
         DBL_MIN))));
  double tangential_damping_constant =
    normal_damping_constant *
    sqrt(tangential_spring_constant / normal_spring_constant);

  // Calculation of normal force
  normal_force =
    ((normal_spring_constant * normal_overlap) * normal_unit_vector) +
    ((normal_damping_constant * normal_relative_velocity_value) *
     normal_unit_vector);

  // Calculation of tangential force. Since we need damping tangential force in
  // the gross sliding again, we define it as a separate variable
  Tensor<1, 3> damping_tangential_force =
    tangential_damping_constant * contact_info.tangential_relative_velocity;
  tangential_force =
    (tangential_spring_constant * contact_info.tangential_overlap) +
    damping_tangential_force;

  double coulomb_threshold =
    this->effective_coefficient_of_friction[particle_one_type]
                                           [particle_two_type] *
    normal_force.norm();


  // Check for gross sliding
  if (tangential_force.norm() > coulomb_threshold)
    {
      // Gross sliding occurs and the tangential overlap and tangential
      // force are limited to Coulumb's criterion
      contact_info.tangential_overlap =
        (coulomb_threshold *
           (tangential_force / (tangential_force.norm() + DBL_MIN)) -
         damping_tangential_force) /
        (tangential_spring_constant + DBL_MIN);

      tangential_force =
        (tangential_spring_constant * contact_info.tangential_overlap) +
        damping_tangential_force;
    }

  // Calculation of torque
  // Torque caused by tangential force (tangential_torque)
  particle_one_tangential_torque =
    cross_product_3d(normal_unit_vector,
                     tangential_force *
                       particle_one_properties[PropertiesIndex::dp] * 0.5);
  particle_two_tangential_torque =
    particle_one_tangential_torque *
    particle_two_properties[PropertiesIndex::dp] /
    particle_one_properties[PropertiesIndex::dp];

  // Rolling resistance torque
  if (rolling_reistance_model ==
      RollingResistanceTorqueModel::no_rolling_resistance)
    rolling_resistance_torque = no_rolling_resistance_torque(
      this->effective_radius,
      particle_one_properties,
      particle_two_properties,
      this->effective_coefficient_of_rolling_friction[particle_one_type]
                                                     [particle_two_type],
      normal_force.norm(),
      normal_unit_vector);
  if (rolling_reistance_model ==
      RollingResistanceTorqueModel::constant_rolling_resistance)
    rolling_resistance_torque = constant_rolling_resistance_torque(
      this->effective_radius,
      particle_one_properties,
      particle_two_properties,
      this->effective_coefficient_of_rolling_friction[particle_one_type]
                                                     [particle_two_type],
      normal_force.norm(),
      normal_unit_vector);
  if (rolling_reistance_model ==
      RollingResistanceTorqueModel::viscous_rolling_resistance)
    rolling_resistance_torque = viscous_rolling_resistance_torque(
      this->effective_radius,
      particle_one_properties,
      particle_two_properties,
      this->effective_coefficient_of_rolling_friction[particle_one_type]
                                                     [particle_two_type],
      normal_force.norm(),
      normal_unit_vector);
}

template class ParticleParticleLinearForce<2>;
template class ParticleParticleLinearForce<3>;
