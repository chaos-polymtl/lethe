#include <dem/pp_nonlinear_force.h>

using namespace DEM;

template <int dim>
PPNonLinearForce<dim>::PPNonLinearForce(
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
            ((youngs_modulus_j * (1 - poisson_ratio_i * poisson_ratio_i)) +
             (youngs_modulus_i * (1 - poisson_ratio_j * poisson_ratio_j)) +
             DBL_MIN);

          this->effective_shear_modulus[i].insert(
            {j,
             (youngs_modulus_i * youngs_modulus_j) /
               (2 * ((youngs_modulus_j * (2 - poisson_ratio_i) *
                      (1 + poisson_ratio_i)) +
                     (youngs_modulus_i * (2 - poisson_ratio_j) *
                      (1 + poisson_ratio_j))) +
                DBL_MIN)});

          this->effective_coefficient_of_restitution[i].insert(
            {j,
             2 * restitution_coefficient_i * restitution_coefficient_j /
               (restitution_coefficient_i + restitution_coefficient_j +
                DBL_MIN)});

          this->effective_coefficient_of_friction[i].insert(
            {j,
             2 * friction_coefficient_i * friction_coefficient_j /
               (friction_coefficient_i + friction_coefficient_j + DBL_MIN)});

          this->effective_coefficient_of_rolling_friction[i].insert(
            {j,
             2 * rolling_friction_coefficient_i *
               rolling_friction_coefficient_j /
               (rolling_friction_coefficient_i +
                rolling_friction_coefficient_j + DBL_MIN)});

          double restitution_coefficient_particle_log =
            std::log(this->effective_coefficient_of_restitution[i][j]);

          model_parameter_beta[i].insert(
            {j,
             restitution_coefficient_particle_log /
               sqrt(restitution_coefficient_particle_log *
                      restitution_coefficient_particle_log +
                    9.8696)});
        }

      // Overriding rolling resistance torque model
      if (dem_parameters.model_parameters.rolling_resistance_method ==
          Parameters::Lagrangian::ModelParameters::RollingResistanceMethod::
            no_resistance)
        {
          calculate_rolling_resistance_torque =
            &PPNonLinearForce<dim>::no_resistance;
        }
      else if (dem_parameters.model_parameters.rolling_resistance_method ==
               Parameters::Lagrangian::ModelParameters::
                 RollingResistanceMethod::constant_resistance)
        {
          calculate_rolling_resistance_torque =
            &PPNonLinearForce<dim>::constant_resistance;
        }
      else if (dem_parameters.model_parameters.rolling_resistance_method ==
               Parameters::Lagrangian::ModelParameters::
                 RollingResistanceMethod::viscous_resistance)
        {
          calculate_rolling_resistance_torque =
            &PPNonLinearForce<dim>::viscous_resistance;
        }
    }
}

template <int dim>
void
PPNonLinearForce<dim>::calculate_pp_contact_force(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
    &                          ghost_adjacent_particles,
  const double &               dt,
  std::vector<Tensor<1, dim>> &momentum,
  std::vector<Tensor<1, dim>> &force)
{
  // Updating contact force of particles for local-local and local-ghost contact
  // pairs are differnet. Consequently, contact forces of local-local and
  // local-ghost particle pairs are performed in separate loops

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
              // Getting information (location and propertis) of particle one
              // and two in contact
              auto             particle_one = contact_info.particle_one;
              auto             particle_two = contact_info.particle_two;
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
                    contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  this->calculate_nonlinear_contact_force_and_torque(
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

                  // Getting particles' momentum and force
#if DEAL_II_VERSION_GTE(10, 0, 0)
                  types::particle_index particle_one_id =
                    particle_one->get_local_index();
                  types::particle_index particle_two_id =
                    particle_two->get_local_index();
#else
                  types::particle_index particle_one_id =
                    particle_one->get_id();
                  types::particle_index particle_two_id =
                    particle_two->get_id();
#endif
                  Tensor<1, dim> &particle_one_momentum =
                    momentum[particle_one_id];
                  Tensor<1, dim> &particle_two_momentum =
                    momentum[particle_two_id];
                  Tensor<1, dim> &particle_one_force = force[particle_one_id];
                  Tensor<1, dim> &particle_two_force = force[particle_two_id];


                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_real(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    particle_two_tangential_torque,
                    rolling_resistance_torque,
                    particle_one_momentum,
                    particle_two_momentum,
                    particle_one_force,
                    particle_two_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d] = 0;
                    }
                }
            }
        }
    }

  // Doing the same calculations for local-ghost particle pairs

  // Looping over ghost_adjacent_particles list with iterator
  // adjacent_particles_list
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      if (!adjacent_particles_list.empty())
        {
          for (auto &&contact_info :
               adjacent_particles_list | boost::adaptors::map_values)

            {
              // Getting information (location and propertis) of particle one
              // and two in contact
              auto             particle_one = contact_info.particle_one;
              auto             particle_two = contact_info.particle_two;
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
                    contact_info,
                    normal_relative_velocity_value,
                    normal_unit_vector,
                    particle_one_properties,
                    particle_two_properties,
                    particle_one_location,
                    particle_two_location,
                    dt);

                  this->calculate_nonlinear_contact_force_and_torque(
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

                  // Getting momentum and force of particle one
#if DEAL_II_VERSION_GTE(10, 0, 0)
                  types::particle_index particle_one_id =
                    particle_one->get_local_index();
#else
                  types::particle_index particle_one_id =
                    particle_one->get_id();
#endif
                  Tensor<1, dim> &particle_one_momentum =
                    momentum[particle_one_id];
                  Tensor<1, dim> &particle_one_force = force[particle_one_id];

                  // Apply the calculated forces and torques on the particle
                  // pair
                  this->apply_force_and_torque_ghost(
                    normal_force,
                    tangential_force,
                    particle_one_tangential_torque,
                    rolling_resistance_torque,
                    particle_one_momentum,
                    particle_one_force);
                }

              else
                {
                  // if the adjacent pair is not in contact anymore, only the
                  // tangential overlap is set to zero
                  for (int d = 0; d < dim; ++d)
                    {
                      contact_info.tangential_overlap[d] = 0;
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
  pp_contact_info_struct<dim> &  contact_info,
  const double &                 normal_relative_velocity_value,
  const Tensor<1, dim> &         normal_unit_vector,
  const double &                 normal_overlap,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  Tensor<1, dim> &               normal_force,
  Tensor<1, dim> &               tangential_force,
  Tensor<1, dim> &               particle_one_tangential_torque,
  Tensor<1, dim> &               particle_two_tangential_torque,
  Tensor<1, dim> &               rolling_resistance_torque)
{
  // Calculation of effective radius and mass
  this->find_effective_radius_and_mass(particle_one_properties,
                                       particle_two_properties);

  const unsigned int particle_one_type =
    particle_one_properties[DEM::PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[DEM::PropertiesIndex::type];

  const double radius_times_overlap_sqrt =
    sqrt(this->effective_radius * normal_overlap);
  const double model_parameter_sn =
    2 * this->effective_youngs_modulus[particle_one_type][particle_two_type] *
    radius_times_overlap_sqrt;
  double model_parameter_st =
    8 * this->effective_shear_modulus[particle_one_type][particle_two_type] *
    radius_times_overlap_sqrt;

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties
  double normal_spring_constant = 0.66665 * model_parameter_sn;
  double normal_damping_constant =
    -1.8257 * model_parameter_beta[particle_one_type][particle_two_type] *
    sqrt(model_parameter_sn * this->effective_mass);
  double tangential_spring_constant =
    8 * this->effective_shear_modulus[particle_one_type][particle_two_type] *
      radius_times_overlap_sqrt +
    DBL_MIN;
  double tangential_damping_constant =
    normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

  // Calculation of normal force
  normal_force =
    ((normal_spring_constant * normal_overlap) * normal_unit_vector) +
    ((normal_damping_constant * normal_relative_velocity_value) *
     normal_unit_vector);

  // Calculation of tangential force. Since we need damping tangential force in
  // the gross sliding again, we define it as a separate variable
  Tensor<1, dim> damping_tangential_force =
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
      // Gross sliding occurs and the tangential overlap and tangnetial
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
  if (dim == 3)
    {
      particle_one_tangential_torque =
        cross_product_3d(normal_unit_vector,
                         tangential_force * particle_one_properties[DEM::dp] *
                           0.5);
      particle_two_tangential_torque =
        particle_one_tangential_torque *
        particle_two_properties[DEM::PropertiesIndex::dp] /
        particle_one_properties[DEM::PropertiesIndex::dp];
    }

  // Rolling resistance torque
  rolling_resistance_torque = (this->*calculate_rolling_resistance_torque)(
    this->effective_radius,
    particle_one_properties,
    particle_two_properties,
    this->effective_coefficient_of_rolling_friction[particle_one_type]
                                                   [particle_two_type],
    normal_force.norm(),
    normal_unit_vector);
}

template class PPNonLinearForce<2>;
template class PPNonLinearForce<3>;
