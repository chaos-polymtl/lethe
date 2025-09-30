// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_particle_contact_force_h
#define lethe_particle_particle_contact_force_h

#include <core/auxiliary_math_functions.h>

#include <dem/contact_info.h>
#include <dem/contact_type.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_heat_transfer.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/rolling_resistance_torque_models.h>

#include <boost/range/adaptor/map.hpp>

#include <vector>

using namespace dealii;

/**
 * @brief Base class for the particle-particle contact force models
 * This class does not implement any of the models, but ensures that
 * an interface without template specialization is available. All the
 * actual implementation of the models are carried out in the
 * ParticleParticleContactForce class which is templated by the contact model
 * type.
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class ParticleParticleContactForceBase
{
public:
  /**
   * @brief Calculate the contact outcomes using the contact pair information
   * obtained in the fine search and physical properties of particles.
   *
   * @param local_adjacent_particles Container of the contact pair candidates
   * information for calculation of the local particle-particle contact forces.
   * @param ghost_adjacent_particles Container of the contact pair candidates
   * information for calculation of the local-ghost particle-particle contact
   * forces.
   * @param local_local_periodic_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local periodic
   * particle-particle contact forces.
   * @param local_ghost_periodic_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local-ghost periodic
   * particle-particle contact forces.
   * @param ghost_local_periodic_adjacent_particles Container of the contact pair
   * candidates information for calculation of the ghost-local periodic
   * particle-particle contact forces.
   * @param dt DEM time step.
   * @param contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_particle_contact(
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_local_periodic_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_ghost_periodic_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
                &ghost_local_periodic_adjacent_particles,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) = 0;

  void
  set_periodic_offset(const Tensor<1, dim> &periodic_offset)
  {
    this->periodic_offset = periodic_offset;
  }

protected:
  Tensor<1, dim> periodic_offset;
};

/**
 * @brief Execute calculation of particle-particle contact forces including
 * linear, non-linear contact models, and non-linear with cohesive forces.
 *
 * Instead of using an inheritance hierarchy to distinguish between
 * the contact model, the class is templated with the type of force model
 * and rolling friction model. Consequently, the code for each
 * combination of force model is generated at compile time.
 *
 * @tparam dim The dimension of the problem.
 * @tparam contact_model The particle-particle contact force model.
 * @tparam rolling_friction_model The rolling resistance model.
 */
template <
  int dim,
  typename PropertiesIndex,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
class ParticleParticleContactForce
  : public ParticleParticleContactForceBase<dim, PropertiesIndex>
{
public:
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters);

  virtual ~ParticleParticleContactForce()
  {}

  /**
   * @brief Calculate the contact outcomes using the contact pair information
   * obtained in the fine search and physical properties of particles.
   *
   * @param local_adjacent_particles Container of the contact pair candidates
   * information for calculation of the local particle-particle contact forces.
   * @param ghost_adjacent_particles Container of the contact pair candidates
   * information for calculation of the local-ghost particle-particle contact
   * forces.
   * @param local_local_periodic_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local periodic
   * particle-particle contact forces.
   * @param local_ghost_periodic_adjacent_particles Container of the contact pair
   * candidates information for calculation of the local-ghost periodic
   * particle-particle contact forces.
   * @param ghost_local_periodic_adjacent_particles Container of the contact pair
   * candidates information for calculation of the ghost-local periodic
   * particle-particle contact forces.
   * @param dt DEM time step.
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_particle_contact(
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_local_periodic_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_ghost_periodic_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
                &ghost_local_periodic_adjacent_particles,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) override;

protected:
  /**
   * @brief Update the contact pair information for all contact force
   * calculations.
   *
   * @param[out] contact_info Contact information of a particle pair in
   * neighborhood.
   * @param[out] tangential_relative_velocity Tangential relative velocity.
   * @param[out] normal_relative_velocity_value Normal relative contact
   * velocity.
   * @param[out] normal_unit_vector Normal vector of the contact.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[in] particle_one_location Location of particle one in contact.
   * @param[in] particle_two_location Location of particle two in contact.
   * @param[in] dt DEM time step.
   */
  inline void
  update_contact_information(
    particle_particle_contact_info<dim> &contact_info,
    Tensor<1, 3>                        &tangential_relative_velocity,
    double                              &normal_relative_velocity_value,
    Tensor<1, 3>                        &normal_unit_vector,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    const Point<3>                      &particle_one_location,
    const Point<3>                      &particle_two_location,
    const double                         dt)
  {
    // Calculation of the contact vector from particle one to particle two
    auto contact_vector = particle_two_location - particle_one_location;

    // Calculation of the normal unit contact vector
    normal_unit_vector = contact_vector / contact_vector.norm();

    // Defining velocities and angular velocities of particles one and
    // two as vectors
    Tensor<1, 3> particle_one_omega, particle_two_omega;

    // Defining relative contact velocity
    Tensor<1, 3> contact_relative_velocity;

    // Assigning velocities and angular velocities of particles
    contact_relative_velocity[0] =
      particle_one_properties[PropertiesIndex::v_x] -
      particle_two_properties[PropertiesIndex::v_x];
    contact_relative_velocity[1] =
      particle_one_properties[PropertiesIndex::v_y] -
      particle_two_properties[PropertiesIndex::v_y];
    contact_relative_velocity[2] =
      particle_one_properties[PropertiesIndex::v_z] -
      particle_two_properties[PropertiesIndex::v_z];

    particle_one_omega[0] = particle_one_properties[PropertiesIndex::omega_x];
    particle_one_omega[1] = particle_one_properties[PropertiesIndex::omega_y];
    particle_one_omega[2] = particle_one_properties[PropertiesIndex::omega_z];

    particle_two_omega[0] = particle_two_properties[PropertiesIndex::omega_x];
    particle_two_omega[1] = particle_two_properties[PropertiesIndex::omega_y];
    particle_two_omega[2] = particle_two_properties[PropertiesIndex::omega_z];

    // Calculation of contact relative velocity
    // v_ij = (v_i - v_j) + (R_i*omega_i + R_j*omega_j) × n_ij
    contact_relative_velocity += (cross_product_3d(
      0.5 * (particle_one_properties[PropertiesIndex::dp] * particle_one_omega +
             particle_two_properties[PropertiesIndex::dp] * particle_two_omega),
      normal_unit_vector));

    // Calculation of normal relative velocity. Note that in the
    // following line the product acts as inner product since both
    // sides are vectors, while in the second line the product is
    // scalar and vector product
    normal_relative_velocity_value =
      contact_relative_velocity * normal_unit_vector;

    // Calculation of tangential relative velocity
    // v_rt = v_ij - (v_ij⋅n_ij)*n_ij
    tangential_relative_velocity =
      contact_relative_velocity -
      (normal_relative_velocity_value * normal_unit_vector);

    // Calculation of new tangential_displacement, since this value is history
    // dependent it needs the value at previous time-step. This variable is the
    // main reason that we have iteration over  two different vectors :
    // tangential_displacement of the particles which were already in contact
    // needs to modified using its history, while the tangential_displacements
    // of new particles are equal to zero delta_t_new = delta_t_old + v_rt*dt
    contact_info.tangential_displacement += tangential_relative_velocity * dt;
    contact_info.tangential_displacement -=
      (contact_info.tangential_displacement * normal_unit_vector) *
      normal_unit_vector;
  }

  /**
   * @brief Get the location of the particle.
   *
   * @param particle The particle to get the location from.
   */
  inline Point<3>
  get_location(const Particles::ParticleIterator<dim> &particle) &
  {
    if constexpr (dim == 3)
      return particle->get_location();

    if constexpr (dim == 2)
      return point_nd_to_3d(particle->get_location());
  }

  /**
   * @brief Get the shifted location of the particle on the periodic boundary.
   *
   * @param particle The particle to get the location from.
   */
  inline Point<3>
  get_periodic_location(const Particles::ParticleIterator<dim> &particle) &
  {
    if constexpr (dim == 3)
      return (particle->get_location() - this->periodic_offset);

    if constexpr (dim == 2)
      return point_nd_to_3d(particle->get_location() - this->periodic_offset);
  }

  /**
   * @brief Calculate the particle-particle contact force and torque
   * according to the contact model.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_contact(particle_particle_contact_info<dim> &contact_info,
                    const Tensor<1, 3> &tangential_relative_velocity,
                    const double        normal_relative_velocity_value,
                    const Tensor<1, 3> &normal_unit_vector,
                    const double        normal_overlap,
                    const double        dt,
                    const ArrayView<const double> &particle_one_properties,
                    const ArrayView<const double> &particle_two_properties,
                    Tensor<1, 3>                  &normal_force,
                    Tensor<1, 3>                  &tangential_force,
                    Tensor<1, 3> &particle_one_tangential_torque,
                    Tensor<1, 3> &particle_two_tangential_torque,
                    Tensor<1, 3> &rolling_resistance_torque)
  {
    using namespace Parameters::Lagrangian;

    if constexpr (contact_model == ParticleParticleContactForceModel::linear)
      {
        calculate_linear_contact(contact_info,
                                 tangential_relative_velocity,
                                 normal_relative_velocity_value,
                                 normal_unit_vector,
                                 normal_overlap,
                                 dt,
                                 particle_one_properties,
                                 particle_two_properties,
                                 normal_force,
                                 tangential_force,
                                 particle_one_tangential_torque,
                                 particle_two_tangential_torque,
                                 rolling_resistance_torque);
      }

    if constexpr (contact_model == ParticleParticleContactForceModel::hertz)
      {
        this->calculate_hertz_contact(contact_info,
                                      tangential_relative_velocity,
                                      normal_relative_velocity_value,
                                      normal_unit_vector,
                                      normal_overlap,
                                      dt,
                                      particle_one_properties,
                                      particle_two_properties,
                                      normal_force,
                                      tangential_force,
                                      particle_one_tangential_torque,
                                      particle_two_tangential_torque,
                                      rolling_resistance_torque);
      }

    if constexpr (contact_model ==
                  ParticleParticleContactForceModel::hertz_mindlin_limit_force)
      {
        calculate_hertz_mindlin_limit_force_contact(
          contact_info,
          tangential_relative_velocity,
          normal_relative_velocity_value,
          normal_unit_vector,
          normal_overlap,
          dt,
          particle_one_properties,
          particle_two_properties,
          normal_force,
          tangential_force,
          particle_one_tangential_torque,
          particle_two_tangential_torque,
          rolling_resistance_torque);
      }

    if constexpr (contact_model == ParticleParticleContactForceModel::
                                     hertz_mindlin_limit_overlap)
      {
        calculate_hertz_mindlin_limit_overlap_contact(
          contact_info,
          tangential_relative_velocity,
          normal_relative_velocity_value,
          normal_unit_vector,
          normal_overlap,
          dt,
          particle_one_properties,
          particle_two_properties,
          normal_force,
          tangential_force,
          particle_one_tangential_torque,
          particle_two_tangential_torque,
          rolling_resistance_torque);
      }


    if constexpr (contact_model == ParticleParticleContactForceModel::hertz_JKR)
      {
        this->calculate_hertz_JKR_contact(contact_info,
                                          tangential_relative_velocity,
                                          normal_relative_velocity_value,
                                          normal_unit_vector,
                                          normal_overlap,
                                          dt,
                                          particle_one_properties,
                                          particle_two_properties,
                                          normal_force,
                                          tangential_force,
                                          particle_one_tangential_torque,
                                          particle_two_tangential_torque,
                                          rolling_resistance_torque);
      }

    if constexpr (contact_model == ParticleParticleContactForceModel::DMT)
      {
        this->calculate_DMT_contact(contact_info,
                                    tangential_relative_velocity,
                                    normal_relative_velocity_value,
                                    normal_unit_vector,
                                    normal_overlap,
                                    dt,
                                    particle_one_properties,
                                    particle_two_properties,
                                    normal_force,
                                    tangential_force,
                                    particle_one_tangential_torque,
                                    particle_two_tangential_torque,
                                    rolling_resistance_torque);
      }
  }

  /**
   * @brief Calculate the minimum overlap at which particle-particle forces are
   * computed.
   *
   * @return minimum overlap for the force calculation.
   */
  inline double
  get_force_calculation_threshold_distance()
  {
    if constexpr (contact_model == Parameters::Lagrangian::
                                     ParticleParticleContactForceModel::DMT)
      {
        // We are looking for the maximum Hamaker constant and minimum surface
        // energy to compute the biggest distance at which force will be
        // computed. In other words, we are maximising the delta_0.
        // This way, force_calculation_threshold_distance can be set
        // to const variable, which should make the code faster.
        const double max_effective_hamaker_constant =
          *(std::max_element(this->effective_hamaker_constant.begin(),
                             this->effective_hamaker_constant.end()));
        const double min_effective_surface_energy =
          *(std::min_element(this->effective_surface_energy.begin(),
                             this->effective_surface_energy.end()));

        // Return the critical delta_0. We put a minus sign in front of since a
        // positive overlap means that particles are in contact.
        return -std::sqrt(
          max_effective_hamaker_constant /
          (12. * M_PI * min_effective_surface_energy * dmt_cut_off_threshold));
      }
    return 0.;
  }

private:
  /**
   * @brief Apply the calculated force and torque on the local-local
   * particle pair in contact.
   *
   * @param[in] normal_force Contact normal force.
   * @param[in] tangential_force Contact tangential force.
   * @param[in] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[in] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[in] rolling_resistance_torque Contact rolling resistance torque.
   * @param[in,out] particle_one_torque
   * @param[in,out] particle_two_torque
   * @param[in,out] particle_one_force Force acting on particle one
   * @param[in,out] particle_two_force Force acting on particle two
   */
  inline void
  apply_force_and_torque_on_local_particles(
    const Tensor<1, 3> &normal_force,
    const Tensor<1, 3> &tangential_force,
    const Tensor<1, 3> &particle_one_tangential_torque,
    const Tensor<1, 3> &particle_two_tangential_torque,
    const Tensor<1, 3> &rolling_resistance_torque,
    Tensor<1, 3>       &particle_one_torque,
    Tensor<1, 3>       &particle_two_torque,
    Tensor<1, 3>       &particle_one_force,
    Tensor<1, 3>       &particle_two_force)
  {
    // Calculation of total force
    Tensor<1, 3> total_force = normal_force + tangential_force;

    // Updating the force and torque of particles in the particle handler
    particle_one_force -= total_force;
    particle_two_force += total_force;
    particle_one_torque +=
      -particle_one_tangential_torque + rolling_resistance_torque;
    particle_two_torque +=
      -particle_two_tangential_torque - rolling_resistance_torque;
  }

  /**
   * @brief Apply the calculated force and torque on the local-ghost particle
   * pair in contact. The contact force is only applied on the local particles
   *
   * @param[in] normal_force Contact normal force.
   * @param[in] tangential_force Contact tangential force.
   * @param[in] particle_one_tangential_torque Contact tangential torque.
   * @param[in] rolling_resistance_torque Contact rolling resistance torque.
   * @param[in,out] particle_one_torque Torque acting on particle one (local).
   * @param[in,out] particle_one_force Force acting on particle one (local).
   */
  inline void
  apply_force_and_torque_on_single_local_particle(
    const Tensor<1, 3> &normal_force,
    const Tensor<1, 3> &tangential_force,
    const Tensor<1, 3> &particle_one_tangential_torque,
    const Tensor<1, 3> &rolling_resistance_torque,
    Tensor<1, 3>       &particle_one_torque,
    Tensor<1, 3>       &particle_one_force)
  {
    // Updating the force and torque acting on particles in the particle handler
    particle_one_force -= normal_force + tangential_force;
    particle_one_torque +=
      -particle_one_tangential_torque + rolling_resistance_torque;
  }

  /**
   * @brief Calculate the effective mass and the effective radius of particle
   * pair.
   *
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   */
  inline std::tuple<double, double>
  find_effective_radius_and_mass(
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties)
  {
    // Calculate the effective radius
    const double diameter_one = particle_one_properties[PropertiesIndex::dp];
    const double diameter_two = particle_two_properties[PropertiesIndex::dp];
    const double sum_diameter = diameter_one + diameter_two;
    double       effective_radius =
      (diameter_one * diameter_two) / (2 * sum_diameter);

    // Calculate the effective mass
    const double mass_one = particle_one_properties[PropertiesIndex::mass];
    const double mass_two = particle_two_properties[PropertiesIndex::mass];
    const double sum_mass = mass_one + mass_two;
    double       effective_mass = (mass_one * mass_two) / sum_mass;

    return std::make_tuple(effective_radius, effective_mass);
  }

  /**
   * @brief Calculate the rolling resistance torque acting on the particles
   * according to the rolling resistance model.
   *
   * @param[in] effective_radius Effective radius of the particle pair.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[in] rolling_friction_coeff Effective rolling friction coefficient.
   * @param[in] rolling_viscous_damping_coeff Effective rolling viscous damping
   * coefficient
   * @param[in] f_coeff Model parameter for the EPSD model.
   * @param[in] dt DEM time step.
   * @param[in] normal_spring_constant normal contact stiffness constant.
   * @param[in] normal_force_norm Norm of the normal force.
   * @param[in] normal_unit_vector Normal unit vector between particles in
   * contact.
   * @param[in,out] cumulative_rolling_resistance_spring_torque Cumulative
   * spring rolling resistance torque applied between particle one and two.
   *
   */
  inline Tensor<1, 3>
  calculate_rolling_resistance_torque(
    [[maybe_unused]] const double                   effective_radius,
    [[maybe_unused]] const ArrayView<const double> &particle_one_properties,
    [[maybe_unused]] const ArrayView<const double> &particle_two_properties,
    [[maybe_unused]] const double                   rolling_friction_coeff,
    [[maybe_unused]] const double        rolling_viscous_damping_coeff,
    [[maybe_unused]] const double        dt,
    [[maybe_unused]] const double        normal_spring_constant,
    [[maybe_unused]] const double        normal_force_norm,
    [[maybe_unused]] const Tensor<1, 3> &normal_unit_vector,
    Tensor<1, 3> &cumulative_rolling_resistance_spring_torque)
  {
    using namespace Parameters::Lagrangian;

    if constexpr (rolling_friction_model == none)
      {
        return no_rolling_resistance_torque();
      }

    if constexpr (rolling_friction_model == constant)
      {
        return constant_rolling_resistance_torque<PropertiesIndex>(
          effective_radius,
          particle_one_properties,
          particle_two_properties,
          rolling_friction_coeff,
          normal_force_norm);
      }

    if constexpr (rolling_friction_model == viscous)
      {
        return viscous_rolling_resistance_torque<PropertiesIndex>(
          effective_radius,
          particle_one_properties,
          particle_two_properties,
          rolling_friction_coeff,
          normal_force_norm,
          normal_unit_vector);
      }
    if constexpr (rolling_friction_model == epsd)
      {
        return epsd_rolling_resistance_torque<dim, PropertiesIndex>(
          effective_radius,
          particle_one_properties,
          particle_two_properties,
          rolling_friction_coeff,
          rolling_viscous_damping_coeff,
          f_coefficient_epsd,
          normal_force_norm,
          dt,
          normal_spring_constant,
          normal_unit_vector,
          cumulative_rolling_resistance_spring_torque);
      }
  }

  /**
   * @brief Calculate the particle-particle linear contact force and torque
   * based on the updated values in contact_info.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_linear_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const double                         dt,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    auto [effective_radius, effective_mass] =
      find_effective_radius_and_mass(particle_one_properties,
                                     particle_two_properties);

    // Get the reference of the effective properties according to the particle
    // types in vectors for easy-to-read equations.
    const unsigned int particle_one_type =
      static_cast<unsigned int>(particle_one_properties[PropertiesIndex::type]);
    const unsigned int particle_two_type =
      static_cast<unsigned int>(particle_two_properties[PropertiesIndex::type]);
    const unsigned int pair_index =
      vec_particle_type_index(particle_one_type, particle_two_type);

    const double youngs_modulus = this->effective_youngs_modulus[pair_index];
    const double beta           = this->model_parameter_beta[pair_index];
    const double rolling_viscous_damping_coeff =
      this->effective_rolling_viscous_damping_coefficient[pair_index];
    const double friction_coeff =
      this->effective_coefficient_of_friction[pair_index];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[pair_index];

    // Get particle diameter references
    const double &diameter_one = particle_one_properties[PropertiesIndex::dp];
    const double &diameter_two = particle_two_properties[PropertiesIndex::dp];

    // Characteristic velocity is set at 1.0 so that the normal and tangential
    // spring constant remain constant throughout a simulation.
    constexpr double characteristic_velocity = 1.0;

    // Pre-calculate common terms to reduce computations
    const double sqrt_effective_radius = sqrt(effective_radius);
    const double effective_mass_times_vel_sq =
      effective_mass * characteristic_velocity * characteristic_velocity;

    // Calculate the normal spring constant using the following formula:
    // kn = 16/15 * sqrt(Re) * Ye * (15/16 * (me * vc^2 / (sqrt(R) * Ye))^0.2
    const double normal_spring_constant =
      1.0667 * sqrt_effective_radius * youngs_modulus *
      pow((0.9375 * effective_mass_times_vel_sq /
           (sqrt_effective_radius * youngs_modulus)),
          0.2);

    // Calculate the tangential spring constant
    // REF :  R. Garg, J. Galvin-Carney, T. Li, and S. Pannala,“Documentation of
    // open-source MFIX–DEM software for gas-solids flows,” Tingwen Li Dr., p.
    // 10, Sep. 2012.
    double tangential_spring_constant = normal_spring_constant * 0.4;

    // Calculate the normal damping constants
    // -2 * beta * sqrt(m * kn)
    double normal_damping_constant =
      -2 * beta * sqrt(effective_mass * normal_spring_constant);

    // Calculate the tangential damping constant
    double tangential_damping_constant =
      normal_damping_constant * 0.6324555320336759; // sqrt(0.4)

    // Calculation of the normal force
    const double normal_force_value =
      normal_spring_constant * normal_overlap +
      normal_damping_constant * normal_relative_velocity_value;
    normal_force = normal_force_value * normal_unit_vector;

    // Calculation of tangential force. Since we need damping tangential force
    // in the gross sliding again, we define it as a separate variable
    Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;

    tangential_force =
      (tangential_spring_constant * contact_info.tangential_displacement) +
      damping_tangential_force;

    double coulomb_threshold = friction_coeff * normal_force_value;

    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement is recalculated
        // from the tangential force limited to Coulomb's criterion
        const Tensor<1, 3> limited_tangential_force =
          coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));

        // Calculate the tangential displacement from the limited tangential
        // force minus the spring tangential force divided by the spring
        // constant.
        contact_info.tangential_displacement =
          (limited_tangential_force - damping_tangential_force) /
          (tangential_spring_constant + DBL_MIN);

        // Recalculate the tangential force using the new tangential
        // displacement.
        tangential_force =
          (tangential_spring_constant * contact_info.tangential_displacement) +
          damping_tangential_force;
      }

    // Calculation of caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force * diameter_one * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque * diameter_two / diameter_one;

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      effective_radius,
      particle_one_properties,
      particle_two_properties,
      rolling_viscous_damping_coeff,
      rolling_friction_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      normal_unit_vector,
      contact_info.rolling_resistance_spring_torque);
  }

  /**
   * @brief Calculate the particle-particle non-linear contact force and torque
   * based on the updated values in contact_info. It uses the Hertz-Mindlin
   * contact model. If gross sliding occurs, the tangential displacement is
   * recalculated from the limited tangential to the Coulomb's criterion.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_hertz_mindlin_limit_overlap_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const double                         dt,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    auto [effective_radius, effective_mass] =
      find_effective_radius_and_mass(particle_one_properties,
                                     particle_two_properties);

    const unsigned int particle_one_type =
      static_cast<unsigned int>(particle_one_properties[PropertiesIndex::type]);
    const unsigned int particle_two_type =
      static_cast<unsigned int>(particle_two_properties[PropertiesIndex::type]);
    const unsigned int pair_index =
      vec_particle_type_index(particle_one_type, particle_two_type);

    const double youngs_modulus = this->effective_youngs_modulus[pair_index];
    const double shear_modulus  = this->effective_shear_modulus[pair_index];
    const double beta           = this->model_parameter_beta[pair_index];
    const double friction_coeff =
      this->effective_coefficient_of_friction[pair_index];
    const double rolling_viscous_damping_coeff =
      this->effective_rolling_viscous_damping_coefficient[pair_index];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[pair_index];

    // Get particle diameter references;
    const double &diameter_one = particle_one_properties[PropertiesIndex::dp];
    const double &diameter_two = particle_two_properties[PropertiesIndex::dp];

    // Calculate intermediate model parameters
    const double radius_times_overlap_sqrt =
      sqrt(effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * youngs_modulus * radius_times_overlap_sqrt;
    const double model_parameter_st =
      8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculation of normal and tangential spring and dashpot constants
    // using particle properties

    // Calculate the normal spring constant using the following formula:
    // kn = 4/3 * Ye * sqrt(Re * delta_n)
    // TODO this is 0.66667
    const double normal_spring_constant = 0.66665 * model_parameter_sn;

    // Calculate the normal damping constants from the following equation:
    // eta_n = -2 * sqrt(5/6) * beta * sqrt(Sn * me)
    const double normal_damping_constant =
      -1.8257 * beta * sqrt(model_parameter_sn * effective_mass);

    // Calculate the tangential spring constant
    // kt = 8 * Ge * sqrt(Re * delta_n)
    const double tangential_spring_constant = model_parameter_st;

    // Calculate the tangential damping constant from ratio with the normal
    // damping constant, but the equation is:
    // eta_t = -2 * sqrt(5/6) * beta * sqrt(St * me)
    const double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of normal force
    const double normal_force_value =
      normal_spring_constant * normal_overlap +
      normal_damping_constant * normal_relative_velocity_value;
    normal_force = normal_force_value * normal_unit_vector;

    // Calculation of tangential force. Since we need damping tangential force
    // in the gross sliding again, we define it as a separate variable
    const Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;
    tangential_force =
      (tangential_spring_constant * contact_info.tangential_displacement) +
      damping_tangential_force;

    const double coulomb_threshold = friction_coeff * normal_force_value;

    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement is recalculated
        // from the tangential force limited to Coulomb's criterion
        const Tensor<1, 3> limited_tangential_force =
          coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));

        // Calculate the tangential displacement from the limited tangential
        // force minus the spring tangential force divided by the spring
        // constant.
        contact_info.tangential_displacement =
          (limited_tangential_force - damping_tangential_force) /
          (tangential_spring_constant + DBL_MIN);

        // Recalculate the tangential force using the new tangential
        // displacement.
        tangential_force =
          (tangential_spring_constant * contact_info.tangential_displacement) +
          damping_tangential_force;
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force * diameter_one * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque * diameter_two / diameter_one;

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      effective_radius,
      particle_one_properties,
      particle_two_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      normal_unit_vector,
      contact_info.rolling_resistance_spring_torque);
  }

  /**
   * @brief Calculate the particle-particle non-linear contact force and torque
   * based on the updated values in contact_info. It uses the Hertz-Mindlin
   * contact model. If gross sliding occurs, the tangential displacement and
   * tangential force are limited to Coulomb's criterion.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_hertz_mindlin_limit_force_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const double                         dt,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    auto [effective_radius, effective_mass] =
      find_effective_radius_and_mass(particle_one_properties,
                                     particle_two_properties);

    const unsigned int particle_one_type =
      static_cast<unsigned int>(particle_one_properties[PropertiesIndex::type]);
    const unsigned int particle_two_type =
      static_cast<unsigned int>(particle_two_properties[PropertiesIndex::type]);
    const unsigned int pair_index =
      vec_particle_type_index(particle_one_type, particle_two_type);

    const double youngs_modulus = this->effective_youngs_modulus[pair_index];
    const double shear_modulus  = this->effective_shear_modulus[pair_index];
    const double beta           = this->model_parameter_beta[pair_index];
    const double friction_coeff =
      this->effective_coefficient_of_friction[pair_index];
    const double rolling_viscous_damping_coeff =
      this->effective_rolling_viscous_damping_coefficient[pair_index];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[pair_index];

    // Get particle diameter references;
    const double &diameter_one = particle_one_properties[PropertiesIndex::dp];
    const double &diameter_two = particle_two_properties[PropertiesIndex::dp];

    // Calculate intermediate model parameters
    const double radius_times_overlap_sqrt =
      sqrt(effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * youngs_modulus * radius_times_overlap_sqrt;
    double model_parameter_st = 8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculate the normal spring constant using the following formula:
    // kn = 4/3 * Ye * sqrt(Re * delta_n)
    // TODO this is 0.66667
    double normal_spring_constant = 0.66665 * model_parameter_sn;

    // Calculate the normal damping constants from the following equation:
    // eta_n = -2 * sqrt(5/6) * beta * sqrt(Sn * me)
    double normal_damping_constant =
      -1.8257 * beta * sqrt(model_parameter_sn * effective_mass);

    // Calculate the tangential spring constant
    // kt = 8 * Ge * sqrt(Re * delta_n)
    double tangential_spring_constant =
      8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculate the tangential damping constant from ratio with the normal
    // damping constant, but the equation is:
    // eta_t = -2 * sqrt(5/6) * beta * sqrt(St * me)
    double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of normal force using spring and dashpot normal forces
    const double normal_force_value =
      normal_spring_constant * normal_overlap +
      normal_damping_constant * normal_relative_velocity_value;
    normal_force = normal_force_value * normal_unit_vector;

    // Calculation of tangential force. Since we need damping tangential force
    // in the gross sliding again, we define it as a separate variable
    Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;
    tangential_force =
      (tangential_spring_constant * contact_info.tangential_displacement) +
      damping_tangential_force;

    double coulomb_threshold = friction_coeff * normal_force_value;

    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement and tangential
        // force are limited to Coulomb's criterion
        tangential_force =
          coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force * diameter_one * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque * diameter_two / diameter_one;

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      effective_radius,
      particle_one_properties,
      particle_two_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      normal_unit_vector,
      contact_info.rolling_resistance_spring_torque);
  }

  /**
   * @brief Calculate the particle-particle non-linear contact force and torque
   * based on the updated values in contact_info. It uses the Hertz contact
   * model, same as the Hertz-Mindlin with limit of the tangential force
   * without the tangential damping force.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_hertz_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3> & /*tangential_relative_velocity*/,
    const double                   normal_relative_velocity_value,
    const Tensor<1, 3>            &normal_unit_vector,
    const double                   normal_overlap,
    const double                   dt,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    Tensor<1, 3>                  &normal_force,
    Tensor<1, 3>                  &tangential_force,
    Tensor<1, 3>                  &particle_one_tangential_torque,
    Tensor<1, 3>                  &particle_two_tangential_torque,
    Tensor<1, 3>                  &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    auto [effective_radius, effective_mass] =
      find_effective_radius_and_mass(particle_one_properties,
                                     particle_two_properties);

    const unsigned int particle_one_type =
      static_cast<unsigned int>(particle_one_properties[PropertiesIndex::type]);
    const unsigned int particle_two_type =
      static_cast<unsigned int>(particle_two_properties[PropertiesIndex::type]);
    const unsigned int pair_index =
      vec_particle_type_index(particle_one_type, particle_two_type);

    const double youngs_modulus = this->effective_youngs_modulus[pair_index];
    const double shear_modulus  = this->effective_shear_modulus[pair_index];
    const double beta           = this->model_parameter_beta[pair_index];
    const double friction_coeff =
      this->effective_coefficient_of_friction[pair_index];
    const double rolling_viscous_damping_coeff =
      this->effective_rolling_viscous_damping_coefficient[pair_index];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[pair_index];

    // Get particle diameter references;
    const double &diameter_one = particle_one_properties[PropertiesIndex::dp];
    const double &diameter_two = particle_two_properties[PropertiesIndex::dp];

    // Calculate intermediate model parameters
    const double radius_times_overlap_sqrt =
      sqrt(effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * youngs_modulus * radius_times_overlap_sqrt;

    // Calculate the normal spring constant using the following formula:
    // kn = 4/3 * Ye * sqrt(Re * delta_n)
    // TODO this is 0.666667
    double normal_spring_constant = 0.66665 * model_parameter_sn;

    // Calculate the normal damping constants from the following equation:
    // eta_n = -2 * sqrt(5/6) * beta * sqrt(Sn * me)
    double normal_damping_constant =
      -1.8257 * beta * sqrt(model_parameter_sn * effective_mass);

    // Calculate the tangential spring constant
    // kt = 8 * Ge * sqrt(Re * delta_n)
    double tangential_spring_constant =
      8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculation of normal force using spring and dashpot normal constants.
    const double normal_force_value =
      normal_spring_constant * normal_overlap +
      normal_damping_constant * normal_relative_velocity_value;
    normal_force = normal_force_value * normal_unit_vector;

    // Calculation of tangential force using the spring tangential constant.
    tangential_force =
      tangential_spring_constant * contact_info.tangential_displacement;

    double coulomb_threshold = friction_coeff * normal_force_value;

    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement and tangential
        // force are limited to Coulomb's criterion
        tangential_force =
          coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force * diameter_one * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque * diameter_two / diameter_one;

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      effective_radius,
      particle_one_properties,
      particle_two_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      normal_unit_vector,
      contact_info.rolling_resistance_spring_torque);
  }
  /**
   * @brief Calculate the particle-particle contact and cohesive forces and
   * contact torque based on the updated values in contact_info. It uses the JKR
   * cohesive force model.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_hertz_JKR_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const double                         dt,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    auto [effective_radius, effective_mass] =
      find_effective_radius_and_mass(particle_one_properties,
                                     particle_two_properties);

    const unsigned int particle_one_type =
      static_cast<unsigned int>(particle_one_properties[PropertiesIndex::type]);
    const unsigned int particle_two_type =
      static_cast<unsigned int>(particle_two_properties[PropertiesIndex::type]);
    const unsigned int pair_index =
      vec_particle_type_index(particle_one_type, particle_two_type);

    const double youngs_modulus = this->effective_youngs_modulus[pair_index];
    const double shear_modulus  = this->effective_shear_modulus[pair_index];
    const double beta           = this->model_parameter_beta[pair_index];
    const double friction_coeff =
      this->effective_coefficient_of_friction[pair_index];
    const double rolling_viscous_damping_coeff =
      this->effective_rolling_viscous_damping_coefficient[pair_index];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[pair_index];
    const double surface_energy = this->effective_surface_energy[pair_index];

    // Get particle diameter references;
    const double &diameter_one = particle_one_properties[PropertiesIndex::dp];
    const double &diameter_two = particle_two_properties[PropertiesIndex::dp];

    // Calculate intermediate model parameters
    const double radius_times_overlap_sqrt =
      sqrt(effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * youngs_modulus * radius_times_overlap_sqrt;
    double model_parameter_st = 8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculation of the  contact path radius using the Ferrari analitycal
    // solution.
    const double c0 =
      Utilities::fixed_power<2>(effective_radius * normal_overlap);
    const double c1 = -2. * Utilities::fixed_power<2>(effective_radius) * M_PI *
                      surface_energy / youngs_modulus;

    const double c2 = -2. * normal_overlap * effective_radius;
    const double P  = -Utilities::fixed_power<2>(c2) / 12. - c0;
    const double Q  = -Utilities::fixed_power<3>(c2) / 108. + c0 * c2 / 3. -
                     Utilities::fixed_power<2>(c1) * 0.125;
    const double root1 =
      0.25 * Utilities::fixed_power<2>(Q) + Utilities::fixed_power<3>(P) / 27.;
    const double U      = std::cbrt(-0.5 * Q + std::sqrt(root1));
    const double s      = -c2 * (5. / 6.) + U - P / (3. * U);
    const double w      = std::sqrt(std::max(1e-16, c2 + 2. * s));
    const double lambda = 0.5 * c1 / w;
    const double root2  = std::max(1e-16, w * w - 4. * (c2 + s + lambda));
    const double a      = 0.5 * (w + std::sqrt(root2));

    // Calculation of the normal damping constant.
    const double normal_damping_constant =
      -1.8257 * beta * sqrt(model_parameter_sn * effective_mass);

    // Calculation of the tangential spring constant
    const double tangential_spring_constant =
      8.0 * radius_times_overlap_sqrt * shear_modulus;
    // Calculation of the tangential damping constant
    const double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of the normal force coefficient (F_n_JKR) # Eq 20
    const double normal_force_coefficient =
      4. * Utilities::fixed_power<3>(a) / (3. * effective_radius) *
        youngs_modulus -
      std::sqrt(8. * M_PI * surface_energy * youngs_modulus *
                Utilities::fixed_power<3>(a));

    // Calculation of the final normal force vector
    normal_force = (normal_force_coefficient +
                    normal_damping_constant * normal_relative_velocity_value) *
                   normal_unit_vector;

    tangential_force =
      tangential_spring_constant * contact_info.tangential_displacement +
      tangential_damping_constant * tangential_relative_velocity;

    // JKR theory says that the coulomb threshold must be modified with the
    // pull-out force.
    const double two_pull_off_force =
      3. * M_PI * surface_energy * effective_radius;

    const double modified_coulomb_threshold =
      (normal_force_coefficient + two_pull_off_force) * friction_coeff;

    if (tangential_force.norm() > modified_coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement and tangential
        // force are limited to Coulomb's criterion
        tangential_force =
          modified_coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force * diameter_one * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque * diameter_two / diameter_one;

    // We need to compute the normal spring constant in case if we use the EPSD
    // rolling resistance model.
    double normal_spring_constant = 0.66665 * model_parameter_sn;

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      effective_radius,
      particle_one_properties,
      particle_two_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      normal_unit_vector,
      contact_info.rolling_resistance_spring_torque);
  }
  /**
   * @brief Calculate the particle-particle contact and cohesive forces and
   * contact torque based on the updated values in contact_info. It uses the DMT
   * cohesive force model.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative contact velocity.
   * @param[in] normal_unit_vector Contact normal unit vector.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_one_properties Properties of particle one in contact.
   * @param[in] particle_two_properties Properties of particle two in contact.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] particle_one_tangential_torque Contact tangential torque on
   * particle one.
   * @param[out] particle_two_tangential_torque Contact tangential torque on
   * particle two.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_DMT_contact(particle_particle_contact_info<dim> &contact_info,
                        const Tensor<1, 3> &tangential_relative_velocity,
                        const double        normal_relative_velocity_value,
                        const Tensor<1, 3> &normal_unit_vector,
                        const double        normal_overlap,
                        const double        dt,
                        const ArrayView<const double> &particle_one_properties,
                        const ArrayView<const double> &particle_two_properties,
                        Tensor<1, 3>                  &normal_force,
                        Tensor<1, 3>                  &tangential_force,
                        Tensor<1, 3> &particle_one_tangential_torque,
                        Tensor<1, 3> &particle_two_tangential_torque,
                        Tensor<1, 3> &rolling_resistance_torque)
  {
    constexpr double M_2PI = 2. * M_PI;

    // Calculation of effective radius and mass
    auto [effective_radius, effective_mass] =
      find_effective_radius_and_mass(particle_one_properties,
                                     particle_two_properties);

    const unsigned int particle_one_type =
      static_cast<unsigned int>(particle_one_properties[PropertiesIndex::type]);
    const unsigned int particle_two_type =
      static_cast<unsigned int>(particle_two_properties[PropertiesIndex::type]);
    const unsigned int pair_index =
      vec_particle_type_index(particle_one_type, particle_two_type);

    const double surface_energy = this->effective_surface_energy[pair_index];
    const double hamaker_constant =
      this->effective_hamaker_constant[pair_index];

    const double F_po = M_2PI * effective_radius * surface_energy;

    const double delta_0 =
      -std::sqrt(hamaker_constant * effective_radius / (6. * F_po));

    // Cohesive force. This will need to be added to the
    // first vector inside the tuple.
    double cohesive_term;

    // Contact between particle + constant cohesive force.
    if (normal_overlap > 0.)
      {
        cohesive_term = -F_po;

        calculate_hertz_mindlin_limit_overlap_contact(
          contact_info,
          tangential_relative_velocity,
          normal_relative_velocity_value,
          normal_unit_vector,
          normal_overlap,
          dt,
          particle_one_properties,
          particle_two_properties,
          normal_force,
          tangential_force,
          particle_one_tangential_torque,
          particle_two_tangential_torque,
          rolling_resistance_torque);
      }
    // No contact, but still in the constant zone for the cohesive force.
    else if (normal_overlap > delta_0)
      {
        cohesive_term = -F_po;
        contact_info.tangential_displacement.clear();
        contact_info.rolling_resistance_spring_torque.clear();
      }
    // No contact. Particle are far from each other. The cohesive force is not
    // constant. It needs to be computed.
    else
      {
        cohesive_term = -hamaker_constant * effective_radius /
                        (6. * Utilities::fixed_power<2>(normal_overlap));
        contact_info.tangential_displacement.clear();
        contact_info.rolling_resistance_spring_torque.clear();
      }
    normal_force += cohesive_term * normal_unit_vector;
  }

  /**
   * @brief Return the index for accessing the properties vectors for a given
   * combinations of particle types.
   *
   * @param[in] i First particle type index.
   * @param[in] j Second particle type index.
   * @return index associated with the combinations of particle types.
   */
  inline unsigned int
  vec_particle_type_index(const unsigned int i, const unsigned int j)
  {
    return i * n_particle_types + j;
  }

  /**
   * @brief Set every containers needed to carry the particle-particle force
   * calculation.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm file.
   */
  void
  set_effective_properties(const DEMSolverParameters<dim> &dem_parameters)
  {
    auto properties = dem_parameters.lagrangian_physical_properties;

    n_particle_types = properties.particle_type_number;
    effective_youngs_modulus.resize(n_particle_types * n_particle_types);
    effective_shear_modulus.resize(n_particle_types * n_particle_types);
    effective_coefficient_of_restitution.resize(n_particle_types *
                                                n_particle_types);
    effective_coefficient_of_friction.resize(n_particle_types *
                                             n_particle_types);
    effective_rolling_viscous_damping_coefficient.resize(n_particle_types *
                                                         n_particle_types);
    effective_coefficient_of_rolling_friction.resize(n_particle_types *
                                                     n_particle_types);
    model_parameter_beta.resize(n_particle_types * n_particle_types);
    effective_surface_energy.resize(n_particle_types * n_particle_types);
    effective_hamaker_constant.resize(n_particle_types * n_particle_types);

    for (unsigned int i = 0; i < n_particle_types; ++i)
      {
        const double youngs_modulus_i =
          properties.youngs_modulus_particle.at(i);
        const double poisson_ratio_i = properties.poisson_ratio_particle.at(i);
        const double restitution_coefficient_i =
          properties.restitution_coefficient_particle.at(i);
        const double friction_coefficient_i =
          properties.friction_coefficient_particle.at(i);
        const double rolling_viscous_damping_coefficient_i =
          properties.rolling_viscous_damping_coefficient_particle.at(i);
        const double rolling_friction_coefficient_i =
          properties.rolling_friction_coefficient_particle.at(i);
        const double surface_energy_i =
          properties.surface_energy_particle.at(i);
        const double hamaker_constant_i =
          properties.hamaker_constant_particle.at(i);

        for (unsigned int j = 0; j < n_particle_types; ++j)
          {
            const unsigned int k = i * n_particle_types + j;

            const double youngs_modulus_j =
              properties.youngs_modulus_particle.at(j);
            const double poisson_ratio_j =
              properties.poisson_ratio_particle.at(j);
            const double restitution_coefficient_j =
              properties.restitution_coefficient_particle.at(j);
            const double friction_coefficient_j =
              properties.friction_coefficient_particle.at(j);
            const double rolling_viscous_damping_coefficient_j =
              properties.rolling_viscous_damping_coefficient_particle.at(j);
            const double rolling_friction_coefficient_j =
              properties.rolling_friction_coefficient_particle.at(j);
            const double surface_energy_j =
              properties.surface_energy_particle.at(j);
            const double hamaker_constant_j =
              properties.hamaker_constant_particle.at(j);

            this->effective_youngs_modulus[k] =
              (youngs_modulus_i * youngs_modulus_j) /
              ((youngs_modulus_j * (1.0 - poisson_ratio_i * poisson_ratio_i)) +
               (youngs_modulus_i * (1.0 - poisson_ratio_j * poisson_ratio_j)) +
               DBL_MIN);

            this->effective_shear_modulus[k] =
              (youngs_modulus_i * youngs_modulus_j) /
              (2.0 * ((youngs_modulus_j * (2.0 - poisson_ratio_i) *
                       (1.0 + poisson_ratio_i)) +
                      (youngs_modulus_i * (2.0 - poisson_ratio_j) *
                       (1.0 + poisson_ratio_j))) +
               DBL_MIN);

            this->effective_coefficient_of_restitution[k] =
              harmonic_mean(restitution_coefficient_i,
                            restitution_coefficient_j);

            this->effective_coefficient_of_friction[k] =
              harmonic_mean(friction_coefficient_i, friction_coefficient_j);

            this->effective_rolling_viscous_damping_coefficient[k] =
              harmonic_mean(rolling_viscous_damping_coefficient_i,
                            rolling_viscous_damping_coefficient_j);

            this->effective_coefficient_of_rolling_friction[k] =
              harmonic_mean(rolling_friction_coefficient_i,
                            rolling_friction_coefficient_j);

            this->effective_surface_energy[k] =
              surface_energy_i + surface_energy_j -
              std::pow(std::sqrt(surface_energy_i) -
                         std::sqrt(surface_energy_j),
                       2);

            this->effective_hamaker_constant[k] =
              0.5 * (hamaker_constant_i + hamaker_constant_j);

            double restitution_coefficient_particle_log =
              std::log(this->effective_coefficient_of_restitution[k]);

            this->model_parameter_beta[k] =
              restitution_coefficient_particle_log /
              sqrt(restitution_coefficient_particle_log *
                     restitution_coefficient_particle_log +
                   9.8696);
          }
      }
  }

  /**
   * @brief Set every containers needed to carry the heat transfer rate
   * calculation.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm
   * file.
   */
  void
  set_multiphysic_properties(const DEMSolverParameters<dim> &dem_parameters)
  {
    auto properties = dem_parameters.lagrangian_physical_properties;

    n_particle_types = properties.particle_type_number;
    effective_real_youngs_modulus.resize(n_particle_types * n_particle_types);
    equivalent_surface_roughness.resize(n_particle_types * n_particle_types);
    equivalent_surface_slope.resize(n_particle_types * n_particle_types);
    effective_microhardness.resize(n_particle_types * n_particle_types);
    thermal_conductivity_particle.resize(n_particle_types);
    gas_parameter_m.resize(n_particle_types * n_particle_types);
    this->thermal_conductivity_gas = properties.thermal_conductivity_gas;

    for (unsigned int i = 0; i < n_particle_types; ++i)
      {
        const double real_youngs_modulus_i =
          properties.real_youngs_modulus_particle.at(i);
        const double poisson_ratio_i = properties.poisson_ratio_particle.at(i);
        const double surface_roughness_i =
          properties.surface_roughness_particle.at(i);
        const double surface_slope_i = properties.surface_slope_particle.at(i);
        const double microhardness_i = properties.microhardness_particle.at(i);
        const double thermal_accommodation_i =
          properties.thermal_accommodation_particle.at(i);


        this->thermal_conductivity_particle[i] =
          properties.thermal_conductivity_particle.at(i);

        for (unsigned int j = 0; j < n_particle_types; ++j)
          {
            const unsigned int k = i * n_particle_types + j;

            const double real_youngs_modulus_j =
              properties.real_youngs_modulus_particle.at(j);
            const double poisson_ratio_j =
              properties.poisson_ratio_particle.at(j);
            const double surface_roughness_j =
              properties.surface_roughness_particle.at(j);
            const double surface_slope_j =
              properties.surface_slope_particle.at(j);
            const double microhardness_j =
              properties.microhardness_particle.at(j);
            const double thermal_accommodation_j =
              properties.thermal_accommodation_particle.at(j);

            this->effective_real_youngs_modulus[k] =
              (real_youngs_modulus_i * real_youngs_modulus_j) /
              ((real_youngs_modulus_j *
                (1.0 - poisson_ratio_i * poisson_ratio_i)) +
               (real_youngs_modulus_i *
                (1.0 - poisson_ratio_j * poisson_ratio_j)) +
               DBL_MIN);
            this->equivalent_surface_roughness[k] =
              sqrt(surface_roughness_i * surface_roughness_i +
                   surface_roughness_j * surface_roughness_j);
            this->equivalent_surface_slope[k] =
              sqrt(surface_slope_i * surface_slope_i +
                   surface_slope_j * surface_slope_j);
            this->effective_microhardness[k] =
              harmonic_mean(microhardness_i, microhardness_j);
            this->gas_parameter_m[k] =
              ((2. - thermal_accommodation_i) / thermal_accommodation_i +
               (2. - thermal_accommodation_j) / thermal_accommodation_j) *
              (2. * properties.specific_heats_ratio_gas) /
              (1. + properties.specific_heats_ratio_gas) *
              properties.molecular_mean_free_path_gas /
              (properties.dynamic_viscosity_gas * properties.specific_heat_gas /
               properties.thermal_conductivity_gas);
          }
      }
  }

  /**
   * @brief Execute the contact calculation step for the particle-particle
   * contact according to the contact type
   *
   * @param[in] adjacent_particles_list Container of the adjacent particles of a
   * particles
   * @param[in] dt DEM time step.
   * @param[out] contact_outcome Interaction outcomes.
   */
  template <ContactType contact_type>
  inline void
  execute_contact_calculation(
    typename DEM::dem_data_structures<dim>::particle_contact_info
                                                 &adjacent_particles_list,
    const double                                  dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
  {
    // No contact calculation if no adjacent particles
    if (adjacent_particles_list.empty())
      return;

    // Define local variables which will be used within the contact calculation
    Tensor<1, 3> normal_unit_vector;
    Tensor<1, 3> normal_force;
    Tensor<1, 3> tangential_force;
    Tensor<1, 3> particle_one_tangential_torque;
    Tensor<1, 3> particle_two_tangential_torque;
    Tensor<1, 3> rolling_resistance_torque;
    double       normal_relative_velocity_value;
    Tensor<1, 3> tangential_relative_velocity;

    // Gather information about particle 1 and set it up.
    auto first_contact_info      = adjacent_particles_list.begin();
    auto particle_one            = first_contact_info->second.particle_one;
    auto particle_one_properties = particle_one->get_properties();

    types::particle_index particle_one_id = particle_one->get_local_index();
    Tensor<1, 3> &particle_one_torque = contact_outcome.torque[particle_one_id];
    Tensor<1, 3> &particle_one_force  = contact_outcome.force[particle_one_id];

    // Fix particle one location for 2d and 3d
    Point<3> particle_one_location = get_location(particle_one);

    for (auto &&contact_info :
         adjacent_particles_list | boost::adaptors::map_values)
      {
        // Getting information (location and properties) of particle 2 in
        // contact with particle 1
        auto particle_two            = contact_info.particle_two;
        auto particle_two_properties = particle_two->get_properties();

        // Get particle 2 location
        Point<3> particle_two_location;
        if constexpr (contact_type == ContactType::local_particle_particle ||
                      contact_type == ContactType::ghost_particle_particle)
          {
            particle_two_location = get_location(particle_two);
          }

        // Get particle 2 location in periodic boundary
        if constexpr (contact_type ==
                        ContactType::local_periodic_particle_particle ||
                      contact_type ==
                        ContactType::ghost_periodic_particle_particle ||
                      contact_type ==
                        ContactType::ghost_local_periodic_particle_particle)
          {
            particle_two_location = get_periodic_location(particle_two);
          }

        // Calculation of normal overlap
        double normal_overlap =
          0.5 * (particle_one_properties[PropertiesIndex::dp] +
                 particle_two_properties[PropertiesIndex::dp]) -
          particle_one_location.distance(particle_two_location);

        // Get the threshold distance for contact force, this is useful for non-
        // contact cohesive force models such as the DMT.
        const double force_calculation_threshold_distance =
          get_force_calculation_threshold_distance();

        if (normal_overlap > force_calculation_threshold_distance)
          {
            // Update of contact information and calculation of contact force
            // are the same for all local-local and local-ghost contact.
            // However, they are based on particle two for ghost-local periodic
            // contact, where particle one is the ghost particle, so the order
            // of particles given to the function are reversed.
            if constexpr (contact_type ==
                            ContactType::local_particle_particle ||
                          contact_type ==
                            ContactType::ghost_particle_particle ||
                          contact_type ==
                            ContactType::local_periodic_particle_particle ||
                          contact_type ==
                            ContactType::ghost_periodic_particle_particle)
              {
                // Update all the information
                this->update_contact_information(contact_info,
                                                 tangential_relative_velocity,
                                                 normal_relative_velocity_value,
                                                 normal_unit_vector,
                                                 particle_one_properties,
                                                 particle_two_properties,
                                                 particle_one_location,
                                                 particle_two_location,
                                                 dt);

                // Calculation the contact force
                this->calculate_contact(contact_info,
                                        tangential_relative_velocity,
                                        normal_relative_velocity_value,
                                        normal_unit_vector,
                                        normal_overlap,
                                        dt,
                                        particle_one_properties,
                                        particle_two_properties,
                                        normal_force,
                                        tangential_force,
                                        particle_one_tangential_torque,
                                        particle_two_tangential_torque,
                                        rolling_resistance_torque);
              }

            if constexpr (contact_type ==
                          ContactType::ghost_local_periodic_particle_particle)
              {
                // In ghost_local periodic contacts, particle one and two are
                // swapped when calling the update_contact_information and the
                // calculate_contact. This is because the first particle
                // received by these functions
                // should always be a local particle. In this case, since this
                // is a ghost_local container, particle_two is the local
                // particle and particle_one is a ghost particle.
                //  Update all the information
                this->update_contact_information(contact_info,
                                                 tangential_relative_velocity,
                                                 normal_relative_velocity_value,
                                                 normal_unit_vector,
                                                 particle_two_properties,
                                                 particle_one_properties,
                                                 particle_two_location,
                                                 particle_one_location,
                                                 dt);

                // Calculation the contact force
                this->calculate_contact(contact_info,
                                        tangential_relative_velocity,
                                        normal_relative_velocity_value,
                                        normal_unit_vector,
                                        normal_overlap,
                                        dt,
                                        particle_two_properties,
                                        particle_one_properties,
                                        normal_force,
                                        tangential_force,
                                        particle_two_tangential_torque,
                                        particle_one_tangential_torque,
                                        rolling_resistance_torque);
              }

            // Apply the calculated forces and torques on both particles
            // of the pair for local-local contacts
            if constexpr (contact_type ==
                            ContactType::local_particle_particle ||
                          contact_type ==
                            ContactType::local_periodic_particle_particle)
              {
                types::particle_index particle_two_id =
                  particle_two->get_local_index();

                Tensor<1, 3> &particle_two_torque =
                  contact_outcome.torque[particle_two_id];
                Tensor<1, 3> &particle_two_force =
                  contact_outcome.force[particle_two_id];

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

            // Apply the calculated forces and torques only on the local
            // particle of the pair for local-ghost contacts
            if constexpr (contact_type ==
                            ContactType::ghost_periodic_particle_particle ||
                          contact_type == ContactType::ghost_particle_particle)
              {
                this->apply_force_and_torque_on_single_local_particle(
                  normal_force,
                  tangential_force,
                  particle_one_tangential_torque,
                  rolling_resistance_torque,
                  particle_one_torque,
                  particle_one_force);
              }

            // Apply the calculated forces and torques only on the local
            // particle of the pair for ghost-local contacts
            if constexpr (contact_type ==
                          ContactType::ghost_local_periodic_particle_particle)
              {
                types::particle_index particle_two_id =
                  particle_two->get_local_index();

                Tensor<1, 3> &particle_two_torque =
                  contact_outcome.torque[particle_two_id];
                Tensor<1, 3> &particle_two_force =
                  contact_outcome.force[particle_two_id];

                this->apply_force_and_torque_on_single_local_particle(
                  normal_force,
                  tangential_force,
                  particle_two_tangential_torque,
                  rolling_resistance_torque,
                  particle_two_torque,
                  particle_two_force);
              }
          }
        else
          {
            // If the adjacent pair is not in contact anymore, only the
            // tangential displacement is set to zero
            contact_info.tangential_displacement.clear();
            contact_info.rolling_resistance_spring_torque.clear();
          }

        if constexpr (std::is_same_v<PropertiesIndex,
                                     DEM::DEMMPProperties::PropertiesIndex>)
          {
            if (normal_overlap > 0)
              {
                const unsigned int particle_one_type =
                  static_cast<unsigned int>(
                    particle_one_properties[PropertiesIndex::type]);
                const unsigned int particle_two_type =
                  static_cast<unsigned int>(
                    particle_two_properties[PropertiesIndex::type]);
                const unsigned int pair_index =
                  vec_particle_type_index(particle_one_type, particle_two_type);
                const double temperature_one =
                  particle_one_properties[PropertiesIndex::T];
                const double temperature_two =
                  particle_two_properties[PropertiesIndex::T];
                types::particle_index particle_two_id =
                  particle_two->get_local_index();
                double &particle_one_heat_transfer_rate =
                  contact_outcome.heat_transfer_rate[particle_one_id];
                double &particle_two_heat_transfer_rate =
                  contact_outcome.heat_transfer_rate[particle_two_id];

                double thermal_conductance;
                calculate_contact_thermal_conductance<contact_type>(
                  0.5 * particle_one_properties[PropertiesIndex::dp],
                  0.5 * particle_two_properties[PropertiesIndex::dp],
                  this->effective_youngs_modulus[pair_index],
                  this->effective_real_youngs_modulus[pair_index],
                  this->equivalent_surface_roughness[pair_index],
                  this->equivalent_surface_slope[pair_index],
                  this->effective_microhardness[pair_index],
                  this->thermal_conductivity_particle[particle_one_type],
                  this->thermal_conductivity_particle[particle_two_type],
                  this->thermal_conductivity_gas,
                  this->gas_parameter_m[pair_index],
                  normal_overlap,
                  normal_force.norm(),
                  thermal_conductance);

                // Apply the heat transfer to both particles
                // of the pair for local-local contacts
                if constexpr (contact_type ==
                                ContactType::local_particle_particle ||
                              contact_type ==
                                ContactType::local_periodic_particle_particle)
                  {
                    apply_heat_transfer_on_local_particles(
                      temperature_one,
                      temperature_two,
                      thermal_conductance,
                      particle_one_heat_transfer_rate,
                      particle_two_heat_transfer_rate);
                  }

                // Apply the heat transfer only to the local
                // particle of the pair for local-ghost contacts
                if constexpr (contact_type ==
                                ContactType::ghost_periodic_particle_particle ||
                              contact_type ==
                                ContactType::ghost_particle_particle)
                  {
                    apply_heat_transfer_on_single_local_particle(
                      temperature_one,
                      temperature_two,
                      thermal_conductance,
                      particle_one_heat_transfer_rate);
                  }

                // Apply the heat transfer only to the local
                // particle of the pair for ghost-local contacts
                if constexpr (contact_type ==
                              ContactType::
                                ghost_local_periodic_particle_particle)
                  {
                    apply_heat_transfer_on_single_local_particle(
                      temperature_two,
                      temperature_one,
                      thermal_conductance,
                      particle_two_heat_transfer_rate);
                  }
              }
          }
      }
  }

  // Members of the class
  // Contact model parameter. It is calculated in the constructor for
  // different combinations of particle types. For different combinations, a
  // map of map is used to store this variable
  unsigned int        n_particle_types;
  std::vector<double> effective_youngs_modulus;
  std::vector<double> effective_real_youngs_modulus;
  std::vector<double> effective_shear_modulus;
  std::vector<double> effective_coefficient_of_restitution;
  std::vector<double> effective_coefficient_of_friction;
  std::vector<double> effective_rolling_viscous_damping_coefficient;
  std::vector<double> effective_coefficient_of_rolling_friction;
  std::vector<double> effective_surface_energy;
  std::vector<double> effective_hamaker_constant;
  std::vector<double> model_parameter_beta;
  const double        dmt_cut_off_threshold;
  const double        f_coefficient_epsd;
  std::vector<double> equivalent_surface_roughness;
  std::vector<double> equivalent_surface_slope;
  std::vector<double> effective_microhardness;
  std::vector<double> thermal_conductivity_particle;
  std::vector<double> gas_parameter_m;
  double              thermal_conductivity_gas;
};

#endif
