// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_contact_force_h
#define lethe_particle_wall_contact_force_h

#include <core/auxiliary_math_functions.h>

#include <dem/contact_info.h>
#include <dem/contact_type.h>
#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/rolling_resistance_torque_models.h>

#include <boost/range/adaptor/map.hpp>

#include <vector>

using namespace dealii;

/**
 * @brief Base class for the particle-wall contact force models
 * This class does not implement any of the models, but ensures that
 * an interface without template specialization is available. All of the
 * actual implementation of the models are carried out in the
 * ParticleWallContactForce class which is templated by the contact model
 * type.
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class ParticleWallContactForceBase
{
public:
  /**
   * @brief Calculate the contact outcomes using the contact pair information
   * obtained in the fine search and physical properties of particles.
   *
   * @param particle_wall_pairs_in_contact Required information for the calculation of the
   * particle-wall contact force.
   * @param dt DEM time step.
   * @param contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                &particle_wall_pairs_in_contact,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) = 0;

  /**
   * @brief Calculate the contact outcomes using the contact pair information
   * obtained in the fine search and physical properties.
   *
   * @param particle_floating_mesh_in_contact A container that stores the information
   * of particle-floating mesh contact
   * @param dt DEM time step.
   * @param solids Floating solids
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                &particle_floating_mesh_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) = 0;
};

/**
 * @brief Execute calculation of particle-wall contact forces.
 *
 * @tparam dim The dimension of the problem.
 * @tparam contact_model The particle-wall contact force model.
 * @tparam rolling_friction_model The rolling resistance model.
 */
template <
  int dim,
  typename PropertiesIndex,
  Parameters::Lagrangian::ParticleWallContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod       rolling_friction_model>
class ParticleWallContactForce
  : public ParticleWallContactForceBase<dim, PropertiesIndex>
{
public:
  ParticleWallContactForce(const DEMSolverParameters<dim> &dem_parameters);

  virtual ~ParticleWallContactForce()
  {}

  /**
   * @brief Calculate the contact outcomes using the contact pair information
   * obtained in the fine search and physical properties.
   *
   * @param particle_wall_pairs_in_contact Required information for the calculation of the
   * particle-wall contact force.
   * @param dt DEM time step.
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                &particle_wall_pairs_in_contact,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) override;

  /**
   * @brief Calculate the contact outcomes using the contact pair information
   * obtained in the fine search and physical properties.
   *
   * @param particle_floating_mesh_in_contact A container that stores the information
   * of particle-floating mesh contact
   * @param dt DEM time step.
   * @param solids Floating solids
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                &particle_floating_mesh_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) override;

}

protected :
  /**
   * @brief Update the contact pair information for both non-linear and
   * linear contact force calculations
   *
   * @param contact_pair_information Contact information of a particle-wall pair
   * in neighborhood
   * @param particle_properties Properties of particle in contact
   * @param dt DEM time step
   */
  inline void
  update_contact_information(particle_wall_contact_info<dim> &contact_info,
                             const Point<3>                  &particle_position,
                             const ArrayView<const double> &particle_properties,
                             const double                   dt)
{
  // i is the particle, j is the wall.
  // we need to put a minus sign infront of the normal_vector to respect the
  // convention (i -> j)
  auto               normal_vector = -contact_info.normal_vector;
  const unsigned int boundary_id   = contact_info.boundary_id;

  // Using velocity and angular velocity of particle as
  // local vectors
  Tensor<1, 3> particle_velocity;
  particle_velocity[0] = particle_properties[PropertiesIndex::v_x];
  particle_velocity[1] = particle_properties[PropertiesIndex::v_y];
  particle_velocity[2] = particle_properties[PropertiesIndex::v_z];

  Tensor<1, 3> particle_angular_velocity;
  particle_angular_velocity[0] = particle_properties[PropertiesIndex::omega_x];
  particle_angular_velocity[1] = particle_properties[PropertiesIndex::omega_y];
  particle_angular_velocity[2] = particle_properties[PropertiesIndex::omega_z];

  // Calculate approximation of the contact point using the normal vector
  Point<3> contact_point =
    particle_position +
    0.5 * particle_properties[PropertiesIndex::dp] * normal_vector;

  // Get vector pointing from the contact point to the origin of the rotation
  // axis
  Tensor<1, 3> vector_to_rotating_axis =
    contact_point - this->point_on_rotation_vector[boundary_id];

  // Remove the rotating axis component of that vector
  vector_to_rotating_axis =
    vector_to_rotating_axis -
    (vector_to_rotating_axis * this->boundary_rotational_vector[boundary_id]) *
      this->boundary_rotational_vector[boundary_id];

  // Tensor<1,3> vector_to_rotation_axis = this->boundary_rotational_speed_map

  // Defining relative contact velocity using the convention
  // v_ij = v_j - v_i
  Tensor<1, 3> contact_relative_velocity =
    this->boundary_translational_velocity_map[boundary_id] - particle_velocity +
    cross_product_3d((-0.5 * particle_properties[PropertiesIndex::dp] *
                      particle_angular_velocity),
                     normal_vector) +
    cross_product_3d(this->boundary_rotational_speed_map[boundary_id] *
                       this->boundary_rotational_vector[boundary_id],
                     vector_to_rotating_axis);

  // Calculation of normal relative velocity
  double normal_relative_velocity_value =
    contact_relative_velocity * normal_vector;
  Tensor<1, 3> normal_relative_velocity =
    normal_relative_velocity_value * normal_vector;

  // Calculation of tangential relative velocity
  Tensor<1, 3> tangential_relative_velocity =
    contact_relative_velocity - normal_relative_velocity;

  // Calculation of new tangential_displacement, since this value is
  // history-dependent it needs the value at previous time-step
  // This variable is the main reason that we have iteration over
  // two different vectors (pairs_in_contact and
  // contact_pair_candidates): tangential_displacement of the particles
  // which were already in contact (pairs_in_contact) needs to be
  // modified using its history, while the tangential_displacements of
  // new particles are equal to zero
  Tensor<1, 3> modified_tangential_displacement =
    contact_info.tangential_displacement + tangential_relative_velocity * dt;

  // Updating the contact_info container based on the new calculated values
  contact_info.normal_relative_velocity     = normal_relative_velocity_value;
  contact_info.tangential_displacement      = modified_tangential_displacement;
  contact_info.tangential_relative_velocity = tangential_relative_velocity;
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
 * @brief Calculate the particle-wall contact force and torque
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
calculate_contact(particle_wall_contact_info<dim> &contact_info,
                  const double                   dt,
                  const ArrayView<const double> &particle_properties,
                  Tensor<1, 3>                  &normal_force,
                  Tensor<1, 3>                  &tangential_force,
                  Tensor<1, 3>                  &tangential_torque,
                  Tensor<1, 3>                  &rolling_resistance_torque)
{
  using namespace Parameters::Lagrangian;

  if constexpr (contact_model == ParticleWallContactForceModel::linear)
    {
      calculate_linear_contact(contact_info,
                               dt,
                               particle_properties,
                               normal_force,
                               tangential_force,
                               tangential_torque,
                               rolling_resistance_torque);
    }

  if constexpr (contact_model == ParticleWallContactForceModel::nonlinear)
    {
      this->calculate_hertz_contact(contact_info,
        dt,
        particle_properties,
        normal_force,
        tangential_force,
        tangential_torque,
        rolling_resistance_torque);
    }

  if constexpr (contact_model ==
                ParticleWallContactForceModel::JRK)
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

  if constexpr (contact_model ==
                ParticleWallContactForceModel::hertz_mindlin_limit_overlap)
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


  if constexpr (contact_model == ParticleWallContactForceModel::hertz_JKR)
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

  if constexpr (contact_model == ParticleWallContactForceModel::DMT)
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
 * @brief Calculate the minimum overlap at which particle-wall forces are
 * computed.
 *
 * @return minimum overlap for the force calculation.
 */
inline double
get_force_calculation_threshold_distance()
{
  if constexpr (contact_model ==
                Parameters::Lagrangian::ParticleWallContactForceModel::DMT)
    {
      // We are looking for the maximum hamaker constant and minimum surface
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
  double       effective_radius =
    (diameter_one * diameter_two) / (2 * (diameter_one + diameter_two));

  // Calculate the effective mass
  const double mass_one       = particle_one_properties[PropertiesIndex::mass];
  const double mass_two       = particle_two_properties[PropertiesIndex::mass];
  double       effective_mass = (mass_one * mass_two) / (mass_one + mass_two);

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
  [[maybe_unused]] const double                   rolling_viscous_damping_coeff,
  [[maybe_unused]] const double                   dt,
  [[maybe_unused]] const double                   normal_spring_constant,
  [[maybe_unused]] const double                   normal_force_norm,
  [[maybe_unused]] const Tensor<1, 3>            &normal_unit_vector,
  Tensor<1, 3> &cumulative_rolling_resistance_spring_torque)
{
  using namespace Parameters::Lagrangian;

  if constexpr (rolling_friction_model ==
                RollingResistanceMethod::no_resistance)
    {
      return no_rolling_resistance_torque();
    }

  if constexpr (rolling_friction_model ==
                RollingResistanceMethod::constant_resistance)
    {
      return constant_rolling_resistance_torque<PropertiesIndex>(
        effective_radius,
        particle_one_properties,
        particle_two_properties,
        rolling_friction_coeff,
        normal_force_norm);
    }

  if constexpr (rolling_friction_model == viscous_resistance)
    {
      return viscous_rolling_resistance_torque<PropertiesIndex>(
        effective_radius,
        particle_one_properties,
        particle_two_properties,
        rolling_friction_coeff,
        normal_force_norm,
        normal_unit_vector);
    }
  if constexpr (rolling_friction_model == epsd_resistance)
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
 * @brief Calculate the particle-wall linear contact force and torque
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
calculate_linear_contact(particle_particle_contact_info<dim> &contact_info,
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
  // Calculation of effective radius and mass
  auto [effective_radius, effective_mass] =
    find_effective_radius_and_mass(particle_one_properties,
                                   particle_two_properties);

  // Get the reference of the effective properties according to the particle
  // types in vectors for easy-to-read equations.
  const unsigned int particle_one_type =
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];
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

  // Calculate the normal spring constant using the following formula:
  // kn = 16/15 * sqrt(Re) * Ye * (15/16 * (me * vc^2 / (sqrt(R) * Ye))^0.2
  double normal_spring_constant =
    1.0667 * sqrt(effective_radius) * youngs_modulus *
    pow((0.9375 * effective_mass * characteristic_velocity *
         characteristic_velocity / (sqrt(effective_radius) * youngs_modulus)),
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
    cross_product_3d(normal_unit_vector, tangential_force * diameter_one * 0.5);
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
 * @brief Calculate the particle-wall non-linear contact force and torque
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
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];
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

  // Calculation of normal and tangential spring and dashpot constants
  // using particle properties

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

  // Calculation of normal force
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

  // Calculation of torque caused by tangential force (tangential_torque)
  particle_one_tangential_torque =
    cross_product_3d(normal_unit_vector, tangential_force * diameter_one * 0.5);
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
 * @brief Calculate the particle-wall non-linear contact force and torque
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
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];
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
    cross_product_3d(normal_unit_vector, tangential_force * diameter_one * 0.5);
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
 * @brief Calculate the particle-wall non-linear contact force and torque
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
calculate_hertz_contact(particle_particle_contact_info<dim> &contact_info,
                        const Tensor<1, 3> & /*tangential_relative_velocity*/,
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
  // Calculation of effective radius and mass
  auto [effective_radius, effective_mass] =
    find_effective_radius_and_mass(particle_one_properties,
                                   particle_two_properties);

  const unsigned int particle_one_type =
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];
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
    cross_product_3d(normal_unit_vector, tangential_force * diameter_one * 0.5);
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
 * @brief Calculate the particle-wall contact and cohesive forces and
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
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];
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
    cross_product_3d(normal_unit_vector, tangential_force * diameter_one * 0.5);
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
 * @brief Calculate the particle-wall contact and cohesive forces and
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
    particle_one_properties[PropertiesIndex::type];
  const unsigned int particle_two_type =
    particle_two_properties[PropertiesIndex::type];
  const unsigned int pair_index =
    vec_particle_type_index(particle_one_type, particle_two_type);

  const double surface_energy   = this->effective_surface_energy[pair_index];
  const double hamaker_constant = this->effective_hamaker_constant[pair_index];

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
 * @brief Set every containers needed to carry the particle-wall force
 * calculation.
 *
 * @param[in] dem_parameters DEM parameters declared in the .prm file.
 */
void
set_effective_wall_properties(const DEMSolverParameters<dim> &dem_parameters)
{
  auto properties = dem_parameters.lagrangian_physical_properties;

  // Wall properties
  const double wall_youngs_modulus = properties.youngs_modulus_wall;
  const double wall_poisson_ratio  = properties.poisson_ratio_wall;
  const double wall_restitution_coefficient =
    properties.restitution_coefficient_wall;
  const double wall_friction_coefficient = properties.friction_coefficient_wall;
  const double wall_rolling_friction_coefficient =
    properties.rolling_friction_wall;
  const double wall_rolling_viscous_damping =
    properties.rolling_viscous_damping_wall;
  const double wall_surface_energy = properties.surface_energy_wall;

  for (unsigned int i = 0; i < properties.particle_type_number; ++i)
    {
      // Particle properties
      const double particle_youngs_modulus =
        properties.youngs_modulus_particle.at(i);
      const double particle_poisson_ratio =
        properties.poisson_ratio_particle.at(i);
      const double particle_restitution_coefficient =
        properties.restitution_coefficient_particle.at(i);
      const double particle_friction_coefficient =
        properties.friction_coefficient_particle.at(i);
      const double particle_rolling_friction_coefficient =
        properties.rolling_friction_coefficient_particle.at(i);
      const double particle_rolling_viscous_damping_coefficient =
        properties.rolling_viscous_damping_coefficient_particle.at(i);
      const double particle_surface_energy =
        properties.surface_energy_particle.at(i);

      // Effective particle-wall properties.
      this->effective_youngs_modulus[i] =
        (particle_youngs_modulus * wall_youngs_modulus) /
        (wall_youngs_modulus *
           (1 - particle_poisson_ratio * particle_poisson_ratio) +
         particle_youngs_modulus *
           (1 - wall_poisson_ratio * wall_poisson_ratio) +
         DBL_MIN);

      this->effective_shear_modulus[i] =
        (particle_youngs_modulus * wall_youngs_modulus) /
        ((2 * wall_youngs_modulus * (2 - particle_poisson_ratio) *
          (1 + particle_poisson_ratio)) +
         (2 * particle_youngs_modulus * (2 - wall_poisson_ratio) *
          (1 + wall_poisson_ratio)) +
         DBL_MIN);

      this->effective_coefficient_of_restitution[i] =
        2 * particle_restitution_coefficient * wall_restitution_coefficient /
        (particle_restitution_coefficient + wall_restitution_coefficient +
         DBL_MIN);

      this->effective_coefficient_of_friction[i] =
        2 * particle_friction_coefficient * wall_friction_coefficient /
        (particle_friction_coefficient + wall_friction_coefficient + DBL_MIN);

      this->effective_coefficient_of_rolling_friction[i] =
        2 * particle_rolling_friction_coefficient *
        wall_rolling_friction_coefficient /
        (particle_rolling_friction_coefficient +
         wall_rolling_friction_coefficient + DBL_MIN);

      this->effective_coefficient_of_rolling_viscous_damping[i] =
        2 * particle_rolling_viscous_damping_coefficient *
        wall_rolling_viscous_damping /
        (particle_rolling_viscous_damping_coefficient +
         wall_rolling_viscous_damping + DBL_MIN);

      this->effective_surface_energy[i] =
        particle_surface_energy + wall_surface_energy -
        std::pow(std::sqrt(particle_surface_energy) -
                   std::sqrt(wall_surface_energy),
                 2);

      const double log_coeff_restitution =
        log(this->effective_coefficient_of_restitution[i]);
      this->model_parameter_beta[i] =
        log_coeff_restitution /
        sqrt((log_coeff_restitution * log_coeff_restitution) + 9.8696);
    }

  /**
   * @brief Execute the contact calculation step for the particle-wall
   * contact according to the contact type
   *
   * @param[in] adjacent_particles_list Container of the adjacent particles of
   * a particles
   * @param[in] dt DEM time step.
   * @param[out] contact_outcome Interaction outcomes.
   */
  template <ContactType contact_type>
  inline void execute_contact_calculation(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact &
      particle_wall_pairs_in_contact,
    const double                                  dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
  {
    ParticleWallContactForce<dim, PropertiesIndex>::force_on_walls =
      ParticleWallContactForce<dim, PropertiesIndex>::initialize();
    ParticleWallContactForce<dim, PropertiesIndex>::torque_on_walls =
      ParticleWallContactForce<dim, PropertiesIndex>::initialize();

    // Looping over particle_wall_pairs_in_contact, which means looping over
    // all the active particles with iterator
    // particle_wall_pairs_in_contact_iterator
    for (auto &&pairs_in_contact_content :
         particle_wall_pairs_in_contact | boost::adaptors::map_values)
      {
        // Now an iterator (particle_wall_contact_info_iterator) on each
        // element of the particle_wall_pairs_in_contact vector is defined.
        // This iterator iterates over a map which contains the required
        // information for calculation of the contact force for each particle
        for (auto &&contact_info :
             pairs_in_contact_content | boost::adaptors::map_values)
          {
            // Defining total force of the contact, properties of particle as
            // local parameters
            auto particle            = contact_info.particle;
            auto particle_properties = particle->get_properties();

            Tensor<1, 3> normal_vector     = contact_info.normal_vector;
            auto         point_on_boundary = contact_info.point_on_boundary;

            // Fix particle one location for 2d and 3d
            Point<3> particle_location_3d = get_location(particle);

            // A vector (point_to_particle_vector) is defined which connects
            // the center of particle to the point_on_boundary. This vector
            // will then be projected on the normal vector of the boundary to
            // obtain the particle-wall distance
            Tensor<1, 3> point_to_particle_vector =
              particle_location_3d - point_on_boundary;

            // Finding the projected vector on the normal vector of the
            // boundary. Here we have used the private function
            // find_projection. Using this projected vector, the particle-wall
            // distance is calculated
            Tensor<1, 3> projected_vector =
              this->find_projection(point_to_particle_vector, normal_vector);
            double normal_overlap =
              ((particle_properties[PropertiesIndex::dp]) * 0.5) -
              (projected_vector.norm());

            // Get the threshold distance for contact force, this is useful
            // for non- contact cohesive force models such as the DMT.
            const double force_calculation_threshold_distance =
              get_force_calculation_threshold_distance();

            if (normal_overlap > force_calculation_threshold_distance)
              {
                // Update all the information
                this->update_contact_information(contact_info,
                                                 particle_location_3d,
                                                 particle_properties,
                                                 dt);

                // Calculation the contact force
                this->calculate_contact(contact_info,
                                        dt,
                                        particle_properties,
                                        normal_force,
                                        tangential_force,
                                        tangential_torque,
                                        rolling_resistance_torque);


                // Apply the calculated forces and torques on the particle
                types::particle_index particle_id = particle->get_local_index();

                Tensor<1, 3> &particle_torque =
                  contact_outcome.torque[particle_id];
                Tensor<1, 3> &particle_force =
                  contact_outcome.force[particle_id];

                this->apply_force_and_torque(forces_and_torques,
                                             particle_torque,
                                             particle_force,
                                             point_on_boundary,
                                             contact_info.boundary_id);
              }
            else
              {
                // If the adjacent pair is not in contact anymore, only the
                // tangential displacement is set to zero
                contact_info.tangential_displacement.clear();
                contact_info.rolling_resistance_spring_torque.clear();
              }
          }
      }
  }

  // Members of the class
  std::unordered_map<unsigned int, Tensor<1, 3>>
                                           boundary_translational_velocity_map;
  std::unordered_map<unsigned int, double> boundary_rotational_speed_map;
  std::unordered_map<unsigned int, Tensor<1, 3>> boundary_rotational_vector;
  std::unordered_map<unsigned int, Point<3>>     point_on_rotation_vector;

  unsigned int n_particle_types;

  // Effective properties
  std::vector<double> effective_youngs_modulus;
  std::vector<double> effective_shear_modulus;
  std::vector<double> effective_coefficient_of_restitution;
  std::vector<double> effective_coefficient_of_friction;
  std::vector<double> effective_coefficient_of_rolling_friction;
  std::vector<double> effective_coefficient_of_rolling_viscous_damping;
  std::vector<double> effective_surface_energy;
  std::vector<double> effective_hamaker_constant;
  std::vector<double> model_parameter_beta;

  std::map<unsigned int, Tensor<1, 3>> force_on_walls;
  std::map<unsigned int, Tensor<1, 3>> torque_on_walls;

  bool                            calculate_force_torque_on_boundary;
  Point<3>                        center_mass_container;
  std::vector<types::boundary_id> boundary_index;
  const unsigned int              vertices_per_triangle = 3;
};

#endif
