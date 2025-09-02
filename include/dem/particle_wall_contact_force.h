// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_contact_force_h
#define lethe_particle_wall_contact_force_h

#include <core/auxiliary_math_functions.h>
#include <core/dem_properties.h>
#include <core/lethe_grid_tools.h>
#include <core/serial_solid.h>

#include <dem/contact_info.h>
#include <dem/contact_type.h>
#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_heat_transfer.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/particle_wall_rolling_resistance_torque.h>

#include <boost/math/special_functions.hpp>
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
   * @brief Calculate the contact outcomes for particle-wall contacts
   * using the contact pair information and physical properties.
   *
   * @param[in] particle_wall_pairs_in_contact Required information for the
   * calculation of the particle-wall contact.
   * @param[in] dt DEM time step.
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_wall_contact(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                &particle_wall_pairs_in_contact,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) = 0;

  /**
   * @brief Calculate the contact outcomes for particle-solid objects contacts
   * using the contact pair information and physical properties.
   *
   * @param[in] particle_floating_mesh_in_contact A container that stores the
   * information of particle-floating mesh contact.
   * @param[in] dt DEM time step.
   * @param[in] solids Floating solids.
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_solid_object_contact(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                &particle_floating_mesh_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) = 0;
};

/**
 * @brief Execute calculation of particle-wall contacts.
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
   * @brief Calculate the contact outcomes for particle-wall contacts
   * using the contact pair information and physical properties.
   *
   * @param[in] particle_wall_pairs_in_contact Required information for the
   * calculation of the particle-wall contact.
   * @param[in] dt DEM time step.
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_wall_contact(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                &particle_wall_pairs_in_contact,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) override;

  /**
   * @brief Calculate the contact outcomes for particle-solid objects contacts
   * using the contact pair information and physical properties.
   *
   * @param[in] particle_floating_mesh_in_contact A container that stores the
   * information of particle-floating mesh contact.
   * @param[in] dt DEM time step.
   * @param[in] solids Floating solids.
   * @param[out] contact_outcome Interaction outcomes.
   */
  virtual void
  calculate_particle_solid_object_contact(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                &particle_floating_mesh_in_contact,
    const double dt,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome) override;

protected:
  /**
   * @brief Update the contact pair information for contact calculations.
   *
   * @param[in,out] contact_info Contact information of a particle-wall pair.
   * @param[out] tangential_relative_velocity Tangential relative velocity.
   * @param[out] normal_relative_velocity_value Normal relative contact
   * velocity value.
   * @param[in] particle_location Location of particle in 3d.
   * @param[in] particle_properties Properties of particle.
   * @param[in] dt DEM time step.
   */
  inline void
  update_contact_information(particle_wall_contact_info<dim> &contact_info,
                             Tensor<1, 3>   &tangential_relative_velocity,
                             double         &normal_relative_velocity_value,
                             const Point<3> &particle_location,
                             const ArrayView<const double> &particle_properties,
                             const double                   dt)
  {
    // i is the particle, j is the wall.
    // There is a minus sign in front of the normal_vector to respect the
    // convention i -> j
    Tensor<1, 3>       normal_vector = -contact_info.normal_vector;
    const unsigned int boundary_id   = contact_info.boundary_id;

    Tensor<1, 3> particle_velocity{{particle_properties[PropertiesIndex::v_x],
                                    particle_properties[PropertiesIndex::v_y],
                                    particle_properties[PropertiesIndex::v_z]}};

    Tensor<1, 3> particle_angular_velocity{
      {particle_properties[PropertiesIndex::omega_x],
       particle_properties[PropertiesIndex::omega_y],
       particle_properties[PropertiesIndex::omega_z]}};

    // Calculating approximation of the contact point using the normal vector
    Point<3> contact_point =
      particle_location +
      0.5 * particle_properties[PropertiesIndex::dp] * normal_vector;

    // Getting vector pointing from the contact point to the origin of the
    // rotation axis
    Tensor<1, 3> vector_to_rotating_axis =
      contact_point - this->point_on_rotation_vector[boundary_id];

    // Removing the rotating axis component of that vector
    vector_to_rotating_axis = vector_to_rotating_axis -
                              (vector_to_rotating_axis *
                               this->boundary_rotational_vector[boundary_id]) *
                                this->boundary_rotational_vector[boundary_id];

    // Defining relative contact velocity using the convention
    // v_ij = v_j - v_i
    Tensor<1, 3> contact_relative_velocity =
      this->boundary_translational_velocity_map[boundary_id] -
      particle_velocity +
      cross_product_3d((-0.5 * particle_properties[PropertiesIndex::dp] *
                        particle_angular_velocity),
                       normal_vector) +
      cross_product_3d(this->boundary_rotational_speed_map[boundary_id] *
                         this->boundary_rotational_vector[boundary_id],
                       vector_to_rotating_axis);

    // Calculating normal relative velocity
    normal_relative_velocity_value = contact_relative_velocity * normal_vector;

    // Calculating tangential relative velocity
    tangential_relative_velocity =
      contact_relative_velocity -
      (normal_relative_velocity_value * normal_vector);

    // Calculating new tangential_displacement
    // Since this value is history-dependent, it needs the value at previous
    // time-step. This variable is the main reason that we have iteration over
    // two different vectors (pairs_in_contact and contact_pair_candidates):
    // tangential_displacement of the particles which were already in contact
    // (pairs_in_contact) needs to be modified using its history, while the
    // tangential_displacements of new particles are equal to zero.
    contact_info.tangential_displacement += tangential_relative_velocity * dt;
  }

  /**
   * @brief Update the contact pair information for particle-solid object contacts.
   *
   * @param[in,out] contact_info Contact information of a particle-wall pair.
   * @param[out] tangential_relative_velocity Tangential relative velocity.
   * @param[out] normal_relative_velocity_value Normal relative contact
   * velocity value.
   * @param[in] particle_properties Properties of particle.
   * @param[in] dt DEM time step.
   * @param[in] cut_cell_translational_velocity Translational velocity of the
   * cut cell.
   * @param[in] cut_cell_rotational_velocity Rotational veclocity of the cut
   * cell.
   * @param[in] center_of_rotation_particle_distance Distance between the
   * particle and the center of rotation of the floating mesh.
   */
  inline void
  update_particle_solid_object_contact_information(
    particle_wall_contact_info<dim> &contact_info,
    Tensor<1, 3>                    &tangential_relative_velocity,
    double                          &normal_relative_velocity_value,
    const ArrayView<const double>   &particle_properties,
    const double                     dt,
    const Tensor<1, 3>              &cut_cell_translational_velocity,
    const Tensor<1, 3>              &cut_cell_rotational_velocity,
    const double                     center_of_rotation_particle_distance)
  {
    // i is the particle, j is the wall.
    // we need to put a minus sign infront of the normal_vector to respect the
    // convention (i -> j)
    const Tensor<1, 3> normal_vector = -contact_info.normal_vector;

    // Using velocity and angular velocity of particle as
    // local vectors
    Tensor<1, 3> particle_velocity{{particle_properties[PropertiesIndex::v_x],
                                    particle_properties[PropertiesIndex::v_y],
                                    particle_properties[PropertiesIndex::v_z]}};

    Tensor<1, 3> particle_angular_velocity{
      {particle_properties[PropertiesIndex::omega_x],
       particle_properties[PropertiesIndex::omega_y],
       particle_properties[PropertiesIndex::omega_z]}};

    // Defining relative contact velocity
    // v_ij = v_j - v_i
    Tensor<1, 3> contact_relative_velocity =
      cut_cell_translational_velocity - particle_velocity +
      cross_product_3d((center_of_rotation_particle_distance *
                          cut_cell_rotational_velocity -
                        0.5 * particle_properties[PropertiesIndex::dp] *
                          particle_angular_velocity),
                       normal_vector);

    // Calculation of normal relative velocity
    normal_relative_velocity_value = contact_relative_velocity * normal_vector;

    // Calculation of tangential relative velocity
    tangential_relative_velocity =
      contact_relative_velocity -
      normal_relative_velocity_value * normal_vector;

    // Calculate the new tangential_displacement
    // Since this value is history-dependent, it needs the value at
    // previous time-step.
    contact_info.tangential_displacement += tangential_relative_velocity * dt;
  }

  /**
   * @brief This function is used to find the projection of vector_a on
   * vector_b
   * @param[in] vector_a A vector which is going to be projected on vector_b
   * @param[in] vector_b The projection vector of vector_a
   * @return The projection of vector_a on vector_b
   */
  inline Tensor<1, 3>
  find_projection(const Tensor<1, 3> &vector_a, const Tensor<1, 3> &vector_b)
  {
    Tensor<1, 3> vector_c;
    vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

    return vector_c;
  }

  /**
   * @brief Get the 3d location of the particle.
   *
   * @param[in] particle The particle to get the location from.
   * @return 3d location of particle.
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
   * @param[in,out] contact_info A container with the required
   * information for calculation of the contact force for a particle-wall pair
   * in contact.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative velocity value.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_properties Properties of particle.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] tangential_torque Contact tangential torque.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_contact(particle_wall_contact_info<dim> &contact_info,
                    const Tensor<1, 3> &tangential_relative_velocity,
                    const double        normal_relative_velocity_value,
                    const double        normal_overlap,
                    const double        dt,
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
                                 tangential_relative_velocity,
                                 normal_relative_velocity_value,
                                 normal_overlap,
                                 dt,
                                 particle_properties,
                                 normal_force,
                                 tangential_force,
                                 tangential_torque,
                                 rolling_resistance_torque);
      }

    if constexpr (contact_model == ParticleWallContactForceModel::nonlinear)
      {
        calculate_nonlinear_contact(contact_info,
                                    tangential_relative_velocity,
                                    normal_relative_velocity_value,
                                    normal_overlap,
                                    dt,
                                    particle_properties,
                                    normal_force,
                                    tangential_force,
                                    tangential_torque,
                                    rolling_resistance_torque);
      }

    if constexpr (contact_model == ParticleWallContactForceModel::JKR)
      {
        calculate_JKR_contact(contact_info,
                              tangential_relative_velocity,
                              normal_relative_velocity_value,
                              normal_overlap,
                              dt,
                              particle_properties,
                              normal_force,
                              tangential_force,
                              tangential_torque,
                              rolling_resistance_torque);
      }

    if constexpr (contact_model == ParticleWallContactForceModel::DMT)
      {
        calculate_DMT_contact(contact_info,
                              tangential_relative_velocity,
                              normal_relative_velocity_value,
                              normal_overlap,
                              dt,
                              particle_properties,
                              normal_force,
                              tangential_force,
                              tangential_torque,
                              rolling_resistance_torque);
      }
  }

  /**
   * @brief Calculate the minimum overlap at which particle-wall forces are
   * computed. This function is useful for non-contact cohesive force models
   * such as the DMT.
   *
   * @return minimum overlap for the force calculation.
   */
  inline double
  get_force_calculation_threshold_distance()
  {
    if constexpr (contact_model ==
                  Parameters::Lagrangian::ParticleWallContactForceModel::DMT)
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
        // positive overlap means that there is contact.
        return -std::sqrt(
          max_effective_hamaker_constant /
          (12. * M_PI * min_effective_surface_energy * dmt_cut_off_threshold));
      }
    return 0.;
  }

private:
  /**
   * @brief Apply the calculated force and torque on the
   * particle in contact.
   *
   * @param[in] normal_force Contact normal force.
   * @param[in] tangential_force Contact tangential force.
   * @param[in] tangential_torque Contact torque induced by the tangential
   * force.
   * @param[in] rolling_resistance_torque Contact rolling resistance torque.
   * @param[in,out] particle_torque Torque acting on particle.
   * @param[in,out] particle_force Force acting on particle.
   */
  inline void
  apply_force_and_torque(const Tensor<1, 3> &normal_force,
                         const Tensor<1, 3> &tangential_force,
                         const Tensor<1, 3> &tangential_torque,
                         const Tensor<1, 3> &rolling_resistance_torque,
                         Tensor<1, 3>       &particle_torque,
                         Tensor<1, 3>       &particle_force)
  {
    // Calculating total force
    Tensor<1, 3> total_force = normal_force + tangential_force;

    // Updating force and torque on particle in the particle handler
    // Since the force was calculated on the wall, we use -= operator
    particle_force -= total_force;
    // Since the torque was direcly calculated on the particle, we use +=
    // operator
    particle_torque += tangential_torque + rolling_resistance_torque;
  }

  /**
   * @brief Calculate the rolling resistance torque acting on the particle
   * according to the rolling resistance model.
   *
   * @param[in] particle_radius Radius of the particle.
   * @param[in] particle_properties Properties of particle.
   * @param[in] rolling_friction_coeff Effective rolling friction coefficient.
   * @param[in] rolling_viscous_damping_coeff Effective rolling viscous damping
   * coefficient
   * @param[in] dt DEM time step.
   * @param[in] normal_spring_constant Normal contact stiffness constant.
   * @param[in] normal_force_norm Norm of the normal force.
   * @param[in] normal_unit_vector Normal unit vector.
   * @param[in,out] cumulative_rolling_resistance_spring_torque Cumulative
   * spring rolling resistance torque.
   *
   * @return Rolling resistance torque.
   *
   */
  inline Tensor<1, 3>
  calculate_rolling_resistance_torque(
    [[maybe_unused]] const double                   particle_radius,
    [[maybe_unused]] const ArrayView<const double> &particle_properties,
    [[maybe_unused]] const double                   rolling_friction_coeff,
    [[maybe_unused]] const double        rolling_viscous_damping_coeff,
    [[maybe_unused]] const double        dt,
    [[maybe_unused]] const double        normal_spring_constant,
    [[maybe_unused]] const double        normal_force_norm,
    [[maybe_unused]] const Tensor<1, 3> &normal_unit_vector,
    [[maybe_unused]] Tensor<1, 3> &cumulative_rolling_resistance_spring_torque)
  {
    using namespace Parameters::Lagrangian;

    if constexpr (rolling_friction_model == none)
      {
        return no_rolling_torque();
      }

    if constexpr (rolling_friction_model == constant)
      {
        return constant_rolling_torque<PropertiesIndex>(particle_radius,
                                                        particle_properties,
                                                        rolling_friction_coeff,
                                                        normal_force_norm);
      }

    if constexpr (rolling_friction_model == viscous)
      {
        return viscous_rolling_torque<PropertiesIndex>(particle_radius,
                                                       particle_properties,
                                                       rolling_friction_coeff,
                                                       normal_force_norm,
                                                       normal_unit_vector);
      }
    if constexpr (rolling_friction_model == epsd)
      {
        return epsd_rolling_torque<dim, PropertiesIndex>(
          particle_radius,
          particle_properties,
          rolling_friction_coeff,
          rolling_viscous_damping_coeff,
          this->f_coefficient_epsd,
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
   * information for calculation of the contact force for a particle-wall pair.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative velocity value.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_properties Properties of particle.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] tangential_torque Contact tangential torque.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_linear_contact(particle_wall_contact_info<dim> &contact_info,
                           const Tensor<1, 3> &tangential_relative_velocity,
                           const double        normal_relative_velocity_value,
                           const double        normal_overlap,
                           const double        dt,
                           const ArrayView<const double> &particle_properties,
                           Tensor<1, 3>                  &normal_force,
                           Tensor<1, 3>                  &tangential_force,
                           Tensor<1, 3>                  &tangential_torque,
                           Tensor<1, 3> &rolling_resistance_torque)
  {
    // i is the particle, j is the wall.
    // there is a minus sign in front of the normal_vector to respect the
    // convention (i -> j), the forces are thus calculated on the wall
    const Tensor<1, 3> normal_vector = -contact_info.normal_vector;
    const unsigned int particle_type =
      particle_properties[PropertiesIndex::type];
    const double particle_radius =
      particle_properties[PropertiesIndex::dp] * 0.5;

    const double youngs_modulus = this->effective_youngs_modulus[particle_type];
    const double beta           = this->model_parameter_beta[particle_type];
    const double friction_coeff =
      this->effective_coefficient_of_friction[particle_type];
    const double rolling_viscous_damping_coeff =
      this->effective_coefficient_of_rolling_viscous_damping[particle_type];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[particle_type];

    // Calculation of normal and tangential spring and dashpot constants
    // using particle properties
    const double rp_sqrt = sqrt(particle_radius);

    const double normal_spring_constant =
      1.0667 * rp_sqrt * youngs_modulus *
      pow((0.9375 * particle_properties[PropertiesIndex::mass] * 1.0 *
           1.0 / // Characteristic velocity is set to 1.0
           (rp_sqrt * youngs_modulus)),
          0.2);

    // There is no minus sign here since model_parameter_beta is negative or
    // equal to zero.
    const double normal_damping_constant =
      2 * beta *
      sqrt(particle_properties[PropertiesIndex::mass] * normal_spring_constant);

    // REF :  R. Garg, J. Galvin-Carney, T. Li, and S. Pannala, “Documentation
    // of open-source MFIX–DEM software for gas-solids flows,” Tingwen Li Dr.,
    // p. 10, Sep. 2012. There is a minus sign since the tangential force is
    // applied in the opposite direction of the tangential_displacement
    const double tangential_spring_constant = -normal_spring_constant * 0.4;

    const double tangential_damping_constant =
      normal_damping_constant * 0.6324555320336759; // sqrt(0.4)

    // Calculation of normal force using spring and dashpot normal forces
    normal_force = (normal_spring_constant * normal_overlap +
                    normal_damping_constant * normal_relative_velocity_value) *
                   normal_vector;

    // Calculation of tangential force
    tangential_force =
      (tangential_spring_constant * contact_info.tangential_displacement +
       tangential_damping_constant * tangential_relative_velocity);

    const double coulomb_threshold = friction_coeff * normal_force.norm();
    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement and tangential
        // force are limited to Coulomb's criterion
        tangential_force =
          coulomb_threshold * (tangential_force / tangential_force.norm());

        contact_info.tangential_displacement =
          tangential_force / (tangential_spring_constant + DBL_MIN);
      }

    // Calculation torque caused by tangential force
    // We add the minus sign here since the tangential_force is applied on the
    // particle in the opposite direction
    tangential_torque =
      cross_product_3d((particle_radius * normal_vector), -tangential_force);

    rolling_resistance_torque = calculate_rolling_resistance_torque(
      particle_radius,
      particle_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      contact_info.normal_vector,
      contact_info.rolling_resistance_spring_torque);
  }

  /**
   * @brief Calculate the particle-wall non-linear contact force and torque
   * based on the updated values in contact_info.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle-wall pair.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative velocity value.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_properties Properties of particle.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] tangential_torque Contact tangential torque..
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_nonlinear_contact(
    particle_wall_contact_info<dim> &contact_info,
    const Tensor<1, 3>              &tangential_relative_velocity,
    const double                     normal_relative_velocity_value,
    const double                     normal_overlap,
    const double                     dt,
    const ArrayView<const double>   &particle_properties,
    Tensor<1, 3>                    &normal_force,
    Tensor<1, 3>                    &tangential_force,
    Tensor<1, 3>                    &tangential_torque,
    Tensor<1, 3>                    &rolling_resistance_torque)
  {
    // i is the particle, j is the wall.
    // there is a minus sign in front of the normal_vector to respect the
    // convention (i -> j), the forces are thus calculated on the wall
    const Tensor<1, 3> normal_vector = -contact_info.normal_vector;
    const double       particle_radius =
      particle_properties[PropertiesIndex::dp] * 0.5;
    const unsigned int particle_type =
      particle_properties[PropertiesIndex::type];

    const double youngs_modulus = this->effective_youngs_modulus[particle_type];
    const double shear_modulus  = this->effective_shear_modulus[particle_type];
    const double beta           = this->model_parameter_beta[particle_type];
    const double friction_coeff =
      this->effective_coefficient_of_friction[particle_type];
    const double rolling_viscous_damping_coeff =
      this->effective_coefficient_of_rolling_viscous_damping[particle_type];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[particle_type];

    // Calculate intermediate model parameters
    // These values are used to consider non-linear relation of the contact
    // force to the normal overlap
    const double radius_times_overlap_sqrt =
      sqrt(particle_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * youngs_modulus * radius_times_overlap_sqrt;
    const double model_parameter_st =
      8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculation of normal and tangential spring and dashpot constants
    // using particle and wall properties
    const double normal_spring_constant =
      1.3333 * youngs_modulus * radius_times_overlap_sqrt;

    // There is no minus sign here since beta is negative or
    // equal to zero.
    const double normal_damping_constant =
      1.8257 * beta *
      sqrt(model_parameter_sn * particle_properties[PropertiesIndex::mass]);

    // There is a minus sign since the tangential force is applied in the
    // opposite direction of the tangential_displacement
    const double tangential_spring_constant =
      -8.0 * shear_modulus * radius_times_overlap_sqrt + DBL_MIN;

    const double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of normal force using spring and dashpot normal forces
    normal_force = (normal_spring_constant * normal_overlap +
                    normal_damping_constant * normal_relative_velocity_value) *
                   normal_vector;

    // Calculation of tangential force. Since we need damping tangential force
    // in the gross sliding again, we define it as a separate variable
    const Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;
    tangential_force =
      tangential_spring_constant * contact_info.tangential_displacement +
      damping_tangential_force;

    const double coulomb_threshold = friction_coeff * normal_force.norm();

    // Check for gross sliding
    const double tangential_force_norm = tangential_force.norm();
    if (tangential_force_norm > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement and tangential
        // force are limited to Coulomb's criterion
        contact_info.tangential_displacement =
          (coulomb_threshold *
             (tangential_force / (tangential_force_norm + DBL_MIN)) -
           damping_tangential_force) /
          (tangential_spring_constant + DBL_MIN);

        tangential_force =
          (tangential_spring_constant * contact_info.tangential_displacement) +
          damping_tangential_force;
      }

    // Calculating torque caused by tangential force
    // We add the minus sign here since the tangential_force applied on the
    // particle is in the opposite direction
    tangential_torque =
      cross_product_3d((particle_radius * normal_vector), -tangential_force);

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      particle_radius,
      particle_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      contact_info.normal_vector,
      contact_info.rolling_resistance_spring_torque);
  }

  /**
   * @brief Calculate the particle-wall JKR contact force and torque
   * based on the updated values in contact_info.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle-wall pair.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative velocity value.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_properties Properties of particle.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] tangential_torque Contact tangential torque.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_JKR_contact(particle_wall_contact_info<dim> &contact_info,
                        const Tensor<1, 3> &tangential_relative_velocity,
                        const double        normal_relative_velocity_value,
                        const double        normal_overlap,
                        const double        dt,
                        const ArrayView<const double> &particle_properties,
                        Tensor<1, 3>                  &normal_force,
                        Tensor<1, 3>                  &tangential_force,
                        Tensor<1, 3>                  &tangential_torque,
                        Tensor<1, 3> &rolling_resistance_torque)
  {
    // i is the particle, j is the wall.
    // there is a minus sign in front of the normal_vector to respect the
    // convention (i -> j), the forces are thus calculated on the wall
    const Tensor<1, 3> normal_vector = -contact_info.normal_vector;
    const double       particle_radius =
      0.5 * particle_properties[PropertiesIndex::dp];
    const unsigned int particle_type =
      particle_properties[PropertiesIndex::type];

    const double youngs_modulus = this->effective_youngs_modulus[particle_type];
    const double shear_modulus  = this->effective_shear_modulus[particle_type];
    const double beta           = this->model_parameter_beta[particle_type];
    const double friction_coeff =
      this->effective_coefficient_of_friction[particle_type];
    const double rolling_viscous_damping_coeff =
      this->effective_coefficient_of_rolling_viscous_damping[particle_type];
    const double rolling_friction_coeff =
      this->effective_coefficient_of_rolling_friction[particle_type];
    const double surface_energy = this->effective_surface_energy[particle_type];

    // Calculate intermediate model parameters
    // These values are used to consider non-linear relation of the contact
    // force to the normal overlap
    const double radius_times_overlap_sqrt =
      sqrt(particle_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * youngs_modulus * radius_times_overlap_sqrt;
    const double model_parameter_st =
      8.0 * shear_modulus * radius_times_overlap_sqrt;

    // Calculation of the contact patch radius (a) using the analytical solution
    // described in the theory guide.
    const double c0 =
      Utilities::fixed_power<2>(particle_radius * normal_overlap);
    const double c1 = -2. * Utilities::fixed_power<2>(particle_radius) * M_PI *
                      surface_energy / youngs_modulus;
    const double c2 = -2. * normal_overlap * particle_radius;
    const double P  = -Utilities::fixed_power<2>(c2) / 12. - c0;
    const double Q  = -Utilities::fixed_power<3>(c2) / 108. + c0 * c2 / 3. -
                     Utilities::fixed_power<2>(c1) * 0.125;
    const double root1  = std::max(0.,
                                  (0.25 * Utilities::fixed_power<2>(Q)) +
                                    (Utilities::fixed_power<3>(P) / 27.));
    const double U      = std::cbrt(-0.5 * Q + std::sqrt(root1));
    const double s      = -c2 * (5. / 6.) + U - P / (3. * U);
    const double w      = std::sqrt(std::max(1e-16, c2 + 2. * s));
    const double lambda = 0.5 * c1 / w;
    const double root2  = std::max(1e-16, w * w - 4. * (c2 + s + lambda));
    const double a      = 0.5 * (w + std::sqrt(root2));

    // Calculation of normal damping and tangential spring and dashpot constants
    // using particle and wall properties.
    // There is no minus sign here since model_parameter_beta is negative or
    // equal to zero.
    const double normal_damping_constant =
      1.8257 * beta * // 2. * sqrt(5./6.)
      sqrt(model_parameter_sn * particle_properties[PropertiesIndex::mass]);

    // Tangential spring constant is set as a negative just like in the other
    // particle-wall models. This must be taken into account for the square root
    // in the tangential_damping_calculation.
    const double tangential_spring_constant =
      -8.0 * shear_modulus * radius_times_overlap_sqrt;

    // There is no minus sign here since model_parameter_beta is negative or
    // equal to zero.
    const double tangential_damping_constant =
      normal_damping_constant *
      sqrt(model_parameter_st / (model_parameter_sn + DBL_MIN));

    // Calculation of the normal force coefficient (F_n_JKR)
    const double normal_force_norm =
      4. * youngs_modulus * Utilities::fixed_power<3>(a) /
        (3. * particle_radius) -
      std::sqrt(8. * M_PI * surface_energy * youngs_modulus *
                Utilities::fixed_power<3>(a)) +
      normal_damping_constant * normal_relative_velocity_value;

    // Calculation of normal force using the normal_force_norm and the
    // normal vector.
    normal_force = normal_force_norm * normal_vector;

    // Calculation of tangential forces.
    Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;
    tangential_force =
      tangential_spring_constant * contact_info.tangential_displacement +
      damping_tangential_force;

    // JKR theory says that the coulomb threshold must be modified with the
    // pull-out force. (Thornton 1991)
    const double modified_coulomb_threshold =
      (normal_force_norm + 3. * M_PI * surface_energy * particle_radius) *
      friction_coeff;

    // Check for gross sliding
    const double tangential_force_norm = tangential_force.norm();
    if (tangential_force_norm > modified_coulomb_threshold)
      {
        // Gross sliding occurs and the tangential displacement and tangential
        // force are limited to Coulomb's criterion
        contact_info.tangential_displacement =
          (modified_coulomb_threshold *
             (tangential_force / (tangential_force_norm + DBL_MIN)) -
           damping_tangential_force) /
          (tangential_spring_constant + DBL_MIN);

        tangential_force =
          (tangential_spring_constant * contact_info.tangential_displacement) +
          damping_tangential_force;
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    // We add the minus sign here since the tangential_force applied on the
    // particle is in the opposite direction
    tangential_torque =
      cross_product_3d((particle_radius * normal_vector), -tangential_force);

    // We need to compute the normal spring constant for case where we use the
    // EPSD rolling resistance model.
    const double normal_spring_constant = 0.66665 * model_parameter_sn;

    // Rolling resistance torque
    rolling_resistance_torque = calculate_rolling_resistance_torque(
      particle_radius,
      particle_properties,
      rolling_friction_coeff,
      rolling_viscous_damping_coeff,
      dt,
      normal_spring_constant,
      normal_force.norm(),
      contact_info.normal_vector,
      contact_info.rolling_resistance_spring_torque);
  }

  /**
   * @brief Calculate the particle-wall cohesive force and
   * contact torque based on the updated values in contact_info. It uses the DMT
   * cohesive force model.
   *
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle-wall.
   * @param[in] tangential_relative_velocity Tangential relative velocity.
   * @param[in] normal_relative_velocity_value Normal relative velocity value.
   * @param[in] normal_overlap Contact normal overlap.
   * @param[in] dt DEM time step.
   * @param[in] particle_properties Properties of particle.
   * @param[out] normal_force Contact normal force.
   * @param[out] tangential_force Contact tangential force.
   * @param[out] tangential_torque Contact tangential torque.
   * @param[out] rolling_resistance_torque Contact rolling resistance torque.
   */
  inline void
  calculate_DMT_contact(particle_wall_contact_info<dim> &contact_info,
                        const Tensor<1, 3> &tangential_relative_velocity,
                        const double        normal_relative_velocity_value,
                        const double        normal_overlap,
                        const double        dt,
                        const ArrayView<const double> &particle_properties,
                        Tensor<1, 3>                  &normal_force,
                        Tensor<1, 3>                  &tangential_force,
                        Tensor<1, 3>                  &tangential_torque,
                        Tensor<1, 3> &rolling_resistance_torque)
  {
    constexpr double M_2PI = 2. * M_PI;

    // i is the particle, j is the wall.
    // there is a minus sign in front of the normal_vector to respect the
    // convention (i -> j), the forces are thus calculated on the wall
    const Tensor<1, 3> normal_vector = -contact_info.normal_vector;
    const double       particle_radius =
      0.5 * particle_properties[PropertiesIndex::dp];
    const unsigned int particle_type =
      particle_properties[PropertiesIndex::type];

    const double surface_energy = this->effective_surface_energy[particle_type];
    const double hamaker_constant =
      this->effective_hamaker_constant[particle_type];

    const double F_po = M_2PI * particle_radius * surface_energy;

    const double delta_0 =
      -std::sqrt(hamaker_constant * particle_radius / (6. * F_po));

    // Cohesive force.
    double cohesive_term;

    // Contact particle-wall + constant cohesive force.
    if (normal_overlap > 0.)
      {
        cohesive_term = -F_po;

        calculate_nonlinear_contact(contact_info,
                                    tangential_relative_velocity,
                                    normal_relative_velocity_value,
                                    normal_overlap,
                                    dt,
                                    particle_properties,
                                    normal_force,
                                    tangential_force,
                                    tangential_torque,
                                    rolling_resistance_torque);
      }
    // No contact, but still in the constant zone for the cohesive force.
    else if (normal_overlap > delta_0)
      {
        cohesive_term = -F_po;
        contact_info.tangential_displacement.clear();
        contact_info.rolling_resistance_spring_torque.clear();
      }
    // No contact. The cohesive force is not constant. It needs to be computed.
    else
      {
        cohesive_term = -hamaker_constant * particle_radius /
                        (6. * Utilities::fixed_power<2>(normal_overlap));
        contact_info.tangential_displacement.clear();
        contact_info.rolling_resistance_spring_torque.clear();
      }
    normal_force += cohesive_term * normal_vector;
  }


  /**
   * @brief Set every containers needed to carry the particle-wall force
   * calculation.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm file.
   */
  void
  set_effective_properties(const DEMSolverParameters<dim> &dem_parameters);

  /**
   * @brief Set every containers needed to carry the heat transfer rate
   * calculation.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm
   * file.
   */
  void
  set_multiphysic_properties(const DEMSolverParameters<dim> &dem_parameters);

  // Members of the class

  unsigned int        n_particle_types;
  std::vector<double> effective_youngs_modulus;
  std::vector<double> effective_real_youngs_modulus;
  std::vector<double> effective_shear_modulus;
  std::vector<double> effective_coefficient_of_restitution;
  std::vector<double> effective_coefficient_of_friction;
  std::vector<double> effective_coefficient_of_rolling_friction;
  std::vector<double> effective_coefficient_of_rolling_viscous_damping;
  std::vector<double> effective_surface_energy;
  std::vector<double> effective_hamaker_constant;
  std::vector<double> model_parameter_beta;
  std::vector<double> equivalent_surface_roughness;
  std::vector<double> equivalent_surface_slope;
  std::vector<double> effective_microhardness;
  std::vector<double> particle_thermal_conductivity;
  std::vector<double> gas_parameter_m;
  double              gas_thermal_conductivity;
  double              wall_thermal_conductivity;
  const double        dmt_cut_off_threshold;
  const double        f_coefficient_epsd;

  std::unordered_map<unsigned int, double> boundary_rotational_speed_map;
  std::unordered_map<unsigned int, Tensor<1, 3>>
    boundary_translational_velocity_map;
  std::unordered_map<unsigned int, Tensor<1, 3>> boundary_rotational_vector;
  std::unordered_map<unsigned int, Point<3>>     point_on_rotation_vector;
  const unsigned int                             vertices_per_triangle = 3;
  Point<3>                                       center_mass_container;
};

#endif
