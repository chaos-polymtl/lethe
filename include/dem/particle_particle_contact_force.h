/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 *
 */

#include <core/auxiliary_math_functions.h>
#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_contact_info.h>
#include <dem/rolling_resistance_torque_models.h>

#include <deal.II/particles/particle_handler.h>

#include <boost/range/adaptor/map.hpp>

#include <vector>

using namespace dealii;

#ifndef particle_particle_contact_force_h
#  define particle_particle_contact_force_h

using namespace DEM;


/**
 * Base class for the particle-particle contact force models
 * This class does not implement any of the models, but ensures that
 * an interface without template specialization is available. All of the
 * actual implementation of the models are carried out in the
 * ParticleParticleContactForce class which is templated by the contact model
 * type.
 */
template <int dim>
class ParticleParticleContactForceBase
{
public:
  /**
   * Carries out the calculation of the contact force using the contact pair
   * information obtained in the fine search and physical properties of
   * particles
   *
   * @param local_adjacent_particles Required information for calculation of the
   * local-local particle-particle contact forces. The information was
   * obtained in the fine search
   * @param ghost_adjacent_particles Required information for calculation of the
   * local-ghost particle-particle contact forces. The information was
   * obtained in the fine search
   * @param dt DEM time step
   * @param torque An unordered_map of torque of particles
   * @param force Force acting on particles
   * @param periodic_offset A tensor of the periodic offset to change the
   * particle location of the particles on the periodic boundary 1 side
   */
  virtual void
  calculate_particle_particle_contact_force(
    DEMContactManager<dim>    &container_manager,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const Tensor<1, dim>       periodic_offset = Tensor<1, dim>()) = 0;
};

/**
 * @brief Class that carries out the calculation of
 * particle-particle contact force including non-linear and linear contact
 * models. Instead of using a inheritance hiearchy to distinguish between
 * the contact model, the class is templated with the type of force model
 * and rolling friction model. Consequently, the code for each
 * combination of force model is generated at compile time.
 *
 * @tparam dim The dimension of the problem
 * @tparam force_model The particle-particle contact force model
 * @tparam rolling_friction_model The rolling resistance model used
 */
template <
  int                                                       dim,
  Parameters::Lagrangian::ParticleParticleContactForceModel force_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
class ParticleParticleContactForce
  : public ParticleParticleContactForceBase<dim>
{
public:
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters);


  virtual ~ParticleParticleContactForce()
  {}

  /**
   * Carries out the calculation of the contact force using the contact pair
   * information
   * obtained in the fine search and physical properties of particles
   *
   * @param local_adjacent_particles Required information for calculation of the
   * loacl-local particle-particle contact force. These information were
   * obtained in the fine search
   * @param ghost_adjacent_particles Required information for calculation of the
   * loacl-ghost particle-particle contact force. These information were
   * obtained in the fine search
   * @param dt DEM time step
   * @param torque An unordered_map of torque of particles
   * @param force Force acting on particles
   * @param periodic_offset A tensor of the periodic offset to change the
   * particle location of the particles on the periodic boundary 1 side
   */
  virtual void
  calculate_particle_particle_contact_force(
    DEMContactManager<dim>    &container_manager,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const Tensor<1, dim>       periodic_offset = Tensor<1, dim>()) override;


protected:
  /**
   * @brief Carries out updating the contact pair information for both non-linear and
   * linear contact force calculations
   *
   * @param contact_info Contact information of a particle pair in
   * neighborhood
   * @param normal_unit_vector Normal vector of the contact. This vector is particle_two_location - particle_one_location.
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param particle_one_location Location of particle one in contact
   * @param particle_two_location Location of particle two in contact
   * @param dt DEM time step
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

    // Calculation of new tangential_overlap, since this value is
    // history-dependent it needs the value at previous time-step
    // This variable is the main reason that we have iteration over
    // two different vectors : tangential_overlap of the particles
    // which were already in contact needs to
    // modified using its history, while the tangential_overlaps of
    // new particles are equal to zero
    // delta_t_new = delta_t_old + v_rt*dt
    contact_info.tangential_overlap += tangential_relative_velocity * dt;
  }

  /**
   * @brief Carries out applying the calculated force and torque on the local-local
   * particle pair in contact, for both non-linear and linear contact force
   * calculations
   *
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   * @param particle_one_torque
   * @param particle_two_torque
   * @param particle_one_force Force acting on particle one
   * @param particle_two_force Force acting on particle two
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
   * Carries out applying the calculated force and torque on the local-ghost
   * particle pair in contact, for both non-linear and linear contact force
   * calculations. The contact force is only applied on the local particles
   *
   * @param normal_force normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   * @param particle_one_torque Torque acting on particle one (local)
   * @param particle_one_force Force acting on particle one
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
   * Carries out the calculation of effective mass and radius of particles i and
   * j in contact.
   *
   * @param particle_one_properties Properties of particle one in
   * contact
   * @param particle_two_properties Properties of particle two in
   * contact
   */
  inline void
  find_effective_radius_and_mass(
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties)
  {
    effective_mass = (particle_one_properties[DEM::PropertiesIndex::mass] *
                      particle_two_properties[DEM::PropertiesIndex::mass]) /
                     (particle_one_properties[DEM::PropertiesIndex::mass] +
                      particle_two_properties[DEM::PropertiesIndex::mass]);
    effective_radius =
      (particle_one_properties[DEM::PropertiesIndex::dp] *
       particle_two_properties[DEM::PropertiesIndex::dp]) /
      (2 * (particle_one_properties[DEM::PropertiesIndex::dp] +
            particle_two_properties[DEM::PropertiesIndex::dp]));
  }

  /**
   * Carries out the calculation of the particle-particle linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param tangential_relative_velocity Tangential relative velocity
   * @param normal_relative_velocity_value Normal relative contact velocity
   * @param normal_unit_vector Contact normal unit vector
   * @param normal_overlap Contact normal overlap
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param particle_one_tangential_torque Contact tangential torque on particle one
   * @param particle_two_tangential_torque Contact tangential torque on particle two
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  inline void
  calculate_linear_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    this->find_effective_radius_and_mass(particle_one_properties,
                                         particle_two_properties);
    const unsigned int particle_one_type =
      particle_one_properties[PropertiesIndex::type];
    const unsigned int particle_two_type =
      particle_two_properties[PropertiesIndex::type];

    const double characteristic_velocity =
      1.0; // Characteristic velocity is set at 1.0 so that the normal and
           // tangential spring_constant remain constant throughout a
           // simulation.

    // Calculation of normal and tangential spring and dashpot constants
    // using particle properties
    double normal_spring_constant =
      1.0667 * sqrt(this->effective_radius) *
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_properties[PropertiesIndex::type])] *
      pow((0.9375 * this->effective_mass * characteristic_velocity *
           characteristic_velocity /
           (sqrt(this->effective_radius) *
            this->effective_youngs_modulus[vec_particle_type_index(
              particle_one_type, particle_two_type)])),
          0.2);

    // REF :  R. Garg, J. Galvin-Carney, T. Li, and S. Pannala,“Documentation of
    // open-source MFIX–DEM software for gas-solids flows,” Tingwen Li Dr., p.
    // 10, Sep. 2012.
    double tangential_spring_constant = normal_spring_constant * 0.4;

    double normal_damping_constant =
      -2 *
      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] *
      sqrt(this->effective_mass * normal_spring_constant);

    double tangential_damping_constant =
      normal_damping_constant * 0.6324555320336759; // sqrt(0.4)
    // Calculation of the normal force
    normal_force = (normal_spring_constant * normal_overlap +
                    normal_damping_constant * normal_relative_velocity_value) *
                   normal_unit_vector;

    // Calculation of tangential force. Since we need damping tangential force
    // in the gross sliding again, we define it as a separate variable
    Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;

    tangential_force =
      (tangential_spring_constant * contact_info.tangential_overlap) +
      damping_tangential_force;

    double coulomb_threshold =
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      normal_force.norm();


    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential overlap and tangential
        // force are limited to Coulomb's criterion
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
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance)
      rolling_resistance_torque = no_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance)
      rolling_resistance_torque = constant_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::viscous_resistance)
      rolling_resistance_torque = viscous_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
  }


  /**
   * @brief Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param tangential_relative_velocity Tangential relative velocity
   * @param normal_relative_velocity_value Normal relative contact velocity
   * @param normal_unit_vector Contact normal unit vector
   * @param normal_overlap Contact normal overlap
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param particle_one_tangential_torque Contact tangential torque on particle one
   * @param particle_two_tangential_torque Contact tangential torque on particle two
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  inline void
  calculate_hertz_mindlin_limit_overlap_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    this->find_effective_radius_and_mass(particle_one_properties,
                                         particle_two_properties);

    const unsigned int particle_one_type =
      particle_one_properties[PropertiesIndex::type];
    const unsigned int particle_two_type =
      particle_two_properties[PropertiesIndex::type];

    const double radius_times_overlap_sqrt =
      sqrt(this->effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 *
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      radius_times_overlap_sqrt;
    double model_parameter_st =
      8.0 *
      this->effective_shear_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      radius_times_overlap_sqrt;

    // Calculation of normal and tangential spring and dashpot constants
    // using particle properties
    double normal_spring_constant = 0.66665 * model_parameter_sn;
    double normal_damping_constant =
      -1.8257 *
      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] *
      sqrt(model_parameter_sn * this->effective_mass);
    double tangential_spring_constant =
      8.0 *
        this->effective_shear_modulus[vec_particle_type_index(
          particle_one_type, particle_two_type)] *
        radius_times_overlap_sqrt +
      DBL_MIN;
    double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of normal force
    const double normal_force_norm =
      normal_spring_constant * normal_overlap +
      normal_damping_constant * normal_relative_velocity_value;
    normal_force = normal_force_norm * normal_unit_vector;

    // Calculation of tangential force. Since we need damping tangential force
    // in the gross sliding again, we define it as a separate variable
    Tensor<1, 3> damping_tangential_force =
      tangential_damping_constant * tangential_relative_velocity;
    tangential_force =
      (tangential_spring_constant * contact_info.tangential_overlap) +
      damping_tangential_force;

    double coulomb_threshold =
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      normal_force_norm;

    // Check for gross sliding
    const double tangential_force_norm = tangential_force.norm();
    if (tangential_force_norm > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential overlap and tangential
        // force are limited to Coulomb's criterion
        contact_info.tangential_overlap =
          (coulomb_threshold *
             (tangential_force / (tangential_force_norm + DBL_MIN)) -
           damping_tangential_force) /
          (tangential_spring_constant + DBL_MIN);

        tangential_force =
          (tangential_spring_constant * contact_info.tangential_overlap) +
          damping_tangential_force;
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force *
                         particle_one_properties[PropertiesIndex::dp] * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque *
      particle_two_properties[PropertiesIndex::dp] /
      particle_one_properties[PropertiesIndex::dp];


    // Rolling resistance torque
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance)
      rolling_resistance_torque = no_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance)
      rolling_resistance_torque = constant_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::viscous_resistance)
      rolling_resistance_torque = viscous_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
  }


  /**
   * @brief Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param tangential_relative_velocity Tangential relative velocity
   * @param normal_relative_velocity_value Normal relative contact velocity
   * @param normal_unit_vector Contact normal unit vector
   * @param normal_overlap Contact normal overlap
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  inline void
  calculate_hertz_mindlin_limit_force_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    this->find_effective_radius_and_mass(particle_one_properties,
                                         particle_two_properties);

    const unsigned int particle_one_type =
      particle_one_properties[PropertiesIndex::type];
    const unsigned int particle_two_type =
      particle_two_properties[PropertiesIndex::type];

    const double radius_times_overlap_sqrt =
      sqrt(this->effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 *
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      radius_times_overlap_sqrt;
    double model_parameter_st =
      8.0 *
      this->effective_shear_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      radius_times_overlap_sqrt;

    // Calculation of normal and tangential spring and dashpot constants
    // using particle properties
    double normal_spring_constant = 0.66665 * model_parameter_sn;
    double normal_damping_constant =
      -1.8257 *
      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] *
      sqrt(model_parameter_sn * this->effective_mass);
    double tangential_spring_constant =
      8.0 *
        this->effective_shear_modulus[vec_particle_type_index(
          particle_one_type, particle_two_type)] *
        radius_times_overlap_sqrt +
      DBL_MIN;
    double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of normal force using spring and dashpot normal forces
    normal_force =
      ((normal_spring_constant * normal_overlap) * normal_unit_vector) +
      ((normal_damping_constant * normal_relative_velocity_value) *
       normal_unit_vector);

    // Calculation of tangential force using spring and dashpot tangential
    // forces. Since we need dashpot tangential force in the gross sliding
    // again, we define it as a separate variable
    tangential_force =
      tangential_spring_constant * contact_info.tangential_overlap +
      tangential_damping_constant * tangential_relative_velocity;

    double coulomb_threshold =
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      normal_force.norm();

    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential overlap and tangential
        // force are limited to Coulomb's criterion
        tangential_force =
          coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force *
                         particle_one_properties[PropertiesIndex::dp] * 0.5);

    particle_two_tangential_torque =
      particle_one_tangential_torque *
      particle_two_properties[PropertiesIndex::dp] /
      particle_one_properties[PropertiesIndex::dp];


    // Rolling resistance torque
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance)
      rolling_resistance_torque = no_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance)
      rolling_resistance_torque = constant_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::viscous_resistance)
      rolling_resistance_torque = viscous_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
  }

  /**
   * @brief Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param tangential_relative_velocity Tangential relative velocity
   * @param normal_relative_velocity_value Normal relative contact velocity
   * @param normal_unit_vector Contact normal unit vector
   * @param normal_overlap Contact normal overlap
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  inline void
  calculate_hertz_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3> & /*tangential_relative_velocity*/,
    const double                   normal_relative_velocity_value,
    const Tensor<1, 3>            &normal_unit_vector,
    const double                   normal_overlap,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    Tensor<1, 3>                  &normal_force,
    Tensor<1, 3>                  &tangential_force,
    Tensor<1, 3>                  &particle_one_tangential_torque,
    Tensor<1, 3>                  &particle_two_tangential_torque,
    Tensor<1, 3>                  &rolling_resistance_torque)
  {
    // Calculation of effective radius and mass
    this->find_effective_radius_and_mass(particle_one_properties,
                                         particle_two_properties);

    const unsigned int particle_one_type =
      particle_one_properties[PropertiesIndex::type];
    const unsigned int particle_two_type =
      particle_two_properties[PropertiesIndex::type];

    const double radius_times_overlap_sqrt =
      sqrt(this->effective_radius * normal_overlap);
    const double model_parameter_sn =
      2 *
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      radius_times_overlap_sqrt;

    // Calculation of normal and tangential spring and dashpot constants
    // using particle properties
    double normal_spring_constant = 0.66665 * model_parameter_sn;
    double normal_damping_constant =
      -1.8257 *
      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)] *
      sqrt(model_parameter_sn * this->effective_mass);
    double tangential_spring_constant =
      8.0 *
        this->effective_shear_modulus[vec_particle_type_index(
          particle_one_type, particle_two_type)] *
        radius_times_overlap_sqrt +
      DBL_MIN;


    // Calculation of normal force using spring and dashpot normal forces
    normal_force =
      ((normal_spring_constant * normal_overlap) * normal_unit_vector) +
      ((normal_damping_constant * normal_relative_velocity_value) *
       normal_unit_vector);

    // Calculation of tangential force using spring and dashpot tangential
    // forces. Since we need dashpot tangential force in the gross sliding
    // again, we define it as a separate variable
    tangential_force =
      tangential_spring_constant * contact_info.tangential_overlap;

    double coulomb_threshold =
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      normal_force.norm();

    // Check for gross sliding
    if (tangential_force.norm() > coulomb_threshold)
      {
        // Gross sliding occurs and the tangential overlap and tangential
        // force are limited to Coulumb's criterion
        tangential_force =
          coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));
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
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance)
      rolling_resistance_torque = no_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance)
      rolling_resistance_torque = constant_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::viscous_resistance)
      rolling_resistance_torque = viscous_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
  }


  /**
   * @brief Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param tangential_relative_velocity Tangential relative velocity
   * @param normal_relative_velocity_value Normal relative contact velocity
   * @param normal_unit_vector Contact normal unit vector
   * @param normal_overlap Contact normal overlap
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param particle_one_tangential_torque Contact tangential torque on particle one
   * @param particle_two_tangential_torque Contact tangential torque on particle two
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  inline void
  calculate_hertz_JKR_contact(
    particle_particle_contact_info<dim> &contact_info,
    const Tensor<1, 3>                  &tangential_relative_velocity,
    const double                         normal_relative_velocity_value,
    const Tensor<1, 3>                  &normal_unit_vector,
    const double                         normal_overlap,
    const ArrayView<const double>       &particle_one_properties,
    const ArrayView<const double>       &particle_two_properties,
    Tensor<1, 3>                        &normal_force,
    Tensor<1, 3>                        &tangential_force,
    Tensor<1, 3>                        &particle_one_tangential_torque,
    Tensor<1, 3>                        &particle_two_tangential_torque,
    Tensor<1, 3>                        &rolling_resistance_torque)
  {
    const unsigned int particle_one_type =
      particle_one_properties[PropertiesIndex::type];
    const unsigned int particle_two_type =
      particle_two_properties[PropertiesIndex::type];

    // Model parameters
    this->find_effective_radius_and_mass(particle_one_properties,
                                         particle_two_properties);
    const double radius_times_overlap_sqrt =
      sqrt(this->effective_radius * normal_overlap);
    const double model_parameter_sn =
      2.0 * radius_times_overlap_sqrt *
      this->effective_youngs_modulus[vec_particle_type_index(
        particle_one_type, particle_two_type)];
    double model_parameter_st =
      8.0 * radius_times_overlap_sqrt *
      this->effective_shear_modulus[vec_particle_type_index(particle_one_type,
                                                            particle_two_type)];

    // Calculation of the  contact path radius using the Ferrari analitycal
    // solution.
    this->find_effective_radius_and_mass(particle_one_properties,
                                         particle_two_properties);
    const double c0 =
      Utilities::fixed_power<2>(this->effective_radius * normal_overlap);
    const double c1 = -2. * Utilities::fixed_power<2>(this->effective_radius) *
                      M_PI *
                      this->effective_surface_energy[vec_particle_type_index(
                        particle_one_type, particle_two_type)] /
                      this->effective_youngs_modulus[vec_particle_type_index(
                        particle_one_type, particle_two_type)];

    const double c2 = -2. * normal_overlap * this->effective_radius;
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
      -1.8257 * sqrt(model_parameter_sn * this->effective_mass) *
      this->model_parameter_beta[vec_particle_type_index(particle_one_type,
                                                         particle_two_type)];
    // Calculation of the tangential spring constant
    const double tangential_spring_constant =
      8.0 * radius_times_overlap_sqrt *
        this->effective_shear_modulus[vec_particle_type_index(
          particle_one_type, particle_two_type)] +
      DBL_MIN;
    // Calculation of the tangential damping constant
    const double tangential_damping_constant =
      normal_damping_constant * sqrt(model_parameter_st / model_parameter_sn);

    // Calculation of the normal force coefficient (F_n_JKR) # Eq 20
    const double normal_force_coefficient =
      4. * Utilities::fixed_power<3>(a) / (3. * this->effective_radius) *
        this->effective_youngs_modulus[vec_particle_type_index(
          particle_one_type, particle_two_type)] -
      std::sqrt(8 * M_PI *
                this->effective_surface_energy[vec_particle_type_index(
                  particle_one_type, particle_two_type)] *
                this->effective_youngs_modulus[vec_particle_type_index(
                  particle_one_type, particle_two_type)] *
                Utilities::fixed_power<3>(a));

    // Calculation of the final normal force vector
    normal_force = (normal_force_coefficient +
                    normal_damping_constant * normal_relative_velocity_value) *
                   normal_unit_vector;

    tangential_force =
      tangential_spring_constant * contact_info.tangential_overlap +
      tangential_damping_constant * tangential_relative_velocity;

    // JKR theory says that the coulomb threshold must be modified with the
    // pull-out force.
    const double pull_off_force =
      3. * M_PI *
      this->effective_surface_energy[vec_particle_type_index(
        particle_one_type, particle_two_type)] *
      this->effective_radius;
    const double modified_coulomb_threshold =
      (normal_force_coefficient + 2. * pull_off_force) *
      this->effective_coefficient_of_friction[vec_particle_type_index(
        particle_one_type, particle_two_type)];

    if (tangential_force.norm() > modified_coulomb_threshold)
      {
        // Gross sliding occurs and the tangential overlap and tangential
        // force are limited to Coulumb's criterion
        tangential_force =
          modified_coulomb_threshold *
          (tangential_force / (tangential_force.norm() + DBL_MIN));
      }

    // Calculation of torque caused by tangential force (tangential_torque)
    particle_one_tangential_torque =
      cross_product_3d(normal_unit_vector,
                       tangential_force *
                         particle_one_properties[PropertiesIndex::dp] * 0.5);
    particle_two_tangential_torque =
      particle_one_tangential_torque *
      particle_two_properties[PropertiesIndex::dp] /
      particle_one_properties[PropertiesIndex::dp];

    // Rolling resistance torque
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    no_resistance)
      rolling_resistance_torque = no_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::RollingResistanceMethod::
                    constant_resistance)
      rolling_resistance_torque = constant_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
    if constexpr (rolling_friction_model ==
                  Parameters::Lagrangian::viscous_resistance)
      rolling_resistance_torque = viscous_rolling_resistance_torque(
        this->effective_radius,
        particle_one_properties,
        particle_two_properties,
        this->effective_coefficient_of_rolling_friction[vec_particle_type_index(
          particle_one_type, particle_two_type)],
        normal_force.norm(),
        normal_unit_vector);
  }

private:
  inline unsigned int
  vec_particle_type_index(const unsigned int i, const unsigned int j)
  {
    return i * n_particle_types + j;
  }

  // Members of the class
  // Contact model parameter. It is calculated in the constructor for different
  // combinations of particle types. For different combinations, a map of map is
  // used to store this variable
  unsigned int        n_particle_types;
  std::vector<double> effective_youngs_modulus;
  std::vector<double> effective_shear_modulus;
  std::vector<double> effective_coefficient_of_restitution;
  std::vector<double> effective_coefficient_of_friction;
  std::vector<double> effective_coefficient_of_rolling_friction;
  std::vector<double> effective_surface_energy;
  std::vector<double> model_parameter_beta;

  double effective_radius;
  double effective_mass;
};

#endif /* particle_particle_contact_force_h */
