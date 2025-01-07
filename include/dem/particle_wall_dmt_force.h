// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_dmt_force_h
#define lethe_particle_wall_dmt_force_h

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_nonlinear_force.h>

#include <deal.II/particles/particle.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace dealii;

/**
 * @brief Calculation of the DMT particle-wall contact force using the
 * information obtained from the fine search and physical properties of
 * particles and walls
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam solve_type Type of solver used for the DEM.
 */
template <int dim, typename PropertiesIndex>
class ParticleWallDMTForce : public ParticleWallNonLinearForce<dim, PropertiesIndex>
{
  using FuncPtrType = Tensor<1, 3> (ParticleWallDMTForce<dim, PropertiesIndex>::*)(
    const ArrayView<const double> &,
    const double,
    const double,
    const Tensor<1, 3> &);
  FuncPtrType calculate_rolling_resistance_torque;

public:
  ParticleWallDMTForce(
    const DEMSolverParameters<dim>        &dem_parameters,
    const std::vector<types::boundary_id> &boundary_index = {});

  /**
   * @brief Carries out the calculation of the particle-wall contact force using
   * DMT model
   *
   * @param particle_wall_pairs_in_contact Required information for the calculation
   * of the particle-wall contact force. These information were obtained in the
   * fine search
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                              &particle_wall_pairs_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force) override;

  /**
   * @brief Carries out the calculation of particle-floating mesh contact force
   * using DMT model
   *
   * @param particle_floating_mesh_in_contact A container that stores the information
   * of particle-floating mesh contact
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param solids Floating solids
   */
  virtual void
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                              &particle_floating_mesh_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids)
    override;

private:
  /**
   * @brief Carries out calculation of the rolling resistance torque using the no
   * resistance model
   *
   * @param particle_properties Particle properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   * @return rolling resistance torque
   */
  inline Tensor<1, 3>
  no_resistance(const ArrayView<const double> & /*particle_properties*/,
                const double /*effective_rolling_friction_coefficient*/,
                const double /*normal_force_norm*/,
                const Tensor<1, 3> & /*normal_contact_vector*/)
  {
    Tensor<1, 3> rolling_resistance({0, 0, 0});
    return rolling_resistance;
  }

  /**
   * @brief Carries out calculation of the rolling resistance torque using the
   * constant model
   *
   * @param particle_properties Particle properties
   * @param effective_rolling_friction_coefficient Effective rolling friction
   * coefficient
   * @param normal_force_norm Normal force norm
   * @return rolling resistance torque
   */
  inline Tensor<1, 3>
  constant_resistance(const ArrayView<const double> &particle_properties,
                      const double effective_rolling_friction_coefficient,
                      const double normal_force_norm,
                      const Tensor<1, 3> & /*normal_contact_vector*/)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] =
          particle_properties[PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    // Calculation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      (particle_properties[PropertiesIndex::dp] * 0.5) *
      normal_force_norm * particle_wall_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief Carries out calculation of the rolling resistance torque using the
   * viscous model
   *
   * @param particle_properties Particle properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   * @return rolling resistance torque
   */
  inline Tensor<1, 3>
  viscous_resistance(const ArrayView<const double> &particle_properties,
                     const double        effective_rolling_friction_coefficient,
                     const double        normal_force_norm,
                     const Tensor<1, 3> &normal_contact_vector)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] =
          particle_properties[PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    Tensor<1, 3> v_omega = cross_product_3d(
      angular_velocity,
      particle_properties[PropertiesIndex::dp] * 0.5 *
        normal_contact_vector);

    // Calculation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      particle_properties[PropertiesIndex::dp] * 0.5 *
      normal_force_norm * v_omega.norm() * particle_wall_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief Initialize every particle-wall properties required
   * for the force calculation. This function is required for the DMT force
   * model, since the dmt_cut_off_distance attribute needs the be set to const.
   * To do so, the effective surface energy and of the effective hamaker contant
   * need to be initialized before hand. For code clarity, we also initialized
   * every other parameter, otherwise the initialization would take place at two
   * different places in the code.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm file
   */
  void
  initialize_particle_wall_properties(
    const DEMSolverParameters<dim> &dem_parameters)
  {
    const double wall_youngs_modulus =
      dem_parameters.lagrangian_physical_properties.youngs_modulus_wall;
    const double wall_poisson_ratio =
      dem_parameters.lagrangian_physical_properties.poisson_ratio_wall;
    const double wall_restitution_coefficient =
      dem_parameters.lagrangian_physical_properties
        .restitution_coefficient_wall;
    const double wall_friction_coefficient =
      dem_parameters.lagrangian_physical_properties.friction_coefficient_wall;
    const double wall_rolling_friction_coefficient =
      dem_parameters.lagrangian_physical_properties.rolling_friction_wall;
    const double wall_surface_energy =
      dem_parameters.lagrangian_physical_properties.surface_energy_wall;
    const double wall_hamaker_constant =
      dem_parameters.lagrangian_physical_properties.hamaker_constant_wall;
    for (unsigned int i = 0;
         i < dem_parameters.lagrangian_physical_properties.particle_type_number;
         ++i)
      {
        const double particle_youngs_modulus =
          dem_parameters.lagrangian_physical_properties.youngs_modulus_particle
            .at(i);
        const double particle_poisson_ratio =
          dem_parameters.lagrangian_physical_properties.poisson_ratio_particle
            .at(i);
        const double particle_restitution_coefficient =
          dem_parameters.lagrangian_physical_properties
            .restitution_coefficient_particle.at(i);
        const double particle_friction_coefficient =
          dem_parameters.lagrangian_physical_properties
            .friction_coefficient_particle.at(i);
        const double particle_rolling_friction_coefficient =
          dem_parameters.lagrangian_physical_properties
            .rolling_friction_coefficient_particle.at(i);
        const double particle_surface_energy =
          dem_parameters.lagrangian_physical_properties.surface_energy_particle
            .at(i);
        const double particle_hamaker_constant =
          dem_parameters.lagrangian_physical_properties
            .hamaker_constant_particle.at(i);

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

        const double log_coeff_restitution =
          log(this->effective_coefficient_of_restitution[i]);
        this->model_parameter_beta[i] =
          log_coeff_restitution /
          sqrt((log_coeff_restitution * log_coeff_restitution) + 9.8696);

        this->effective_coefficient_of_friction[i] =
          2 * particle_friction_coefficient * wall_friction_coefficient /
          (particle_friction_coefficient + wall_friction_coefficient + DBL_MIN);

        this->effective_coefficient_of_rolling_friction[i] =
          2 * particle_rolling_friction_coefficient *
          wall_rolling_friction_coefficient /
          (particle_rolling_friction_coefficient +
           wall_rolling_friction_coefficient + DBL_MIN);

        this->effective_surface_energy[i] =
          particle_surface_energy + wall_surface_energy -
          std::pow(std::sqrt(particle_surface_energy) -
                     std::sqrt(wall_surface_energy),
                   2);
        this->effective_hamaker_constant[i] =
          0.5 * (particle_hamaker_constant + wall_hamaker_constant);
      }
  }


  /**
   * @brief This function uses the effective surface energy and effective hamaker
   * constant to compute the cut off distance for the force calculation.
   *
   * @param[in] dem_parameters DEM parameters declared in the .prm file
   *
   * @return Cut off distance for the force calculation.
   */
  double
  set_dmt_cut_off_distance()
  {
    // We are searching for the large delta_o, thus the smallest effective
    // surface energy and largest Hamaker constant.
    double max_effective_hamaker_constant =
      *(std::max_element(this->effective_hamaker_constant.begin(),
                         this->effective_hamaker_constant.end()));
    double min_effective_surface_energy =
      *(std::min_element(this->effective_surface_energy.begin(),
                         this->effective_surface_energy.end()));

    const double delta_0 =
      -std::sqrt(max_effective_hamaker_constant /
                 (12. * M_PI * min_effective_surface_energy));

    return delta_0 / std::sqrt(dmt_cut_off_threshold);
  };

  const double dmt_cut_off_threshold;
};

#endif
