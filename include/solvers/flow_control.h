// SPDX-FileCopyrightText: Copyright (c) 2020-2021, 2023-2024, 2026 The Lethe
// Authors SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR
// LGPL-2.1-or-later

#ifndef lethe_flow_control_h
#define lethe_flow_control_h

#include <core/parameters.h>

#include <deal.II/base/tensor.h>

#include <string>

using namespace dealii;

/**
 * @brief Dynamic flow controller that drives the average velocity of the flow
 * towards a target value.
 *
 * At each time step, the controller calculates a force coefficient beta,
 * applied as a source term in the momentum equation along the flow direction,
 * from the average velocity of the previous time steps. The main controller is
 * based on Wang, 2009, "A non-body conformal grid method for simulations of
 * laminar and turbulent flows with a compressible large eddy simulation
 * solver".
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class FlowControl
{
public:
  /**
   * @brief Constructor of the flow controller.
   *
   * @param[in] flow_control Parameters of the dynamic flow control. When the
   * flow control is enabled, its flow direction must be smaller than @p dim.
   */
  explicit FlowControl(const Parameters::DynamicFlowControl &flow_control);

  /**
   * @brief Calculate the beta force coefficient applied during the current
   * time step.
   *
   * The coefficient adjusts the flow rate towards the target average
   * velocity. Once the target is reached, the controller calculates the beta
   * force that keeps the flow rate constant. The new coefficient is stored and
   * can be obtained with get_beta().
   *
   * @param[in] average_velocity_n Average velocity of the last time step.
   *
   * @param[in] dt Current time step.
   *
   * @param[in] step_number Current time step number.
   */
  void
  calculate_beta(const double       &average_velocity_n,
                 const double       &dt,
                 const unsigned int &step_number);

  /**
   * @brief Get the beta force of the current time step.
   *
   * @return Beta force tensor, whose only non-zero component is along the flow
   * direction.
   */
  Tensor<1, dim>
  get_beta() const
  {
    return beta;
  }

  /**
   * @brief Get the beta force of the current time step scaled for the
   * particles.
   *
   * The particles are subjected to the same pressure gradient as the fluid,
   * which gives beta_p * rho_p = beta_f * rho_f, in units of a pressure
   * gradient M/(T²L²).
   *
   * @param[in] fluid_density Density of the fluid.
   *
   * @param[in] particle_density Density of the particles.
   *
   * @return Beta force for the particles. It is a tensor of dimension 3 since
   * it is added to the gravity of the DEM parameters, which is always a tensor
   * of dimension 3.
   */
  inline Tensor<1, 3>
  get_beta_particles(const double &fluid_density,
                     const double &particle_density) const
  {
    Tensor<1, 3> beta_particle;
    beta_particle[0] = beta[0];
    beta_particle[1] = beta[1];

    if constexpr (dim == 3)
      beta_particle[2] = beta[2];

    return beta_particle * fluid_density / particle_density;
  }

  /**
   * @brief Save the state of the flow controller to a checkpoint file.
   *
   * @param[in] prefix Prefix of the checkpoint files. The state is written in
   * the file named @p prefix followed by the extension .flowcontrol.
   */
  void
  save(const std::string &prefix) const;

  /**
   * @brief Read the state of the flow controller from a checkpoint file.
   *
   * @param[in] prefix Prefix of the checkpoint files. The state is read from
   * the file named @p prefix followed by the extension .flowcontrol.
   */
  void
  read(const std::string &prefix);

private:
  /**
   * @brief Calculate a beta force from the last average velocity only. It is
   * used for the first and second time steps, when no initial beta is given or
   * when the history of the main controller is not available yet.
   *
   * @param[in] average_velocity_n Average velocity of the last time step.
   *
   * @param[in] dt Current time step.
   *
   * @return Beta force calculated by the proportional controller.
   */
  double
  proportional_flow_controller(const double &average_velocity_n,
                               const double &dt) const;

  /**
   * @brief Calculate a beta force from the previous beta force and the
   * average velocities of the last two time steps, following Wang (2009).
   *
   * Once the target average velocity has been reached, a beta force opposed to
   * the target flow is set to zero.
   *
   * @param[in] average_velocity_n Average velocity of the last time step.
   *
   * @param[in] dt Current time step.
   *
   * @return Beta force calculated by the main controller.
   */
  double
  main_flow_controller(const double &average_velocity_n,
                       const double &dt) const;

  /// Below this magnitude, a beta force is considered to be zero.
  static constexpr double zero_beta_tolerance = 1e-6;

  /// Initial ratio of the upper limit of the force-free velocity band to the
  /// target average velocity.
  static constexpr double initial_force_free_velocity_band_factor = 1.01;

  /// Smallest ratio of the upper limit of the force-free velocity band to the
  /// target average velocity. The band stops shrinking below this value.
  static constexpr double minimum_force_free_velocity_band_factor = 1.00125;

  /// Beta force of the current time step. Only its component along the flow
  /// direction is non-zero.
  Tensor<1, dim> beta;

  /// Target average velocity.
  double target_average_velocity;

  /// Initial beta force given by the user. It is only used when its magnitude
  /// is larger than zero_beta_tolerance.
  double initial_beta;

  /// Relaxation coefficient of the controllers.
  double alpha;

  /// Relative change of the beta force below which the previous beta force is
  /// kept.
  double beta_threshold;

  /// Beta force applied during the previous time step.
  double previous_beta = 0.0;

  /// Average velocity of the previous time step.
  double previous_average_velocity = 0.0;

  /// Direction of the flow (0 for x, 1 for y and 2 for z).
  unsigned int flow_direction;

  /// True while the target average velocity has not been reached once without
  /// any force applied. In that state, an average velocity slightly above the
  /// target (within the force-free velocity band) is allowed to decrease
  /// without any force, through the pressure drop only.
  bool force_free_deceleration_is_allowed = true;

  /// Ratio of the upper limit of the force-free velocity band to the target
  /// average velocity. Its distance to 1 is halved when the pressure drop alone
  /// is too slow to bring the average velocity back to the target.
  double force_free_velocity_band_factor =
    initial_force_free_velocity_band_factor;
};

#endif
