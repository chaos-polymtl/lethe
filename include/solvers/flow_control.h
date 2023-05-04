/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_flow_control_h
#define lethe_flow_control_h

// Lethe includes
#include <core/parameters.h>

#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>



using namespace dealii;

/**
 * @brief FlowControl. The FlowControl class allows to dynamically
 * control the flow with a beta coefficient calculated at each step time.
 */

template <int dim>
class FlowControl
{
public:
  FlowControl(const Parameters::DynamicFlowControl &flow_control);
  /**
   * @brief calculate_beta. This function calculates a beta coefficient which
   * applies a force to the flow in order to adjust the flow rate to a desired
   * average velocity through the previous flow rate. Once the flow rate is
   * reached, the algorithm calculates a new beta to keep a constant flow rate.
   *
   * @param average_velocity. The last step average_velocity
   *
   * @param dt. The current time step
   *
   * @param step_number. The current step
   */
  void
  calculate_beta(const double &      average_velocity,
                 const double &      dt,
                 const unsigned int &step_number);

  /**
   * @brief proportional_flow_controller. This function calculates a beta only
   * with the last step average velocity. It is used to initialize the beta
   * if there's no initial beta value at 1st time step and the 2nd time step.
   */
  inline double
  proportional_flow_controller(const double &average_velocity_n,
                               const double &dt)
  {
    return 0.5 * alpha * (average_velocity_0 - average_velocity_n) / dt;
  }

  /**
   * @brief main_flow_controller. This function calculates a beta coefficient
   * based on a formula described by Wang, 2009 in "A non-body conformal grid
   * method for simulations
   * of laminar and turbulent flows with a compressible
   * large eddy simulation solver"
   */
  inline double
  main_flow_controller(const double &average_velocity_n, const double &dt)
  {
    double beta_n1 = beta_n + alpha *
                                (average_velocity_0 - 2 * average_velocity_n +
                                 average_velocity_1n) /
                                dt;

    // If desired average velocity is reached, new beta only maintains the force
    // to keep the flow at the desired value. If calculated beta is
    // negative it is set to 0 to avoided +/- force.
    if (average_velocity_0 * beta_n1 < 0 && no_force == false)
      return 0.0;
    else
      return beta_n1;
  }

  /**
   * @brief get_beta. This function gives the beta force of the step time
   */
  Tensor<1, dim>
  get_beta()
  {
    return beta;
  }

  /**
   * @brief get_beta_particle. This function gives the beta force of the
   * CFD time step scaled for particles.
   */
  inline Tensor<1, 3>
  get_beta_particles(const double &fluid_density,
                     const double &particle_density)
  {
    // The beta force for particles must have a tensor of dim = 3 since the body
    // force g parameters is a tensor of dim 3 (those forces will add up)
    Tensor<1, 3> beta_particle;
    beta_particle[0] = beta[0];
    beta_particle[1] = beta[1];

    if constexpr (dim == 3)
      beta_particle[2] = beta[2];

    // From beta_p * rho_p = beta_f * rho_f [=] M/(T²L²) units of ΔP/L
    return beta_particle * fluid_density / particle_density;
  }

  void
  save(std::string prefix);

  void
  read(std::string prefix);

private:
  // The coefficients are stored in the following fashion :
  // 0 - target average velocity, n - n, 1n - n-1, n1 - n+1

  // Beta value as tensor for force application
  Tensor<1, dim> beta;

  // Target average velocity
  double average_velocity_0;

  // User defined parameters to improve convergence
  double beta_0;
  double alpha;
  double beta_threshold;

  // Average velocity and beta history
  double beta_n;
  double average_velocity_1n;

  // Flow direction
  unsigned int flow_direction;

  // Parameters of class to improve convergence
  bool   no_force;
  double threshold_factor;
};

#endif
