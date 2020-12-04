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

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>

#include <core/parameters.h>

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
   * value through the previous flow rate. Once the flow rate is reached, the
   * algorithm calculates a new beta to keep a constant flow rate.
   *
   * @param flow_rate. The last step flow rate
   *
   * @param dt. The current time step
   *
   * @param step_number. The current step
   */
  void
  calculate_beta(const std::pair<double, double> &flow_rate,
                 const double &                   dt,
                 const unsigned int &             step_number);

  /**
   * @brief get_beta. This function gives the beta coefficient of the
   * step time
   */
  Tensor<1, dim>
  get_beta()
  {
    return beta;
  }

  void
  save(std::string prefix);

  void
  read(std::string prefix);

private:
  // The coefficients are stored in the following fashion :
  // 0 - flow rate intended, n - n, 1n - n-1, n1 - n+1
  Tensor<1, dim> beta;
  double         beta_0;
  double         beta_n;
  double         beta_n1;
  double         flow_rate_0;
  double         flow_rate_1n;
  double         flow_rate_n;
  double         area;
  unsigned int   flow_direction;

  // Variables used to improve convergence
  bool   no_force;
  double threshold_factor;
};

#endif
