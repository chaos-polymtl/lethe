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
 * Author: Lucka Barbeau, Bruno Blais, Polytechnique Montreal, 2019 -
 */

#include <deal.II/base/parsed_function.h>
#include <deal.II/base/point.h>

#ifndef lethe_ib_particle_h
#  define lethe_ib_particle_h

using namespace dealii;

template <int dim>
class IBParticle
{
public:
  // Function to initialise the value associated with each particle.
  enum PropertiesIndex : int
  {
    id           = 0,
    dp           = 1,
    vx           = 2,
    vy           = 3,
    vz           = 4,
    fx           = 5,
    fy           = 6,
    fz           = 7,
    ox           = 8,
    oy           = 9,
    oz           = 10,
    tx           = 11,
    ty           = 12,
    tz           = 13,
    n_properties = 14,
  };
  /**
   * @brief
   * initialised the particle
   *
   */
  void
  initialise_all();
  /**
   * @brief
   * initialize the value of the last state of the particle
   *
   */
  void
  initialise_last();
  /**
   * @brief
   * Return the names of properties of the IB_particle for visualisation.
   *
   */
  std::vector<std::pair<std::string, int>>
  get_properties_name();
  /**
   * @brief
   * Return the value of the properties of the particle for visualisation.
   *
   */
  std::vector<double>
  get_properties();
  /**
   * @brief
   * Return the number of properties of the particle for visualisation.
   *
   */
  unsigned int
  get_number_properties();

  // This class defines values related to a particle used in the sharp interface
  // IB. Each particle defined will have these value used in the solver.
  Point<dim>     position;
  Point<dim>     last_position;
  Tensor<1, dim> forces;
  Tensor<1, dim> last_forces;
  unsigned int   particle_id;

  Tensor<1, 3> torques;
  Tensor<1, 3> last_torques;
  double       mass;
  Tensor<2, 3> inertia;
  // Translational velocity
  Tensor<1, dim> velocity;
  // Store the last velocity of the fix point iteration.
  Tensor<1, dim> last_velocity;
  Tensor<1, dim> velocity_iter;
  // Angular velocity

  // By default the angular position is always 0 on every axis.
  Tensor<1, 3> angular_position;
  // Store the last angular position of the particle for integration.
  Tensor<1, 3> last_angular_position;

  // Angular velocity
  Tensor<1, 3> omega;
  // Store the last angular velocity of the particle for integration.
  Tensor<1, 3> last_omega;
  // Store the last  angular velocity of the of the particle for the fix point
  // iteration.
  Tensor<1, 3> omega_iter;




  Tensor<1, dim> impulsion;
  Tensor<1, dim> impulsion_iter;

  std::shared_ptr<Functions::ParsedFunction<dim>> f_velocity;
  std::shared_ptr<Functions::ParsedFunction<dim>> f_position;
  std::shared_ptr<Functions::ParsedFunction<dim>> f_omega;


  Tensor<1, 3> omega_impulsion;
  Tensor<1, 3> omega_impulsion_iter;

  // Allow the definition of a local relaxation parameter for each particle in
  // the integration process.
  double local_alpha_torque;
  double local_alpha_force;

  double radius;

  double youngs_modulus;
  double restitution_coefficient;
  double friction_coefficient;
  double poisson_ratio;
  double rolling_friction_coefficient;


  // Pressure imposition location
  Point<dim> pressure_location;
};

#endif
