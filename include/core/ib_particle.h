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

#include <core/shape.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/point.h>

#include <deal.II/physics/transformations.h>

#include <vector>

#ifndef lethe_ib_particle_h
#  define lethe_ib_particle_h

using namespace dealii;

template <int dim>
class IBParticle
{
public:
  // Function to initialise the value associated with each particle.

  // Properties Index is necessary to interface with some DEM functions and
  // simplify the creation of vtu output file.
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
    m            = 11,
    tx           = 12,
    ty           = 13,
    tz           = 14,
    n_properties = 15,
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

  /**
   * @brief
   * Return the evaluation of the signed distance function of this solid at the
   * given point p
   * Most levelset functions come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   */
  double
  get_levelset(const Point<dim> &p);

  /**
   * @brief
   * Return the point on the surface of the solid which has the minimal distance
   * from the given point p
   */
  void
  closest_surface_point(const Point<dim> &p, Point<dim> &closest_point);

  /**
   * @brief
   * Return true if the given point is inside the crown whose limits are defined
   * by inner and outer radius factors. An effective radius is used for
   * non spherical shapes.
   */
  bool
  is_inside_crown(const Point<dim> &evaluation_pt,
                  const double      outer_radius,
                  const double      inside_radius);

  /**
   * @brief
   * Sets up a default shape
   */
  void
  initialize_shape();

  /**
   * @brief
   * Sets up a shape in accordance with the given type and arguments
   */
  void
  initialize_shape(const std::string type, Tensor<1, 3> solid_arguments);

  /**
   * @brief
   * Sets the position of the particle and dependant members to the argument
   */
  void
  set_position(const Point<dim> position);

  /**
   * @brief
   * Increments one component of the position and updates the relevant members
   */
  void
  move(const double position_update, const unsigned int component = 0);

  /**
   * @brief
   * Sets the orientation of the particle and dependant members to the argument
   */
  void
  set_orientation(const Tensor<1, 3> orientation);

  /**
   * @brief
   * Sets the center of rotation offset of the particle and dependant members to
   * the argument
   */
  void
  set_center_of_rotation_offset(const Point<dim> cor_offset);


  // This class defines values related to a particle used in the sharp interface
  // IB. Each particle defined will have these value used in the solver.

  // The unique identification number of the particle.
  unsigned int particle_id;

  // The geometrical information regarding the particle
  std::shared_ptr<Shape<dim>> shape;
  // The function from which the solid arguments are parsed
  std::shared_ptr<Functions::ParsedFunction<3>> f_solid_arguments;
  // The arguments defining the shape. Their meaning depend of the shape type.
  Tensor<1, 3> solid_arguments;

  // The particle effective radius. It is the actual radius for spheres.
  double radius;
  // The particle Young's modulus. Used in case of contact.
  double youngs_modulus;
  // The particle restitution coefficient. Used in case of contact.
  double restitution_coefficient;
  // The particle Poisson's ratio. Used in case of contact.
  double poisson_ratio;
  // The particle friction coefficient. Used in case of contact.
  double friction_coefficient;
  // The particle rolling friction coefficient. Used in case of contact.
  double rolling_friction_coefficient;
  // The particle mass.
  double mass;
  // The particle rotational inertia. It is uniform in all directions.
  Tensor<2, 3> inertia;
  // The particle position.
  Point<dim> position;
  // The vector of particle positions at the end of the last n time steps.
  std::vector<Point<dim>> previous_positions;
  // The center of mass offset of the solid.
  Point<dim> center_of_mass_offset;
  // The vector of solid center of mass offset at the end of the last n time
  // steps.
  std::vector<Point<dim>> previous_center_of_mass_offset;
  // The fluid force applied on the particle.
  Tensor<1, 3> fluid_forces;
  // The fluid force applied on the particle at the end of the last time step.
  Tensor<1, 3> previous_fluid_forces;
  // The fluid torque applied on the particle.
  Tensor<1, 3> fluid_torque;
  // The fluid torque is applied on the particle at the end of the last time
  // step.
  Tensor<1, 3> previous_fluid_torque;
  // The translational velocity
  Tensor<1, 3> velocity;
  // The vector of particle translational velocity at the end of the last n time
  // steps.
  std::vector<Tensor<1, 3>> previous_velocity;
  // The last non-linear iteration of the velocity vector.
  Tensor<1, 3> velocity_iter;
  // The last correction vector of the velocity value without any relaxation.
  Tensor<1, 3> previous_d_velocity;

  // Angular velocity
  // By default, the angular position is always 0 on every axis.
  Tensor<1, 3> orientation;
  //  The vector of the particle angular position of the last n time steps.
  std::vector<Tensor<1, 3>> previous_orientation;
  // Angular velocity
  Tensor<1, 3> omega;
  // The angular velocity at the end of the last time step.
  std::vector<Tensor<1, 3>> previous_omega;
  // The last iteration of the omega vector.
  Tensor<1, 3> omega_iter;
  // The last correction vector of the velocity value without any relaxation.
  Tensor<1, 3> previous_d_omega;

  // The total impulsion that the particle feels during the current time step.
  Tensor<1, 3> impulsion;
  // The impulsion from contact that the particle feels.
  Tensor<1, 3> contact_impulsion;
  // The last non-linear iteration of the total impulsion felt by the particle
  Tensor<1, 3> impulsion_iter;
  // The total angular impulsion that the particle as felt during the current
  // time step.
  Tensor<1, 3> omega_impulsion;
  // The angular impulsion from contact that the particle feels.
  Tensor<1, 3> omega_contact_impulsion;
  // The last non-linear iteration of the total angular impulsion felt by the
  // particle.
  Tensor<1, 3> omega_impulsion_iter;

  // The function from which the initial particle velocity is determined.
  // If the dynamic is not resolved, this function determines the velocity at
  // every time step.
  std::shared_ptr<Functions::ParsedFunction<dim>> f_velocity;
  // The function from which the initial particle position is determined.
  // If the dynamic is not resolved, this function determines the position at
  // every time step.
  std::shared_ptr<Functions::ParsedFunction<dim>> f_position;
  // The function from which the particle center of mass offset is
  // determined.
  std::shared_ptr<Functions::ParsedFunction<dim>> f_center_of_mass_offset;
  // The function from which the particle initial angular velocity is
  // determined. If the dynamic is not resolved, this function determines the
  // angular velocity at every time step.
  std::shared_ptr<Functions::ParsedFunction<dim>> f_omega;
  // The function from which the initial particle angular position is
  // determined. If the dynamic is not resolved, this function determines the
  // angular position at every time step.
  std::shared_ptr<Functions::ParsedFunction<dim>> f_orientation;

  // Allow the definition of a local relaxation parameter for each particle in
  // the integration process.

  // Current residual of the particle velocity.
  double residual_velocity;
  // Current residual of the particle angular velocity.
  double residual_omega;
  // Last relaxation parameter used for this particle translational velocity
  // iteration.
  double previous_local_alpha_velocity;
  // Last relaxation parameter used for this particle angular velocity
  // iteration.
  double previous_local_alpha_omega;

  // Location of the pressure reference point relative to the center of the
  // particle. This point is used to define a constant on the pressure.
  Point<dim> pressure_location;
};

#endif
