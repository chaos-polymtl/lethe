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
 */

#include <core/dem_properties.h>
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
  /**
   * @brief
   * initialised the particle
   *
   */
  void
  initialize_all();

  /**
   * @brief
   * initialize the value of the last state of the particle
   */
  void
  initialize_previous_solution();

  /**
   * @brief
   * Return the names of properties of the IB_particle for visualisation.
   */
  std::vector<std::pair<std::string, int>>
  get_properties_name();

  /**
   * @brief
   * Return the value of the properties of the particle for visualisation.
   */
  std::vector<double>
  get_properties();

  /**
   * @brief
   * Return the number of properties of the particle for visualisation.
   */
  inline unsigned int
  get_number_properties()
  {
    return DEM::PropertiesIndex::n_properties;
  }

  /**
   * @brief
   * Returns the evaluation of the signed distance function of this shape
   * Most levelset functions come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   *
   * @param p The point at which the evaluation is performed
   * @param cell_guess A guess of the cell containing the evaluation point, which
   * is useful to reduce computation time
   */
  inline double
  get_levelset(const Point<dim>                                     &p,
               const typename DoFHandler<dim>::active_cell_iterator &cell_guess)
  {
    return shape->value_with_cell_guess(p, cell_guess);
  }

  /**
   * See overloaded function
   */
  inline double
  get_levelset(const Point<dim> &p)
  {
    return shape->value(p);
  }

  /**
   * @brief
   * Sets the closest_point parameter to be the point on the surface of the
   * shape which has the minimal distance from the given point p
   *
   * @param p The point at which the evaluation is performed
   * @param closest_point The reference to the closest point. This point will be modified by the function.
   * @param cell_guess A guess of the cell containing the evaluation point, which
   * is useful to reduce computation time
   */
  void
  closest_surface_point(
    const Point<dim>                                     &p,
    Point<dim>                                           &closest_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * See overloaded function
   */
  void
  closest_surface_point(const Point<dim> &p, Point<dim> &closest_point);

  /**
   * @brief
   * Returns true if the given point is inside the crown for which the limits
   * are defined by inner and outer radius factors. An effective radius is used
   * for non spherical shapes.
   *
   *  @param evaluation_point The point at which the evaluation is performed
   *  @param outer_radius The factor to be multiplied by the effective radius to check if the evaluation point is inside the outer limits
   *  @param inside_radius The factor to be multiplied by the effective radius to check if the evaluation point is outside the inner limits
   *  @param absolute_distance Whether the distance is defined as absolute or relative to radius
   *  @param cell_guess A guess of the cell containing the evaluation point, which
   *  is useful to reduce computation time
   */
  bool
  is_inside_crown(
    const Point<dim>                                     &evaluation_point,
    const double                                          outer_radius,
    const double                                          inside_radius,
    const bool                                            absolute_distance,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * See overloaded function
   */
  bool
  is_inside_crown(const Point<dim> &evaluation_point,
                  const double      outer_radius,
                  const double      inside_radius,
                  const bool        absolute_distance);

  /**
   * @brief
   * Sets the position of the particle and dependent members to the argument
   *
   * @param position The new position to set the particle at
   */
  inline void
  set_position(const Point<dim> position)
  {
    this->position = position;
    this->shape->set_position(this->position);
  }


  /**
   * @brief Clear the cache used to evaluate the value and vector defining the sign distance function of the shape
   */
  void
  clear_shape_cache();

  /**
   * @brief
   * Sets up a shape in accordance with the given type and arguments
   *
   * @param type The type of shape to be initialized: sphere, cone, ellipsoid,
   * rectangle, death star, cut hollow sphere, torus, or rbf
   * @param shape_arguments The dimensions to be used for shape initialization
   */
  void
  initialize_shape(const std::string         type,
                   const std::vector<double> shape_arguments);
  /**
   * @brief
   * Sets up a shape in accordance with the given type and file name
   *
   * @param type The type of shape to be initialized: sphere, cone, ellipsoid,
   * rectangle, death star, cut hollow sphere, torus, or rbf
   * @param raw_arguments The untreated arguments of the shape
   */
  void
  initialize_shape(const std::string type, const std::string raw_arguments);

  /**
   * @brief
   * Sets the position of the particle and dependent members to the argument for
   * one component
   *
   * @param position_component The component of the new position to set the particle at
   * @param component The component index for which the position will be updated
   */
  inline void
  set_position(const double       position_component,
               const unsigned int component = 0)
  {
    this->position[component] = position_component;
    this->shape->set_position(this->position);
  }

  /**
   * @brief
   * Sets the orientation of the particle and dependent members to the argument
   *
   * @param orientation The new orientation to set the particle at
   */
  void
  set_orientation(const Tensor<1, 3> orientation);

  /**
   * @brief
   * Sets the layer thickening value (positive or negative) of the particle's
   * shape
   *
   * @param layer_thickening Thickness to be artificially added to the particle.
   * A negative value will decrease the particle's thickness by subtracting a
   * layer of specified width.
   */
  void
  set_layer_thickening(const double layer_thickening)
  {
    shape->set_layer_thickening(layer_thickening);
  }

  /**
   * @brief Sets the proper dof handler, then computes/updates the map of cells
   * and their likely non-null nodes
   * @param updated_dof_handler the reference to the new dof_handler
   * @param mesh_based_precalculations mesh-based precalculations that can lead to slight shape misrepresentation (if type=RBF)
   */
  void
  update_precalculations(DoFHandler<dim> &updated_dof_handler,
                         const bool       mesh_based_precalculations);

  /**
   * @brief Updates precalculations if needed, then computes and removes superfluous data
   * @param updated_dof_handler the reference to the new dof_handler
   * @param mesh_based_precalculations mesh-based precalculations that can lead to slight shape misrepresentation (if type=RBF)
   */
  void
  remove_superfluous_data(DoFHandler<dim> &updated_dof_handler,
                          const bool       mesh_based_precalculations);

  /**
   * @brief Loads data from the files for file-based Shapes (RBF at the moment)
   */
  void
  load_data_from_file();

  void
  compute_local_inertia(Tensor<2, 3> global_inertia);




  // This class defines values related to a particle used in the sharp interface
  // IB. Each particle defined will have these value used in the solver.

  // The unique identification number of the particle.
  unsigned int particle_id;

  // The geometrical information regarding the particle
  std::shared_ptr<Shape<dim>> shape;

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
  // The surface energy of the particle
  double surface_energy;
  // The particle tensor of inertia in the global reference frame
  Tensor<2, 3> inertia;
  // The particle position.
  Point<dim> position;
  // The vector of particle positions at the end of the last n time steps.
  std::vector<Point<dim>> previous_positions;
  // The fluid force applied on the particle.
  Tensor<1, 3> fluid_forces;
  Tensor<1, 3> fluid_viscous_forces;
  Tensor<1, 3> fluid_pressure_forces;
  // The fluid force applied on the particle at the end of the last time step.
  Tensor<1, 3> previous_fluid_forces;
  Tensor<1, 3> previous_fluid_viscous_forces;
  Tensor<1, 3> previous_fluid_pressure_forces;
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
  Tensor<1, 3> previous_velocity_residual;

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
  Tensor<1, 3> previous_omega_residual;

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

  // Bool that indicates if the motion of this particle must be integrated. If
  // it is false, the position and velocity are defined by the function.
  bool integrate_motion;
  // Bool that indicates if mesh-based precalculations should be performed.
  // For RBF shapes with nodes outside the background mesh, slight deformations
  // can happen near the boundary when mesh-based precalculations is used (due
  // to the RBF nodes partitioning algorithm).
  bool mesh_based_precalculations;

  // Location of the pressure reference point relative to the center of the
  // particle. This point is used to define a constant on the pressure.
  Point<dim> pressure_location;

  // Location of the center of mass relative to the frame of reference of the
  // particule.
  Point<dim> center_of_mass_location;

  // Rotation matrix of the particle in the global space
  Tensor<2, 3> rotation_matrix;
};

#endif
