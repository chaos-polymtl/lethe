// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ib_particle_h
#define lethe_ib_particle_h

#include <core/dem_properties.h>
#include <core/shape.h>

#include <deal.II/base/parsed_function.h>
#include <deal.II/base/point.h>

#include <vector>

using namespace dealii;

/**
 * @brief Define the behavior of a specific particle for the sharp interface IB solver.
 *
 * Each particle defined will have these values used in the solver.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 */
template <int dim>
class IBParticle
{
public:
  /**
   * @brief Initialize all particles.
   */
  void
  initialize_all();

  /**
   * @brief Initialize the value of the last state of the particle.
   */
  void
  initialize_previous_solution();

  /**
   * @brief Return the names of the DEM properties, for visualization.
   */
  std::vector<std::pair<std::string, int>>
  get_properties_name();

  /**
   * @brief Return the value of the DEM properties, for visualization.
   */
  std::vector<double>
  get_properties();

  /**
   * @brief Return the number of DEM properties, for visualization.
   */
  inline unsigned int
  get_number_properties()
  {
    return DEM::CFDDEMProperties::PropertiesIndex::n_properties;
  }

  /**
   * @brief Compute the evaluation of the signed distance function of this shape.
   *
   * Most level set functions used here come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   *
   * @param[in] p Point at which the evaluation is performed.
   *
   * @param[in] cell_guess Guess of the cell containing the evaluation point,
   * which is useful to reduce computation time.
   *
   * @return Signed distance function at evaluation point.
   */
  inline double
  get_levelset(const Point<dim>                                     &p,
               const typename DoFHandler<dim>::active_cell_iterator &cell_guess)
  {
    return shape->value_with_cell_guess(p, cell_guess);
  }

  /**
   * @brief See overloaded function.
   */
  inline double
  get_levelset(const Point<dim> &p)
  {
    return shape->value(p);
  }

  /**
   * @brief Find the closest surface point on the shape.
   *
   * This point has the minimal absolute distance from the given point p.
   *
   * @param[in] p Point at which the evaluation is performed.
   *
   * @param[out] closest_point Reference to the closest point. This point
   * will be modified by the function.
   *
   * @param[in] cell_guess A guess of the cell containing the evaluation point,
   * which is useful to reduce computation time.
   */
  void
  closest_surface_point(
    const Point<dim>                                     &p,
    Point<dim>                                           &closest_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * @brief See overloaded function.
   */
  void
  closest_surface_point(const Point<dim> &p, Point<dim> &closest_point);

  /**
   * @brief Determine if the given point is inside the crown for which the limits
   * are defined by inner and outer radius factors. An effective radius is used
   * for non spherical shapes.
   *
   *  @param[in] evaluation_point Point at which the evaluation is performed.
   *
   *  @param[in] outer_radius Relative or absolute distance to check if the
   *  evaluation point is inside the outer limits.
   *
   *  @param[in] inside_radius Relative or absolute distance to check if the
   *  evaluation point is outside the inner limits.
   *
   *  @param[in] absolute_distance Boolean for whether the radii are defined as
   * absolute or relative to radius. If the value is false, the outer and inner
   * radii are multiplied by the effective radius to construct the absolute
   * distances to be used.
   *
   *  @param[in] cell_guess Guess of the cell containing the evaluation point,
   * which is useful to reduce computation time.
   *
   *  @return true if the evaluation point is inside the crown.
   */
  bool
  is_inside_crown(
    const Point<dim>                                     &evaluation_point,
    const double                                          outer_radius,
    const double                                          inside_radius,
    const bool                                            absolute_distance,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * @brief See overloaded function.
   */
  bool
  is_inside_crown(const Point<dim> &evaluation_point,
                  const double      outer_radius,
                  const double      inside_radius,
                  const bool        absolute_distance);

  /**
   * @brief Set the position and dependent members to the argument.
   *
   * @param[in] position New position to use.
   */
  inline void
  set_position(const Point<dim> position)
  {
    this->position = position;
    this->shape->set_position(this->position);
  }


  /**
   * @brief Clear the cache used to evaluate the value and vector defining
   * the signed distance function.
   */
  void
  clear_shape_cache();

  /**
   * @brief Set up a shape in accordance with the given type and arguments.
   *
   * @param[in] type Type of the shape to be initialized.
   *
   * @param[in] shape_arguments Arguments to be used for shape initialization.
   */
  void
  initialize_shape(const std::string         &type,
                   const std::vector<double> &shape_arguments);

  /**
   * @brief Set up a shape in accordance with the given type and raw text argument.
   *
   * The raw text argument is typically a file name.
   *
   * @param[in] type Type of shape to be initialized.
   *
   * @param[in] raw_arguments Untreated arguments of the shape.
   */
  void
  initialize_shape(const std::string &type, const std::string &raw_arguments);

  /**
   * @brief Set the position and dependent members to the position for
   * one component only.
   *
   * @param[in] position_component Component of the new position to set the
   * particle at.
   *
   * @param[in] component Component index for which the position will be
   * updated.
   */
  inline void
  set_position(const double       position_component,
               const unsigned int component = 0)
  {
    this->position[component] = position_component;
    this->shape->set_position(this->position);
  }

  /**
   * @brief Set the orientation and dependent members to the given orientation.
   *
   * @param[in] orientation New orientation to set the particle at.
   */
  void
  set_orientation(const Tensor<1, 3> &orientation);

  /**
   * @brief Set the layer thickening value (positive or negative) of the shape.
   *
   * @param[in] layer_thickening Thickness to be artificially added to the
   * particle. A negative value will decrease the particle's thickness by
   * subtracting a layer of specified width.
   */
  void
  set_layer_thickening(const double layer_thickening)
  {
    shape->set_layer_thickening(layer_thickening);
  }

  /**
   * @brief Set the proper DOF handler, then compute/update the map of cells
   * and their likely non-null nodes.
   *
   * @param[in,out] updated_dof_handler Reference to the new dof_handler.
   *
   * @param[in] mesh_based_precalculations Define if mesh-based precalculations
   * are used (if type=RBF).
   */
  void
  update_precalculations(DoFHandler<dim> &updated_dof_handler,
                         const bool       mesh_based_precalculations);

  /**
   * @brief Update precalculations if needed, then compute and remove superfluous data.
   *
   * @param[in,out] updated_dof_handler Reference to the new dof_handler.
   *
   * @param[in] mesh_based_precalculations Define if mesh-based precalculations
   * are used (if type=RBF).
   */
  void
  remove_superfluous_data(DoFHandler<dim> &updated_dof_handler,
                          const bool       mesh_based_precalculations);

  /**
   * @brief Load data from the files (for file-based Shapes).
   */
  void
  load_data_from_file();

  /**
   * @brief Unique identification number.
   */
  unsigned int particle_id;

  /**
   * @brief Geometrical information for IB application.
   */
  std::shared_ptr<Shape<dim>> shape;

  /**
   * @brief Effective radius, which is the actual radius for spheres.
   */
  double radius;

  /**
   * @brief Young's modulus used for contact force calculations.
   */
  double youngs_modulus;

  /**
   * @brief Restitution coefficient used for contact force calculations.
   */
  double restitution_coefficient;

  /**
   * @brief Poisson's ratio used for contact force calculations.
   */
  double poisson_ratio;

  /**
   * @brief Friction coefficient used for contact force calculations.
   */
  double friction_coefficient;

  /**
   * @brief Rolling friction coefficient used for contact force calculations.
   */
  double rolling_friction_coefficient;

  /**
   * @brief Particle's mass.
   */
  double mass;

  /**
   * @brief Particle's volume
   */
  double volume;

  /**
   * @brief Surface energy
   */
  double surface_energy;

  /**
   * @brief Inertia matrix in the particle frame of reference.
   */
  Tensor<2, 3> inertia;

  /**
   * @brief Position, which corresponds to the centroid for most implemented shapes.
   */
  Point<dim> position;

  /**
   * @brief Positions at the last n time steps. n depends on the time integration scheme.
   */
  std::vector<Point<dim>> previous_positions;

  /**
   * @brief Total fluid force applied on the particle.
   */
  Tensor<1, 3> fluid_forces;

  /**
   * @brief Viscous contribution to the fluid force applied on the particle.
   */
  Tensor<1, 3> fluid_viscous_forces;

  /**
   * @brief Pressure contribution to the fluid force applied on the particle.
   */
  Tensor<1, 3> fluid_pressure_forces;

  /**
   * @brief Previous time step's total fluid force applied on the particle.
   */
  Tensor<1, 3> previous_fluid_forces;

  /**
   * @brief Previous time step's viscous contribution to the fluid force
   * applied on the particle.
   */
  Tensor<1, 3> previous_fluid_viscous_forces;

  /**
   * @brief Previous time step's pressure contribution to the fluid force.
   */
  Tensor<1, 3> previous_fluid_pressure_forces;

  /**
   * @brief Fluid torque applied on the particle.
   */
  Tensor<1, 3> fluid_torque;

  /**
   * @brief Previous time step's fluid torque applied on the particle.
   */
  Tensor<1, 3> previous_fluid_torque;
  /**
   * @brief Translational velocity.
   */
  Tensor<1, 3> velocity;

  /**
   * @brief Previous time step's translational velocity.
   */
  std::vector<Tensor<1, 3>> previous_velocity;

  /**
   * @brief Last non-linear iteration's velocity for CFD-DEM coupling.
   */
  Tensor<1, 3> velocity_iter;

  /**
   * @brief Last non-linear correction (without relaxation) to the velocity.
   */
  Tensor<1, 3> previous_velocity_residual;

  /**
   * @brief Orientation, with 0 angular position by default on every axis.
   */
  Tensor<1, 3> orientation;

  /**
   * @brief Previous time step's orientation.
   */
  std::vector<Tensor<1, 3>> previous_orientation;

  /**
   * @brief Angular velocity.
   */
  Tensor<1, 3> omega;

  /**
   * @brief Previous time step's angular velocity.
   */
  std::vector<Tensor<1, 3>> previous_omega;

  /**
   * @brief Last non-linear iteration's angular velocity for CFD-DEM coupling.
   */
  Tensor<1, 3> omega_iter;

  /**
   * @brief Last non-linear correction (without relaxation) to the angular velocity.
   */
  Tensor<1, 3> previous_omega_residual;

  /**
   * @brief Total impulsion felt by the particle during the current time step.
   */
  Tensor<1, 3> impulsion;

  /**
   * @brief Contact-generated impulsion felt by the particle.
   */
  Tensor<1, 3> contact_impulsion;

  /**
   * @brief Last non-linear iteration's total impulsion on the particle.
   */
  Tensor<1, 3> impulsion_iter;

  /**
   * @brief Total angular impulsion felt by the particle during the current time step.
   */
  Tensor<1, 3> omega_impulsion;

  /**
   * @brief Contact-generated angular impulsion felt by the particle.
   */
  Tensor<1, 3> omega_contact_impulsion;

  /**
   * @brief Last non-linear iteration's total angular impulsion on the particle.
   */
  Tensor<1, 3> omega_impulsion_iter;

  /**
   * @brief Function that determines the velocity.
   *
   * If the particle dynamics is solved, this function determines the
   * value only at the initial state.
   */
  std::shared_ptr<Functions::ParsedFunction<dim>> f_velocity;

  /**
   * @brief Function that determines the position.
   *
   * If the particle dynamics is solved, this function determines the
   * value only at the initial state.
   */
  std::shared_ptr<Functions::ParsedFunction<dim>> f_position;

  /**
   * @brief Function that determines the angular velocity.
   *
   * If the particle dynamics is solved, this function determines the
   * value only at the initial state.
   */
  std::shared_ptr<Functions::ParsedFunction<dim>> f_omega;

  /**
   * @brief Function that determines the orientation.
   *
   * If the particle dynamics is solved, this function determines the
   * value only at the initial state.
   */
  std::shared_ptr<Functions::ParsedFunction<dim>> f_orientation;

  /**
   * @brief Indicate if the particle dynamics are solved.
   *
   * If not, functions are used to define position, velocity, orientation and
   * angular velocity.
   */
  bool integrate_motion;

  /**
   * @brief Indicate if mesh-based precalculations are performed.
   *
   * For RBF shapes with nodes outside the background mesh, slight deformations
   * can happen near the boundary when mesh-based precalculations are enabled
   * (due to the RBF nodes partitioning algorithm).
   */
  bool mesh_based_precalculations;

  /**
   * @brief Current residual of the particle velocity.
   */
  double residual_velocity;

  /**
   * @brief Current residual of the particle angular velocity.
   */
  double residual_omega;

  /**
   * @brief Previous iteration's relaxation parameter used for the translational velocity
   * iteration.
   */
  double previous_local_alpha_velocity;

  /**
   * @brief Previous iteration's relaxation parameter used for the angular velocity
   * iteration.
   */
  double previous_local_alpha_omega;

  /**
   * @brief Position of the pressure reference point in the particle's
   * referential. This point is used to define a constant on the pressure.
   */
  Point<dim> pressure_location;

  /**
   * @brief Rotation matrix in the global space.
   */
  Tensor<2, 3> rotation_matrix;
};

#endif
