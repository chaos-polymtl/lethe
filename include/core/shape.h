// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_shape_h
#define lethe_shape_h

#include <core/tensors_and_points_dimension_manipulation.h>
#include <core/utilities.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#include <boost/range/adaptor/map.hpp>

#ifdef DEAL_II_WITH_OPENCASCADE
#  include <deal.II/opencascade/manifold_lib.h>
#  include <deal.II/opencascade/utilities.h>

#  include <BRepBuilderAPI_MakeVertex.hxx>
#  include <BRepClass3d_SolidClassifier.hxx>
#  include <BRepExtrema_DistShapeShape.hxx>
#  include <BRepGProp.hxx>
#  include <GProp_GProps.hxx>
#endif

#include <deal.II/base/function_signed_distance.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/manifold_lib.h>

#include <deal.II/physics/transformations.h>

#include <cfloat>
#include <memory>

using namespace dealii;

/**
 * @brief A base class used to represent geometrical entities. Its main uses
 * are to return signed distance and its gradient. It inherits
 * AutoDerivativeFunction so that it can evaluate gradients even when they are
 * not implemented analytically.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 */
template <int dim>
class Shape : public AutoDerivativeFunction<dim>
{
public:
  /**
   * @brief enum class that associate an integer index tp each type of shape
   */
  enum ShapeType : int
  {
    sphere,
    hyper_rectangle,
    ellipsoid,
    torus,
    cone,
    cylinder,
    cylindrical_tube,
    cylindrical_helix,
    cut_hollow_sphere,
    death_star,
    composite_shape,
    rbf_shape,
    opencascade_shape,
  } type;

  /**
   * @brief A general constructor for the Shapes
   *
   * @param radius The effective radius to be set for the shape. It's necessary for some calculations since not all shapes are spheres.
   * @param position The position to set the shape at
   * @param orientation The orientation to set the shape at
   */
  Shape(double              radius,
        const Point<dim>   &position,
        const Tensor<1, 3> &orientation)

    : AutoDerivativeFunction<dim>(1e-8)
    , effective_radius(radius)
    , position(position)
    , orientation(orientation)
    , part_of_a_composite(false)
    , layer_thickening(0.)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point
   * Most level set functions implemented come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  virtual double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override = 0;

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point with a guess for the cell containing
   * the evaluation point
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  virtual double
  value_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int                                   component = 0);

  /**
   * @brief Return the smoothed maximum of two variables used for shape contact calculation.
   * @param a first variable
   * @param b second variable
   */
  double
  smooth_max(double a, double b, double smooth_factor = 10)
  {
    return (a * std::exp(a * smooth_factor) + b * std::exp(b * smooth_factor)) /
           ((std::exp(a * smooth_factor) + std::exp(b * smooth_factor)));
  }

  /**
   * @brief Return the distance, the center point, and the unit normal vector between the current shape and the shape given in the argument. The center point is the point where both shapes are at the same distance from each other. The unit normal vector is defined using the closest surface point of the two shapes. By default, the function does not calculate the distance between the two shapes if their bounding boxes are not in contact. This behavior can be modified by setting exact_distance_outside_of_contact to true.
   * @param shape The shape with which the distance is evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param candidate_points This is the initial guess points used in the calculation.
   * @param precision This is the precision of the distance between the two shapes.
   * @param exact_distance_outside_of_contact This is a boolean to force the exact distance evaluation if the shapes are not in contact.
   */
  virtual std::tuple<double, Tensor<1, dim>, Point<dim>>
  distance_to_shape_with_cell_guess(
    Shape<dim>                                           &shape,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    std::vector<Point<dim>>                              &candidate_points,
    double                                                precision = 1e-6,
    bool exact_distance_outside_of_contact                          = false)
  {
    // Initialized the default values.
    double         distance = DBL_MAX;
    Tensor<1, dim> normal;
    Point<dim>     contact_point;
    bool           bounding_box_contact = true;

    // Check the bounding boxes. First we check the spherical bounding boxes
    // then the rectangular boxes.
    if (exact_distance_outside_of_contact == false)
      {
        Point<dim> bounding_box_center_one;
        Point<dim> bounding_box_center_two;
        double     radius_one = this->bounding_box_half_length.norm();
        double     radius_two = shape.bounding_box_half_length.norm();
        // In the case of the plane, the radius of the bounding box is zero, so
        // we replace it with the radius of the other object.
        if (radius_one == 0)
          {
            radius_one = radius_two;
          }
        if (radius_two == 0)
          {
            radius_two = radius_one;
          }
        if constexpr (dim == 3)
          {
            bounding_box_center_one =
              this->position +
              this->rotation_matrix * this->bounding_box_center;
            bounding_box_center_two =
              shape.position +
              shape.rotation_matrix * shape.bounding_box_center;
          }
        else
          {
            bounding_box_center_one =
              this->position +
              tensor_nd_to_2d(this->rotation_matrix *
                              point_nd_to_3d(this->bounding_box_center));
            bounding_box_center_two =
              shape.position +
              tensor_nd_to_2d(shape.rotation_matrix *
                              point_nd_to_3d(shape.bounding_box_center));
          }
        if (((bounding_box_center_one - bounding_box_center_two).norm() -
             radius_one - radius_two) <= 0)
          {
            bounding_box_contact = this->bounding_box_contact(shape);
          }
        else
          {
            bounding_box_contact = false;
          }
      }

    // The following algorithm does a minimization of the level set obtained by
    // the intersection of two shapes. A usual gradient descent does not work as
    // such. This algorithm uses a cartesian search combined with an additional
    // guess based on the gradient of the level set and a numerical gradient
    // descent.
    if (bounding_box_contact)
      {
        std::vector<Tensor<1, dim>> search_direction(2 * dim);
        // Define the cartesian direction search directions.
        search_direction[0][0] = 1.;
        search_direction[1][0] = -1.;
        search_direction[2][1] = 1.;
        search_direction[3][1] = -1.;
        if constexpr (dim == 3)
          {
            search_direction[4][2] = 1.;
            search_direction[5][2] = -1.;
          }
        // Loop over all the initial candidate points.
        for (unsigned int i = 0; i < candidate_points.size(); ++i)
          {
            // Initialize variables used in the calculation.
            Point<dim> current_point = candidate_points[i];
            Point<dim> dx{}, distance_gradient{}, previous_position{},
              previous_gradient{};
            // Initialize the iteration counter.
            unsigned int       iteration     = 0;
            const unsigned int iteration_max = 2e2;

            // The initial step size is set to 25% of the effective radius of
            // the shape. Tests showed that this initial step size is generally
            // adequate.
            double max_step           = shape.effective_radius * 0.25;
            double previous_step_size = max_step;

            // Initialize the value. In this minimization of the intersection of
            // two level set we use the smooth max function for the union as it
            // significantly smooths the minimization problem as we get close to
            // the local minimum of the intersection.
            double value_first_component =
              this->value_with_cell_guess(current_point, cell);
            double value_second_component =
              shape.value_with_cell_guess(current_point, cell);
            double current_distance =
              smooth_max(value_first_component, value_second_component);

            previous_position = current_point;
            // The initial previous value and position
            double       previous_value     = DBL_MAX;
            unsigned int consecutive_center = 0;

            // Iterate to find the minimum of the intersection.
            while (iteration < iteration_max && previous_step_size > precision)
              {
                // Check the local gradient direction.
                value_first_component =
                  this->value_with_cell_guess(current_point, cell);
                value_second_component =
                  shape.value_with_cell_guess(current_point, cell);
                if (value_first_component > value_second_component)
                  {
                    distance_gradient =
                      this->gradient_with_cell_guess(current_point, cell);
                  }
                else
                  {
                    distance_gradient =
                      shape.gradient_with_cell_guess(current_point, cell);
                  }
                Tensor<1, dim> direction = distance_gradient;
                dx = -max_step * direction / direction.norm();
                if (distance_gradient.norm() == 0)
                  {
                    dx = -max_step * direction;
                  }
                // Check the first candidate point in the direction of the
                // gradient.
                Point<dim> new_point = current_point + dx;
                value_first_component =
                  this->value_with_cell_guess(new_point, cell);
                value_second_component =
                  shape.value_with_cell_guess(new_point, cell);
                double new_distance =
                  smooth_max(value_first_component, value_second_component);
                // Check if the guess is better than the previous guess. If it
                // is not better, we do the cartesian directional search.
                if (new_distance > current_distance - precision * precision)
                  {
                    // Initialize the container for the value.
                    std::vector<double> diff_results;
                    diff_results.resize(search_direction.size());
                    Point<dim> best_point;
                    double     best_dist = current_distance;
                    for (unsigned int d = 0; d < search_direction.size(); ++d)
                      {
                        Tensor<1, dim> perturbation =
                          search_direction[d] * max_step;
                        value_first_component = this->value_with_cell_guess(
                          current_point + perturbation, cell);
                        value_second_component = shape.value_with_cell_guess(
                          current_point + perturbation, cell);
                        new_distance = smooth_max(value_first_component,
                                                  value_second_component);
                        diff_results[d] =
                          (new_distance - current_distance) / max_step;
                        // Check if the Cartesian search is better than the last
                        // candidate.
                        if (new_distance < best_dist - precision * precision)
                          {
                            best_dist  = new_distance;
                            best_point = current_point + perturbation;
                          }
                      }
                    // Store the current results.
                    if (best_dist < current_distance - precision * precision)
                      {
                        current_distance = best_dist;
                        current_point    = best_point;
                      }
                    Point<dim> extra_guess;
                    // If none of the cartesian candidates are better, we use
                    // the value to find the point that minimizes the numerical
                    // gradient we evaluated with the cartesian guess.
                    if (current_distance >
                        previous_value - precision * precision)
                      {
                        // Define the guess point.
                        new_point[0] =
                          previous_position[0] -
                          ((diff_results[0] - diff_results[1]) / 2) /
                            ((diff_results[0] + diff_results[1]) / max_step);
                        new_point[1] =
                          previous_position[1] -
                          ((diff_results[2] - diff_results[3]) / 2) /
                            ((diff_results[2] + diff_results[3]) / max_step);
                        if constexpr (dim == 3)
                          {
                            new_point[2] =
                              previous_position[2] -
                              ((diff_results[4] - diff_results[5]) / 2) /
                                ((diff_results[4] + diff_results[5]) /
                                 max_step);
                          }
                        value_first_component =
                          this->value_with_cell_guess(new_point, cell);
                        value_second_component =
                          shape.value_with_cell_guess(new_point, cell);
                        new_distance = smooth_max(value_first_component,
                                                  value_second_component);
                        // Check if it is better.
                        if (new_distance <
                            previous_value - precision * precision)
                          {
                            // If it is, we reduce the step size to the size of
                            // the step made by this guess.
                            max_step = (new_point - previous_position).norm();
                            current_distance = new_distance;
                            current_point    = new_point;
                          }
                      }
                    // If the guess is not better, we reduce the search size and
                    // store the value. If one of the points is better, we reset
                    // the consecutive_center counter. This variable helps to
                    // converge faster when one of the guesses is very close to
                    // the optimal point.
                    if (current_distance >
                        previous_value - precision * precision)
                      {
                        consecutive_center += 1;
                        previous_step_size = max_step;
                        max_step *= 1.0 / std::pow(2.0, consecutive_center);
                        current_point    = previous_position;
                        current_distance = previous_value;
                      }
                    else
                      {
                        max_step *= 1;
                        consecutive_center = 0;
                      }
                  }
                else
                  {
                    current_point    = new_point;
                    current_distance = new_distance;
                    max_step *= 1;
                    consecutive_center = 0;
                  }

                // Store the previous values.
                previous_value    = current_distance;
                previous_position = current_point;
                iteration++;
              }
            // Store the minimum obtained with this initial guess and compare it
            // with the results for other initial points.
            if (distance > current_distance)
              {
                distance      = current_distance;
                contact_point = current_point;
                // The normal is defined using the closest surface point on the
                // two shapes.
                Point<dim> closest_surface_point_on_this_shape;
                Point<dim> closest_surface_point_on_the_input_shape;
                this->closest_surface_point(current_point,
                                            closest_surface_point_on_this_shape,
                                            cell);
                shape.closest_surface_point(
                  current_point,
                  closest_surface_point_on_the_input_shape,
                  cell);
                if (distance > 0)
                  {
                    normal = (closest_surface_point_on_the_input_shape -
                              closest_surface_point_on_this_shape) /
                             ((closest_surface_point_on_the_input_shape -
                               closest_surface_point_on_this_shape)
                                .norm() +
                              DBL_MIN);
                  }
                else
                  {
                    normal = (closest_surface_point_on_this_shape -
                              closest_surface_point_on_the_input_shape) /
                             ((closest_surface_point_on_this_shape -
                               closest_surface_point_on_the_input_shape)
                                .norm() +
                              DBL_MIN);
                  }
              }
          }
      }
    return std::make_tuple(distance, normal, contact_point);
  }

  /**
   * @brief Return the distance, the center point, and the normal between the current shape and the shape given in the argument. The center point is the point where both shapes are at the same distance from each other. The normal is defined using the closest surface point on the two shapes. By default, the function does not calculate the distance between the two shapes if their bounding boxes are not in contact. This behavior can be modified using the appropriate parameter.
   * @param shape The shape with which the distance is evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param candidate_points This is the initial guess points used in the calculation.
   * @param precision This is the precision of the distance between the two shapes.
   * @param exact_distance_outside_of_contact This is a boolean to force the exact distance evaluation if the shapes are not in contact.
   */
  virtual std::tuple<double, Tensor<1, dim>, Point<dim>>
  distance_to_shape(Shape<dim>              &shape,
                    std::vector<Point<dim>> &candidate_points,
                    double                   precision     = 1e-6,
                    bool exact_distance_outside_of_contact = false)
  {
    // Initialized the default values.
    double         distance = DBL_MAX;
    Tensor<1, dim> normal;
    Point<dim>     contact_point;
    bool           bounding_box_contact = true;

    // Check the bounding boxes. First we check the spherical bounding boxes
    // then the rectangular boxes.
    if (exact_distance_outside_of_contact == false)
      {
        Point<dim> bounding_box_center_one;
        Point<dim> bounding_box_center_two;
        double     radius_one = this->bounding_box_half_length.norm();
        double     radius_two = shape.bounding_box_half_length.norm();
        // In the case of the plane, the radius of the bounding box is zero, so
        // we replace it with the radius of the other object.
        if (radius_one == 0)
          {
            radius_one = radius_two;
          }
        if (radius_two == 0)
          {
            radius_two = radius_one;
          }
        if constexpr (dim == 3)
          {
            bounding_box_center_one =
              this->position +
              this->rotation_matrix * this->bounding_box_center;
            bounding_box_center_two =
              shape.position +
              shape.rotation_matrix * shape.bounding_box_center;
          }
        else
          {
            bounding_box_center_one =
              this->position +
              tensor_nd_to_2d(this->rotation_matrix *
                              point_nd_to_3d(this->bounding_box_center));
            bounding_box_center_two =
              shape.position +
              tensor_nd_to_2d(shape.rotation_matrix *
                              point_nd_to_3d(shape.bounding_box_center));
          }
        if (((bounding_box_center_one - bounding_box_center_two).norm() -
             radius_one - radius_two) <= 0)
          {
            bounding_box_contact = this->bounding_box_contact(shape);
          }
        else
          {
            bounding_box_contact = false;
          }
      }

    // The following algorithm does a minimization of the level_set obtained by
    // the intersection of two shapes. A usual gradient descent does not work as
    // such. This algorithm uses a cartesian search combined with an additional
    // guess based on the gradient of the level set and a numerical gradient
    // descent.
    if (bounding_box_contact)
      {
        std::vector<Tensor<1, dim>> search_direction(2 * dim);
        // Define the cartesian direction search directions.
        search_direction[0][0] = 1.;
        search_direction[1][0] = -1.;
        search_direction[2][1] = 1.;
        search_direction[3][1] = -1.;
        if constexpr (dim == 3)
          {
            search_direction[4][2] = 1.;
            search_direction[5][2] = -1.;
          }
        // Loop over all the initial candidate points.
        for (unsigned int i = 0; i < candidate_points.size(); ++i)
          {
            // Initialize variable used in the calculation.
            Point<dim> current_point = candidate_points[i];
            Point<dim> dx{}, distance_gradient{}, previous_position{},
              previous_gradient{};
            // Initialize the iteration counter. We limit the number of
            // iterations to 200.
            unsigned int       iteration     = 0;
            const unsigned int iteration_max = 2e2;

            // The initial step size is set to 25% of the effective radius of
            // the shape. Tests showed that this initial step size is generally
            // adequate.
            double max_step           = shape.effective_radius * 0.25;
            double previous_step_size = max_step;

            // Initialize the value. In this minimisation of the intersection of
            // two level set we use the smooth max function for the union as it
            // significantly smooths the minimization problem as we get close to
            // the local minimum of the intersection.
            double value_first_component  = this->value(current_point);
            double value_second_component = shape.value(current_point);
            double current_distance =
              smooth_max(value_first_component, value_second_component);

            previous_position = current_point;
            // The initial previous value and position
            double       previous_value     = DBL_MAX;
            unsigned int consecutive_center = 0;

            // Iterate to find the minimum of the intersection.
            while (iteration < iteration_max &&
                   (previous_step_size) > precision)
              {
                // Check the local gradient direction.
                value_first_component  = this->value(current_point);
                value_second_component = shape.value(current_point);
                if (value_first_component > value_second_component)
                  {
                    distance_gradient = this->gradient(current_point);
                  }
                else
                  {
                    distance_gradient = shape.gradient(current_point);
                  }
                Tensor<1, dim> direction = distance_gradient;
                dx = -max_step * direction / direction.norm();
                if (distance_gradient.norm() == 0)
                  {
                    dx = -max_step * direction;
                  }
                // Check the first candidate point in the direction of the
                // gradient.
                Point<dim> new_point   = current_point + dx;
                value_first_component  = this->value(new_point);
                value_second_component = shape.value(new_point);
                double new_distance =
                  smooth_max(value_first_component, value_second_component);
                // Check if the guess is better than the previous guess. If it
                // is not better, we do the cartesian directional search.
                if (new_distance > current_distance - precision * precision)
                  {
                    // Initialize the container for the value.
                    std::vector<double> diff_results;
                    diff_results.resize(search_direction.size());
                    Point<dim> best_point;
                    double     best_dist = current_distance;
                    for (unsigned int d = 0; d < search_direction.size(); ++d)
                      {
                        Tensor<1, dim> perturbation =
                          search_direction[d] * max_step;
                        value_first_component =
                          this->value(current_point + perturbation);
                        value_second_component =
                          shape.value(current_point + perturbation);
                        new_distance = smooth_max(value_first_component,
                                                  value_second_component);
                        diff_results[d] =
                          (new_distance - current_distance) / max_step;
                        // Check if the Cartesian search is better than the last
                        // candidate.
                        if (new_distance < best_dist - precision * precision)
                          {
                            best_dist  = new_distance;
                            best_point = current_point + perturbation;
                          }
                      }
                    // Store the current results.
                    if (best_dist < current_distance - precision * precision)
                      {
                        current_distance = best_dist;
                        current_point    = best_point;
                      }
                    Point<dim> extra_guess;
                    // If none of the cartesian candidates are better, we use
                    // the value to find the point that minimizes the numerical
                    // gradient we evaluated with the cartesian guess.
                    if (current_distance >
                        previous_value - precision * precision)
                      {
                        // Define the guess point.
                        new_point[0] =
                          previous_position[0] -
                          ((diff_results[0] - diff_results[1]) / 2) /
                            ((diff_results[0] + diff_results[1]) / max_step);
                        new_point[1] =
                          previous_position[1] -
                          ((diff_results[2] - diff_results[3]) / 2) /
                            ((diff_results[2] + diff_results[3]) / max_step);
                        if constexpr (dim == 3)
                          {
                            new_point[2] =
                              previous_position[2] -
                              ((diff_results[4] - diff_results[5]) / 2) /
                                ((diff_results[4] + diff_results[5]) /
                                 max_step);
                          }
                        value_first_component  = this->value(new_point);
                        value_second_component = shape.value(new_point);
                        new_distance = smooth_max(value_first_component,
                                                  value_second_component);
                        // Check if it is better.
                        if (new_distance <
                            previous_value - precision * precision)
                          {
                            // If it is, we reduce the step size to the size of
                            // the step made by this guess.
                            max_step = (new_point - previous_position).norm();
                            current_distance = new_distance;
                            current_point    = new_point;
                          }
                      }
                    // If the guess is not better, we reduce the search size and
                    // store the value. If one of the points is better, we reset
                    // the consecutive_center counter. This variable helps to
                    // converge faster when one of the guesses is very close to
                    // the optimal point.
                    if (current_distance >
                        previous_value - precision * precision)
                      {
                        consecutive_center += 1;
                        previous_step_size = max_step;
                        max_step *= 1.0 / std::pow(2.0, consecutive_center);
                        current_point    = previous_position;
                        current_distance = previous_value;
                      }
                    else
                      {
                        max_step *= 1;
                        consecutive_center = 0;
                      }
                  }
                else
                  {
                    current_point    = new_point;
                    current_distance = new_distance;
                    max_step *= 1;
                    consecutive_center = 0;
                  }

                // Store the previous values.
                previous_value    = current_distance;
                previous_position = current_point;
                iteration++;
              }
            // Store the minimum obtained with this initial guess and compare it
            // with the results for other initial points.
            if (distance > current_distance)
              {
                distance      = current_distance;
                contact_point = current_point;
                Point<dim> closest_surface_point_on_this_shape;
                Point<dim> closest_surface_point_on_the_input_shape;
                // The normal is defined using the closest surface point on the
                // two shapes.
                this->closest_surface_point(
                  current_point, closest_surface_point_on_this_shape);
                shape.closest_surface_point(
                  current_point, closest_surface_point_on_the_input_shape);
                if (distance > 0)
                  {
                    normal = (closest_surface_point_on_the_input_shape -
                              closest_surface_point_on_this_shape) /
                             ((closest_surface_point_on_the_input_shape -
                               closest_surface_point_on_this_shape)
                                .norm() +
                              DBL_MIN);
                  }
                else
                  {
                    normal = (closest_surface_point_on_this_shape -
                              closest_surface_point_on_the_input_shape) /
                             ((closest_surface_point_on_this_shape -
                               closest_surface_point_on_the_input_shape)
                                .norm() +
                              DBL_MIN);
                  }
              }
          }
      }
    return std::make_tuple(distance, normal, contact_point);
  }


  /**
   * @brief Check if the bounding box of this shape and another one are overlapping.
   * This function uses the Separating Axis Theorem (SAT) to perform the
   * bounding box check.
   * @param shape The shape with which the distance is evaluated
   * @return This function returns false if the shapes are not overlapping.
   */
  bool
  bounding_box_contact(Shape<dim> &shape)
  {
    if constexpr (dim == 3)
      {
        // First, if one of the bounding boxes is a plane, then the size of the
        // bounding box is set to zero. As such, we test the other box vertices
        // directly.
        if (this->bounding_box_half_length.norm() == 0)
          {
            // loop over the vertices;
            for (int i = -1; i < 2; i += 2)
              {
                for (int j = -1; j < 2; j += 2)
                  {
                    for (int k = -1; k < 2; k += 2)
                      {
                        // create the vertex to test and check the value
                        Tensor<1, dim> vertex_position(
                          {(shape.bounding_box_half_length[0] +
                            shape.layer_thickening) *
                             i,
                           (shape.bounding_box_half_length[1] +
                            shape.layer_thickening) *
                             j,
                           (shape.bounding_box_half_length[2] +
                            shape.layer_thickening) *
                             k});
                        Point<3> vertex =
                          point_nd_to_3d(shape.position) +
                          shape.rotation_matrix *
                            point_nd_to_3d(shape.bounding_box_center +
                                           vertex_position);
                        if (this->value(vertex) <= 0)
                          {
                            return true;
                          }
                      }
                  }
              }
            return false;
          }
        else if (shape.bounding_box_half_length.norm() == 0)
          {
            // loop over the vertices;
            for (int i = -1; i < 2; i += 2)
              {
                for (int j = -1; j < 2; j += 2)
                  {
                    for (int k = -1; k < 2; k += 2)
                      {
                        // create the vertex to test and check the value
                        Tensor<1, dim> vertex_position(
                          {(this->bounding_box_half_length[0] +
                            this->layer_thickening) *
                             i,
                           (this->bounding_box_half_length[1] +
                            this->layer_thickening) *
                             j,
                           (this->bounding_box_half_length[2] +
                            this->layer_thickening) *
                             k});
                        Point<3> vertex =
                          point_nd_to_3d(this->position) +
                          this->rotation_matrix *
                            point_nd_to_3d(this->bounding_box_center +
                                           vertex_position);
                        if (shape.value(vertex) <= 0)
                          {
                            return true;
                          }
                      }
                  }
              }
            return false;
          }

        double         ra, rb;
        Tensor<2, dim> R, AbsR;
        // Compute rotation matrix expressing each box in the other's coordinate
        // frame
        for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
              {
                R[i][j] = scalar_product(this->rotation_matrix[i],
                                         shape.rotation_matrix[j]);
              }
          }
        // Compute translation vector t in this shape frame
        Tensor<1, dim> t =
          shape.position + shape.rotation_matrix * shape.bounding_box_center -
          (this->position + this->rotation_matrix * this->bounding_box_center);
        t = this->rotation_matrix * t;
        // Compute common subexpressions. Add in an epsilon term to counteract
        // arithmetic errors
        for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
              {
                AbsR[i][j] =
                  std::abs(R[i][j]) + std::numeric_limits<double>::epsilon();
              }
          }
        // Test axes L = A0, L = A1, L = A2
        for (int i = 0; i < 3; i++)
          {
            ra = (this->bounding_box_half_length[i] + this->layer_thickening);
            rb = (shape.bounding_box_half_length[0] + shape.layer_thickening) *
                   AbsR[i][0] +
                 (shape.bounding_box_half_length[1] + shape.layer_thickening) *
                   AbsR[i][1] +
                 (shape.bounding_box_half_length[2] + shape.layer_thickening) *
                   AbsR[i][2];
            if (std::abs(t[i]) > ra + rb)
              return false;
          }
        // Test axes L = B0, L = B1, L = B2
        for (int i = 0; i < 3; i++)
          {
            ra = (this->bounding_box_half_length[0] + this->layer_thickening) *
                   AbsR[0][i] +
                 (this->bounding_box_half_length[1] + this->layer_thickening) *
                   AbsR[1][i] +
                 (this->bounding_box_half_length[2] + this->layer_thickening) *
                   AbsR[2][i];
            rb = shape.bounding_box_half_length[i] + shape.layer_thickening;
            if (std::abs(t[0] * R[0][i] + t[1] * R[1][i] + t[2] * R[2][i]) >
                ra + rb)
              return false;
          }
        // Check the axis obtained from the cross product of all the bases axis.
        Tensor<1, dim> L;
        for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
              {
                L  = cross_product_3d(this->rotation_matrix[i],
                                     shape.rotation_matrix[j]);
                ra = (this->bounding_box_half_length[(i + 1) % 3] +
                      this->layer_thickening) *
                       AbsR[(i + 2) % 3][j] +
                     (this->bounding_box_half_length[(i + 2) % 3] +
                      this->layer_thickening) *
                       AbsR[(i + 1) % 3][j];
                rb = (shape.bounding_box_half_length[(j + 1) % 3] +
                      shape.layer_thickening) *
                       AbsR[i][(j + 2) % 3] +
                     (shape.bounding_box_half_length[(j + 2) % 3] +
                      shape.layer_thickening) *
                       AbsR[i][(j + 1) % 3];
                if (std::abs(scalar_product(t, L)) > ra + rb)
                  return false;
              }
          }
        // No separating axis found, the boxes must be intersecting
        return true;
      }
    else
      {
        // First, if one of the bounding boxes is a plane, then the size of the
        // bounding box is zero. As such, we test the other box vertices
        // directly.
        if (this->bounding_box_half_length.norm() == 0)
          {
            // loop over the vertices;
            for (int i = -1; i < 2; i += 2)
              {
                for (int j = -1; j < 2; j += 2)
                  {
                    // create the vertex to test and check the value
                    Tensor<1, dim> vertex_position(
                      {(shape.bounding_box_half_length[0] +
                        shape.layer_thickening) *
                         i,
                       (shape.bounding_box_half_length[1] +
                        shape.layer_thickening) *
                         j});
                    Point<3> vertex =
                      point_nd_to_3d(shape.position) +
                      shape.rotation_matrix *
                        point_nd_to_3d(shape.bounding_box_center +
                                       vertex_position);
                    if (this->value(point_nd_to_2d(vertex)) <= 0)
                      {
                        return false;
                      }
                  }
              }
            return true;
          }
        else if (shape.bounding_box_half_length.norm() == 0)
          {
            // loop over the vertices;
            for (int i = -1; i < 2; i += 2)
              {
                for (int j = -1; j < 2; j += 2)
                  {
                    // create the vertex to test and check the value
                    Tensor<1, 3> vertex_position(
                      {(this->bounding_box_half_length[0] +
                        this->layer_thickening) *
                         i,
                       (this->bounding_box_half_length[1] +
                        this->layer_thickening) *
                         j,
                       0});
                    Point<3> vertex =
                      point_nd_to_3d(this->position) +
                      this->rotation_matrix *
                        (point_nd_to_3d(this->bounding_box_center) +
                         vertex_position);
                    if (shape.value(point_nd_to_2d(vertex)) <= 0)
                      {
                        return false;
                      }
                  }
              }
            return true;
          }
        double       ra, rb;
        Tensor<2, 3> R, AbsR;

        // Compute rotation matrix expressing each box in the other's coordinate
        // frame
        for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
              {
                R[i][j] = scalar_product(this->rotation_matrix[i],
                                         shape.rotation_matrix[j]);
              }
          }

        // Compute translation vector t in this shape frame
        Tensor<1, 3> t =
          point_nd_to_3d(shape.position) +
          shape.rotation_matrix * point_nd_to_3d(shape.bounding_box_center) -
          (point_nd_to_3d(this->position) +
           this->rotation_matrix * point_nd_to_3d(this->bounding_box_center));
        t = this->rotation_matrix * t;

        // Compute common subexpressions. Add in an epsilon term to counteract
        // arithmetic errors
        for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
              {
                AbsR[i][j] =
                  std::abs(R[i][j]) + std::numeric_limits<double>::epsilon();
              }
          }

        // Test axes L = A0, L = A1, L = A2
        for (int i = 0; i < 2; i++)
          {
            ra = this->bounding_box_half_length[i] + this->layer_thickening;
            rb = (shape.bounding_box_half_length[0] + shape.layer_thickening) *
                   AbsR[i][0] +
                 (shape.bounding_box_half_length[1] + shape.layer_thickening) *
                   AbsR[i][1];
            if (std::abs(t[i]) > ra + rb)
              return false;
          }

        // Test axes L = B0, L = B1, L = B2
        for (int i = 0; i < 2; i++)
          {
            ra = (this->bounding_box_half_length[0] + this->layer_thickening) *
                   AbsR[0][i] +
                 (this->bounding_box_half_length[1] + this->layer_thickening) *
                   AbsR[1][i];
            rb = shape.bounding_box_half_length[i] + shape.layer_thickening;
            if (std::abs(t[0] * R[0][i] + t[1] * R[1][i]) > ra + rb)
              return false;
          }

        // No separating axis found, the boxes must be intersecting
        return true;
      }
  }

  /**
   * @brief Return the analytical gradient of the distance
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  virtual Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int                                   component = 0);


  /**
   * @brief Return a pointer to a copy of the Shape
   */
  virtual std::shared_ptr<Shape<dim>>
  static_copy() const = 0;


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
  virtual void
  closest_surface_point(
    const Point<dim>                                     &p,
    Point<dim>                                           &closest_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);
  virtual void
  closest_surface_point(const Point<dim> &p, Point<dim> &closest_point) const;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  virtual double
  displaced_volume();


  /**
   * @brief
   * Return the curvature of the level set at a given point.
   */
  virtual double
  local_curvature_radius(Point<dim> p)
  {
    double curvature = 0;
    double epsilone  = 1e-6 * this->effective_radius;
    for (unsigned int d = 0; d < dim; ++d)
      {
        Tensor<1, dim> perturbation;
        perturbation[d]           = epsilone;
        Tensor<1, dim> gradient_p = this->gradient(p + perturbation);
        gradient_p                = gradient_p / (gradient_p.norm() + DBL_MIN);
        Tensor<1, dim> gradient_m = this->gradient(p - perturbation);
        gradient_m                = gradient_m / (gradient_m.norm() + DBL_MIN);
        curvature += (gradient_p - gradient_m)[d] / (2 * epsilone);
      }
    // Here we bound the curvature radius obtained
    double radius_of_curvature = 1.0 / abs(curvature);
    return std::min(std::max(radius_of_curvature,
                             this->effective_radius * 1e-2),
                    this->effective_radius * 1e2);
  }
  /**
   * @brief
   * Return the curvature of the level at a given point using the cell guess
   * functions.
   * @param cell The cell that is likely to contain the evaluation point
   * @param p The point where the curvature is evaluated.
   */
  virtual double
  local_curvature_radius_with_cell_guess(
    const Point<dim>                                     p,
    const typename DoFHandler<dim>::active_cell_iterator cell)
  {
    double curvature = 0;
    double epsilone  = 1e-6 * this->effective_radius;
    // Calculate the curvature using the divergence of the normal and define the
    // local curvature radius as the absolute value of the inverse of the
    // curvature.
    for (unsigned int d = 0; d < dim; ++d)
      {
        Tensor<1, dim> perturbation;
        perturbation[d] = epsilone;
        Tensor<1, dim> gradient_p =
          this->gradient_with_cell_guess(p + perturbation, cell);
        gradient_p = gradient_p / (gradient_p.norm() + DBL_MIN);
        Tensor<1, dim> gradient_m =
          this->gradient_with_cell_guess(p - perturbation, cell);
        gradient_m = gradient_m / (gradient_m.norm() + DBL_MIN);
        curvature += (gradient_p - gradient_m)[d] / (2 * epsilone);
      }
    // Here we bound the curvature radius obtained
    double radius_of_curvature = 1.0 / abs(curvature);
    return std::min(std::max(radius_of_curvature,
                             this->effective_radius * 1e-2),
                    this->effective_radius * 1e2);
  }

  /**
   * @brief
   * Sets a new position for the shape
   *
   * @param The new position the shape will be placed at
   */
  inline virtual void
  set_position(const Point<dim> &new_position)
  {
    clear_cache();
    position = new_position;
  }

  /**
   * @brief
   * Sets a new orientation for the shape
   *
   * @param The new orientation the shape will be set at
   */
  inline virtual void
  set_orientation(const Tensor<1, 3> &new_orientation)
  {
    clear_cache();
    this->rotation_matrix       = 0;
    this->rotation_matrix[0][0] = 1.0;
    this->rotation_matrix[1][1] = 1.0;
    this->rotation_matrix[2][2] = 1.0;
    if constexpr (dim == 2)
      {
        Tensor<1, 3> axis;
        axis[2] = 1.0;
        this->rotation_matrix =
          Physics::Transformations::Rotations::rotation_matrix_3d(
            axis, new_orientation[2]);
      }
    if constexpr (dim == 3)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            Tensor<1, 3> axis;
            axis[2 - i] = 1.0;
            this->rotation_matrix =
              Physics::Transformations::Rotations::rotation_matrix_3d(
                axis, new_orientation[2 - i]) *
              this->rotation_matrix;
          }
      }
    orientation = new_orientation;
  }
  /**
   * @brief
   * Returns the position of the shape
   *
   */
  inline virtual Point<dim>
  get_position()
  {
    return position;
  }

  /**
   * @brief
   * Returns the orientation of the shape
   *
   */
  inline virtual Tensor<1, 3>
  get_orientation()
  {
    return orientation;
  }

  /**
   * @brief
   * Returns the orientation of the shape
   *
   */
  inline virtual Tensor<2, 3>
  get_rotation_matrix()
  {
    return rotation_matrix;
  }

  /**
   * @brief
   * Returns the bounding box half length
   *
   */
  inline virtual Tensor<1, dim>
  get_bounding_box_half_length()
  {
    return bounding_box_half_length;
  }

  /**
   * @brief
   * Returns the bounding box position
   *
   */
  inline virtual Point<dim>
  get_bounding_box_center()
  {
    return bounding_box_center;
  }

  /**
   * @brief
   * Returns the default manifold of the shape. If not redefined, it is a flat
   * manifold.
   */
  virtual std::shared_ptr<Manifold<dim - 1, dim>>
  get_shape_manifold();


  /**
   * @brief
   * Clear the cache of the shape
   *
   */
  virtual void
  clear_cache();

  /**
   * @brief
   * Most value functions assume that the particle's position is at the origin
   * and that the shape is aligned with one of the main axes. This function
   * returns a point that is rotated and translated, in accordance with the
   * current shape position and orientation, so that subsequent calculations
   * for the value function are made more easily; it abstracts a step that is
   * required in the value function for most shapes.
   *
   * Returns the centered and aligned point used on the level set evaluation.
   *
   * @param evaluation_point The point that will be recentered and realigned
   * @return The aligned and centered point
   */
  Point<dim>
  align_and_center(const Point<dim> &evaluation_point) const;

  /**
   * @brief
   * This function applies the inverse operation of align_and_center
   *
   * Returns the centered and aligned point used on the level set evaluation in
   * the global reference frame.
   *
   * @param evaluation_point The point that will be centered and aligned in the global reference frame.
   * @return The aligned and centered point
   */
  Point<dim>
  reverse_align_and_center(const Point<dim> &evaluation_point) const;

  /**
   * @brief
   * This function returns a point in a string of text. This is used in the
   * cache of the shape.
   *
   * @param evaluation_point is the point that is transformed to its text form.
   */
  std::string
  point_to_string(const Point<dim> &evaluation_point) const;

  /**
   * @brief Defines if this shape is part of a composite.
   * If true, cache management is deactivated and delegated to the upper level
   * shape, to avoid cache duplication.
   * @param part_of_a_composite is true if this shape is a constituent of a composite shape
   */
  void
  set_part_of_a_composite(const bool part_of_a_composite)
  {
    this->part_of_a_composite = part_of_a_composite;
  }

  /**
   * @brief
   * Sets the layer thickening value (positive or negative) of the particle's
   * shape
   *
   * @param layer_thickening Thickness to be artificially added to the particle.
   * A negative value will decrease the particle's thickness by subtracting a
   * layer of specified width.
   */
  virtual void
  set_layer_thickening(const double layer_thickening)
  {
    this->layer_thickening = layer_thickening;
  }

  /**
   * @brief
   * Function that return the XYZ rotation angles from a rotation matrix
   *
   * @param rotation_matrix_representation The rotation matrix.
   */
  Tensor<1, 3>
  rotation_matrix_to_xyz_angles(
    Tensor<2, 3> &rotation_matrix_representation) const
  {
    Tensor<1, 3> xyz_rotation;
    xyz_rotation[0] = std::atan2(-rotation_matrix_representation[1][2],
                                 rotation_matrix_representation[2][2]);
    xyz_rotation[1] = std::asin(rotation_matrix_representation[0][2]);
    xyz_rotation[2] = std::atan2(-rotation_matrix_representation[0][1],
                                 rotation_matrix_representation[0][0]);
    return xyz_rotation;
  };

  // Effective radius used for crown refinement
  double effective_radius;

  // The string contains additional information on the shape. This may refer
  // to the file type used to define the shape or any other information
  // relative to how the shape was defined.
  std::string additional_info_on_shape;

  // The half lengths of bounding box of the shape
  Tensor<1, dim> bounding_box_half_length;
  // The center of the bounding box of the shape in the shape reference frame
  Point<dim> bounding_box_center;

protected:
  // Position of the center of the Shape. It doesn't always correspond to the
  // center of mass
  Point<dim> position;
  // The solid orientation, which is defined as the sequential rotation around
  // the axes x->y->z by each of the tensor components, in radian
  Tensor<1, 3> orientation;

  // The cache of the evaluation of the shape. This is used to avoid costly
  // reevaluation of the shape.
  std::unordered_map<std::string, double>         value_cache;
  std::unordered_map<std::string, Tensor<1, dim>> gradient_cache;
  std::unordered_map<std::string, Point<dim>>     closest_point_cache;
  bool                                            part_of_a_composite;


  // The rotation matrix that describes the solid orientation.
  Tensor<2, 3> rotation_matrix;
  // Layer thickening: used to artificially inflate/deflate the shape
  double layer_thickening;
};


template <int dim>
class Sphere : public Shape<dim>
{
public:
  /**
   * @brief Constructor for a sphere
   * @param radius The sphere radius
   * @param position The sphere center
   * @param orientation The sphere orientation
   */
  Sphere(double              radius,
         const Point<dim>   &position,
         const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
  {
    sphere_function =
      std::make_shared<Functions::SignedDistance::Sphere<dim>>(position,
                                                               radius);
    for (unsigned int d = 0; d < dim; ++d)
      {
        this->bounding_box_half_length[d] = radius;
        this->bounding_box_center[d]      = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;


  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief Return the analytical gradient of the distance
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  Tensor<1, dim>
  gradient(const Point<dim>  &evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

  void
  set_position(const Point<dim> &position) override;

  /**
   * @brief Return the manifold of the sphere as a spherical manifold center at the position of the shape.
   */
  virtual std::shared_ptr<Manifold<dim - 1, dim>>
  get_shape_manifold() override;


  /**
   * @brief Return the distance, the center point, and the normal between the current shape and the shape given in the argument. The center point is the point that minimized the level set obtained from the intersection of the two shapes. The normal is defined using the closest surface point on the two shapes. This function is an optimized version of the general shape distance calculation for the case of a sphere.
   * @param shape The shape with which the distance is evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param candidate_points This is the initial guess points used in the calculation. (In this implementation this parameter is void)
   * @param precision This is the precision of the distance between the two shapes. (In this implementation this parameter is void)
   * @param exact_distance_outside_of_contact This is a boolean to force the exact distance evaluation if the shapes are not in contact. (In this implementation this parameter is void)
   */
  std::tuple<double, Tensor<1, dim>, Point<dim>>
  distance_to_shape_with_cell_guess(
    Shape<dim>                                           &shape,
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    std::vector<Point<dim>>                              &candidate_points,
    double                                                precision = 1e-6,
    bool exact_distance_outside_of_contact = false) override;

  /**
   * @brief Return the distance, the center point, and the normal between the current shape and the shape given in the argument. The center point is the point that minimized the level set obtained from the intersection of the two shapes. The normal is defined using the closest surface point on the two shapes. This function is an optimized version of the general shape distance calculation for the case of a sphere.
   * @param shape The shape with which the distance is evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param candidate_points This is the initial guess points used in the calculation. (In this implementation, this parameter is void)
   * @param precision This is the precision of the distance between the two shapes. (In this implementation, this parameter is void)
   * @param exact_distance_outside_of_contact This is a boolean to force the exact distance evaluation if the shapes are not in contact. (In this implementation, this parameter is void)
   */
  std::tuple<double, Tensor<1, dim>, Point<dim>>
  distance_to_shape(Shape<dim>              &shape,
                    std::vector<Point<dim>> &candidate_points,
                    double                   precision     = 1e-6,
                    bool exact_distance_outside_of_contact = false) override;

  /**
   * @brief
   * Return the curvature of the level set at a given point
   *
   * @param p The point where the curvature is evaluated.
   */
  virtual double
  local_curvature_radius(Point<dim> p) override;

  /**
   * @brief
   * Return the curvature of the level set at a given point using the
   * @param cell The cell that is likely to contain the evaluation point
   * @param p The point where the curvature is evaluated.
   */
  virtual double
  local_curvature_radius_with_cell_guess(
    const Point<dim>                                     p,
    const typename DoFHandler<dim>::active_cell_iterator cell) override;


private:
  std::shared_ptr<Functions::SignedDistance::Sphere<dim>> sphere_function;
};



template <int dim>
class Plane : public Shape<dim>
{
public:
  /**
   * @brief Constructor for a infinite plane. The plane is normal to the Z axis in 3D and Y axis in 2D in the reference frame of the particle.
   * @param position The point on the plane.
   * @param orientation The orientation.
   */
  Plane(const Point<dim> &position, const Tensor<1, 3> &orientation)
    : Shape<dim>(1.0, position, orientation)
  {
    if constexpr (dim == 3)
      {
        normal = Tensor<1, dim>({0, 0, 1});
      }
    else
      {
        normal = Tensor<1, dim>({0, 1});
      }
    // This is a special case since the plane has no bounding box so here we
    // fill it with zeros.
    for (unsigned int d = 0; d < dim; ++d)
      {
        this->bounding_box_half_length[d] = 0;
        this->bounding_box_center[d]      = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;


  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief Return the analytical gradient of the distance
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  Tensor<1, dim>
  gradient(const Point<dim>  &evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;



private:
  Tensor<1, dim> normal;
};


/**
 * @brief This class defines superquadric shapes. Their signed distance
 * function is: \f$\left|\frac{x}{a}\right|^r + \left|\frac{y}{b}\right|^s +
 * \left|\frac{z}{c}\right|^t - 1 = 0\f$
 * @tparam dim Dimension of the shape
 */
template <int dim>
class Superquadric : public Shape<dim>
{
public:
  /**
   * @brief Constructor for a superquadric shape
   * @param half_lengths The half-lengths of each direction
   * @param exponents The blockiness in each direction
   * @param epsilon The tolerance for surface representation
   * @param position The superquadric center
   * @param orientation The superquadric orientation
   */
  Superquadric(const Tensor<1, dim> half_lengths,
               const Tensor<1, dim> exponents,
               const double         epsilon,
               const Point<dim>    &position,
               const Tensor<1, 3>  &orientation)
    : Shape<dim>(half_lengths.norm(), position, orientation)
    , half_lengths(half_lengths)
    , exponents(exponents)
    , epsilon(epsilon)
  {
    // Initialize the bounding box.
    this->bounding_box_half_length = half_lengths;
    for (unsigned int d = 0; d < dim; ++d)
      {
        this->bounding_box_center[d] = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point with a guess for the cell containing
   * the evaluation point
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int /*component = 0*/) override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief Return the analytical gradient of the distance
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  Tensor<1, dim>
  gradient(const Point<dim>  &evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief
   * Sets the closest_point parameter to be the point on the surface of the
   * shape which has the minimal distance from the given point p. Since this
   * function is at the core of the calculations for superquadrics, it is
   * redefined here.
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
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess) override;
  void
  closest_surface_point(const Point<dim> &p,
                        Point<dim>       &closest_point) const override;

  /**
   * @brief Return the sign of the parameter. To be removed once PR#794 is merged.
   * @param prm parameter
   */
  inline double
  sign(const double prm) const
  {
    if (prm > 0)
      return 1;
    else if (prm < 0)
      return -1;
    else
      return 0;
  }

  /**
   * @brief Computes the value of the superquadric from its equation
   * @param centered_point point at which we make the evaluation, in the shape referential
   */
  inline double
  superquadric(const Point<dim> &centered_point) const
  {
    return pow(abs(centered_point[0] / half_lengths[0]), exponents[0]) +
           pow(abs(centered_point[1] / half_lengths[1]), exponents[1]) +
           pow(abs(centered_point[2] / half_lengths[2]), exponents[2]) - 1.0;
  }

  /**
   * @brief Computes the gradient of the superquadric from its equation
   * @param centered_point point at which we make the evaluation, in the shape referential
   */
  inline Point<dim>
  superquadric_gradient(const Point<dim> &centered_point) const
  {
    Point<dim> gradient{};
    for (unsigned int d = 0; d < dim; d++)
      {
        // For cases where the coordinate is of value 0, we avoid computing the
        // gradient. That is an issue when the exponent is lower than or equal
        // to 1
        if (abs(centered_point[d]) > epsilon)
          gradient[d] = exponents[d] *
                        pow(abs(half_lengths[d]), -exponents[d]) *
                        pow(abs(centered_point[d]), exponents[d]) /
                        (centered_point[d] + DBL_MIN);
        else
          gradient[d] = exponents[d] *
                        pow(abs(half_lengths[d]), -exponents[d]) *
                        pow(abs(epsilon), exponents[d]) / (epsilon + DBL_MIN);
      }
    return gradient;
  }

private:
  Tensor<1, dim> half_lengths;
  Tensor<1, dim> exponents;
  double         epsilon;
};

template <int dim>
class HyperRectangle : public Shape<dim>
{
public:
  /**
   * @brief Constructs a box with the given parameters
   * @param half_lengths The half lengths of each direction
   * @param position The hyper rectangle center
   * @param orientation The hyper rectangle orientation
   */
  HyperRectangle(const Tensor<1, dim> &half_lengths,
                 const Point<dim>     &position,
                 const Tensor<1, 3>   &orientation)
    : Shape<dim>(half_lengths.norm(), position, orientation)
    , half_lengths(half_lengths)
  {
    this->bounding_box_half_length = half_lengths;
    for (unsigned int d = 0; d < dim; ++d)
      {
        this->bounding_box_center[d] = 0;
      }
  }


  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

  // Half-lengths of every side of the box
  Tensor<1, dim> half_lengths;
};

template <int dim>
class Ellipsoid : public Shape<dim>
{
public:
  /**
   * @brief Constructs an ellipsoid with the given arguments
   * @param radii The radii of each direction
   * @param position The ellipsoid center
   * @param orientation The ellipsoid orientation
   */
  Ellipsoid(const Tensor<1, dim> &radii,
            const Point<dim>     &position,
            const Tensor<1, 3>   &orientation)
    : Shape<dim>(radii.norm(), position, orientation)
    , radii(radii)
  {
    this->bounding_box_half_length = radii;
    for (unsigned int d = 0; d < dim; ++d)
      {
        this->bounding_box_center[d] = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  // The radii of all directions in which the ellipsoid is defined
  Tensor<1, dim> radii;
};

template <int dim>
class Torus : public Shape<dim>
{
public:
  /**
   * @brief Constructs a torus with a given thickness and radius
   * @param ring_radius The ring radius
   * @param ring_thickness The ring thickness radius/half-thickness
   * @param position The torus center
   * @param orientation The orientation of the axis at the center of the torus
   */
  Torus(double              ring_radius,
        double              ring_thickness,
        const Point<dim>   &position,
        const Tensor<1, 3> &orientation)
    : Shape<dim>(ring_thickness, position, orientation)
    , ring_radius(ring_radius)
    , ring_thickness(ring_thickness)
  {
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[0] = ring_radius + ring_thickness;
        this->bounding_box_half_length[1] = ring_thickness;
        this->bounding_box_half_length[2] = ring_radius + ring_thickness;
      };
    for (unsigned int d = 0; d < dim; ++d)
      {
        this->bounding_box_center[d] = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double ring_radius;
  double ring_thickness;
};

template <int dim>
class Cone : public Shape<dim>
{
public:
  /**
   * @brief Constructs a cone
   * @param tan_base_angle The tangent of the angle between the base of the cone and its curve side
   * @param height The height of the cone
   * @param position The position of the center of cone's base
   * @param orientation The orientation of the cone axis, from its base to its tip
   */
  Cone(double              tan_base_angle,
       double              height,
       const Point<dim>   &position,
       const Tensor<1, 3> &orientation)
    : Shape<dim>(height, position, orientation)
    , tan_base_angle(tan_base_angle)
    , height(height)
    , base_radius(height / tan_base_angle)
    , intermediate_q({height * tan_base_angle, -height})
  {
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[0] = base_radius;
        this->bounding_box_half_length[1] = height / 2;
        this->bounding_box_half_length[2] = base_radius;
        this->bounding_box_center[0]      = 0;
        this->bounding_box_center[1]      = -height / 2;
        this->bounding_box_center[2]      = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double tan_base_angle;
  double height;
  double base_radius;

  Tensor<1, 2> intermediate_q;
};

template <int dim>
class CutHollowSphere : public Shape<dim>
{
public:
  /**
   * @brief Constructs a hollow sphere that has a wall thickness and that is cut
   * by a given depth
   * @param radius The radius of the smallest sphere containing the cut hollow sphere
   * @param cut_depth The height of the slice removed from the sphere
   * @param shell_thickness The thickness of the hollow sphere shell
   * @param position The center of the sphere
   * @param orientation The orientation of the sphere, from it's center to the cut's center
   */
  CutHollowSphere(double              radius,
                  double              cut_depth,
                  double              shell_thickness,
                  const Point<dim>   &position,
                  const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
    , radius(radius)
    , cut_depth(cut_depth)
    , shell_thickness(shell_thickness)
    , intermediate_w(sqrt(radius * radius - cut_depth * cut_depth))
  {
    this->bounding_box_half_length[0] = radius;
    this->bounding_box_half_length[1] = radius / 2 + cut_depth / 2;
    this->bounding_box_center[0]      = 0;
    this->bounding_box_center[1]      = -radius / 2 + cut_depth / 2;
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[2] = radius;
        this->bounding_box_center[2]      = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double radius;
  double cut_depth;
  double shell_thickness;

  double intermediate_w;
};

template <int dim>
class DeathStar : public Shape<dim>
{
public:
  /**
   * @brief The Death Star is the result of a boolean substraction of one sphere from
   * another
   * @param radius The main sphere radius
   * @param hole_radius The removed sphere radius
   * @param spheres_distance The distance between the centers of the sphere
   * @param position The main sphere's center
   * @param orientation The orientation from the main sphere's center to the removed sphere's center
   */
  DeathStar(double              radius,
            double              hole_radius,
            double              spheres_distance,
            const Point<dim>   &position,
            const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
    , radius(radius)
    , hole_radius(hole_radius)
    , spheres_distance(spheres_distance)
    , intermediate_a((radius * radius - hole_radius * hole_radius +
                      spheres_distance * spheres_distance) /
                     (2. * spheres_distance))
    , intermediate_b(
        sqrt(std::max(radius * radius - intermediate_a * intermediate_a, 0.)))
  {
    this->bounding_box_half_length[0] = radius;
    this->bounding_box_half_length[0] = 0;
    this->bounding_box_half_length[1] = radius;
    this->bounding_box_half_length[1] = 0;
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[2] = radius;
        this->bounding_box_half_length[2] = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double radius;
  double hole_radius;
  double spheres_distance;

  double intermediate_a;
  double intermediate_b;
};

/**
 * @brief This class was designed so that specific composite shapes could be
 * defined through the parameter file. Boolean operations such as union,
 * difference, and intersection are allowed.
 * @tparam dim Dimension of the shape
 */
template <int dim>
class CompositeShape : public Shape<dim>
{
public:
  enum class BooleanOperation : int
  {
    Union,
    Difference,
    Intersection,
  };

  /**
   * @brief Constructs an assembly of shapes into a composite shape
   * @param constituents The shapes from which this composite shape will be composed
   * @param operations The list of operations to perform to construct the composite
   */
  CompositeShape(
    std::map<unsigned int, std::shared_ptr<Shape<dim>>> constituents,
    std::map<unsigned int,
             std::tuple<BooleanOperation, unsigned int, unsigned int>>
                        operations,
    const Point<dim>   &position,
    const Tensor<1, 3> &orientation)
    : Shape<dim>(0., position, orientation)
    , constituents(constituents)
    , operations(operations)
  {
    Point<dim> point_max;
    Point<dim> point_min;
    bool       xyz_min_max_is_initialized = false;

    // Calculation of the effective radius and setting of constituents' status
    for (auto const &constituent : constituents | boost::adaptors::map_values)
      {
        this->effective_radius =
          std::max(this->effective_radius, constituent->effective_radius);
        constituent->set_part_of_a_composite(true);
        std::vector<Point<3>> component_bounding_box_vertex(std::pow(2, dim));
        // For each of the components update the bounding box.
        if constexpr (dim == 2)
          {
            std::vector<Point<3>> component_bounding_box_vertex(4);
            unsigned int          index = 0;
            for (int i = -1; i < 2; i += 2)
              {
                for (int j = -1; j < 2; j += 2)
                  {
                    component_bounding_box_vertex[index][0] =
                      constituent->get_bounding_box_center()[0] +
                      constituent->get_bounding_box_half_length()[0] * i;
                    component_bounding_box_vertex[index][1] =
                      constituent->get_bounding_box_center()[1] +
                      constituent->get_bounding_box_half_length()[1] * j;
                    index += 1;
                  }
              }
            for (unsigned int i = 0; i < 4; ++i)
              {
                // Position the vertex point in the reference frame of the
                // composite
                constituent->set_position(constituent->get_position());
                constituent->set_orientation(constituent->get_orientation());
                auto new_point = constituent->get_rotation_matrix() *
                                   component_bounding_box_vertex[i] +
                                 point_nd_to_3d(constituent->get_position());
                if (xyz_min_max_is_initialized == false)
                  {
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        point_min[d] = new_point[d];
                        point_max[d] = new_point[d];
                      }
                    xyz_min_max_is_initialized = true;
                  }
                else
                  {
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        point_min[d] = std::min(point_min[d], new_point[d]);
                        point_max[d] = std::max(point_max[d], new_point[d]);
                      }
                  }
              }
          }
        else
          {
            std::vector<Point<3>> component_bounding_box_vertex(8);
            unsigned int          index = 0;
            for (int i = -1; i < 2; i += 2)
              {
                for (int j = -1; j < 2; j += 2)
                  {
                    for (int k = -1; k < 2; k += 2)
                      {
                        component_bounding_box_vertex[index][0] =
                          constituent->get_bounding_box_center()[0] +
                          constituent->get_bounding_box_half_length()[0] * i;
                        component_bounding_box_vertex[index][1] =
                          constituent->get_bounding_box_center()[1] +
                          constituent->get_bounding_box_half_length()[1] * j;
                        component_bounding_box_vertex[index][2] =
                          constituent->get_bounding_box_center()[2] +
                          constituent->get_bounding_box_half_length()[2] * k;
                        index += 1;
                      }
                  }
              }
            for (unsigned int i = 0; i < 8; ++i)
              {
                // Position the vertex point in the reference frame of the
                // composite
                constituent->set_position(constituent->get_position());
                constituent->set_orientation(constituent->get_orientation());
                auto new_point = constituent->get_rotation_matrix() *
                                   component_bounding_box_vertex[i] +
                                 point_nd_to_3d(constituent->get_position());

                if (xyz_min_max_is_initialized == false)
                  {
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        point_min[d] = new_point[d];
                        point_max[d] = new_point[d];
                      }
                    xyz_min_max_is_initialized = true;
                  }
                else
                  {
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        point_min[d] = std::min(point_min[d], new_point[d]);
                        point_max[d] = std::max(point_max[d], new_point[d]);
                      }
                  }
              }
          }
      }
    this->bounding_box_half_length = (point_max - point_min) / 2;
    this->bounding_box_center      = (point_max + point_min) / 2;
  }

  /**
   * @brief Constructs an assembly of shapes into a composite shape from a vector of shapes.
   * This constructor is mainly used for outputting multiple shapes with a
   * global level set function defined as a union.
   * @param constituents_vector The shapes from which this composite sphere will be composed
   */
  CompositeShape(std::vector<std::shared_ptr<Shape<dim>>> constituents_vector,
                 const Point<dim>                        &position,
                 const Tensor<1, 3>                      &orientation)
    : Shape<dim>(0., position, orientation)
  {
    size_t number_of_constituents = constituents_vector.size();

    for (size_t i = 0; i < number_of_constituents; i++)
      constituents[i] = constituents_vector[i];
    if (number_of_constituents > 1)
      {
        // If there are at least two components, the first operation should
        // always be a union of 0 and 1
        operations[number_of_constituents] =
          std::make_tuple(BooleanOperation::Union, 0, 1);
        // We make the union until the before last component
        for (size_t i = 1; i < number_of_constituents - 1; i++)
          {
            operations[i + number_of_constituents] =
              std::make_tuple(BooleanOperation::Union,
                              i + 1,
                              i + number_of_constituents - 1);
          }
      }
    // Calculation of the effective radius
    for (auto const &constituent : constituents | boost::adaptors::map_values)
      {
        this->effective_radius =
          std::max(this->effective_radius, constituent->effective_radius);
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point with a guess for the cell containing
   * the evaluation point
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int /*component = 0*/) override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient(const Point<dim>  &evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief Sets the proper dof handler, then computes/updates the map of cells
   * and their likely non-null nodes
   * @param updated_dof_handler the reference to the new dof_handler
   * @param mesh_based_precalculations mesh-based precalculations that can lead to slight shape misrepresentation (if RBF typed)
   */
  void
  update_precalculations(DoFHandler<dim> &updated_dof_handler,
                         const bool       mesh_based_precalculations);

  /**
   * @brief Load data from file. To be called at initialization, after repartitioning or when shape has moved.
   * This function is used only for RBF shapes and its composites at the
   * moment.
   */
  void
  load_data_from_file();

  /**
   * @brief Remove data that affects only artificial cells (not locally owned and not ghost).
   * The data is removed if it would never be accessed by the local process.
   * @param dof_handler the reference to the new dof_handler
   * @param mesh_based_precalculations mesh-based precalculations that can lead to slight shape misrepresentation (if RBF typed)
   */
  void
  remove_superfluous_data(DoFHandler<dim> &updated_dof_handler,
                          const bool       mesh_based_precalculations);

  /**
   * @brief Computes the assigned boolean operations
   * @param constituent_shapes_values map containing the computed values for the component shapes
   * @param constituent_shapes_gradients map containing the computed gradients for the component shapes
   */
  inline std::pair<double, Tensor<1, dim>>
  apply_boolean_operations(
    std::map<unsigned int, double>         constituent_shapes_values,
    std::map<unsigned int, Tensor<1, dim>> constituent_shapes_gradients) const;

  /**
   * @brief
   * Clear the cache of the shape
   */
  virtual void
  clear_cache() override;

  /**
   * @brief
   * See base
   */
  void
  set_layer_thickening(const double layer_thickening) override;

private:
  // The members of this class are all the constituent and operations that are
  // to be performed to construct the composite shape
  // This map link all primitive constituents of the composite shape to an id
  std::map<unsigned int, std::shared_ptr<Shape<dim>>> constituents;
  // This map links all operations between primitive constituents or
  // intermediate constituents (resulting from each operation) to an id. The
  // unsigned integers correspond to the first and second ids of the shapes
  // used for an operation
  std::map<unsigned int,
           std::tuple<BooleanOperation, unsigned int, unsigned int>>
    operations;
};

template <int dim>
class OpenCascadeShape : public Shape<dim>
{
public:
  /**
   * @brief Constructor for an OpenCascade shape
   * @param file_name The name of the file describing the shape
   * @param position The shape center
   * @param orientation The shape orientation
   */
  OpenCascadeShape(const std::string   file_name,
                   const Point<dim>   &position,
                   const Tensor<1, 3> &orientation)
    : Shape<dim>(0.1, position, orientation)
  {
    // First, we read the shape file name
    local_file_name = file_name;
#ifdef DEAL_II_WITH_OPENCASCADE
    // Checks the file name extension to identify which type of OpenCascade
    // shape we are working with.
    std::vector<std::string> file_name_and_extension(
      Utilities::split_string_list(local_file_name, "."));

    // Load the shape with the appropriate tool.
    if (file_name_and_extension[1] == "step" ||
        file_name_and_extension[1] == "stp")
      {
        shape = OpenCASCADE::read_STEP(local_file_name);
        this->additional_info_on_shape = "step";
      }
    else if (file_name_and_extension[1] == "iges" ||
             file_name_and_extension[1] == "igs")
      {
        shape = OpenCASCADE::read_IGES(local_file_name);
        this->additional_info_on_shape = "iges";
      }
    else if (file_name_and_extension[1] == "stl")
      {
        shape                          = OpenCASCADE::read_STL(local_file_name);
        this->additional_info_on_shape = "stl";
      }
    else
      {
        throw std::runtime_error(
          "Wrong file extension for an opencascade shape. The possible file type are: step, stp, iges, igs, stl.");
      }

    // used this local variable as the shape tolerance in the calculations.
    shape_tol = OpenCASCADE::get_shape_tolerance(shape);

    // Initialize some variables and the OpenCascade distance tool.
    OpenCASCADE::extract_compound_shapes(
      shape, compounds, compsolids, solids, shells, wires);
    vertex_position = OpenCASCADE::point(Point<dim>());
    vertex          = BRepBuilderAPI_MakeVertex(vertex_position);
    distancetool    = BRepExtrema_DistShapeShape(shape, vertex);

    // Check if the shape has a shell. If it has a shell, we initialize a
    // distance tool with just the shell.
    if (shells.size() > 0)
      {
        // Check if the number of solids is precisely 1. If it is, we redefine
        // the shape as only the solid. If it is not the case and there are
        // multiple shells, we throw an error since we won't be able to
        // represent the shape correctly.
        if (solids.size() == 1)
          {
            // Extract the solid
            shape = solids[0];
            // Extract the shell 0.
            OpenCASCADE::extract_compound_shapes(
              shape, compounds, compsolids, solids, shells, wires);
            // Load the tools
            distancetool = BRepExtrema_DistShapeShape(shells[0], vertex);
            point_classifier.Load(shape);
          }
        else if (shells.size() > 1)
          {
            throw std::runtime_error(
              "Error!: The shape has more than one shell. The code does not support shapes with multiple shells or solids. If your shape has more than one shell or solid, it is usually possible to recombine them into one. Otherwise, it is possible to split the shape into sub-shells and sub-solids and then define one particle for each of them.");
          }
        distancetool = BRepExtrema_DistShapeShape(shells[0], vertex);
        point_classifier.Load(shape);
      }
    else
      {
        distancetool = BRepExtrema_DistShapeShape(shape, vertex);
        point_classifier.Load(shape);
      }

    // Define the effective radius as the raidus of the sphere with the same
    // volume as the shape.
    GProp_GProps system;
    BRepGProp::LinearProperties(shape, system);
    BRepGProp::SurfaceProperties(shape, system);
    BRepGProp::VolumeProperties(shape, system);
    this->effective_radius =
      std::pow(system.Mass() * 3.0 / (4 * numbers::PI), 1.0 / dim);

    // Define an "infinite" bounding box size
    this->bounding_box_half_length[0] = DBL_MAX;
    this->bounding_box_half_length[1] = DBL_MAX;
    this->bounding_box_center[0]      = 0;
    this->bounding_box_center[1]      = 0;
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[2] = DBL_MAX;
        this->bounding_box_center[2]      = 0;
      }

#endif
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluate.
   * @param component Not applicable
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point. Use
   * @param component Not applicable
   */
  double
  value_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient(const Point<dim>  &evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

  /**
   * @brief
   * Sets a new position for the shape
   *
   * @param The new position the shape will be placed at
   */
  void
  set_position(const Point<dim> &position) override;

private:
  // Keep in memory the file name of the shape.
  std::string local_file_name;
#ifdef DEAL_II_WITH_OPENCASCADE
  // The shape define by the opencascade file.
  TopoDS_Shape shape;

  // We split the shape into its components we store them in these containers.
  std::vector<TopoDS_Compound>  compounds;
  std::vector<TopoDS_CompSolid> compsolids;
  std::vector<TopoDS_Solid>     solids;
  std::vector<TopoDS_Shell>     shells;
  std::vector<TopoDS_Wire>      wires;

  // The point at which we are going to evaluate the shape.
  gp_Pnt vertex_position;

  // The shape object (vertex) used in the distance evaluation.
  TopoDS_Vertex vertex;

  // The tool used for the distance evaluation.
  BRepClass3d_SolidClassifier point_classifier;
  BRepExtrema_DistShapeShape  distancetool;
  BRepExtrema_DistShapeShape  distancetool_shell;
  double                      shape_tol;
#endif
};


/**
 * @tparam dim Dimension of the shape
 * @brief RBF Shapes express the signed distance function as a linear
 * combination of Radial Basis Functions (RBF), which have a defined support
 * radius and basis function. A collection of nodes and weights compose the
 * object. Outside of the domain covered by the nodes, the distance is
 * computed by using the distance to a bounding box instead.
 */
template <int dim>
class RBFShape : public Shape<dim>
{
public:
  /**
   * Class taken from Optimad Bitpit. https://github.com/optimad/bitpit
   * @enum RBFBasisFunction
   * @brief Enum class defining types of RBF kernel functions that could be used
   * in the class
   */
  enum class RBFBasisFunction : int
  {
    CUSTOM,
    WENDLANDC2,
    LINEAR,
    GAUSS90,
    GAUSS95,
    GAUSS99,
    C1C0,
    C2C0,
    C0C1,
    C1C1,
    C2C1,
    C0C2,
    C1C2,
    C2C2,
    COS,
  };

  /**
   * @brief An RBFShape represents a physical object by describing its signed
   * distance field with a linear combination of radial basis functions. Each
   * radial basis function has a location and properties that are used in the
   * sum.
   * @param shape_arguments_str the name of the file used to load the data
   * @param position the location of the RBF shape
   * @param orientation the orientation of the shape with respect to each main
   * axis
   */
  RBFShape(const std::string  &shape_arguments_str,
           const Point<dim>   &position,
           const Tensor<1, 3> &orientation);

  /**
   * @brief Load RBF data from file. To be called at initialization, after repartitioning or when shape has moved.
   * This function is used only for RBF shapes and its composites at the
   * moment.
   */
  void
  load_data_from_file();

  /**
   * @brief Remove data that affects only artificial cells (not locally owned and not ghost).
   * The data is removed if it would never be accessed by the local process.
   * @param dof_handler the reference to the new dof_handler
   * @param mesh_based_precalculations mesh-based precalculations that can lead to slight shape misrepresentation (if RBF typed)
   */
  void
  remove_superfluous_data(DoFHandler<dim> &dof_handler,
                          const bool       mesh_based_precalculations);

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point. The computation and addition of the
   * bounding box distance are necessary since the RBF nodes may not cover the
   * whole simulation domain. In that case, it is assumed that the distance
   * from the RBF object is approximately the same as the distance from the
   * corresponding bounding box.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point with a guess for the cell containing
   * the evaluation point
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value_with_cell_guess(
    const Point<dim> &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
    const unsigned int /*component = 0*/) override;

  /**
   * @brief Return the analytical gradient of the distance
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim>                                    &evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief Return the analytical gradient of the distance for the current RBF
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  Tensor<1, dim>
  gradient(const Point<dim>  &evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

  /**
   * A bounding box is constructed around the collection of nodes defining the
   * RBF. It solves an issue where the collection of nodes is located only
   * around the object itself, which would result in an undefined distance
   * when the value is evaluated outside of all support radii. The hyper
   * rectangle shape doesn't have this limitation, as its distance can be
   * evaluated anywhere. The distance computed by an RBF object will therefore
   * use an approximated distance when the evaluation point is too far.
   *
   * @brief Initializes the bounding box around the nodes which enables distance
   * calculation even if the RBF nodes don't cover the whole domain
   */
  void
  initialize_bounding_box();

  /**
   * @brief Returns the value of the basis function for a given distance.
   * Inspired by Optimad Bitpit. https://github.com/optimad/bitpit
   * @param basis_function basis function to be used for calculation
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  evaluate_basis_function(const RBFBasisFunction basis_function,
                          const double           distance) const;

  /**
   * @brief Returns the derivative of the basis function for a given distance.
   * Inspired by Optimad Bitpit. https://github.com/optimad/bitpit
   * @param basis_function basis function to be used for calculation
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  evaluate_basis_function_derivative(const RBFBasisFunction basis_function,
                                     const double           distance) const;

  /**
   * @brief Establishes which nodes bring a non null contribution to the RBF
   * @param cell the cell for which the likely nodes are to be found
   * @param support_point one point that is located inside the cell
   */
  void
  determine_likely_nodes_for_one_cell(
    const typename DoFHandler<dim>::cell_iterator &cell,
    const Point<dim>                               support_point);

  /**
   * @brief Sets the proper dof handler, then computes/updates the map of cells
   * and their likely non-null nodes
   * @param dof_handler the reference to the new dof_handler
   * @param mesh_based_precalculations mesh-based precalculations that can lead to slight shape misrepresentation (if RBF typed)
   * */
  void
  update_precalculations(DoFHandler<dim> &updated_dof_handler,
                         const bool       mesh_based_precalculations);

  /**
   * @brief Rotate RBF nodes in the global reference frame (the reference frame of the triangulation).
   */
  void
  rotate_nodes();

  /**
   * @brief Compact Wendland C2 function defined from 0 to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  wendlandc2(const double distance) const
  {
    return distance > 1.0 ?
             0.0 :
             std::pow(1. - distance, 4.0) * (4.0 * distance + 1.0);
  }

  /**
   * @brief Compact linear function defined from 0 to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  linear(const double distance) const
  {
    return distance > 1.0 ? 0.0 : (1.0 - distance);
  }

  /**
   * @brief Non-compact Gaussian function with 0.1 value at distance equal to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  gauss90(const double distance) const
  {
    double eps = std::pow(-1.0 * std::log(0.1), 0.5);
    return std::exp(-1.0 * std::pow(distance * eps, 2.0));
  }

  /**
   * @brief Non-compact Gaussian function with 0.05 value at distance equal to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  gauss95(const double distance) const
  {
    double eps = std::pow(-1.0 * std::log(0.05), 0.5);
    return std::exp(-1.0 * std::pow(distance * eps, 2.0));
  }

  /**
   * @brief Non-compact Gaussian function with 0.01 value at distance equal to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  gauss99(const double distance) const
  {
    double eps = std::pow(-1.0 * std::log(0.01), 0.5);
    return std::exp(-1.0 * std::pow(distance * eps, 2.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C0 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c0(const double distance) const
  {
    return distance > 1.0 ? 0.0 : (1.0 - std::pow(distance, 2.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C0 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c0(const double distance) const
  {
    return distance > 1.0 ? 0.0 : (1.0 - std::pow(distance, 3.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C0 continuity at distance=0, and C1 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c0c1(const double distance) const
  {
    return distance > 1.0 ? 0.0 :
                            (1.0 - 2.0 * distance + std::pow(distance, 2.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C1 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c1(const double distance) const
  {
    return distance > 1.0 ? 0.0 :
                            (1.0 - 3.0 * std::pow(distance, 2.0) +
                             2.0 * std::pow(distance, 3.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C1 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c1(const double distance) const
  {
    return distance > 1.0 ? 0.0 :
                            (1.0 - 4.0 * std::pow(distance, 3.0) +
                             3.0 * std::pow(distance, 4.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C0 continuity at distance=0, and C2 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c0c2(const double distance) const
  {
    return distance > 1.0 ?
             0.0 :
             (1.0 - 3.0 * distance + 3.0 * std::pow(distance, 2.0) -
              std::pow(distance, 3.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C2 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c2(const double distance) const
  {
    return distance > 1.0 ?
             0.0 :
             (1.0 - 6.0 * std::pow(distance, 2.0) +
              8.0 * std::pow(distance, 3.0) - 3.0 * std::pow(distance, 4.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C2 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c2(const double distance) const
  {
    return distance > 1.0 ?
             0.0 :
             (1.0 - 10.0 * std::pow(distance, 3.0) +
              15.0 * std::pow(distance, 4.0) - 6.0 * std::pow(distance, 5.0));
  }

  /**
   * @brief Compact cosinusoidal basis function. It is null when r>1
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  cos(const double distance) const
  {
    return distance > 1.0 ? 0.0 : 0.5 + 0.5 * std::cos(distance * M_PI);
  }


  /**
   * @brief Derivative of a compact Wendland C2 function defined from 0 to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  wendlandc2_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 :
                            -20.0 * distance * std::pow(1.0 - distance, 3.0);
  }

  /**
   * @brief Derivative of a compact linear function defined from 0 to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  linear_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -1.0;
  }

  /**
   * @brief Derivative of a non-compact Gaussian function with 0.1 value at distance equal to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  gauss90_derivative(const double distance) const
  {
    double eps = std::pow(-1.0 * std::log(0.1), 0.5);
    return -2.0 * std::pow(eps, 2.0) * distance *
           std::exp(-1.0 * std::pow(distance * eps, 2.0));
  }

  /**
   * @brief Derivative of a non-compact Gaussian function with 0.05 value at distance equal to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  gauss95_derivative(const double distance) const
  {
    double eps = std::pow(-1.0 * std::log(0.05), 0.5);
    return -2.0 * std::pow(eps, 2.0) * distance *
           std::exp(-1.0 * std::pow(distance * eps, 2.0));
  }

  /**
   * @brief Derivative of a non-compact Gaussian function with 0.01 value at distance equal to 1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  gauss99_derivative(const double distance) const
  {
    double eps = std::pow(-1.0 * std::log(0.01), 0.5);
    return -2.0 * std::pow(eps, 2.0) * distance *
           std::exp(-1.0 * std::pow(distance * eps, 2.0));
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C0 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c0_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -2.0 * distance;
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C0 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c0_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -3.0 * std::pow(distance, 2.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C0 continuity at distance=0, and C1 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c0c1_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : 2.0 * (distance - 1.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C1 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c1_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : 6.0 * (distance - 1.0) * distance;
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C1 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c1_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 :
                            12 * (distance - 1.0) * std::pow(distance, 2.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C0 continuity at distance=0, and C2 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c0c2_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -3.0 * std::pow(distance - 1.0, 2.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C2 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c2_derivative(const double distance) const
  {
    return distance > 1.0 ?
             0.0 :
             (1.0 - 6.0 * std::pow(distance, 2.0) +
              8.0 * std::pow(distance, 3.0) - 3.0 * std::pow(distance, 4.0));
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C2 continuity at
   * distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c2_derivative(const double distance) const
  {
    return distance > 1.0 ?
             0.0 :
             -30.0 * std::pow(distance - 1.0, 2.0) * std::pow(distance, 2.0);
  }

  /**
   * @brief Derivative of a compact cosinusoidal basis function.
   * It preserves continuity at every point
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  cosinus_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -M_PI_2 * std::sin(M_PI * distance);
  }

  /**
   * @brief Swap the vector of all nodes with a likely node vector
   * @param cell A likely one where the evaluation point is located
   */
  void
  swap_iterable_nodes(
    const typename DoFHandler<dim>::active_cell_iterator cell);

private:
  std::string                          filename;
  size_t                               number_of_nodes;
  std::shared_ptr<HyperRectangle<dim>> bounding_box;

  // Elements of this vector are tuples containing: the cell barycenter, the
  // cell diameter, and the RBF nodes located inside the active cell
  std::vector<
    std::tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>>
    iterable_nodes;
  // Entries of this map contain vectors of tuples containing the cell
  // barycenter, diameter and likely nodes that are in that cell. The various
  // elements of the vector are the tuples which contain information on the
  // RBF nodes that affect the level set evaluation in this cell. Note: the
  // division of RBF nodes is made by separating nodes by which active cells
  // they are in, and then these RBF nodes groups are referenced by the cells
  // in a subsequent cell (by using a vector of pointers). This is required to
  // avoid repeating the nodes IDs multiple times (to keep memory requirements
  // low); they only appear once in a vector, then this vector is used by
  // passing its shared pointer.
  std::map<
    const typename DoFHandler<dim>::cell_iterator,
    std::shared_ptr<std::vector<
      std::tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>>>>
                   likely_nodes_map;
  size_t           max_number_of_inside_nodes;
  DoFHandler<dim> *dof_handler;

  double maximal_support_radius;

public:
  std::vector<double>     weights;
  std::vector<Point<dim>> nodes_positions;
  std::vector<Point<dim>> rotated_nodes_positions;
  std::vector<double>     support_radii;
  std::vector<double>     basis_functions;
  std::vector<bool>       useful_rbf_nodes;
};


template <int dim>
class Cylinder : public Shape<dim>
{
public:
  /**
   * @brief Constructs a cylinder aligned with the z axis
   * @param radius the radius of the cylinder
   * @param half_length the half-length of the cylinder
   * @param position position of the barycenter of the cylinder
   * @param orientation orientation of the cylinder
   */
  Cylinder(double              radius,
           double              half_length,
           const Point<dim>   &position,
           const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
    , radius(radius)
    , half_length(half_length)
  {
    this->bounding_box_half_length[0] = radius;
    this->bounding_box_half_length[1] = radius;
    this->bounding_box_half_length[2] = half_length;
    this->bounding_box_center[0]      = 0;
    this->bounding_box_center[1]      = 0;
    this->bounding_box_center[2]      = 0;
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double radius;
  double half_length;
};

template <int dim>
class CylindricalTube : public Shape<dim>
{
public:
  /**
   * @brief Constructs a tube by boolean difference of two tubes aligned with
   * the z axis when orientation is set to 0;0;0
   * @param radius_inside the radius of the negative (inside) cylinder
   * @param radius_outside the radius of the positive (outside) cylinder
   * @param half_length the half-length of the cylinders
   * @param position position of the barycenter of the cylinder
   * @param orientation orientation of the cylinder
   */
  CylindricalTube(double              radius_inside,
                  double              radius_outside,
                  double              half_length,
                  const Point<dim>   &position,
                  const Tensor<1, 3> &orientation)
    : Shape<dim>((radius_outside + radius_inside) / 2., position, orientation)
    , radius((radius_outside + radius_inside) / 2.)
    , height(half_length * 2.)
    , rectangular_base(radius_outside - radius_inside)
  {
    this->bounding_box_half_length[0] = radius_outside;
    this->bounding_box_half_length[1] = radius_outside;
    this->bounding_box_center[0]      = 0;
    this->bounding_box_center[1]      = 0;
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[2] = half_length;
        this->bounding_box_center[2]      = 0;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double radius;
  double height;
  double rectangular_base;
};

template <int dim>
class CylindricalHelix : public Shape<dim>
{
public:
  /**
   * @brief Constructs a cylindrical helix by extruding a disk through a helicoidal path
   * aligned with the z axis when orientation is set to 0;0;0
   * @param radius_helix the radius of the helicoidal path
   * @param radius_disk the radius of the disk that is extruded along the helicoidal
   * path
   * @param height the total height of the helicoidal path
   * @param pitch the height difference between each helix loop around its axis
   * @param position the position of the helix base
   * @param orientation the orientation of the helix axis compared to the z axis
   */
  CylindricalHelix(double              radius_helix,
                   double              radius_disk,
                   double              height,
                   double              pitch,
                   const Point<dim>   &position,
                   const Tensor<1, 3> &orientation)
    : Shape<dim>(radius_disk, position, orientation)
    , radius(radius_helix)
    , height(height)
    , pitch(pitch)
    , radius_disk(radius_disk)
  {
    this->bounding_box_half_length[0] = radius_helix + radius_disk;
    this->bounding_box_half_length[1] = radius_helix + radius_disk;
    this->bounding_box_center[0]      = 0;
    this->bounding_box_center[1]      = 0;
    if constexpr (dim == 3)
      {
        this->bounding_box_half_length[2] = height / 2 + radius_disk;
        this->bounding_box_center[2]      = height / 2;
      }
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim>  &evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid.
   */
  double
  displaced_volume() override;

private:
  double radius;
  double height;
  double pitch;
  double radius_disk;
};

#endif // lethe_shape_h
