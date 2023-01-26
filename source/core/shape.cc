/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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
 */
#include <core/shape.h>

#include <deal.II/grid/manifold_lib.h>

template <int dim>
double
Shape<dim>::displaced_volume(const double /*fluid_density*/)
{
  StandardExceptions::ExcNotImplemented();
  return 1.0;
}

template <int dim>
Point<dim>
Shape<dim>::align_and_center(const Point<dim> &evaluation_point) const
{
  // Translation and rotation to standard position and orientation for
  // distance calculations
  Point<dim> center_of_rotation = position;

  Point<dim> rotated_point;
  Point<dim> translated_point;

  // Rotation from the solid orientation
  // Angular position around x, y and z axis
  Tensor<1, 3> theta = orientation;

  // The centralized point is the one to be rotated, and it is updated after
  // each rotation around one axis The centralized rotated point is the result
  // of each rotation, and it is initialized in case no rotation is performed
  Point<dim> centralized_point;
  centralized_point              = evaluation_point - center_of_rotation;
  Point<dim> centralized_rotated = centralized_point;

  // Selection of the first axis around which to rotate:
  // x -> 0, y -> 1, z -> 2
  // In 2D, only rotation around the z axis is possible
  if constexpr (dim == 2)
    {
      if (std::abs(theta[2]) > 1e-10)
        {
          Tensor<2, 2> rotation_matrix =
            Physics::Transformations::Rotations::rotation_matrix_2d(-theta[2]);

          // Multiplication
          centralized_rotated.clear();
          for (unsigned int j = 0; j < dim; ++j)
            {
              for (unsigned int k = 0; k < dim; ++k)
                {
                  centralized_rotated[j] +=
                    rotation_matrix[j][k] * centralized_point[k];
                }
            }
          centralized_point = centralized_rotated;
        }
    }
  else // (dim == 3)
    {
      for (unsigned int i = 0; i < 3; ++i)
        {
          if (std::abs(theta[i]) > 1e-10)
            {
              Tensor<1, 3> axis;
              axis[i] = 1.0;
              Tensor<2, 3> rotation_matrix =
                Physics::Transformations::Rotations::rotation_matrix_3d(
                  axis, -theta[i]);

              // Multiplication
              centralized_rotated.clear();
              for (unsigned int j = 0; j < dim; ++j)
                {
                  for (unsigned int k = 0; k < dim; ++k)
                    {
                      centralized_rotated[j] +=
                        rotation_matrix[j][k] * centralized_point[k];
                    }
                }
              centralized_point = centralized_rotated;
            }
        }
    }
  rotated_point = centralized_rotated + center_of_rotation;

  // Translation from the solid position
  translated_point = rotated_point - position;

  return translated_point;
}


template <int dim>
double
Shape<dim>::value_with_cell_guess(
  const Point<dim> &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
  const unsigned int /*component*/)
{
  return this->value(evaluation_point);
}

template <int dim>
Tensor<1, dim>
Shape<dim>::gradient_with_cell_guess(
  const Point<dim> &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
  const unsigned int /*component*/)
{
  return this->gradient(evaluation_point);
}

template <int dim>
double
Sphere<dim>::value(const Point<dim> &evaluation_point,
                   const unsigned int /*component*/) const
{
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  return evaluation_point.distance(this->position) - this->effective_radius;
#else
  return sphere_function->value(evaluation_point);
#endif
}

template <int dim>
std::shared_ptr<Manifold<dim - 1, dim>>
Shape<dim>::get_shape_manifold()
{
  return std::make_shared<FlatManifold<dim - 1, dim>>();
}

template <int dim>
std::shared_ptr<Shape<dim>>
Sphere<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<Sphere<dim>>(this->effective_radius,
                                  this->position,
                                  this->orientation);
  return copy;
}

template <int dim>
Tensor<1, dim>
Sphere<dim>::gradient(const Point<dim> &evaluation_point,
                      const unsigned int /*component*/) const
{
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  const Tensor<1, dim> center_to_point = evaluation_point - this->position;
  const Tensor<1, dim> grad = center_to_point / center_to_point.norm();
  return grad;
#else
  return sphere_function->gradient(evaluation_point);
#endif
}

template <int dim>
std::shared_ptr<Manifold<dim - 1, dim>>
Sphere<dim>::get_shape_manifold()
{
  return std::make_shared<SphericalManifold<dim - 1, dim>>(this->position);
}

template <int dim>
double
Sphere<dim>::displaced_volume(const double fluid_density)
{
  double solid_volume;
  using numbers::PI;
  if (dim == 2)
    solid_volume =
      this->effective_radius * this->effective_radius * PI * fluid_density;

  else if (dim == 3)
    solid_volume = 4.0 / 3.0 * this->effective_radius * this->effective_radius *
                   this->effective_radius * PI;
  return solid_volume;
}

template <int dim>
void
Sphere<dim>::set_position(const Point<dim> &position)
{
  this->Shape<dim>::set_position(position);
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
#else
  sphere_function = std::make_shared<Functions::SignedDistance::Sphere<dim>>(
    position, this->effective_radius);
#endif
}

template <int dim>
double
Rectangle<dim>::value(const Point<dim> &evaluation_point,
                      const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  Point<dim> abs_p;
  Point<dim> half_lengths_dim;
  for (unsigned int i = 0; i < dim; ++i)
    {
      abs_p[i]            = std::abs(centered_point[i]);
      half_lengths_dim[i] = half_lengths[i];
    }
  Point<dim> q;
  q = abs_p - half_lengths_dim;
  Point<dim> max_q_0;
  for (unsigned int i = 0; i < dim; ++i)
    {
      max_q_0[i] = std::max(q[i], 0.0);
    }
  double max_q = std::max(q[0], std::max(q[1], q[dim - 1]));
  return max_q_0.norm() + std::min(max_q, 0.0);
}

template <int dim>
std::shared_ptr<Shape<dim>>
Rectangle<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<Rectangle<dim>>(this->half_lengths,
                                     this->position,
                                     this->orientation);
  return copy;
}

template <int dim>
double
Rectangle<dim>::displaced_volume(const double /*fluid_density*/)
{
  double solid_volume = 1.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      solid_volume = solid_volume * 2.0 * half_lengths[dim];
    }
  return solid_volume;
}

template <int dim>
double
Ellipsoid<dim>::value(const Point<dim> &evaluation_point,
                      const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  Point<dim> v_k0;
  Point<dim> v_k1;
  for (unsigned int i = 0; i < dim; ++i)
    {
      // Ellipsoid parameters[0->2]:= radii x, y, z
      v_k0[i] = centered_point[i] / radii[i];
      v_k1[i] = centered_point[i] / (radii[i] * radii[i]);
    }
  double k0 = v_k0.norm();
  double k1 = v_k1.norm();
  return k0 * (k0 - 1.) / k1;
}

template <int dim>
std::shared_ptr<Shape<dim>>
Ellipsoid<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<Ellipsoid<dim>>(this->radii,
                                     this->position,
                                     this->orientation);
  return copy;
}

template <int dim>
double
Ellipsoid<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  double solid_volume = PI * 4.0 / 3.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      solid_volume = solid_volume * radii[dim];
    }
  return solid_volume;
}

template <int dim>
double
Torus<dim>::value(const Point<dim> &evaluation_point,
                  const unsigned int /*component*/) const
{
  AssertDimension(dim, 3);

  Point<dim> centered_point = this->align_and_center(evaluation_point);

  Point<2> p_xz({centered_point[0], centered_point[2]});
  Point<2> q({p_xz.norm() - ring_radius, centered_point[1]});
  return q.norm() - ring_thickness;
}

template <int dim>
std::shared_ptr<Shape<dim>>
Torus<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy = std::make_shared<Torus<dim>>(
    this->ring_radius, this->ring_thickness, this->position, this->orientation);
  return copy;
}

template <int dim>
double
Torus<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  return 2.0 * PI * PI * ring_radius * ring_thickness * ring_thickness;
}

template <int dim>
double
Cone<dim>::value(const Point<dim> &evaluation_point,
                 const unsigned int /*component*/) const
{
  AssertDimension(dim, 3);

  Point<dim> centered_point = this->align_and_center(evaluation_point);

  // For a cone, the parameters are tan(base angle) and height
  Point<2> p_xz({centered_point[0], centered_point[2]});
  Point<2> w({p_xz.norm(), centered_point[1]});
  double   dot_w_q = scalar_product<1, 2, double>(w, intermediate_q);
  double dot_q_q = scalar_product<1, 2, double>(intermediate_q, intermediate_q);
  Point<2> a;
  a = w - intermediate_q * std::clamp(dot_w_q / dot_q_q, 0., 1.0);
  Point<2> b_intermediate1({std::clamp(w[0] / intermediate_q[0], 0.0, 1.), 1.});
  Point<2> b_intermediate2({intermediate_q[0] * b_intermediate1[0],
                            intermediate_q[1] * b_intermediate1[1]});
  Point<2> b;
  b        = w - b_intermediate2;
  double k = (intermediate_q[1] > 0) ? 1 : ((intermediate_q[1] < 0) ? -1 : 0);
  double d = std::min(scalar_product<1, 2, double>(a, a),
                      scalar_product<1, 2, double>(b, b));
  double s = std::max(k * (w[0] * intermediate_q[1] - w[1] * intermediate_q[0]),
                      k * (w[1] - intermediate_q[1]));

  return sqrt(d) * ((s > 0) ? 1 : ((s < 0) ? -1 : 0));
}

template <int dim>
std::shared_ptr<Shape<dim>>
Cone<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy = std::make_shared<Cone<dim>>(
    this->tan_base_angle, this->height, this->position, this->orientation);
  return copy;
}

template <int dim>
double
Cone<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  return PI / 3.0 * base_radius * base_radius * height;
}

template <int dim>
double
CutHollowSphere<dim>::value(const Point<dim> &evaluation_point,
                            const unsigned int /*component*/) const
{
  AssertDimension(dim, 3);

  Point<dim> centered_point = this->align_and_center(evaluation_point);

  Point<2> p_xz({centered_point[0], centered_point[2]});
  Point<2> q({p_xz.norm(), centered_point[1]});

  if (cut_depth * q[0] < intermediate_w * q[1])
    {
      Point<2> wh({intermediate_w, cut_depth});
      return (q - wh).norm();
    }
  else
    {
      return std::abs(q.norm() - radius) - shell_thickness;
    }
}

template <int dim>
std::shared_ptr<Shape<dim>>
CutHollowSphere<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<CutHollowSphere<dim>>(this->radius,
                                           this->cut_depth,
                                           this->shell_thickness,
                                           this->position,
                                           this->orientation);
  return copy;
}

template <int dim>
double
CutHollowSphere<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  double small_radius = radius - shell_thickness;
  return 4.0 * PI / 3.0 *
         (radius * radius * radius -
          small_radius * small_radius * small_radius);
}

template <int dim>
double
DeathStar<dim>::value(const Point<dim> &evaluation_point,
                      const unsigned int /*component*/) const
{
  AssertDimension(dim, 3);

  Point<dim> centered_point = this->align_and_center(evaluation_point);

  Point<2> p_yz({centered_point[1], centered_point[2]});
  Point<2> corrected_p_2d({centered_point[0], p_yz.norm()});
  if (corrected_p_2d[0] * intermediate_b - corrected_p_2d[1] * intermediate_a >
      spheres_distance * std::max(intermediate_b - corrected_p_2d[1], 0.0))
    {
      Point<2> ab({intermediate_a, intermediate_b});
      return (corrected_p_2d - ab).norm();
    }
  else
    {
      Point<2> d0({spheres_distance, 0.});
      return std::max(corrected_p_2d.norm() - radius,
                      -((corrected_p_2d - d0).norm() - hole_radius));
    }
}

template <int dim>
std::shared_ptr<Shape<dim>>
DeathStar<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<DeathStar<dim>>(this->radius,
                                     this->hole_radius,
                                     this->spheres_distance,
                                     this->position,
                                     this->orientation);
  return copy;
}

template <int dim>
double
DeathStar<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  return 4. * PI / 3. * radius * radius * radius;
}

template <int dim>
double
CompositeShape<dim>::value(const Point<dim> &evaluation_point,
                           const unsigned int /*component*/) const
{
  // We align and center the evaluation point according to the shape referential
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  // The levelset value of all constituent shapes is computed
  std::map<unsigned int, double> components_value;
  for (auto const &[component_id, component] : constituents)
    {
      components_value[component_id] = component->value(centered_point);
    }

  // The boolean operations between the shapes are applied in order
  // The last computed levelset value is considered to be the right value
  double levelset = components_value[0];
  for (auto const &[operation_id, op_triplet] : operations)
    {
      BooleanOperation operation;
      unsigned int     first_id;
      unsigned int     second_id;
      std::tie(operation, first_id, second_id) = op_triplet;

      double value_first_component, value_second_component;
      value_first_component  = components_value.at(first_id);
      value_second_component = components_value.at(second_id);
      switch (operation)
        {
          case BooleanOperation::Union:
            components_value[operation_id] =
              std::min(value_first_component, value_second_component);
            break;
          case BooleanOperation::Difference:
            // We substract the first component to the second
            components_value[operation_id] =
              std::max(-value_first_component, value_second_component);
            break;
          case BooleanOperation::Intersection:
            components_value[operation_id] =
              std::max(value_first_component, value_second_component);
            break;
          default:
            throw std::logic_error(
              "The BooleanOperation isn't supported. Either it is not supported "
              "yet or it is simply not valid.");
        }
      levelset = components_value[operation_id];
    }
  return levelset;
}

template <int dim>
double
CompositeShape<dim>::value_with_cell_guess(
  const Point<dim> &                                   evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  // We align and center the evaluation point according to the shape referential
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  // The levelset value of all component shapes is computed
  std::map<unsigned int, double> components_value;
  for (auto const &[component_id, component] : constituents)
    {
      components_value[component_id] =
        component->value_with_cell_guess(centered_point, cell);
    }

  // The boolean operations between the shapes are applied in order
  // The last computed levelset value is considered to be the right value
  double levelset = components_value[0];
  for (auto const &[operation_id, op_triplet] : operations)
    {
      BooleanOperation operation;
      unsigned int     first_id;
      unsigned int     second_id;
      std::tie(operation, first_id, second_id) = op_triplet;

      double value_first_component, value_second_component;
      value_first_component  = components_value.at(first_id);
      value_second_component = components_value.at(second_id);
      switch (operation)
        {
          case BooleanOperation::Union:
            components_value[operation_id] =
              std::min(value_first_component, value_second_component);
            break;
          case BooleanOperation::Difference:
            components_value[operation_id] =
              std::max(-value_first_component, value_second_component);
            break;
          default: // BooleanOperation::Intersection
            components_value[operation_id] =
              std::max(value_first_component, value_second_component);
        }
      levelset = components_value[operation_id];
    }
  return levelset;
}

template <int dim>
std::shared_ptr<Shape<dim>>
CompositeShape<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy = std::make_shared<CompositeShape<dim>>(
    constituents, operations, this->position, this->orientation);
  return copy;
}

template <int dim>
void
CompositeShape<dim>::update_precalculations(
  DoFHandler<dim> &             updated_dof_handler,
  std::shared_ptr<Mapping<dim>> mapping)
{
  for (auto const &[component_id, component] : constituents)
    {
      if (typeid(*component) == typeid(RBFShape<dim>))
        {
          std::static_pointer_cast<RBFShape<dim>>(component)
            ->update_precalculations(updated_dof_handler, mapping);
        }
      else if (typeid(*component) == typeid(CompositeShape<dim>))
        {
          std::static_pointer_cast<CompositeShape<dim>>(component)
            ->update_precalculations(updated_dof_handler, mapping);
        }
    }
}

template <int dim>
RBFShape<dim>::RBFShape(const std::vector<double> &          support_radii,
                        const std::vector<RBFBasisFunction> &basis_functions,
                        const std::vector<double> &          weights,
                        const std::vector<Point<dim>> &      nodes,
                        const Point<dim> &                   position,
                        const Tensor<1, 3> &                 orientation)
  : Shape<dim>(support_radii[0], position, orientation)
  , number_of_nodes(weights.size())
  , iterable_nodes(weights.size())
  , likely_nodes_map()
  , max_number_of_nodes(1)
  , minimal_mesh_level(std::numeric_limits<int>::max())
  , highest_level_searched(-1)
  , max_cell_diameter(0.)
  , number_of_ignored_levels(2)
  , nodes_id(weights.size())
  , weights(weights)
  , nodes_positions(nodes)
  , support_radii(support_radii)
  , basis_functions(basis_functions)
{
  std::iota(std::begin(nodes_id), std::end(nodes_id), 0);
  iterable_nodes = nodes_id;
  initialize_bounding_box();
  this->effective_radius = bounding_box->half_lengths.norm();
}

template <int dim>
RBFShape<dim>::RBFShape(const std::vector<double> &shape_arguments,
                        const Point<dim> &         position,
                        const Tensor<1, 3> &       orientation)
  : Shape<dim>(shape_arguments[shape_arguments.size() / (dim + 3)],
               position,
               orientation)
// The effective radius is extracted at the proper index from shape_arguments
{
  number_of_nodes = shape_arguments.size() / (dim + 3);
  weights.resize(number_of_nodes);
  support_radii.resize(number_of_nodes);
  basis_functions.resize(number_of_nodes);
  nodes_positions.resize(number_of_nodes);
  nodes_id.resize(number_of_nodes);
  std::iota(std::begin(nodes_id), std::end(nodes_id), 0);
  iterable_nodes         = nodes_id;
  max_number_of_nodes    = 1;
  minimal_mesh_level     = std::numeric_limits<int>::max();
  highest_level_searched = -1;
  max_cell_diameter      = 0.;
  number_of_ignored_levels = 2;

  for (size_t n_i = 0; n_i < number_of_nodes; n_i++)
    {
      weights[n_i]         = shape_arguments[0 * number_of_nodes + n_i];
      support_radii[n_i]   = shape_arguments[1 * number_of_nodes + n_i];
      basis_functions[n_i] = static_cast<enum RBFShape<dim>::RBFBasisFunction>(
        round(shape_arguments[2 * number_of_nodes + n_i]));
      nodes_positions[n_i][0] = shape_arguments[3 * number_of_nodes + n_i];
      nodes_positions[n_i][1] = shape_arguments[4 * number_of_nodes + n_i];
      nodes_positions[n_i][2] = shape_arguments[5 * number_of_nodes + n_i];
    }
  initialize_bounding_box();
  this->effective_radius = bounding_box->half_lengths.norm();
}

template <int dim>
double
RBFShape<dim>::value_with_cell_guess(
  const Point<dim> &                                   evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  prepare_iterable_nodes(cell);
  double value = this->value(evaluation_point);
  reset_iterable_nodes(cell);
  return value;
}

template <int dim>
Tensor<1, dim>
RBFShape<dim>::gradient_with_cell_guess(
  const Point<dim> &                                   evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  prepare_iterable_nodes(cell);
  Tensor<1, dim> gradient = this->gradient(evaluation_point);
  reset_iterable_nodes(cell);
  return gradient;
}

template <int dim>
double
RBFShape<dim>::value(const Point<dim> &evaluation_point,
                     const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  double bounding_box_distance = bounding_box->value(centered_point);
  double value                 = std::max(bounding_box_distance, 0.0);

  double normalized_distance, basis;

  // Algorithm inspired by Optimad Bitpit. https://github.com/optimad/bitpit
  for (const size_t &node_id : iterable_nodes)
    {
      normalized_distance = (centered_point - nodes_positions[node_id]).norm() /
                            support_radii[node_id];
      basis =
        evaluate_basis_function(basis_functions[node_id], normalized_distance);
      value += basis * weights[node_id];
    }
  return value;
}

template <int dim>
Tensor<1, dim>
RBFShape<dim>::gradient(const Point<dim> &evaluation_point,
                        const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  double         bounding_box_distance = bounding_box->value(centered_point);
  Tensor<1, dim> bounding_box_gradient = bounding_box->gradient(centered_point);
  if (bounding_box_distance >
      *std::max_element(support_radii.begin(), support_radii.end()))
    return bounding_box_gradient;

  double     distance, normalized_distance;
  Point<dim> relative_position{};
  // We use chain derivation to express the gradient of the basis in regards to
  // the position d(basis)/dx = d(basis)              /d(normalized_distance)
  //             * d(normalized_distance)/d(distance)
  //             * d(distance)           /dx
  double         dbasis_drnorm_derivative;
  double         drnorm_dr_derivative;
  Tensor<1, dim> dr_dx_derivative{};
  Tensor<1, dim> gradient{};

  if (iterable_nodes.size() > 0)
    {
      for (const size_t &node_id : iterable_nodes)
        {
          // Calculation of the dr/dx
          relative_position   = centered_point - nodes_positions[node_id];
          distance            = (relative_position).norm();
          normalized_distance = distance / support_radii[node_id];
          if (distance > 0.0)
            dr_dx_derivative = relative_position / distance;
          else
            for (int d = 0; d < dim; d++)
              // Can be proved by taking the limit (definition of a derivative)
              dr_dx_derivative[d] = 1.0;
          // Calculation of the dr_norm/dr
          drnorm_dr_derivative = 1.0 / support_radii[node_id];
          // Calculation of the d(basis)/dr
          dbasis_drnorm_derivative =
            evaluate_basis_function_derivative(basis_functions[node_id],
                                               normalized_distance);
          // Sum
          gradient += dbasis_drnorm_derivative * drnorm_dr_derivative *
                      dr_dx_derivative * weights[node_id];
        }
    }
  else
    gradient = bounding_box_gradient;
  return gradient;
}

template <int dim>
std::shared_ptr<Shape<dim>>
RBFShape<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<RBFShape<dim>>(this->support_radii,
                                    this->basis_functions,
                                    this->weights,
                                    this->nodes_positions,
                                    this->position,
                                    this->orientation);
  return copy;
}

template <int dim>
double
RBFShape<dim>::displaced_volume(const double fluid_density)
{
  return bounding_box->displaced_volume(fluid_density);
}

template <int dim>
void
RBFShape<dim>::initialize_bounding_box()
{
  Point<dim>     high_bounding_point{};
  Point<dim>     low_bounding_point{};
  Point<dim>     bounding_box_center{};
  Tensor<1, dim> half_lengths = Tensor<1, dim>();
  for (int d = 0; d < dim; d++)
    {
      high_bounding_point[d] = std::numeric_limits<double>::lowest();
      low_bounding_point[d]  = std::numeric_limits<double>::max();
      for (size_t i = 0; i < number_of_nodes; i++)
        {
          low_bounding_point[d] =
            std::min(low_bounding_point[d], nodes_positions[i][d]);
          high_bounding_point[d] =
            std::max(high_bounding_point[d], nodes_positions[i][d]);
        }
      bounding_box_center[d] =
        0.5 * (low_bounding_point[d] + high_bounding_point[d]);
      half_lengths[d] = 0.5 * (high_bounding_point[d] - low_bounding_point[d]);
    }
  bounding_box = std::make_shared<Rectangle<dim>>(half_lengths,
                                                  bounding_box_center,
                                                  Tensor<1, 3>());
}

template <int dim>
void
RBFShape<dim>::determine_likely_nodes_for_one_cell(
  const typename DoFHandler<dim>::cell_iterator &cell,
  const Point<dim>                               support_point)
{
  // We exit the function immediately if the cell is already in the map
  if (likely_nodes_map.find(cell) != likely_nodes_map.end())
    return;

  bool parent_found  = false;
  minimal_mesh_level = std::min(minimal_mesh_level, cell->level());
  std::vector<size_t> temporary_iterable_nodes;
  if (cell->level() > minimal_mesh_level)
    try
      {
        const auto cell_parent     = cell->parent();
        const auto parent_iterator = likely_nodes_map.find(cell_parent);
        if (parent_iterator != likely_nodes_map.end())
          {
            parent_found = true;
            temporary_iterable_nodes.swap(iterable_nodes);
            iterable_nodes.swap(*parent_iterator->second);
          }
      }
    catch (TriaAccessorExceptions::ExcCellHasNoParent())
      {}

  double     distance, max_distance;
  double     cell_diameter = cell->diameter();
  Point<dim> centered_support_point;

  likely_nodes_map[cell] = std::make_shared<std::vector<size_t>>();
  likely_nodes_map[cell]->reserve(max_number_of_nodes);
  for (const auto &node_id : iterable_nodes)
    {
      centered_support_point = this->align_and_center(support_point);
      distance = (centered_support_point - nodes_positions[node_id]).norm();
      // We check if the distance is lower than 1 cell diagonal, since we
      // only check the distance with 1 support point, added to the support
      // radius
      max_distance = 2.*max_cell_diameter + support_radii[node_id];
      if (distance < max_distance)
        likely_nodes_map[cell]->push_back(node_id);
    }
  max_number_of_nodes =
    std::max(max_number_of_nodes, likely_nodes_map[cell]->size());
  if (cell->level() > minimal_mesh_level)
    if (parent_found)
      {
        const auto cell_parent     = cell->parent();
        const auto parent_iterator = likely_nodes_map.find(cell_parent);
        iterable_nodes.swap(*parent_iterator->second);
        temporary_iterable_nodes.swap(iterable_nodes);
      }
}

template <int dim>
void
RBFShape<dim>::update_precalculations(DoFHandler<dim> &dof_handler,
                                      std::shared_ptr<Mapping<dim>> /*mapping*/)
{
  int            maximal_level    = dof_handler.get_triangulation().n_levels();
  const MPI_Comm mpi_communicator = dof_handler.get_communicator();

  for (int level = highest_level_searched + 1; level < maximal_level; level++)
    {
      double local_max_cell_diameter = 0.;

      const auto &active_cell_iterator = dof_handler.active_cell_iterators();
      for (const auto &cell : active_cell_iterator)
        local_max_cell_diameter =
          std::max(local_max_cell_diameter, cell->diameter());
      const auto &prep_cell_iterator =
        dof_handler.cell_iterators_on_level(level);
      for (const auto &cell : prep_cell_iterator)
        local_max_cell_diameter =
          std::max(local_max_cell_diameter, cell->diameter());

      max_cell_diameter =
        Utilities::MPI::max(local_max_cell_diameter, mpi_communicator);

     // bool skip = (level > maximal_level - 1 - TEST_UNSEEKED_LEVELS);
      //std::cout<<"level "<<level<< " is skipped : " << skip << std::endl;

      const auto &cell_iterator = dof_handler.cell_iterators_on_level(level);
      for (const auto &cell : cell_iterator)
        {
          // We first check if we are in the zone where we simply assume that
          // children have the same likely nodes as their parents
          if (level > maximal_level - 1 - number_of_ignored_levels)
            {
              likely_nodes_map[cell] = likely_nodes_map[cell->parent()];
              continue;
            }

          if (level == maximal_level)
            if (!cell->is_locally_owned())
              continue;
          if (!cell->is_artificial_on_level())
            determine_likely_nodes_for_one_cell(cell, cell->barycenter());
        }
      highest_level_searched = level;

      // We remove unnecessary cells
      for (auto it = likely_nodes_map.cbegin(); it != likely_nodes_map.cend();)
        {
          auto cell = it->first;
          bool cell_still_needed =
            cell->is_active() ||
            (cell->level() > level - 1 - number_of_ignored_levels);
          if (!cell_still_needed || cell->is_artificial_on_level())
            {
              likely_nodes_map.erase(it++);
            }
          else
            ++it;
        }
    }
}


template <int dim>
double
RBFShape<dim>::evaluate_basis_function(const RBFBasisFunction basis_function,
                                       const double           distance) const
{
  switch (basis_function)
    {
      case RBFBasisFunction::WENDLANDC2:
        return RBFShape<dim>::wendlandc2(distance);
      case RBFBasisFunction::LINEAR:
        return RBFShape<dim>::linear(distance);
      case RBFBasisFunction::GAUSS90:
        return RBFShape<dim>::gauss90(distance);
      case RBFBasisFunction::GAUSS95:
        return RBFShape<dim>::gauss95(distance);
      case RBFBasisFunction::GAUSS99:
        return RBFShape<dim>::gauss99(distance);
      case RBFBasisFunction::C1C0:
        return RBFShape<dim>::c1c0(distance);
      case RBFBasisFunction::C2C0:
        return RBFShape<dim>::c2c0(distance);
      case RBFBasisFunction::C0C1:
        return RBFShape<dim>::c0c1(distance);
      case RBFBasisFunction::C1C1:
        return RBFShape<dim>::c1c1(distance);
      case RBFBasisFunction::C2C1:
        return RBFShape<dim>::c2c1(distance);
      case RBFBasisFunction::C0C2:
        return RBFShape<dim>::c0c2(distance);
      case RBFBasisFunction::C1C2:
        return RBFShape<dim>::c1c2(distance);
      case RBFBasisFunction::C2C2:
        return RBFShape<dim>::c2c2(distance);
      case RBFBasisFunction::COS:
        return RBFShape<dim>::cos(distance);
      default:
        return RBFShape<dim>::linear(distance);
    }
}

template <int dim>
double
RBFShape<dim>::evaluate_basis_function_derivative(
  const RBFBasisFunction basis_function,
  const double           distance) const
{
  switch (basis_function)
    {
      case RBFBasisFunction::WENDLANDC2:
        return RBFShape<dim>::wendlandc2_derivative(distance);
      case RBFBasisFunction::LINEAR:
        return RBFShape<dim>::linear_derivative(distance);
      case RBFBasisFunction::GAUSS90:
        return RBFShape<dim>::gauss90_derivative(distance);
      case RBFBasisFunction::GAUSS95:
        return RBFShape<dim>::gauss95_derivative(distance);
      case RBFBasisFunction::GAUSS99:
        return RBFShape<dim>::gauss99_derivative(distance);
      case RBFBasisFunction::C1C0:
        return RBFShape<dim>::c1c0_derivative(distance);
      case RBFBasisFunction::C2C0:
        return RBFShape<dim>::c2c0_derivative(distance);
      case RBFBasisFunction::C0C1:
        return RBFShape<dim>::c0c1_derivative(distance);
      case RBFBasisFunction::C1C1:
        return RBFShape<dim>::c1c1_derivative(distance);
      case RBFBasisFunction::C2C1:
        return RBFShape<dim>::c2c1_derivative(distance);
      case RBFBasisFunction::C0C2:
        return RBFShape<dim>::c0c2_derivative(distance);
      case RBFBasisFunction::C1C2:
        return RBFShape<dim>::c1c2_derivative(distance);
      case RBFBasisFunction::C2C2:
        return RBFShape<dim>::c2c2_derivative(distance);
      case RBFBasisFunction::COS:
        return RBFShape<dim>::cosinus_derivative(distance);
      default:
        return RBFShape<dim>::linear_derivative(distance);
    }
}

template <int dim>
void
RBFShape<dim>::prepare_iterable_nodes(
  const typename DoFHandler<dim>::active_cell_iterator cell)
{
  // Here we check if the likely nodes have been identified
  auto iterator = likely_nodes_map.find(cell);
  if (iterator != likely_nodes_map.end())
    iterable_nodes.swap(*likely_nodes_map[cell]);
  else
    iterable_nodes.swap(nodes_id);
}

template <int dim>
void
RBFShape<dim>::reset_iterable_nodes(
  const typename DoFHandler<dim>::active_cell_iterator cell)
{
  // Here we check if the likely nodes have been identified
  auto iterator = likely_nodes_map.find(cell);
  if (iterator != likely_nodes_map.end())
    iterable_nodes.swap(*likely_nodes_map[cell]);
  else
    iterable_nodes.swap(nodes_id);
}

template <int dim>
double
Cylinder<dim>::value(const Point<dim> &evaluation_point,
                     const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);


  double p_radius    = std::pow(centered_point[0] * centered_point[0] +
                               centered_point[1] * centered_point[1],
                             0.5);
  double radius_diff = p_radius - radius;
  double h_diff      = abs(centered_point[2]) - half_length;

  if (radius_diff > 0 && h_diff > 0)
    return std::pow(radius_diff * radius_diff + h_diff * h_diff, 0.5);
  else if (radius_diff <= 0 && h_diff > 0)
    return h_diff;
  else if (radius_diff > 0 && h_diff <= 0)
    return radius_diff;

  return std::max(radius_diff, h_diff);
}

template <int dim>
std::shared_ptr<Shape<dim>>
Cylinder<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy = std::make_shared<Cylinder<dim>>(
    this->radius, this->half_length, this->position, this->orientation);
  return copy;
}

template <int dim>
double
Cylinder<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  double solid_volume = PI * radius * radius * half_length * 2;

  return solid_volume;
}

template <int dim>
double
CylindricalTube<dim>::value(const Point<dim> &evaluation_point,
                            const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);


  // external cylinder
  double level_set_of_cylinder_hallow = 0;
  double p_radius      = std::pow(centered_point[0] * centered_point[0] +
                               centered_point[1] * centered_point[1],
                             0.5);
  double radius_diff_o = p_radius - (radius + rectangular_base / 2);
  double radius_diff_i = p_radius - (radius - rectangular_base / 2);
  double h_diff_o      = abs(centered_point[2] - height / 2) - height / 2;

  if (radius_diff_o > 0 && h_diff_o > 0)
    level_set_of_cylinder_hallow =
      std::pow(radius_diff_o * radius_diff_o + h_diff_o * h_diff_o, 0.5);
  else if (radius_diff_o > 0 && h_diff_o <= 0)
    level_set_of_cylinder_hallow = radius_diff_o;
  else if (radius_diff_o <= 0 && radius_diff_i > 0 && h_diff_o > 0)
    level_set_of_cylinder_hallow = h_diff_o;
  else if (radius_diff_i <= 0 && h_diff_o > 0)
    level_set_of_cylinder_hallow =
      std::pow(radius_diff_i * radius_diff_i + h_diff_o * h_diff_o, 0.5);
  else if (radius_diff_i <= 0 && h_diff_o <= 0)
    level_set_of_cylinder_hallow = -radius_diff_i;
  else
    level_set_of_cylinder_hallow =
      std::max(std::max(radius_diff_o, h_diff_o), -radius_diff_i);

  return level_set_of_cylinder_hallow;
}

template <int dim>
std::shared_ptr<Shape<dim>>
CylindricalTube<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<CylindricalTube<dim>>(this->radius + rectangular_base / 2,
                                           this->radius - rectangular_base / 2,
                                           this->height,
                                           this->position,
                                           this->orientation);
  return copy;
}

template <int dim>
double
CylindricalTube<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  double solid_volume = height * PI *
                        ((this->radius + rectangular_base / 2) *
                           (this->radius + rectangular_base / 2) -
                         (this->radius - rectangular_base / 2) *
                           (this->radius - rectangular_base / 2));

  return solid_volume;
}


template <int dim>
double
CylindricalHelix<dim>::value(const Point<dim> &evaluation_point,
                             const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);

  // The evaluation proposed in this function is exact close to the helix but
  // far away above or bellow the helix some error may happen since the wrong
  // helix loop may be evaluated.

  // Distance to the center of the helix
  // First we find a good initial guess for the non linear resolution.
  double phase = std::atan2(centered_point[1], centered_point[0]);
  if (phase != phase)
    phase = 0;
  if (phase < 0)
    phase = phase + 2 * numbers::PI;

  // Some calculation to find a good initial guess
  double nb_turns     = height / pitch;
  double phase_at_top = (nb_turns - std::floor(nb_turns)) * 2 * numbers::PI;
  double x            = centered_point[2] - phase * pitch / (2 * numbers::PI);

  // Initialize the parametric variable t for the helix
  double t               = 0;
  double t_initial_guess = 0;

  // We select two initial guess for which we do the full computation. Both
  // point have the same phase as the point we are evaluating one is above on
  // the helix the other bellow. We keep the closest one.
  t_initial_guess = std::floor(x / pitch) * 2 * numbers::PI + phase;

  // Solve the non-linear equation to find the point on the helix with the first
  // guess
  double level_set_tube = 0;
  t                     = t_initial_guess;
  double residual = -(radius * cos(t) - centered_point[0]) * sin(t) * radius +
                    (radius * sin(t) - centered_point[1]) * cos(t) * radius +
                    (pitch / (2 * numbers::PI) * t - centered_point[2]) *
                      pitch / (2 * numbers::PI);

  unsigned int i = 0;

  // Newton iterations.
  while (abs(residual) > 1e-8 && i < 20)
    {
      residual = -(radius * cos(t) - centered_point[0]) * sin(t) * radius +
                 (radius * sin(t) - centered_point[1]) * cos(t) * radius +
                 (pitch / (2 * numbers::PI) * t - centered_point[2]) * pitch /
                   (2 * numbers::PI);
      double dr_dt = centered_point[0] * cos(t) * radius +
                     centered_point[1] * sin(t) * radius +
                     (pitch / (2 * numbers::PI)) * (pitch / (2 * numbers::PI));
      t = t - residual / dr_dt;
      i += 1;
    }
  // Cap the value of the parametric variable t and find the point.
  if (t > nb_turns * numbers::PI * 2)
    {
      t = nb_turns * numbers::PI * 2;
      Point<dim> point_on_helix;
      point_on_helix[0] = radius * cos(t);
      point_on_helix[1] = radius * sin(t);
      point_on_helix[2] = t * pitch / (2 * numbers::PI);
      level_set_tube    = (centered_point - point_on_helix).norm();
    }
  else if (t < 0)
    {
      t = 0;
      Point<dim> point_on_helix;
      point_on_helix[0] = radius * cos(t);
      point_on_helix[1] = radius * sin(t);
      point_on_helix[2] = t * pitch / (2 * numbers::PI);

      level_set_tube = (centered_point - point_on_helix).norm();
    }
  else
    {
      Point<dim> point_on_helix;
      point_on_helix[0] = radius * cos(t);
      point_on_helix[1] = radius * sin(t);
      point_on_helix[2] = t * pitch / (2 * numbers::PI);
      level_set_tube = (centered_point - point_on_helix).norm() - radius_disk;
    }

  // repeat for second guess
  double level_set_tube_2 = 0;
  t                       = t_initial_guess + 2 * numbers::PI;
  residual = -(radius * cos(t) - centered_point[0]) * sin(t) * radius +
             (radius * sin(t) - centered_point[1]) * cos(t) * radius +
             (pitch / (2 * numbers::PI) * t - centered_point[2]) * pitch /
               (2 * numbers::PI);

  i = 0;
  while (abs(residual) > 1e-8 && i < 20)
    {
      residual = -(radius * cos(t) - centered_point[0]) * sin(t) * radius +
                 (radius * sin(t) - centered_point[1]) * cos(t) * radius +
                 (pitch / (2 * numbers::PI) * t - centered_point[2]) * pitch /
                   (2 * numbers::PI);
      double dr_dt = centered_point[0] * cos(t) * radius +
                     centered_point[1] * sin(t) * radius +
                     (pitch / (2 * numbers::PI)) * (pitch / (2 * numbers::PI));
      t = t - residual / dr_dt;
      i += 1;
    }

  if (t > nb_turns * numbers::PI * 2)
    {
      t = nb_turns * numbers::PI * 2;
      Point<dim> point_on_helix;
      point_on_helix[0] = radius * cos(t);
      point_on_helix[1] = radius * sin(t);
      point_on_helix[2] = t * pitch / (2 * numbers::PI);
      level_set_tube_2  = (centered_point - point_on_helix).norm();
    }
  else if (t < 0)
    {
      t = 0;
      Point<dim> point_on_helix;
      point_on_helix[0] = radius * cos(t);
      point_on_helix[1] = radius * sin(t);
      point_on_helix[2] = t * pitch / (2 * numbers::PI);

      level_set_tube_2 = (centered_point - point_on_helix).norm();
    }
  else
    {
      Point<dim> point_on_helix;
      point_on_helix[0] = radius * cos(t);
      point_on_helix[1] = radius * sin(t);
      point_on_helix[2] = t * pitch / (2 * numbers::PI);
      level_set_tube_2 = (centered_point - point_on_helix).norm() - radius_disk;
    }

  // Keep the best guess
  level_set_tube = std::min(level_set_tube, level_set_tube_2);

  // Cap the helix with a plane at each end. Cap at the base.
  Point<dim>     point_at_base;
  Tensor<1, dim> vector_at_base;
  point_at_base[0]  = radius;
  point_at_base[1]  = 0;
  point_at_base[2]  = 0;
  vector_at_base[0] = 0;
  vector_at_base[1] = -radius;
  vector_at_base[2] = -pitch / (2 * numbers::PI);

  double p_radius_cap_base =
    ((centered_point - point_at_base) -
     scalar_product((centered_point - point_at_base), vector_at_base) /
       vector_at_base.norm_square() * vector_at_base)
      .norm();
  double radial_distance_cap_base = p_radius_cap_base - radius_disk;
  double h_dist_from_cap_0 =
    (scalar_product((centered_point - point_at_base), vector_at_base) /
     vector_at_base.norm_square() * vector_at_base)
      .norm();

  double dist_from_cap = 0;
  if (radial_distance_cap_base > 0)
    dist_from_cap =
      std::pow(h_dist_from_cap_0 * h_dist_from_cap_0 +
                 radial_distance_cap_base * radial_distance_cap_base,
               0.5);
  else
    dist_from_cap = h_dist_from_cap_0;

  // Cap the helix with a plane at each end. The cap at the top.
  Point<dim>     point_at_top;
  Tensor<1, dim> vector_at_top;
  point_at_top[0]  = std::cos(phase_at_top) * radius;
  point_at_top[1]  = std::sin(phase_at_top) * radius;
  point_at_top[2]  = height;
  vector_at_top[0] = -point_at_top[1];
  vector_at_top[1] = point_at_top[0];
  vector_at_top[2] = pitch / (2 * numbers::PI);

  double p_radius_cap_top =
    ((centered_point - point_at_top) -
     scalar_product((centered_point - point_at_top), vector_at_top) /
       vector_at_top.norm_square() * vector_at_top)
      .norm();
  double radial_distance_cap_top = p_radius_cap_top - radius_disk;
  double h_dist_from_cap_top =
    (scalar_product((centered_point - point_at_top), vector_at_top) /
     vector_at_top.norm_square() * vector_at_top)
      .norm();

  double dist_from_cap_top = 0;
  if (radial_distance_cap_top > 0)
    dist_from_cap_top =
      std::pow(h_dist_from_cap_top * h_dist_from_cap_top +
                 radial_distance_cap_top * radial_distance_cap_top,
               0.5);
  else
    dist_from_cap_top = h_dist_from_cap_top;

  // Do the union of the helix and the cap.
  double level_set = 0;
  if (level_set_tube > 0)
    level_set =
      std::min(std::min(level_set_tube, dist_from_cap), dist_from_cap_top);
  else
    level_set =
      std::max(std::max(level_set_tube, -dist_from_cap_top), -dist_from_cap);

  return level_set;
}

template <int dim>
std::shared_ptr<Shape<dim>>
CylindricalHelix<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<CylindricalHelix<dim>>(this->radius,
                                            this->radius_disk,
                                            this->height,
                                            this->pitch,
                                            this->position,
                                            this->orientation);
  return copy;
}

template <int dim>
double
CylindricalHelix<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  double solid_volume = 1;

  return solid_volume;
}

template class Sphere<2>;
template class Sphere<3>;
template class Rectangle<2>;
template class Rectangle<3>;
template class Ellipsoid<2>;
template class Ellipsoid<3>;
template class Torus<3>;
template class Cone<3>;
template class Cylinder<3>;
template class CylindricalTube<3>;
template class CylindricalHelix<3>;
template class CutHollowSphere<3>;
template class DeathStar<3>;
template class CompositeShape<2>;
template class CompositeShape<3>;
template class RBFShape<2>;
template class RBFShape<3>;