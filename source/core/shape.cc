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

#include <cfloat>

template <int dim>
double
Shape<dim>::displaced_volume(const double /*fluid_density*/)
{
  StandardExceptions::ExcNotImplemented();
  return 1.;
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
            Physics::Transformations::Rotations::rotation_matrix_2d(theta[2]);

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
              axis[i] = 1.;
              Tensor<2, 3> rotation_matrix =
                Physics::Transformations::Rotations::rotation_matrix_3d(
                  axis, theta[i]);

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
void
Shape<dim>::set_position(const Point<dim> &position)
{
  this->position = position;
}

template <int dim>
void
Shape<dim>::set_orientation(const Tensor<1, 3> &orientation)
{
  this->orientation = orientation;
}

template <int dim>
Point<dim>
Shape<dim>::get_position()
{
  return position;
}

template <int dim>
Tensor<1, 3>
Shape<dim>::get_orientation()
{
  return orientation;
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
      max_q_0[i] = std::max(q[i], 0.);
    }
  double max_q = std::max(q[0], std::max(q[1], q[dim - 1]));
  return max_q_0.norm() + std::min(max_q, 0.);
}

template <int dim>
std::shared_ptr<Shape<dim>>
Rectangle<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<Rectangle<dim>>(half_lengths,
                                     this->position,
                                     this->orientation);
  return copy;
}

template <int dim>
double
Rectangle<dim>::displaced_volume(const double /*fluid_density*/)
{
  double solid_volume = 1.;
  for (unsigned int i = 0; i < dim; i++)
    {
      solid_volume = solid_volume * 2. * half_lengths[dim];
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
    std::make_shared<Ellipsoid<dim>>(radii, this->position, this->orientation);
  return copy;
}

template <int dim>
double
Ellipsoid<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  double solid_volume = PI * 4. / 3.;
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
    ring_radius, ring_thickness, this->position, this->orientation);
  return copy;
}

template <int dim>
double
Torus<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  return 2. * PI * PI * ring_radius * ring_thickness * ring_thickness;
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
  a = w - intermediate_q * std::clamp(dot_w_q / dot_q_q, 0., 1.);
  Point<2> b_intermediate1({std::clamp(w[0] / intermediate_q[0], 0., 1.), 1.});
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
    tan_base_angle, height, this->position, this->orientation);
  return copy;
}

template <int dim>
double
Cone<dim>::displaced_volume(const double /*fluid_density*/)
{
  using numbers::PI;
  return PI / 3. * base_radius * base_radius * height;
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
  std::shared_ptr<Shape<dim>> copy = std::make_shared<CutHollowSphere<dim>>(
    radius, cut_depth, shell_thickness, this->position, this->orientation);
  return copy;
}

template <int dim>
double
CutHollowSphere<dim>::displaced_volume(const double /*fluid_density*/)
{
  std::cout
    << "Warning: For a cut hollow sphere, the real volume will be lower than "
       "output."
    << std::endl;
  using numbers::PI;
  double small_radius = radius - shell_thickness;
  return 4. * PI / 3. *
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
      spheres_distance * std::max(intermediate_b - corrected_p_2d[1], 0.))
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
  std::shared_ptr<Shape<dim>> copy = std::make_shared<DeathStar<dim>>(
    radius, hole_radius, spheres_distance, this->position, this->orientation);
  return copy;
}

template <int dim>
double
DeathStar<dim>::displaced_volume(const double /*fluid_density*/)
{
  std::cout
    << "Warning: For a death star, the real volume will be lower than output."
    << std::endl;
  using numbers::PI;
  return 4. * PI / 3. * radius * radius * radius;
}

template <int dim>
double
CompositeShape<dim>::value(const Point<dim> &evaluation_point,
                           const unsigned int /*component*/) const
{
  double levelset = DBL_MAX;
  for (const std::shared_ptr<Shape<dim>> &elem : components)
    {
      levelset = std::min(elem->value(evaluation_point), levelset);
    }
  return levelset;
}

template <int dim>
std::shared_ptr<Shape<dim>>
CompositeShape<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<CompositeShape<dim>>(this->components);
  return copy;
}

template <int dim>
double
CompositeShape<dim>::displaced_volume(const double fluid_density)
{
  double solid_volume = 0;
  std::cout
    << "Warning: For composite shapes, the real volume may be bigger than output "
       "since intersections aren't considered in the calculation."
    << std::endl;
  for (const std::shared_ptr<Shape<dim>> &elem : components)
    {
      solid_volume += elem->displaced_volume(fluid_density);
    }
  return solid_volume;
}


template <int dim>
double
RBFShape<dim>::value(const Point<dim> &evaluation_point,
                     const unsigned int /*component*/) const
{
  Point<dim> centered_point = this->align_and_center(evaluation_point);


  double bounding_box_distance = bounding_box->value(centered_point);
  double value = bounding_box_distance > 0. ? bounding_box_distance : 0.;

  double       dist, basis;
  unsigned int number_of_nodes = weight.size();
  // Algorithm inspired by Optimad Bitpit. https://github.com/optimad/bitpit
  for (unsigned int i = 0; i < number_of_nodes; ++i)
    {
      dist  = (centered_point - nodes[i]).norm() / support_radius[i];
      basis = evaluate_basis_function(basis_function[i], dist);
      value += basis * weight[i];
    }
  return value;
}

template <int dim>
std::shared_ptr<Shape<dim>>
RBFShape<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<RBFShape<dim>>(support_radius,
                                    basis_function,
                                    weight,
                                    nodes,
                                    this->position,
                                    this->orientation);
  return copy;
}

template <int dim>
double
RBFShape<dim>::displaced_volume(const double fluid_density)
{
  std::cout
    << "Warning: For a RBF shape, the real volume will be lower than output."
    << std::endl;
  return bounding_box->displaced_volume(fluid_density);
}

/*!
 * Wendland C2 function
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::wendlandc2(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return std::pow(1. - dist, 4.) * (4. * dist + 1.);
    }
}

/*!
 * Linear function
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::linear(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1. - dist);
    }
}

/*!
 * Non compact gaussian function with 0.1 value at dist equal to 1
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::gauss90(double dist) const
{
  double eps = std::pow(-1.0 * std::log(0.1), 0.5);

  return std::exp(-1.0 * std::pow(dist * eps, 2.));
}

/*!
 * Non compact gaussian function with 0.05 value at dist equal to 1
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::gauss95(double dist) const
{
  double eps = std::pow(-1.0 * std::log(0.05), 0.5);

  return std::exp(-1.0 * std::pow(dist * eps, 2.));
}

/*!
 * Non compact gaussian function with 0.01 value at dist equal to 1
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::gauss99(double dist) const
{
  double eps = std::pow(-1.0 * std::log(0.01), 0.5);

  return std::exp(-1.0 * std::pow(dist * eps, 2.));
}

/*!
 * Polynomial function defined between 0,1. Preserve C1 continuity at dist=0, C0
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c1c0(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - std::pow(dist, 2.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C2 continuity at dist=0, C0
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c2c0(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - std::pow(dist, 3.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C0 continuity at dist=0, C1
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c0c1(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - 2.0 * dist + std::pow(dist, 2.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C1 continuity at dist=0, C1
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c1c1(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - 3.0 * std::pow(dist, 2.) + 2.0 * std::pow(dist, 3.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C2 continuity at dist=0, C1
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c2c1(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - 4.0 * std::pow(dist, 3.) + 3.0 * std::pow(dist, 4.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C0 continuity at dist=0, C2
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c0c2(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - 3.0 * dist + 3.0 * std::pow(dist, 2.) - std::pow(dist, 3.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C1 continuity at dist=0, C2
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c1c2(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - 6.0 * std::pow(dist, 2.) + 8.0 * std::pow(dist, 3.) -
              3.0 * std::pow(dist, 4.));
    }
}

/*!
 * Polynomial function defined between 0,1. Preserve C2 continuity at dist=0, C2
 * continuity at dist=1. At dist > 1 is 0.
 * @param[in] dist distance normalized with respect to support radius
 * @return rbf value
 */
template <int dim>
double
RBFShape<dim>::c2c2(double dist) const
{
  if (dist > 1.)
    {
      return 0.;
    }
  else
    {
      return (1.0 - 10.0 * std::pow(dist, 3.) + 15.0 * std::pow(dist, 4.) -
              6.0 * std::pow(dist, 5.));
    }
}

template <int dim>
double
RBFShape<dim>::evaluate_basis_function(const unsigned int basis_function_id,
                                       const double       distance) const
{
  double value;
  switch (basis_function_id)
    {
      case 1:
        value = RBFShape<dim>::wendlandc2(distance);
        break;
      case 2:
        value = RBFShape<dim>::linear(distance);
        break;
      case 3:
        value = RBFShape<dim>::gauss90(distance);
        break;
      case 4:
        value = RBFShape<dim>::gauss95(distance);
        break;
      case 5:
        value = RBFShape<dim>::gauss99(distance);
        break;
      case 6:
        value = RBFShape<dim>::c1c0(distance);
        break;
      case 7:
        value = RBFShape<dim>::c2c0(distance);
        break;
      case 8:
        value = RBFShape<dim>::c0c1(distance);
        break;
      case 9:
        value = RBFShape<dim>::c1c1(distance);
        break;
      case 10:
        value = RBFShape<dim>::c2c1(distance);
        break;
      case 11:
        value = RBFShape<dim>::c0c2(distance);
        break;
      case 12:
        value = RBFShape<dim>::c1c2(distance);
        break;
      case 13:
        value = RBFShape<dim>::c2c2(distance);
        break;
      default:
        value = RBFShape<dim>::linear(distance);
        break;
    }
  return value;
}

template class Sphere<2>;
template class Sphere<3>;
template class Rectangle<2>;
template class Rectangle<3>;
template class Ellipsoid<2>;
template class Ellipsoid<3>;
template class Torus<3>;
template class Cone<3>;
template class CutHollowSphere<3>;
template class DeathStar<3>;
template class CompositeShape<2>;
template class CompositeShape<3>;
template class RBFShape<2>;
template class RBFShape<3>;
