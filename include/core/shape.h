/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * ---------------------------------------------------------------------*/

#ifndef lethe_shape_h
#define lethe_shape_h

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function.h>

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
#else
#  include <deal.II/base/function_signed_distance.h>
#endif


#include <deal.II/physics/transformations.h>

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
   * @brief A general constructor for the Shapes
   */
  Shape(double radius)
    : AutoDerivativeFunction<dim>(1e-8)
  {
    this->effective_radius = radius;
  }

  virtual ~Shape()
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point p
   * Most levelset functions implemented come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   */
  virtual double
  value(const Point<dim> & p,
        const unsigned int component = 0) const override = 0;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  virtual std::shared_ptr<Shape<dim>>
  static_copy() const = 0;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   */
  virtual double
  volume(const double fluid_density)
  {
    StandardExceptions::ExcNotImplemented();
    return 1.;
  }

  /**
   * @brief
   * Returns the centered and aligned point to evaluate properly the levelset
   * afterwards
   */
  Point<dim>
  real_to_centered(const Point<dim> &evaluation_pt) const
  {
    Point<dim> corrected_pt;

    Point<dim> cor = position + cor_offset;
    // Translation and rotation to standard position and orientation for
    // distance calculations
    Point<dim> rotated_pt;
    Point<dim> translated_pt;
    // Rotation from the solid orientation
    // Angular position around x, y and z axis
    Tensor<1, 3> theta = orientation;
    Point<dim>   centralized_point;
    centralized_point              = evaluation_pt - cor;
    Point<dim> centralized_rotated = centralized_point;

    // Selection of the first axis around which to rotate:
    // x -> 0, y -> 1, z -> 2
    // In 2D, rotation is only possible around the z axis
    unsigned int first_rotation_axis = 0;
    if (dim == 2)
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

    else if (dim == 3)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            if (std::abs(theta[i]) > 1e-10)
              {
                Point<3> axis;
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
    rotated_pt = centralized_rotated + cor;
    // Translate
    translated_pt = rotated_pt - position;

    return translated_pt;
  }

  // Position of the center of the Shape. It doesn't always correspond to the
  // center of mass
  Point<dim> position;
  // The offset of the center of rotation in relation to the position
  Point<dim> cor_offset;
  // The solid orientation, which is defined as the sequential rotation around
  // the axes x->y->z by each of the tensor components, in radian
  Tensor<1, 3> orientation;
  // Effective radius used for crown refinement
  double effective_radius;
};


template <int dim>
class Sphere : public Shape<dim>
{
public:
  Sphere<dim>(double radius)
    : Shape<dim>(radius)
  {}

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
    return p.distance(this->position) - this->effective_radius;
#else
    Functions::SignedDistance::Sphere<dim> sphere_function(
      this->position, this->effective_radius);
    return sphere_function.value(p);
#endif
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<Sphere<dim>>(this->effective_radius);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

  Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const override
  {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
    const Tensor<1, dim> center_to_point = p - this->position;
    const Tensor<1, dim> grad = center_to_point / center_to_point.norm();
    return grad;
#else
    Functions::SignedDistance::Sphere<dim> sphere_function(
      this->position, this->effective_radius);
    return sphere_function.gradient(p);
#endif
  }

  double
  volume(const double fluid_density) override
  {
    double solid_volume;
    using numbers::PI;
    if (dim == 2)
      solid_volume =
        this->effective_radius * this->effective_radius * PI * fluid_density;

    else if (dim == 3)
      solid_volume = 4.0 / 3.0 * this->effective_radius *
                     this->effective_radius * this->effective_radius * PI;
    return solid_volume;
  }
};

template <int dim>
class Rectangle : public Shape<dim>
{
public:
  Rectangle<dim>(Tensor<1, 3> half_lengths)
    : Shape<dim>(half_lengths.norm())
  {
    this->half_lengths = half_lengths;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    double     levelset;
    Point<dim> centered_pt = this->real_to_centered(p);

    Point<dim> abs_p;
    Point<dim> half_lengths_dim;
    for (unsigned int i = 0; i < dim; ++i)
      {
        abs_p[i]            = std::abs(centered_pt[i]);
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
    levelset     = max_q_0.norm() + std::min(max_q, 0.);

    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<Rectangle<dim>>(this->half_lengths);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  Tensor<1, 3> half_lengths;
};

template <int dim>
class Ellipsoid : public Shape<dim>
{
public:
  Ellipsoid<dim>(Tensor<1, 3> radii)
    : Shape<dim>(radii.norm())
  {
    this->radii = radii;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    double     levelset;
    Point<dim> centered_pt = this->real_to_centered(p);

    Point<dim> v_k0;
    Point<dim> v_k1;
    for (unsigned int i = 0; i < dim; ++i)
      {
        // Ellipsoid parameters[0->2]:= radii x, y, z
        v_k0[i] = centered_pt[i] / radii[i];
        v_k1[i] = centered_pt[i] / (radii[i] * radii[i]);
      }
    double k0 = v_k0.norm();
    double k1 = v_k1.norm();
    levelset  = k0 * (k0 - 1.) / k1;

    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<Ellipsoid<dim>>(this->radii);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  Tensor<1, 3> radii;
};

template <int dim>
class Torus : public Shape<dim>
{
public:
  Torus<dim>(double ring_radius, double ring_thickness)
    : Shape<dim>(ring_thickness)
  {
    this->ring_radius    = ring_radius;
    this->ring_thickness = ring_thickness;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    AssertDimension(dim, 3);

    double     levelset;
    Point<dim> centered_pt = this->real_to_centered(p);

    Point<2> p_xz({centered_pt[0], centered_pt[2]});
    Point<2> q({p_xz.norm() - ring_radius, centered_pt[1]});
    levelset = q.norm() - ring_thickness;

    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<Torus<dim>>(this->ring_radius, ring_thickness);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  double ring_radius;
  double ring_thickness;
};

template <int dim>
class Cone : public Shape<dim>
{
public:
  Cone<dim>(double tan_theta, double height)
    : Shape<dim>(height)
  {
    this->tan_theta = tan_theta;
    this->height    = height;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    AssertDimension(dim, 3);

    double     levelset;
    Point<dim> centered_pt = this->real_to_centered(p);

    // For a cone, the parameters are tan(base angle) and height
    Point<2> q({height * tan_theta, -height});
    Point<2> p_xz({centered_pt[0], centered_pt[2]});
    Point<2> w({p_xz.norm(), centered_pt[1]});
    double   dot_w_q = scalar_product<1, 2, double>(w, q);
    double   dot_q_q = scalar_product<1, 2, double>(q, q);
    Point<2> a;
    a = w - q * std::clamp(dot_w_q / dot_q_q, 0., 1.);
    Point<2> b_intermediate1({std::clamp(w[0] / q[0], 0., 1.), 1.});
    Point<2> b_intermediate2(
      {q[0] * b_intermediate1[0], q[1] * b_intermediate1[1]});
    Point<2> b;
    b        = w - b_intermediate2;
    double k = (q[1] > 0) ? 1 : ((q[1] < 0) ? -1 : 0);
    double d = std::min(scalar_product<1, 2, double>(a, a),
                        scalar_product<1, 2, double>(b, b));
    double s = std::max(k * (w[0] * q[1] - w[1] * q[0]), k * (w[1] - q[1]));

    levelset = sqrt(d) * ((s > 0) ? 1 : ((s < 0) ? -1 : 0));

    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<Cone<dim>>(this->tan_theta, this->height);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  double tan_theta;
  double height;
};

template <int dim>
class CutHollowSphere : public Shape<dim>
{
public:
  CutHollowSphere<dim>(double r, double h, double t)
    : Shape<dim>(r)
  {
    this->r = r;
    this->h = h;
    this->t = t;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    AssertDimension(dim, 3);

    double     levelset;
    Point<dim> centered_pt = this->real_to_centered(p);

    double   w = sqrt(r * r - h * h);
    Point<2> p_xz({centered_pt[0], centered_pt[2]});
    Point<2> q({p_xz.norm(), centered_pt[1]});

    if (h * q[0] < w * q[1])
      {
        Point<2> wh({w, h});
        levelset = (q - wh).norm();
      }
    else
      {
        levelset = std::abs(q.norm() - r) - t;
      }

    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<CutHollowSphere<dim>>(this->r, this->h, this->t);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  double r;
  double h;
  double t;
};

template <int dim>
class DeathStar : public Shape<dim>
{
public:
  DeathStar<dim>(double ra, double rb, double d)
    : Shape<dim>(ra)
  {
    this->ra = ra;
    this->rb = rb;
    this->d  = d;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    AssertDimension(dim, 3);

    double     levelset;
    Point<dim> centered_pt = this->real_to_centered(p);

    double a = (ra * ra - rb * rb + d * d) / (2. * d);
    double b = sqrt(std::max(ra * ra - a * a, 0.));

    Point<2> p_yz({centered_pt[1], centered_pt[2]});
    Point<2> corrected_p_2d({centered_pt[0], p_yz.norm()});
    if (corrected_p_2d[0] * b - corrected_p_2d[1] * a >
        d * std::max(b - corrected_p_2d[1], 0.))
      {
        Point<2> ab({a, b});
        levelset = (corrected_p_2d - ab).norm();
      }
    else
      {
        Point<2> d0({d, 0.});
        levelset = std::max(corrected_p_2d.norm() - ra,
                            -((corrected_p_2d - d0).norm() - rb));
      }

    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<DeathStar<dim>>(this->ra, this->rb, this->d);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  double ra;
  double rb;
  double d;
};

template <int dim>
class CompositeShape : public Shape<dim>
{
public:
  CompositeShape<dim>(std::vector<std::shared_ptr<Shape<dim>>> components)
    : Shape<dim>(1.)
  {
    this->components = components;
    double radius    = 0.;
    for (const std::shared_ptr<Shape<dim>> &elem : components)
      {
        radius = std::max(radius, elem->effective_radius);
      }
    this->effective_radius = radius;
  }

  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    double levelset = DBL_MAX;
    for (const std::shared_ptr<Shape<dim>> &elem : components)
      {
        levelset = std::min(elem->value(p), levelset);
      }
    return levelset;
  }

  std::shared_ptr<Shape<dim>>
  static_copy() const override
  {
    std::shared_ptr<Shape<dim>> copy =
      std::make_shared<CompositeShape<dim>>(this->components);
    copy->position    = this->position;
    copy->cor_offset  = this->cor_offset;
    copy->orientation = this->orientation;
    return copy;
  }

private:
  std::vector<std::shared_ptr<Shape<dim>>> components;
};

#endif // lethe_shape_h