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

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#ifdef DEAL_II_WITH_OPENCASCADE
#  include <deal.II/opencascade/manifold_lib.h>
#  include <deal.II/opencascade/utilities.h>

#  include <BRepBuilderAPI_MakeVertex.hxx>
#  include <BRepClass3d_SolidClassifier.hxx>
#  include <BRepExtrema_DistShapeShape.hxx>
#  include <BRepGProp.hxx>
#  include <GProp_GProps.hxx>
#endif

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
#else
#  include <deal.II/base/function_signed_distance.h>
#endif

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
        const Point<dim> &  position,
        const Tensor<1, 3> &orientation)

    : AutoDerivativeFunction<dim>(1e-8)
    , effective_radius(radius)
    , position(position)
    , orientation(orientation)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point
   * Most levelset functions implemented come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  virtual double
  value(const Point<dim> & evaluation_point,
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
    const Point<dim> &                                   evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int                                   component = 0);

  /**
   * @brief Return the analytical gradient of the distance
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  virtual Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim> &                                   evaluation_point,
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
    const Point<dim> &                                    p,
    Point<dim> &                                          closest_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);
  virtual void
  closest_surface_point(const Point<dim> &p, Point<dim> &closest_point) const;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  virtual double
  displaced_volume(const double fluid_density);

  /**
   * @brief
   * Sets a new position for the shape
   *
   * @param The new position the shape will be placed at
   */
  inline virtual void
  set_position(const Point<dim> &new_position)
  {
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
   * current shape position and orientation, so that subsequent calculations for
   * the value function are made more easily; it abstracts a step that is
   * required in the value function for most shapes.
   *
   * Returns the centered and aligned point used on the levelset evaluation.
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
   * Returns the centered and aligned point used on the levelset evaluation in
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

  // Effective radius used for crown refinement
  double effective_radius;

  // The string contains additional information on the shape. This may refer to
  // the file type used to define the shape or any other information relative to
  // how the shape was defined.
  std::string additional_info_on_shape;

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
  Sphere<dim>(double              radius,
              const Point<dim> &  position,
              const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
  {
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
#else
    sphere_function =
      std::make_shared<Functions::SignedDistance::Sphere<dim>>(position,
                                                               radius);
#endif
  }

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
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
  gradient(const Point<dim> & evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

  void
  set_position(const Point<dim> &position) override;

  /**
   * @brief Return the manifold of the sphere as a spherical manifold center at the position of the shape.
   */
  virtual std::shared_ptr<Manifold<dim - 1, dim>>
  get_shape_manifold() override;



private:
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
#else
  std::shared_ptr<Functions::SignedDistance::Sphere<dim>> sphere_function;
#endif
};

/**
 * @class This class defines superquadric shapes. Their signed distance function
 * is:
 * \left|\frac{x}{a}\right|^r + \left|\frac{y}{b}\right|^s +
 * \left|\frac{z}{c}\right|^t - 1 = 0
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
  Superquadric<dim>(const Tensor<1, dim> half_lengths,
                    const Tensor<1, dim> exponents,
                    const double         epsilon,
                    const Point<dim> &   position,
                    const Tensor<1, 3> & orientation)
    : Shape<dim>(half_lengths.norm(), position, orientation)
    , half_lengths(half_lengths)
    , exponents(exponents)
    , epsilon(epsilon)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
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
    const Point<dim> &                                   evaluation_point,
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
  gradient(const Point<dim> & evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim> &                                   evaluation_point,
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
    const Point<dim> &                                    p,
    Point<dim> &                                          closest_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess) override;
  void
  closest_surface_point(const Point<dim> &p,
                        Point<dim> &      closest_point) const override;

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
  HyperRectangle<dim>(const Tensor<1, dim> &half_lengths,
                      const Point<dim> &    position,
                      const Tensor<1, 3> &  orientation)
    : Shape<dim>(half_lengths.norm(), position, orientation)
    , half_lengths(half_lengths)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  Ellipsoid<dim>(const Tensor<1, dim> &radii,
                 const Point<dim> &    position,
                 const Tensor<1, 3> &  orientation)
    : Shape<dim>(radii.norm(), position, orientation)
    , radii(radii)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  Torus<dim>(double              ring_radius,
             double              ring_thickness,
             const Point<dim> &  position,
             const Tensor<1, 3> &orientation)
    : Shape<dim>(ring_thickness, position, orientation)
    , ring_radius(ring_radius)
    , ring_thickness(ring_thickness)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  Cone<dim>(double              tan_base_angle,
            double              height,
            const Point<dim> &  position,
            const Tensor<1, 3> &orientation)
    : Shape<dim>(height, position, orientation)
    , tan_base_angle(tan_base_angle)
    , height(height)
    , base_radius(height / tan_base_angle)
    , intermediate_q({height * tan_base_angle, -height})
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  CutHollowSphere<dim>(double              radius,
                       double              cut_depth,
                       double              shell_thickness,
                       const Point<dim> &  position,
                       const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
    , radius(radius)
    , cut_depth(cut_depth)
    , shell_thickness(shell_thickness)
    , intermediate_w(sqrt(radius * radius - cut_depth * cut_depth))
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  DeathStar<dim>(double              radius,
                 double              hole_radius,
                 double              spheres_distance,
                 const Point<dim> &  position,
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
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

private:
  double radius;
  double hole_radius;
  double spheres_distance;

  double intermediate_a;
  double intermediate_b;
};

/**
 * @class This class was designed so that specific composite shapes could be
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
  CompositeShape<dim>(
    std::map<unsigned int, std::shared_ptr<Shape<dim>>> constituents,
    std::map<unsigned int,
             std::tuple<BooleanOperation, unsigned int, unsigned int>>
                        operations,
    const Point<dim> &  position,
    const Tensor<1, 3> &orientation)
    : Shape<dim>(0., position, orientation)
    , constituents(constituents)
    , operations(operations)
  {
    // Calculation of the effective radius
    for (auto const &[component_id, component] : constituents)
      {
        this->effective_radius =
          std::max(this->effective_radius, component->effective_radius);
      }
  }

  /**
   * @brief Constructs an assembly of shapes into a composite shape from a vector of shapes.
   * This constructor is mainly used for outputting multiple shapes with a
   * global levelset function defined as a union.
   * @param constituents_vector The shapes from which this composite sphere will be composed
   */
  CompositeShape<dim>(
    std::vector<std::shared_ptr<Shape<dim>>> constituents_vector,
    const Point<dim> &                       position,
    const Tensor<1, 3> &                     orientation)
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
    for (auto const &[constituent_id, constituent] : constituents)
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
  value(const Point<dim> & evaluation_point,
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
    const Point<dim> &                                   evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int /*component = 0*/) override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient(const Point<dim> & evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim> &                                   evaluation_point,
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
   * @param levels_not_precalculated the number of finer levels that won't be
   * precalculated
   */
  void
  update_precalculations(DoFHandler<dim> &  updated_dof_handler,
                         const unsigned int levels_not_precalculated);

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

private:
  // The members of this class are all the constituent and operations that are
  // to be performed to construct the composite shape
  // This map link all primitive constituents of the composite shape to an id
  std::map<unsigned int, std::shared_ptr<Shape<dim>>> constituents;
  // This map links all operations between primitive constituents or
  // intermediate constituents (resulting from each operation) to an id. The
  // unsigned integers correspond to the first and second ids of the shapes used
  // for an operation
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
  OpenCascadeShape<dim>(const std::string   file_name,
                        const Point<dim> &  position,
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
  value(const Point<dim> & evaluation_point,
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
    const Point<dim> &                                   evaluation_point,
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
  gradient(const Point<dim> & evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return the gradient of the distance function
   * @param evaluation_point The point at which the function will be evaluated
   * @param cell The cell that is likely to contain the evaluation point
   * @param component Not applicable
   */
  Tensor<1, dim>
  gradient_with_cell_guess(
    const Point<dim> &                                   evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
 * @class RBF Shapes express the signed distance function as a linear
 * combination of Radial Basis Functions (RBF), which have a defined support
 * radius and basis function. A collection of nodes and weights compose the
 * object. Outside of the domain covered by the nodes, the distance is computed
 * by using the distance to a bounding box instead.
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
   * @param support_radius the scaling of the reach of the nodes
   * @param basis_function the basis function that is used to parametrize the RBF object
   * @param weight the weighting associated to each node for the sum operation
   * @param nodes the center of each basis function
   * @param position the location of the RBF shape
   * @param orientation the orientation of the shape in relation to each main
   * axis
   */
  RBFShape<dim>(const std::vector<double> &          support_radii,
                const std::vector<RBFBasisFunction> &basis_functions,
                const std::vector<double> &          weights,
                const std::vector<Point<dim>> &      nodes,
                const Point<dim> &                   position,
                const Tensor<1, 3> &                 orientation);

  /**
   * @brief An RBFShape represents a physical object by describing its signed
   * distance field with a linear combination of radial basis functions. Each
   * radial basis function has a location and properties that are used in the
   * sum.
   * @param shape_arguments the concatenated vector of all shape arguments for
   * an RBF in the order: weights, support_radii, basis_functions, nodes_x,
   * nodes_y, nodes_z
   * @param position the location of the RBF shape
   * @param orientation the orientation of the shape in relation to each main
   * axis
   */
  RBFShape<dim>(const std::vector<double> &shape_arguments,
                const Point<dim> &         position,
                const Tensor<1, 3> &       orientation);

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point. The computation and addition of the
   * bounding box distance are necessary since the RBF nodes may not cover the
   * whole simulation domain. In that case, it is assumed that the distance from
   * the RBF object is approximately the same as the distance from the
   * corresponding bounding box.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
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
    const Point<dim> &                                   evaluation_point,
    const typename DoFHandler<dim>::active_cell_iterator cell,
    const unsigned int component = 0) override;

  /**
   * @brief Return the analytical gradient of the distance for the current RBF
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  Tensor<1, dim>
  gradient(const Point<dim> & evaluation_point,
           const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
   * @param levels_not_precalculated the number of finer levels that won't be
   * precalculated
   */
  void
  update_precalculations(DoFHandler<dim> &  dof_handler,
                         const unsigned int levels_not_precalculated);

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
   * It preserves C1 continuity at distance=0, and C0 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c0(const double distance) const
  {
    return distance > 1.0 ? 0.0 : (1.0 - std::pow(distance, 2.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C0 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c0(const double distance) const
  {
    return distance > 1.0 ? 0.0 : (1.0 - std::pow(distance, 3.0));
  }

  /**
   * @brief Compact polynomial function defined from 0 to 1.
   * It preserves C0 continuity at distance=0, and C1 continuity at distance=1.
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
   * It preserves C1 continuity at distance=0, and C1 continuity at distance=1.
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
   * It preserves C2 continuity at distance=0, and C1 continuity at distance=1.
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
   * It preserves C0 continuity at distance=0, and C2 continuity at distance=1.
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
   * It preserves C1 continuity at distance=0, and C2 continuity at distance=1.
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
   * It preserves C2 continuity at distance=0, and C2 continuity at distance=1.
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
   * It preserves C1 continuity at distance=0, and C0 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c0_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -2.0 * distance;
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C0 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c2c0_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -3.0 * std::pow(distance, 2.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C0 continuity at distance=0, and C1 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c0c1_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : 2.0 * (distance - 1.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C1 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c1c1_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : 6.0 * (distance - 1.0) * distance;
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C2 continuity at distance=0, and C1 continuity at distance=1.
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
   * It preserves C0 continuity at distance=0, and C2 continuity at distance=1.
   * @param distance distance to the node normalized by the support radius
   */
  inline double
  c0c2_derivative(const double distance) const
  {
    return distance > 1.0 ? 0.0 : -3.0 * std::pow(distance - 1.0, 2.0);
  }

  /**
   * @brief Derivative of a compact polynomial function defined from 0 to 1.
   * It preserves C1 continuity at distance=0, and C2 continuity at distance=1.
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
   * It preserves C2 continuity at distance=0, and C2 continuity at distance=1.
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
   * @brief Checks if the particle moved since the last precalculations
   */
  bool
  has_shape_moved();

  /**
   * @brief Checks if possible nodes affecting the current cell have been identified, and returns the proper vector to use for iteration
   * @param cell A likely one where the evaluation point is located
   */
  void
  prepare_iterable_nodes(
    const typename DoFHandler<dim>::active_cell_iterator cell);

  /**
   * @brief Resets the iterable nodes to all nodes
   * @param cell A likely one where the evaluation point is located
   */
  void
  reset_iterable_nodes(
    const typename DoFHandler<dim>::active_cell_iterator cell);

private:
  size_t                               number_of_nodes;
  std::shared_ptr<HyperRectangle<dim>> bounding_box;
  std::vector<size_t>                  iterable_nodes;

  std::map<const typename DoFHandler<dim>::cell_iterator,
           std::shared_ptr<std::vector<size_t>>>
                   likely_nodes_map;
  size_t           max_number_of_nodes;
  DoFHandler<dim> *dof_handler;
  Point<dim>       position_precalculated;
  Tensor<1, 3>     orientation_precalculated;

  // levels_not_precalculated is used in order to not precompute and store the
  // likely RBF nodes for the finer levels of cells in the grid. Setting this
  // parameter is a tradeoff between having faster distance evaluations and
  // reducing the memory footprint: increasing levels_not_precalculated
  // increases evaluation time.
  int    levels_not_precalculated;
  double maximal_support_radius;

public:
  std::vector<size_t>           nodes_id;
  std::vector<double>           weights;
  std::vector<Point<dim>>       nodes_positions;
  std::vector<Point<dim>>       rotated_nodes_positions;
  std::vector<double>           support_radii;
  std::vector<RBFBasisFunction> basis_functions;
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
  Cylinder<dim>(double              radius,
                double              half_length,
                const Point<dim> &  position,
                const Tensor<1, 3> &orientation)
    : Shape<dim>(radius, position, orientation)
    , radius(radius)
    , half_length(half_length)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  CylindricalTube<dim>(double              radius_inside,
                       double              radius_outside,
                       double              half_length,
                       const Point<dim> &  position,
                       const Tensor<1, 3> &orientation)
    : Shape<dim>((radius_outside + radius_inside) / 2., position, orientation)
    , radius((radius_outside + radius_inside) / 2.)
    , height(half_length * 2.)
    , rectangular_base(radius_outside - radius_inside)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

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
  CylindricalHelix<dim>(double              radius_helix,
                        double              radius_disk,
                        double              height,
                        double              pitch,
                        const Point<dim> &  position,
                        const Tensor<1, 3> &orientation)
    : Shape<dim>(radius_disk, position, orientation)
    , radius(radius_helix)
    , height(height)
    , pitch(pitch)
    , radius_disk(radius_disk)
  {}

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point.
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component This parameter is not used, but it is necessary because Shapes inherit from the Function class of deal.II.
   */
  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  /**
   * @brief Return a pointer to a copy of the Shape
   */
  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  /**
   * @brief
   * Return the volume displaced by the solid
   *
   * @param fluid_density The density of the fluid that is displaced
   */
  double
  displaced_volume(const double fluid_density) override;

private:
  double radius;
  double height;
  double pitch;
  double radius_disk;
};

#endif // lethe_shape_h
