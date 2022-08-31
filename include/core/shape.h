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
   * @brief enum class that associate an integer index tp each type of shape
   */
  enum ShapeType : int
  {
    sphere            = 0,
    rectangle         = 1,
    ellipsoid         = 2,
    torus             = 3,
    cone              = 4,
    cut_hollow_sphere = 5,
    death_star        = 6,
    composite_shape   = 7,
    rbf_shape         = 8,
  } type;


  virtual std::pair<std::string, int>
  get_shape_name() = 0;

  /**
   * @brief Return the evaluation of the signed distance function of this solid
   * at the given point evaluation point
   * Most levelset functions implemented come from Inigo Quilez:
   * iquilezles.org/articles/distfunctions
   *
   * @param evaluation_point The point at which the function will be evaluated
   * @param component Not applicable
   */
  virtual double
  value(const Point<dim> & evaluation_point,
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
   * @param fluid_density The density of the fluid that is displaced
   */
  virtual double
  displaced_volume(const double fluid_density) = 0;

  /**
   * @brief
   * Sets a new position for the shape
   *
   * @param The new position the shape will be placed at
   */
  virtual void
  set_position(const Point<dim> &position);

  /**
   * @brief
   * Sets a new orientation for the shape
   *
   * @param The new orientation the shape will be set at
   */
  virtual void
  set_orientation(const Tensor<1, 3> &orientation);

  /**
   * @brief
   * Returns the position of the shape
   *
   */
  virtual Point<dim>
  get_position();

  /**
   * @brief
   * Returns the orientation of the shape
   *
   */
  virtual Tensor<1, 3>
  get_orientation();

  /**
   * @brief
   * Most value functions assume that the particle's position is at the origin
   * and that the shape is aligned with one of the main axes. This function
   * returns a point that is rotated and translated, in accordance with the
   * current shape position and orientation, so that subsequent calculations for
   * the value function are made more easily; it abstract a step that is
   * required in the value function for most shapes.
   *
   * Returns the centered and aligned point used on the levelset evaluation.
   *
   * @param evaluation_point The point that will be recentered and realigned
   * @return The aligned and centered point
   */
  Point<dim>
  align_and_center(const Point<dim> &evaluation_point) const;
  // Effective radius used for crown refinement
  double effective_radius;


protected:
  // Position of the center of the Shape. It doesn't always correspond to the
  // center of mass
  Point<dim> position;
  // The solid orientation, which is defined as the sequential rotation around
  // the axes x->y->z by each of the tensor components, in radian
  Tensor<1, 3> orientation;
};


template <int dim>
class Sphere : public Shape<dim>
{
public:
  /**
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

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  Tensor<1, dim>
  gradient(const Point<dim> & evaluation_point,
           const unsigned int component = 0) const override;

  double
  displaced_volume(const double fluid_density) override;


  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("sphere", Shape<dim>::ShapeType::sphere);
  };

  void
  set_position(const Point<dim> &position) override;


private:
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
#else
  std::shared_ptr<Functions::SignedDistance::Sphere<dim>> sphere_function;
#endif
};

template <int dim>
class Rectangle : public Shape<dim>
{
public:
  /**
   * @param half_lengths The half lengths of each direction
   * @param position The rectangle center
   * @param orientation The rectangle orientation
   */
  Rectangle<dim>(const Tensor<1, dim> &half_lengths,
                 const Point<dim> &    position,
                 const Tensor<1, 3> &  orientation)
    : Shape<dim>(half_lengths.norm(), position, orientation)
    , half_lengths(half_lengths)
  {}

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("rectangle", Shape<dim>::ShapeType::rectangle);
  };

private:
  Tensor<1, dim> half_lengths;
};

template <int dim>
class Ellipsoid : public Shape<dim>
{
public:
  /**
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

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("ellipsoid", Shape<dim>::ShapeType::ellipsoid);
  }

private:
  Tensor<1, dim> radii;
};

template <int dim>
class Torus : public Shape<dim>
{
public:
  /**
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

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("torus", Shape<dim>::ShapeType::torus);
  }

private:
  double ring_radius;
  double ring_thickness;
};

template <int dim>
class Cone : public Shape<dim>
{
public:
  /**
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

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("cone", Shape<dim>::ShapeType::cone);
  }

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

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("cut_hollow_sphere",
                          Shape<dim>::ShapeType::cut_hollow_sphere);
  }

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
   * The Death Star is the result of a boolean substraction of one sphere from
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

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("death_star", Shape<dim>::ShapeType::death_star);
  }

private:
  double radius;
  double hole_radius;
  double spheres_distance;

  double intermediate_a;
  double intermediate_b;
};

// Composite Shapes are currently used only to output the signed distance of
// particles in the GLS Sharp Navier Stokes solver. The class was however
// designed so that specific composite shapes could be defined through the
// parameter file, although this functionality has not been implemented yet.
template <int dim>
class CompositeShape : public Shape<dim>
{
public:
  /**
   * @param components The shapes from which this composite sphere will be composed
   */
  CompositeShape<dim>(std::vector<std::shared_ptr<Shape<dim>>> components)
    : Shape<dim>(0.,
                 components[0]->get_position(),
                 components[0]->get_orientation())
    , components(components)
  {
    // Calculation of the effective radius
    for (const std::shared_ptr<Shape<dim>> &elem : components)
      {
        this->effective_radius =
          std::max(this->effective_radius, elem->effective_radius);
      }
  }

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("composite_shape",
                          Shape<dim>::ShapeType::composite_shape);
  }

private:
  std::vector<std::shared_ptr<Shape<dim>>> components;
};


// RBF Shapes express the signed distance function as a linear combination of
// Radial Basis Functions, which have a defined support radius and basis
// function. A collection of nodes and weights composes the object.
template <int dim>
class RBFShape : public Shape<dim>
{
public:
  /**
   * @param support_radius the scaling of the extent of the nodes
   * @param basis_function the basis function that was used to parametrize the RBF object
   * @param weight the weighting associated to each node for the sum operation
   * @param nodes the center of each basis function
   */
  RBFShape<dim>(std::vector<double>         support_radius,
                std::vector<unsigned int>   basis_function,
                std::vector<double>         weight,
                std::vector<Tensor<1, dim>> nodes,
                const Point<dim> &          position,
                const Tensor<1, 3> &        orientation)
    : Shape<dim>(support_radius[0], position, orientation)
    , support_radius(support_radius)
    , basis_function(basis_function)
    , nodes(nodes)
    , weight(weight)
  {
    unsigned int number_of_nodes = weight.size();

    Point<dim>     high_bounding_point = Point<dim>();
    Point<dim>     low_bounding_point  = Point<dim>();
    Point<dim>     bounding_box_center = Point<dim>();
    Tensor<1, dim> half_lengths        = Tensor<1, dim>();
    for (int d = 0; d < dim; d++)
      {
        high_bounding_point[d] = DBL_MIN;
        low_bounding_point[d]  = DBL_MAX;
        for (int i = 0; i < number_of_nodes; i++)
          {
            if (low_bounding_point[d] > nodes[i][d])
              low_bounding_point[d] = nodes[i][d];
            if (high_bounding_point[d] < nodes[i][d])
              high_bounding_point[d] = nodes[i][d];
          }
        bounding_box_center[d] =
          0.5 *
          (low_bounding_point[d] + high_bounding_point[d]); // + position[d];
        half_lengths[d] =
          0.5 * (high_bounding_point[d] - low_bounding_point[d]);
      }
    bounding_box = std::make_shared<Rectangle<dim>>(half_lengths,
                                                    bounding_box_center,
                                                    orientation);
  }

  double
  value(const Point<dim> & evaluation_point,
        const unsigned int component = 0) const override;

  std::shared_ptr<Shape<dim>>
  static_copy() const override;

  double
  displaced_volume(const double fluid_density) override;

  double
  // TODO How to properly cite bitpit
  evaluate_basis_function(const unsigned int basis_function_id,
                          const double       distance) const;

  double
  wendlandc2(double) const;
  double
  linear(double) const;
  double
  gauss90(double) const;
  double
  gauss95(double) const;
  double
  gauss99(double) const;
  double
  c1c0(double) const;
  double
  c2c0(double) const;
  double
  c0c1(double) const;
  double
  c1c1(double) const;
  double
  c2c1(double) const;
  double
  c0c2(double) const;
  double
  c1c2(double) const;
  double
  c2c2(double) const;


  std::pair<std::string, int>
  get_shape_name() override
  {
    return std::make_pair("rbf_shape", Shape<dim>::ShapeType::rbf_shape);
  }

  /**
   * Class taken from Optimad Bitpit. https://github.com/optimad/bitpit
   * // TODO How to properly introduce the code/citation
   * @enum RBFBasisFunction
   * @ingroup RBF
   * @brief Enum class defining types of RBF kernel functions that could be used in bitpit::RBF class
   */
  enum class RBFBasisFunction
  {
    CUSTOM     = 0,
    WENDLANDC2 = 1,
    LINEAR     = 2,
    GAUSS90    = 3,
    GAUSS95    = 4,
    GAUSS99    = 5,
    C1C0       = 6,
    C2C0       = 7,
    C0C1       = 8,
    C1C1       = 9,
    C2C1       = 10,
    C0C2       = 11,
    C1C2       = 12,
    C2C2       = 13,
  };

private:
  std::vector<double>             weight;
  std::vector<Tensor<1, dim>>     nodes;
  std::vector<double>             support_radius;
  std::vector<unsigned int>       basis_function;
  std::shared_ptr<Rectangle<dim>> bounding_box;
};

#endif // lethe_shape_h
