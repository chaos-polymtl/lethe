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

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

#ifdef DEAL_II_WITH_OPENCASCADE
#  include <deal.II/opencascade/manifold_lib.h>
#  include <deal.II/opencascade/utilities.h>

#  include <BRepBuilderAPI_MakeVertex.hxx>
#  include <BRepClass3d_SolidClassifier.hxx>
#  include <BRepExtrema_DistShapeShape.hxx>
#endif



template <int dim>
double
Shape<dim>::displaced_volume(const double /*fluid_density*/)
{
  StandardExceptions::ExcNotImplemented();
  return 0;
}


template <int dim>
void
Shape<dim>::clear_cache()
{
  value_cache.clear();
  gradient_cache.clear();
  closest_point_cache.clear();
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
      if (std::abs(theta[2]) > 1e-10 * this->effective_radius)
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
          if (std::abs(theta[i]) > 1e-10 * this->effective_radius)
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
Point<dim>
Shape<dim>::reverse_align_and_center(const Point<dim> &evaluation_point) const
{
  // Reverse translation and rotation from standard position and orientation to
  // global referential
  Point<dim> center_of_rotation = position;
  Point<dim> rotated_point;
  Point<dim> translated_point;

  // Rotation from the solid orientation
  // Angular position around x, y and z axis
  Tensor<1, 3> theta = orientation;

  // The centralized point is the one to be rotated, and it is updated after
  // each rotation around one axis. The centralized rotated point is the result
  // of each rotation, and it is initialized in case no rotation is performed
  Point<dim> centralized_point;
  centralized_point              = evaluation_point;
  Point<dim> centralized_rotated = centralized_point;

  // Selection of the first axis around which to rotate:
  // x -> 0, y -> 1, z -> 2
  // In 2D, only rotation around the z axis is possible

  if constexpr (dim == 2)
    {
      if (std::abs(theta[2]) > 1e-10 * this->effective_radius)
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
      for (unsigned int i = 0; i < dim; ++i)
        {
          if (std::abs(theta[2 - i]) > 1e-10 * this->effective_radius)
            {
              Tensor<1, 3> axis;
              axis[2 - i] = 1.0;
              Tensor<2, 3> rotation_matrix =
                Physics::Transformations::Rotations::rotation_matrix_3d(
                  axis, theta[2 - i]);

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
  translated_point = rotated_point;

  return translated_point;
}

template <int dim>
std::shared_ptr<Manifold<dim - 1, dim>>
Shape<dim>::get_shape_manifold()
{
  return std::make_shared<FlatManifold<dim - 1, dim>>();
}

template <int dim>
std::string
Shape<dim>::point_to_string(const Point<dim> &evaluation_point) const
{
  // This function transforms a point into a string.
  // The point precision is conserve up to a precision of 1e-12.
  std::string point_in_string = "";
  for (unsigned int d = 0; d < dim; ++d)
    {
      point_in_string =
        point_in_string + ";" +
        std::to_string(std::round(evaluation_point[d] * 1e12) / 1e12);
    }
  return point_in_string;
}


template <int dim>
void
Shape<dim>::closest_surface_point(
  const Point<dim>                                     &p,
  Point<dim>                                           &closest_point,
  const typename DoFHandler<dim>::active_cell_iterator &cell_guess)
{
  Tensor<1, dim> actual_gradient;
  double         distance_from_surface;
  actual_gradient       = this->gradient_with_cell_guess(p, cell_guess);
  distance_from_surface = this->value_with_cell_guess(p, cell_guess);
  // Check if the gradient is well-defined. If the point is on the surface,
  // the gradient can be badly defined for some shapes. We return the point
  // directly in these cases since it would be on the surface.
  closest_point = p - (actual_gradient / (actual_gradient.norm() +
                                          1e-16 * this->effective_radius)) *
                        distance_from_surface;
}


template <int dim>
void
Shape<dim>::closest_surface_point(const Point<dim> &p,
                                  Point<dim>       &closest_point) const
{
  Tensor<1, dim> actual_gradient;
  double         distance_from_surface;
  actual_gradient       = this->gradient(p);
  distance_from_surface = this->value(p);

  // Check if the gradient is well-defined. If the point is on the surface,
  // the gradient can be badly defined for some shapes. We return the point
  // directly in these cases since it would be on the surface.
  closest_point = p - (actual_gradient / (actual_gradient.norm() +
                                          1e-16 * this->effective_radius)) *
                        distance_from_surface;
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
  return evaluation_point.distance(this->position) - this->effective_radius -
         this->layer_thickening;
#else
  return sphere_function->value(evaluation_point) - this->layer_thickening;
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
  // We make sure that the evaluation point and the sphere center are different,
  // because if they are the same the analytical gradient is not defined: the
  // function returns a NaN. We use the numerical gradient if the points are the
  // same.
  if ((evaluation_point - this->position).norm() <
      1e-12 * this->effective_radius)
    return AutoDerivativeFunction<dim>::gradient(evaluation_point);
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
std::tuple<double, Tensor<1, dim>, Point<dim>>
Sphere<dim>::distance_to_shape_with_cell_guess(
  Shape<dim>                                           &shape,
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  std::vector<Point<dim>>                              &candidate_points,
  double                                                precision,
  bool exact_distance_outside_of_contact) {
  (void) candidate_points;
  (void) precision;
  double                      distance = DBL_MAX;
  Tensor<1, dim>              normal;
  Point<dim>                  contact_point;
  if (typeid(shape) == typeid(Sphere<dim>))
    {
      Tensor<1,dim> center_to_center_vector=shape.get_position()-this->position;
      distance=(center_to_center_vector.norm()-this->effective_radius-shape.effective_radius)/2;
      normal=center_to_center_vector/center_to_center_vector.norm();
      contact_point=this->position+normal*(this->effective_radius+distance);
    }
  else
    {
      distance=(shape.value_with_cell_guess(this->position,cell)-this->effective_radius)/2;
      normal=-shape.gradient_with_cell_guess(this->position,cell);
      contact_point=this->position+normal*(this->effective_radius+distance);
    }

  return std::make_tuple(distance, normal, contact_point);

}

template <int dim>
std::tuple<double, Tensor<1, dim>, Point<dim>>
Sphere<dim>::distance_to_shape(
  Shape<dim>                                           &shape,
  std::vector<Point<dim>>                              &candidate_points,
  double                                                precision,
  bool exact_distance_outside_of_contact ) {

  (void) candidate_points;
  (void) precision;
  double                      distance = DBL_MAX;
  Tensor<1, dim>              normal;
  Point<dim>                  contact_point;
  if (typeid(shape) == typeid(Sphere<dim>))
    {
      Tensor<1,dim> center_to_center_vector=shape.get_position()-this->position;
      distance=(center_to_center_vector.norm()-this->effective_radius-shape.effective_radius)/2;
      normal=center_to_center_vector/center_to_center_vector.norm();
      contact_point=this->position+normal*(this->effective_radius+distance);
    }
  else
    {
      distance=(shape.value(this->position)-this->effective_radius)/2;
      normal=-shape.gradient(this->position);
      contact_point=this->position+normal*(this->effective_radius+distance);
    }

  return std::make_tuple(distance, normal, contact_point);
}


template <int dim>
double
Plane<dim>::value(const Point<dim> &evaluation_point,
                  const unsigned int /*component*/) const
{
  Point<dim> current_point = this->align_and_center(evaluation_point);
  double     dot_product   = scalar_product((current_point), normal);
  Point<dim> projected_point =
    current_point - dot_product / normal.norm_square() * normal;

  auto rotate_in_globalpoint = this->reverse_align_and_center(projected_point);
  if (dot_product > 0)
    return (rotate_in_globalpoint - evaluation_point).norm();
  else
    return -(rotate_in_globalpoint - evaluation_point).norm();
}


template <int dim>
std::shared_ptr<Shape<dim>>
Plane<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<Plane<dim>>(this->position, this->orientation);
  return copy;
}

template <int dim>
Tensor<1, dim>
Plane<dim>::gradient(const Point<dim> &evaluation_point,
                     const unsigned int /*component*/) const
{
  // We make sure that the evaluation point and the sphere center are different,
  // because if they are the same the analytical gradient is not defined: the
  // function returns a NaN. We use the numerical gradient if the points are the
  // same.
  Point<dim> current_point = this->align_and_center(evaluation_point);
  double     dot_product   = scalar_product((current_point), normal);
  Point<dim> projected_point =
    current_point - dot_product / normal.norm_square() * normal;

  auto rotate_in_globalpoint = this->reverse_align_and_center(projected_point);
  if (dot_product > 0)
    return (rotate_in_globalpoint - evaluation_point) /
           (rotate_in_globalpoint - evaluation_point).norm();
  else
    return -(rotate_in_globalpoint - evaluation_point) /
           ((rotate_in_globalpoint - evaluation_point).norm()+DBL_MIN);
}

template <int dim>
double
Plane<dim>::displaced_volume(const double /*fluid_density*/)
{
  return 0;
}


template <int dim>
void
Superquadric<dim>::closest_surface_point(
  const Point<dim> &p,
  Point<dim>       &closest_point,
  const typename DoFHandler<dim>::active_cell_iterator & /*cell_guess*/)
{
  auto point_in_string = this->point_to_string(p);
  auto iterator        = this->closest_point_cache.find(point_in_string);
  if (iterator == this->closest_point_cache.end())
    {
      this->closest_surface_point(p, closest_point);

      Point<dim> copy_closest_point{};
      for (unsigned int d = 0; d < dim; d++)
        copy_closest_point[d] = closest_point[d];
      if (!this->part_of_a_composite)
        {
          this->closest_point_cache[point_in_string] = copy_closest_point;
          this->value_cache[point_in_string]         = this->value(p);
          this->gradient_cache[point_in_string]      = this->gradient(p);
        }
    }
  else
    closest_point = iterator->second;
}

template <int dim>
void
Superquadric<dim>::closest_surface_point(const Point<dim> &p,
                                         Point<dim>       &closest_point) const
{
  auto point_in_string = this->point_to_string(p);
  auto iterator        = this->closest_point_cache.find(point_in_string);
  if (iterator == this->closest_point_cache.end())
    {
      // The initial guess is chosen as the centered point (i.e. evaluation
      // point, in the superquadric referential). It should already be somewhat
      // close to the closest point, and it is already in the right octant.
      Point<dim> current_point = this->align_and_center(p);

      unsigned int       iteration     = 0;
      const unsigned int iteration_max = 1e2;

      Point<dim> dx{}, distance_gradient{};
      double     current_distance = superquadric(current_point);

      double relaxation = 1;
      while (iteration < iteration_max && abs(current_distance) > epsilon)
        {
          distance_gradient = superquadric_gradient(current_point);
          // limit the step size

          if (distance_gradient.norm() < epsilon)
            // Gradient can be null if the evaluation point is exactly on the
            // centroid of the shape. In this case it is also the closest point.
            break;
          // This is a modify newton method that limit the step size when the gradient norm is smaller then 1.
          dx = -relaxation * (current_distance * distance_gradient) /
               distance_gradient.norm()/std::max(1.0,distance_gradient.norm());

          current_point    = current_point + dx;
          current_distance = superquadric(current_point);
          //std::cout<<"iteration i "<< iteration <<" current distance "<< current_distance<<std::endl;

          iteration++;
        }
      if(iteration==iteration_max){
          std::cout<<"hi superquadric did not converge after 100 iteration  point"<< p <<std::endl;
        }

      closest_point = this->reverse_align_and_center(current_point);
    }
  else
    closest_point = iterator->second;
}

template <int dim>
double
Superquadric<dim>::value(const Point<dim> &evaluation_point,
                         const unsigned int /*component*/) const
{
  using numbers::PI;

  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator == this->value_cache.end())
    {
      // The closest point has to be found first, because it is used in value
      // calculation.
      Point<dim> closest_point{};
      this->closest_surface_point(evaluation_point, closest_point);

      Point<dim> centered_point = this->align_and_center(evaluation_point);
      if (superquadric(centered_point) > 0)
        return (closest_point - evaluation_point).norm() -
               this->layer_thickening;
      else
        return -(closest_point - evaluation_point).norm() -
               this->layer_thickening;
    }
  else
    return iterator->second;
}

template <int dim>
double
Superquadric<dim>::value_with_cell_guess(
  const Point<dim> &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
  const unsigned int /*component = 0*/)
{
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator == this->value_cache.end())
    {
      // The closest point has to be found first, because it is used in value
      // calculation.
      Point<dim> closest_point{};
      this->closest_surface_point(evaluation_point, closest_point);
      Point<dim> copy_closest_point{};
      for (unsigned int d = 0; d < dim; d++)
        copy_closest_point[d] = closest_point[d];

      double levelset = this->value(evaluation_point);
      if (!this->part_of_a_composite)
        {
          this->closest_point_cache[point_in_string] = copy_closest_point;
          this->value_cache[point_in_string]         = levelset;
          this->gradient_cache[point_in_string] =
            this->gradient(evaluation_point);
        }
      return levelset;
    }
  else
    return iterator->second;
}

template <int dim>
std::shared_ptr<Shape<dim>>
Superquadric<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<Superquadric<dim>>(this->half_lengths,
                                        this->exponents,
                                        this->epsilon,
                                        this->position,
                                        this->orientation);
  return copy;
}

template <int dim>
Tensor<1, dim>
Superquadric<dim>::gradient(const Point<dim> &evaluation_point,
                            const unsigned int /*component*/) const
{
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->gradient_cache.find(point_in_string);
  if (iterator == this->gradient_cache.end())
    {
      Point<dim> closest_point{};
      this->closest_surface_point(evaluation_point, closest_point);

      Tensor<1, dim> gradient;
      if (superquadric(evaluation_point) > 0)
        gradient = -(evaluation_point - closest_point) /
                   ((evaluation_point - closest_point).norm() +
                    1e-16 * this->effective_radius);
      else
        gradient = (evaluation_point - closest_point) /
                   ((evaluation_point - closest_point).norm() +
                    1e-16 * this->effective_radius);
      ;
      return gradient;
    }
  else
    return iterator->second;
}

template <int dim>
Tensor<1, dim>
Superquadric<dim>::gradient_with_cell_guess(
  const Point<dim> &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
  const unsigned int /*component*/)
{
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->gradient_cache.find(point_in_string);
  if (iterator == this->gradient_cache.end())
    {
      // The closest point has to be found first, because it is used in value
      // calculation
      Point<dim> closest_point{};
      this->closest_surface_point(evaluation_point, closest_point);
      Point<dim> copy_closest_point{};
      for (unsigned int d = 0; d < dim; d++)
        copy_closest_point[d] = closest_point[d];

      Tensor<1, dim> gradient = this->gradient(evaluation_point);
      if (!this->part_of_a_composite)
        {
          this->closest_point_cache[point_in_string] = copy_closest_point;
          this->value_cache[point_in_string]    = this->value(evaluation_point);
          this->gradient_cache[point_in_string] = gradient;
        }
      return gradient;
    }
  else
    return iterator->second;
}

template <int dim>
double
OpenCascadeShape<dim>::value(const Point<dim> &evaluation_point,
                             const unsigned int /*component*/) const
{
#ifdef DEAL_II_WITH_OPENCASCADE
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator != this->value_cache.end())
    return iterator->second;

  Point<dim>    centered_point = this->align_and_center(evaluation_point);
  Point<dim>    projected_point;
  auto          pt     = OpenCASCADE::point(centered_point);
  TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(pt);
  BRepExtrema_DistShapeShape distancetool(shape, vertex);
  distancetool.Perform();
  gp_Pnt pt_on_surface = distancetool.PointOnShape1(1);
  projected_point[0]   = pt_on_surface.X();
  projected_point[1]   = pt_on_surface.Y();
  if constexpr (dim == 3)
    {
      projected_point[2] = pt_on_surface.Z();
    }

  BRepClass3d_SolidClassifier point_classifier(shape);
  point_classifier.Perform(pt, shape_tol);
  TopAbs_State point_state = point_classifier.State();
  // Check if the evaluation point is inside the shape (This check can return
  // true only if the shape is a solid). By default, step file format will
  // define a solid. This is necessary to sign the distance function evaluation.
  if (point_state == TopAbs_State::TopAbs_IN)
    {
      // If the evaluation point is inside the shape and the shape is a solid,
      // the distance to the shape will be 0. This is why we need to evaluate
      // the distance of the point with the shell of the shape.
      // Rotate the solution found to the global reference frame.
      auto rotate_in_globalpoint =
        this->reverse_align_and_center(projected_point);
      return -(rotate_in_globalpoint - evaluation_point).norm() -
             this->layer_thickening;
    }
  else
    {
      auto rotate_in_globalpoint =
        this->reverse_align_and_center(projected_point);
      return (rotate_in_globalpoint - evaluation_point).norm() -
             this->layer_thickening;
    }

#else
  (void)evaluation_point;
  return 0;
#endif
}
template <int dim>
double
OpenCascadeShape<dim>::value_with_cell_guess(
  const Point<dim> &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
  const unsigned int /*component*/)
{
#ifdef DEAL_II_WITH_OPENCASCADE
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);

  if (iterator == this->value_cache.end())
    {
      // Transform the point to an OpenCascade point
      Point<dim> centered_point = this->align_and_center(evaluation_point);
      Point<dim> projected_point;
      vertex_position = OpenCASCADE::point(centered_point);
      // Perform the evaluation
      vertex = BRepBuilderAPI_MakeVertex(vertex_position);
      distancetool.LoadS2(vertex);
      distancetool.Perform();

      // Convert OpenCascade solution point to dealii Point.
      gp_Pnt pt_on_surface = distancetool.PointOnShape1(1);
      projected_point[0]   = pt_on_surface.X();
      projected_point[1]   = pt_on_surface.Y();
      if constexpr (dim == 3)
        {
          projected_point[2] = pt_on_surface.Z();
        }

      point_classifier.Perform(vertex_position, shape_tol);
      TopAbs_State point_state = point_classifier.State();
      // Check if the evaluation point is inside the shape (This check can
      // return true only if the shape is a solid). By default, step file format
      // will define a solid. This is necessary to sign the distance function
      // evaluation.
      if (point_state == TopAbs_State::TopAbs_IN)
        {
          // If the evaluation point is inside the shape and the shape is a
          // solid, the distance to the shape will be 0. This is why we need to
          // evaluate the distance of the point with the shell of the shape.
          // Rotate the solution found to the global reference frame and cache
          // the solution.
          if (!this->part_of_a_composite)
            {
              this->value_cache[point_in_string] =
                -(centered_point - projected_point).norm();
              auto rotate_in_globalpoint =
                this->reverse_align_and_center(projected_point);
              this->gradient_cache[point_in_string] =
                (rotate_in_globalpoint - evaluation_point) /
                ((rotate_in_globalpoint - evaluation_point).norm() + 1.0e-16);
            }
          return -(centered_point - projected_point).norm() -
                 this->layer_thickening;
        }
      else
        {
          if (!this->part_of_a_composite)
            {
              this->value_cache[point_in_string] =
                (centered_point - projected_point).norm();
              auto rotate_in_globalpoint =
                this->reverse_align_and_center(projected_point);
              this->gradient_cache[point_in_string] =
                -(rotate_in_globalpoint - evaluation_point) /
                ((rotate_in_globalpoint - evaluation_point).norm() + 1.0e-16);
            }
          return (centered_point - projected_point).norm() -
                 this->layer_thickening;
        }
    }
  else
    {
      return this->value_cache[point_in_string];
    }
#else
  (void)evaluation_point;
  return 0;
#endif
}

template <int dim>
std::shared_ptr<Shape<dim>>
OpenCascadeShape<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<OpenCascadeShape<dim>>(local_file_name,
                                            this->position,
                                            this->orientation);
  return copy;
}

template <int dim>
Tensor<1, dim>
OpenCascadeShape<dim>::gradient(const Point<dim> &evaluation_point,
                                const unsigned int /*component*/) const
{
#ifdef DEAL_II_WITH_OPENCASCADE
  Point<dim> centered_point = this->align_and_center(evaluation_point);
  Point<dim> projected_point;

  // Transform the point to an OpenCascade point
  auto          pt     = OpenCASCADE::point(centered_point);
  TopoDS_Vertex vertex = BRepBuilderAPI_MakeVertex(pt);
  // Perform the evaluation
  BRepExtrema_DistShapeShape distancetool(shape, vertex);
  distancetool.Perform();
  gp_Pnt pt_on_surface = distancetool.PointOnShape1(1);

  // Convert OpenCascade solution point to dealii Point.
  projected_point[0] = pt_on_surface.X();
  projected_point[1] = pt_on_surface.Y();
  if constexpr (dim == 3)
    {
      projected_point[2] = pt_on_surface.Z();
    }
  BRepClass3d_SolidClassifier point_classifier(shape);
  point_classifier.Perform(pt, shape_tol);
  TopAbs_State point_state = point_classifier.State();
  // Check the evaluation point is inside the shape (This check can return true
  // only if the shape is a solid). By default, step file format will define a
  // solid. This is necessary to sign the distance function evaluation.
  if (point_state == TopAbs_State::TopAbs_IN)
    {
      // If the evaluation point is inside the shape and the shape is a solid,
      // the distance to the shape will be 0. This is why we need to evaluate
      // the distance of the point with the shell of the shape.
      // Rotate the solution found to the global reference frame.
      auto rotate_in_globalpoint =
        this->reverse_align_and_center(projected_point);
      return (rotate_in_globalpoint - evaluation_point) /
             ((rotate_in_globalpoint - evaluation_point).norm() + 1.0e-16);
    }
  else
    {
      auto rotate_in_globalpoint =
        this->reverse_align_and_center(projected_point);
      return -(rotate_in_globalpoint - evaluation_point) /
             ((rotate_in_globalpoint - evaluation_point).norm() + 1.0e-16);
    }
#else
  // Empty return if Lethe is not compile with OpenCascade
  (void)evaluation_point;
  return Point<dim>();
#endif
}


template <int dim>
Tensor<1, dim>
OpenCascadeShape<dim>::gradient_with_cell_guess(
  const Point<dim> &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator /*cell*/,
  const unsigned int /*component*/)
{
#ifdef DEAL_II_WITH_OPENCASCADE
  // This function is identical to the Gradient function but has the possibility
  // to use the cache of the shape.

  // First step is to convert the point into a string.
  auto point_in_string = this->point_to_string(evaluation_point);
  // We check the unordered map
  auto iterator = this->gradient_cache.find(point_in_string);
  // If we did not find the solution in the map we perform the calculation.
  if (iterator == this->gradient_cache.end())
    {
      // Transform the point to an OpenCascade point
      Point<dim> centered_point = this->align_and_center(evaluation_point);
      Point<dim> projected_point;
      vertex_position = OpenCASCADE::point(centered_point);
      // Perform the evaluation
      vertex = BRepBuilderAPI_MakeVertex(vertex_position);
      distancetool.LoadS2(vertex);
      distancetool.Perform();

      // Convert OpenCascade solution point to dealii Point.
      gp_Pnt pt_on_surface = distancetool.PointOnShape1(1);
      projected_point[0]   = pt_on_surface.X();
      projected_point[1]   = pt_on_surface.Y();
      if constexpr (dim == 3)
        {
          projected_point[2] = pt_on_surface.Z();
        }
      point_classifier.Perform(vertex_position, shape_tol);
      TopAbs_State point_state = point_classifier.State();
      // Check the evaluation point is inside the shape (This check can return
      // true only if the shape is a solid). By default, step file format will
      // define a solid. This is necessary to sign the distance function
      // evaluation.
      if (point_state == TopAbs_State::TopAbs_IN)
        {
          // Rotate the solution found to the global reference frame and cache
          // the solution.
          auto rotate_in_globalpoint =
            this->reverse_align_and_center(projected_point);
          if (!this->part_of_a_composite)
            {
              this->value_cache[point_in_string] =
                -(centered_point - projected_point).norm();
              this->gradient_cache[point_in_string] =
                (rotate_in_globalpoint - evaluation_point) /
                ((rotate_in_globalpoint - evaluation_point).norm() + 1.0e-16);
              return this->gradient_cache[point_in_string];
            }
          else
            return (rotate_in_globalpoint - evaluation_point) /
                   ((rotate_in_globalpoint - evaluation_point).norm() +
                    1.0e-16);
        }
      else
        {
          auto rotate_in_globalpoint =
            this->reverse_align_and_center(projected_point);
          if (!this->part_of_a_composite)
            {
              this->value_cache[point_in_string] =
                (centered_point - projected_point).norm();
              this->gradient_cache[point_in_string] =
                -(rotate_in_globalpoint - evaluation_point) /
                ((rotate_in_globalpoint - evaluation_point).norm() + 1.0e-16);
              return this->gradient_cache[point_in_string];
            }
          else
            return -(rotate_in_globalpoint - evaluation_point) /
                   ((rotate_in_globalpoint - evaluation_point).norm() +
                    1.0e-16);
        }
    }
  else
    {
      // If we are here it is that this point was already evaluated and the
      // previous evaluation was cache.
      return this->gradient_cache[point_in_string];
    }
#else

  (void)evaluation_point;
  return Tensor<1, dim>();
#endif
}


template <int dim>
double
OpenCascadeShape<dim>::displaced_volume(const double /*fluid_density*/)
{
  return 0;
}

template <int dim>
void
OpenCascadeShape<dim>::set_position(const Point<dim> &position)
{
  this->Shape<dim>::set_position(position);
}


template <int dim>
double
HyperRectangle<dim>::value(const Point<dim> &evaluation_point,
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
  return max_q_0.norm() + std::min(max_q, 0.0) - this->layer_thickening;
}

template <int dim>
std::shared_ptr<Shape<dim>>
HyperRectangle<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<HyperRectangle<dim>>(this->half_lengths,
                                          this->position,
                                          this->orientation);
  return copy;
}

template <int dim>
double
HyperRectangle<dim>::displaced_volume(const double /*fluid_density*/)
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
  return k0 * (k0 - 1.) / k1 - this->layer_thickening;
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
  return q.norm() - ring_thickness - this->layer_thickening;
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

  return sqrt(d) * ((s > 0) ? 1 : ((s < 0) ? -1 : 0)) - this->layer_thickening;
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
      return (q - wh).norm() - this->layer_thickening;
    }
  else
    {
      return std::abs(q.norm() - radius) - shell_thickness -
             this->layer_thickening;
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
      return (corrected_p_2d - ab).norm() - this->layer_thickening;
    }
  else
    {
      Point<2> d0({spheres_distance, 0.});
      return std::max(corrected_p_2d.norm() - radius,
                      -((corrected_p_2d - d0).norm() - hole_radius)) -
             this->layer_thickening;
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
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator == this->value_cache.end())
    {
      // We align and center the evaluation point according to the shape
      // referential
      Point<dim> centered_point = this->align_and_center(evaluation_point);

      // The levelset value of all component shapes is computed
      std::map<unsigned int, double>         constituent_shapes_values;
      std::map<unsigned int, Tensor<1, dim>> constituent_shapes_gradients;
      for (auto const &[component_id, component] : constituents)
        {
          constituent_shapes_values[component_id] =
            component->value(centered_point);
          // A dummy gradient is used here because apply_boolean_operations
          // requires a gradient map as an argument. This design choice of not
          // duplicating apply_boolean_operations was made for brevity of the
          // code, at a negligible additional computing cost.
          constituent_shapes_gradients[component_id] = Tensor<1, dim>{};
        }

      double levelset;
      std::tie(levelset, std::ignore) =
        apply_boolean_operations(constituent_shapes_values,
                                 constituent_shapes_gradients);
      return levelset;
    }
  else
    return iterator->second;
}

template <int dim>
double
CompositeShape<dim>::value_with_cell_guess(
  const Point<dim>                                    &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator == this->value_cache.end())
    {
      // We align and center the evaluation point according to the shape
      // referential
      Point<dim> centered_point = this->align_and_center(evaluation_point);

      // The levelset value of all component shapes is computed
      std::map<unsigned int, double>         constituent_shapes_values;
      std::map<unsigned int, Tensor<1, dim>> constituent_shapes_gradients;
      for (auto const &[component_id, component] : constituents)
        {
          constituent_shapes_values[component_id] =
            component->value_with_cell_guess(centered_point, cell);
          // A dummy gradient is used here because apply_boolean_operations
          // requires a gradient map as an argument. This design choice of not
          // duplicating apply_boolean_operations was made for brevity of the
          // code, at a negligible additional computing cost.
          constituent_shapes_gradients[component_id] = Tensor<1, dim>{};
        }

      double levelset;
      std::tie(levelset, std::ignore) =
        apply_boolean_operations(constituent_shapes_values,
                                 constituent_shapes_gradients);
      if (!this->part_of_a_composite)
        {
          this->value_cache[point_in_string] = levelset;
        }
      return levelset;
    }
  else
    return this->value_cache[point_in_string];
}

template <int dim>
Tensor<1, dim>
CompositeShape<dim>::gradient(const Point<dim> &evaluation_point,
                              const unsigned int /*component*/) const
{
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->gradient_cache.find(point_in_string);
  if (iterator == this->gradient_cache.end())
    {
      // We align and center the evaluation point according to the shape
      // referential
      Point<dim> centered_point = this->align_and_center(evaluation_point);
      // The levelset value and gradient of all component shapes is computed
      std::map<unsigned int, double>         constituent_shapes_values;
      std::map<unsigned int, Tensor<1, dim>> constituent_shapes_gradients;
      for (auto const &[component_id, component] : constituents)
        {
          constituent_shapes_values[component_id] =
            component->value(centered_point);
          constituent_shapes_gradients[component_id] =
            component->gradient(centered_point);
        }

      double         levelset;
      Tensor<1, dim> gradient;
      std::tie(levelset, gradient) =
        apply_boolean_operations(constituent_shapes_values,
                                 constituent_shapes_gradients);

      // The gradient obtained is in the shape referential. It needs to be
      // returned to the global referential. We find the closest point (in the
      // shape referential), we return that point to its global referential
      // equivalent, then compute the gradient in the global referential.
      Point<dim> closest_point_shape_referential =
        centered_point -
        (gradient / (gradient.norm() + 1e-16 * this->effective_radius)) *
          levelset;
      Point<dim> closest_point_global_referential =
        this->reverse_align_and_center(closest_point_shape_referential);
      gradient = -(closest_point_global_referential - evaluation_point) /
                 ((closest_point_global_referential - evaluation_point).norm() +
                  1.0e-16);

      return gradient;
    }
  else
    return iterator->second;
}

template <int dim>
Tensor<1, dim>
CompositeShape<dim>::gradient_with_cell_guess(
  const Point<dim>                                    &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->gradient_cache.find(point_in_string);
  if (iterator == this->gradient_cache.end())
    {
      // We align and center the evaluation point according to the shape
      // referential
      Point<dim> centered_point = this->align_and_center(evaluation_point);
      // The levelset value and gradient of all component shapes is computed
      std::map<unsigned int, double>         constituent_shapes_values;
      std::map<unsigned int, Tensor<1, dim>> constituent_shapes_gradients;
      for (auto const &[component_id, component] : constituents)
        {
          constituent_shapes_values[component_id] =
            component->value_with_cell_guess(centered_point, cell);
          constituent_shapes_gradients[component_id] =
            component->gradient_with_cell_guess(centered_point, cell);
        }

      double         levelset;
      Tensor<1, dim> gradient;
      std::tie(levelset, gradient) =
        apply_boolean_operations(constituent_shapes_values,
                                 constituent_shapes_gradients);

      // The gradient obtained is in the shape referential. It needs to be
      // returned to the global referential. We find the closest point (in the
      // shape referential), we return that point to its global referential
      // equivalent, then compute the gradient in the global referential.
      Point<dim> closest_point_shape_referential =
        centered_point -
        (gradient / (gradient.norm() + 1e-16 * this->effective_radius)) *
          levelset;
      Point<dim> closest_point_global_referential =
        this->reverse_align_and_center(closest_point_shape_referential);
      gradient = -(closest_point_global_referential - evaluation_point) /
                 ((closest_point_global_referential - evaluation_point).norm() +
                  1.0e-16);
      if (!this->part_of_a_composite)
        {
          this->gradient_cache[point_in_string] = gradient;
        }
      return gradient;
    }
  else
    return this->gradient_cache[point_in_string];
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
  DoFHandler<dim> &updated_dof_handler,
  const bool       mesh_based_precalculations)
{
  if (!mesh_based_precalculations)
    return;
  for (auto const &constituent : constituents | boost::adaptors::map_values)
    if (typeid(*constituent) == typeid(RBFShape<dim>))
      std::static_pointer_cast<RBFShape<dim>>(constituent)
        ->update_precalculations(updated_dof_handler,
                                 mesh_based_precalculations);
    else if (typeid(*constituent) == typeid(CompositeShape<dim>))
      std::static_pointer_cast<CompositeShape<dim>>(constituent)
        ->update_precalculations(updated_dof_handler,
                                 mesh_based_precalculations);
}

template <int dim>
void
CompositeShape<dim>::load_data_from_file()
{
  for (auto const &constituent : constituents | boost::adaptors::map_values)
    if (typeid(*constituent) == typeid(RBFShape<dim>))
      std::static_pointer_cast<RBFShape<dim>>(constituent)
        ->load_data_from_file();
    else if (typeid(*constituent) == typeid(CompositeShape<dim>))
      std::static_pointer_cast<CompositeShape<dim>>(constituent)
        ->load_data_from_file();
}

template <int dim>
void
CompositeShape<dim>::remove_superfluous_data(
  DoFHandler<dim> &updated_dof_handler,
  const bool       mesh_based_precalculations)
{
  for (auto const &constituent : constituents | boost::adaptors::map_values)
    if (typeid(*constituent) == typeid(RBFShape<dim>))
      std::static_pointer_cast<RBFShape<dim>>(constituent)
        ->remove_superfluous_data(updated_dof_handler,
                                  mesh_based_precalculations);
    else if (typeid(*constituent) == typeid(CompositeShape<dim>))
      std::static_pointer_cast<CompositeShape<dim>>(constituent)
        ->remove_superfluous_data(updated_dof_handler,
                                  mesh_based_precalculations);
}

template <int dim>
std::pair<double, Tensor<1, dim>>
CompositeShape<dim>::apply_boolean_operations(
  std::map<unsigned int, double>         constituent_shapes_values,
  std::map<unsigned int, Tensor<1, dim>> constituent_shapes_gradients) const
{
  // The boolean operations between the shapes are applied in order
  // The last computed values are considered to be the correct values to return
  double         levelset = constituent_shapes_values[0];
  Tensor<1, dim> gradient = constituent_shapes_gradients[0];
  for (auto const &[operation_id, op_triplet] : operations)
    {
      BooleanOperation operation;
      unsigned int     first_id;
      unsigned int     second_id;
      std::tie(operation, first_id, second_id) = op_triplet;

      double         value_first_component, value_second_component;
      Tensor<1, dim> gradient_first_component{};
      Tensor<1, dim> gradient_second_component{};
      value_first_component     = constituent_shapes_values.at(first_id);
      value_second_component    = constituent_shapes_values.at(second_id);
      gradient_first_component  = constituent_shapes_gradients.at(first_id);
      gradient_second_component = constituent_shapes_gradients.at(second_id);
      switch (operation)
        {
          case BooleanOperation::Union:
            if (value_first_component < value_second_component)
              {
                constituent_shapes_values[operation_id] = value_first_component;
                constituent_shapes_gradients[operation_id] =
                  gradient_first_component;
              }
            else
              {
                constituent_shapes_values[operation_id] =
                  value_second_component;
                constituent_shapes_gradients[operation_id] =
                  gradient_second_component;
              }
            break;
          case BooleanOperation::Difference:
            if (-value_first_component > value_second_component)
              {
                constituent_shapes_values[operation_id] =
                  -value_first_component;
                constituent_shapes_gradients[operation_id] =
                  -gradient_first_component;
              }
            else
              {
                constituent_shapes_values[operation_id] =
                  value_second_component;
                constituent_shapes_gradients[operation_id] =
                  gradient_second_component;
              }
            break;
          default: // BooleanOperation::Intersection
            if (value_first_component > value_second_component)
              {
                constituent_shapes_values[operation_id] = value_first_component;
                constituent_shapes_gradients[operation_id] =
                  gradient_first_component;
              }
            else
              {
                constituent_shapes_values[operation_id] =
                  value_second_component;
                constituent_shapes_gradients[operation_id] =
                  gradient_second_component;
              }
        }
      levelset = constituent_shapes_values[operation_id];
      gradient = constituent_shapes_gradients[operation_id];
    }
  return {levelset, gradient};
}

template <int dim>
void
CompositeShape<dim>::clear_cache()
{
  Shape<dim>::clear_cache();
  // The constituents themselves don't need to have their cache cleared, because
  // everytime the value or gradient functions are called the evaluation point
  // is modified (the relationships between the constituents don't change).
}

template <int dim>
void
CompositeShape<dim>::set_layer_thickening(const double layer_thickening)
{
  this->Shape<dim>::set_layer_thickening(layer_thickening);
  for (auto const &constituent : constituents | boost::adaptors::map_values)
    constituent->set_layer_thickening(layer_thickening);
}

template <int dim>
RBFShape<dim>::RBFShape(const std::string   shape_arguments_str,
                        const Point<dim>   &position,
                        const Tensor<1, 3> &orientation)
  : Shape<dim>(1, position, orientation)
  , number_of_nodes(1)
  , iterable_nodes(1)
  , likely_nodes_map()
  , max_number_of_inside_nodes(1)
  , maximal_support_radius(1)
  , weights(1)
  , nodes_positions(1)
  , support_radii(1)
  , basis_functions(1)
  , useful_rbf_nodes(1)
{
  this->filename = shape_arguments_str;
  this->load_data_from_file();
}

template <int dim>
double
RBFShape<dim>::value_with_cell_guess(
  const Point<dim>                                    &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  double bounding_box_distance = bounding_box->value(evaluation_point);
  if (bounding_box_distance >= 0)
    return bounding_box_distance + this->effective_radius;

  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator == this->value_cache.end())
    {
      swap_iterable_nodes(cell);
      double value = this->value(evaluation_point);
      swap_iterable_nodes(cell);
      if (!this->part_of_a_composite)
        {
          this->value_cache[point_in_string] = value;
        }
      return value;
    }
  else
    {
      return this->value_cache[point_in_string];
    }
}

template <int dim>
Tensor<1, dim>
RBFShape<dim>::gradient_with_cell_guess(
  const Point<dim>                                    &evaluation_point,
  const typename DoFHandler<dim>::active_cell_iterator cell,
  const unsigned int /*component*/)
{
  Point<dim>     centered_point        = evaluation_point;
  double         bounding_box_distance = bounding_box->value(centered_point);
  Tensor<1, dim> bounding_box_gradient = bounding_box->gradient(centered_point);
  if (bounding_box_distance > 0.)
    return bounding_box_gradient;

  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->gradient_cache.find(point_in_string);
  if (iterator == this->gradient_cache.end())
    {
      swap_iterable_nodes(cell);
      Tensor<1, dim> gradient = this->gradient(evaluation_point);
      swap_iterable_nodes(cell);
      if (!this->part_of_a_composite)
        {
          this->gradient_cache[point_in_string] = gradient;
        }
      return gradient;
    }
  else
    {
      return this->gradient_cache[point_in_string];
    }
}

template <int dim>
double
RBFShape<dim>::value(const Point<dim> &evaluation_point,
                     const unsigned int /*component*/) const
{
  double bounding_box_distance = bounding_box->value(evaluation_point);
  if (bounding_box_distance >= 0)
    return bounding_box_distance + this->effective_radius -
           this->layer_thickening;

  auto point_in_string = this->point_to_string(evaluation_point);
  auto iterator        = this->value_cache.find(point_in_string);
  if (iterator != this->value_cache.end())
    return iterator->second;

  double value = std::max(bounding_box_distance, 0.0);
  double normalized_distance, basis;
  // Algorithm inspired by Optimad Bitpit. https://github.com/optimad/bitpit
  // Here we loop on every portion (of RBF nodes located in active cells close
  // to the evaluation point) and on every RBF node located in these active
  // cells.
  for (size_t portion_id = 0; portion_id < iterable_nodes.size(); portion_id++)
    {
      for (const size_t &node_id : *(std::get<2>(iterable_nodes[portion_id])))
        {
          normalized_distance =
            (evaluation_point - rotated_nodes_positions[node_id]).norm() /
            support_radii[node_id];
          // We cast the basis function value to the proper RBFBasisFunction
          basis = evaluate_basis_function(
            static_cast<enum RBFShape<dim>::RBFBasisFunction>(
              round(basis_functions[node_id])),
            normalized_distance);
          value += basis * weights[node_id];
        }
    }
  return value - this->layer_thickening;
}

template <int dim>
Tensor<1, dim>
RBFShape<dim>::gradient(const Point<dim> &evaluation_point,
                        const unsigned int /*component*/) const
{
  Point<dim>     centered_point        = evaluation_point;
  double         bounding_box_distance = bounding_box->value(centered_point);
  Tensor<1, dim> bounding_box_gradient = bounding_box->gradient(centered_point);
  if (bounding_box_distance > 0.)
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

  if (iterable_nodes.size() == 0)
    {
      throw std::logic_error(
        "Every location inside the bounding box should be covered by the radius "
        "of at least one RBF node. It this isn't the case, an error has been "
        "introduced in the code. ");
    }
  // Here we loop on every portion (of RBF nodes located in active cells close
  // to the evaluation point) and on every RBF node located in these active
  // cells.
  for (size_t portion_id = 0; portion_id < iterable_nodes.size(); portion_id++)
    {
      for (const size_t &node_id : *(std::get<2>(iterable_nodes[portion_id])))
        {
          // Calculation of dr/dx
          relative_position =
            (evaluation_point - rotated_nodes_positions[node_id]);
          distance            = (relative_position).norm();
          normalized_distance = distance / support_radii[node_id];
          if (distance > 1e-16 * this->effective_radius)
            dr_dx_derivative = relative_position / distance;
          else
            {
              // If the evaluation point overlaps with the node position, we
              // assume that the contribution of this node is 0. This assumption
              // can be made because radial basis functions are symmetrical. If
              // can be made because radial basis functions are symmetrical. If
              // the basis function is not differentiable at its node (e.g.
              // linear function), this approximation will still hold since the
              // approximated distance field is already imperfect and
              // neighboring nodes will add their contribution at the evaluation
              // point.
              for (int d = 0; d < dim; d++)
                dr_dx_derivative[d] = 0;
            }
          // Calculation of dr_norm/dr
          drnorm_dr_derivative = 1.0 / support_radii[node_id];
          // Calculation of d(basis)/dr
          // We cast the basis function value to the proper RBFBasisFunction
          dbasis_drnorm_derivative = evaluate_basis_function_derivative(
            static_cast<enum RBFShape<dim>::RBFBasisFunction>(
              round(basis_functions[node_id])),
            normalized_distance);
          // Sum
          gradient += dbasis_drnorm_derivative * drnorm_dr_derivative *
                      dr_dx_derivative * weights[node_id];
        }
    }
  return gradient;
}

template <int dim>
std::shared_ptr<Shape<dim>>
RBFShape<dim>::static_copy() const
{
  std::shared_ptr<Shape<dim>> copy =
    std::make_shared<RBFShape<dim>>(this->filename,
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
  bounding_box =
    std::make_shared<HyperRectangle<dim>>(half_lengths,
                                          this->reverse_align_and_center(
                                            bounding_box_center),
                                          this->orientation);
  this->bounding_box_half_length=half_lengths;
  this->bounding_box_center=this->reverse_align_and_center(
    bounding_box_center);
}

template <int dim>
void
RBFShape<dim>::update_precalculations(DoFHandler<dim> &dof_handler,
                                      const bool mesh_based_precalculations)
{
  if (!mesh_based_precalculations)
    return;
  // We first reset the mapping, since the grid partitioning may change between
  // calls of this function. The precalculation cost is low enough that this
  // reset does not have a significant impact on global computational cost.
  likely_nodes_map.clear();

  this->dof_handler                = &dof_handler;
  const unsigned int maximal_level = dof_handler.get_triangulation().n_levels();

  // We start by dividing the list of RBF nodes into manageable portions
  number_of_nodes = std::get<2>(iterable_nodes[0])->size();
  std::shared_ptr<std::vector<size_t>> temporary_node_numbers =
    std::make_shared<std::vector<size_t>>(number_of_nodes);
  std::iota(std::begin(*temporary_node_numbers),
            std::end(*temporary_node_numbers),
            0);
  std::map<const typename DoFHandler<dim>::cell_iterator,
           std::tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>>
    temporary_nodes_portions_map;
  temporary_nodes_portions_map.clear();
  for (unsigned int level = 0; level < maximal_level; level++)
    {
      max_number_of_inside_nodes = 1;
      const auto &cell_iterator  = dof_handler.cell_iterators_on_level(level);
      for (const auto &cell : cell_iterator)
        {
          if (cell->is_artificial_on_level())
            continue;
          if (cell->level() > 0)
            try
              {
                const auto cell_parent = cell->parent();
                const auto parent_iterator =
                  temporary_nodes_portions_map.find(cell_parent);
                if (parent_iterator != temporary_nodes_portions_map.end())
                  {
                    temporary_node_numbers =
                      std::get<2>(parent_iterator->second);
                  }
                else
                  // If the parent doesn't contain any node, the child will not
                  // either
                  continue;
              }
            catch (TriaAccessorExceptions::ExcCellHasNoParent())
              {}

          std::shared_ptr<std::vector<size_t>> current_cell_nodes =
            std::make_shared<std::vector<size_t>>();
          current_cell_nodes->reserve(max_number_of_inside_nodes);
          for (unsigned int i = 0; i < temporary_node_numbers->size(); i++)
            {
              if (cell->point_inside(
                    rotated_nodes_positions[(*temporary_node_numbers)[i]]))
                current_cell_nodes->push_back(temporary_node_numbers->at(i));
            }
          max_number_of_inside_nodes =
            std::max(max_number_of_inside_nodes, current_cell_nodes->size());
          current_cell_nodes->shrink_to_fit();
          if (!current_cell_nodes->empty())
            temporary_nodes_portions_map[cell] =
              std::make_tuple(cell->barycenter(),
                              cell->diameter(),
                              current_cell_nodes);

          // We currently ignore RBF nodes that are in no cells
          // A fix for this will be implemented shortly
        }
      // We need to remove all cells that are not needed anymore
      // from the map
      for (auto it = temporary_nodes_portions_map.cbegin();
           it != temporary_nodes_portions_map.cend();)
        {
          auto               cell       = it->first;
          const unsigned int cell_level = cell->level();
          bool cell_still_needed = cell->is_active() || cell_level == level;

          auto previous_it = it;
          ++it;
          if (!cell_still_needed)
            temporary_nodes_portions_map.erase(previous_it);
        }
    }

  // We give all the subsets to the 0 level, as an initial partitioning
  const auto &cell_iterator = dof_handler.cell_iterators_on_level(0);
  typename DoFHandler<dim>::cell_iterator temp_cell;
  std::tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>
    temp_cell_tuple{};
  for (const auto &cell : cell_iterator)
    {
      if (cell->is_artificial_on_level())
        continue;

      likely_nodes_map[cell] = std::make_shared<std::vector<
        std::
          tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>>>();
      for (auto it = temporary_nodes_portions_map.cbegin();
           it != temporary_nodes_portions_map.cend();
           it++)
        {
          temp_cell = it->first;
          if (temp_cell->is_active())
            {
              temp_cell_tuple = std::make_tuple(temp_cell->barycenter(),
                                                temp_cell->diameter(),
                                                std::get<2>(it->second));
              likely_nodes_map[cell]->push_back(temp_cell_tuple);
            }
        }
    }

  // We then treat all other levels
  for (unsigned int level = 1; level < maximal_level; level++)
    {
      const auto &cell_iterator_on_level =
        dof_handler.cell_iterators_on_level(level);
      for (const auto &cell : cell_iterator_on_level)
        {
          if (cell->is_artificial_on_level())
            continue;

          // We also check if the cell diameter is lower than the minimal
          // support radius. In that case, the likely nodes should stay more
          // or less the same.
          const bool cell_smaller_than_rbf_radius =
            (cell->diameter() < maximal_support_radius);
          if (cell_smaller_than_rbf_radius)
            likely_nodes_map[cell] = likely_nodes_map[cell->parent()];
          else
            determine_likely_nodes_for_one_cell(cell, cell->barycenter());
        }
    }

  // We need to remove all cells that are not needed anymore
  // from the map
  for (auto it = likely_nodes_map.cbegin(); it != likely_nodes_map.cend();)
    {
      auto cell              = it->first;
      bool cell_still_needed = cell->is_active() && !cell->is_artificial();

      auto previous_it = it;
      ++it;
      if (!cell_still_needed)
        likely_nodes_map.erase(previous_it);
    }

  // Here we loop on all (local or ghost) and active cells to define which RBF
  // nodes are useful to keep. We loop over the likely nodes map and take note
  // of the RBF nodes that are required in at least one active cell.
  number_of_nodes = weights.size();
  useful_rbf_nodes.clear();
  useful_rbf_nodes.resize(number_of_nodes);
  for (auto it = temporary_nodes_portions_map.cbegin();
       it != temporary_nodes_portions_map.cend();
       it++)
    {
      auto cell = it->first;
      if (!cell->is_active())
        continue;
      if (std::get<2>(it->second).use_count() <= 1)
        continue;
      for (const size_t &node_id : *(std::get<2>(it->second)))
        useful_rbf_nodes[node_id] = true;
    }
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

  bool parent_found = false;
  std::vector<
    std::tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>>
    temporary_iterable_nodes{
      {Point<dim>{}, 0, std::make_shared<std::vector<size_t>>()}};
  if (cell->level() > 0)
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
  Point<dim> temp_cell_barycenter;
  double     temp_cell_diameter;

  likely_nodes_map[cell]          = std::make_shared<std::vector<
    std::tuple<Point<dim>, double, std::shared_ptr<std::vector<size_t>>>>>();
  const size_t number_of_portions = iterable_nodes.size();
  for (size_t portion_id = 0; portion_id < number_of_portions; portion_id++)
    {
      std::tie(temp_cell_barycenter, temp_cell_diameter, std::ignore) =
        (iterable_nodes[portion_id]);
      // We only check for one support point, but we use a high security
      // factor. This allows not to loop over all support points.
      distance = (support_point - temp_cell_barycenter).norm();
      // We check if the distance is lower than 1 cell diagonal, since we
      // only check the distance with 1 support point, added to the support
      // radius
      max_distance =
        0.5 * cell_diameter + 0.5 * temp_cell_diameter + maximal_support_radius;
      if (distance < max_distance)
        {
          likely_nodes_map[cell]->push_back(iterable_nodes[portion_id]);
        }
    }

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
RBFShape<dim>::load_data_from_file()
{
  weights.clear();
  support_radii.clear();
  basis_functions.clear();
  useful_rbf_nodes.clear();
  nodes_positions.clear();
  rotated_nodes_positions.clear();

  // The following lines retrieve information regarding an RBF
  // with a given file name. Then, it converts the information
  // into a vector which is used to initialize the RBF shape.
  // All the information is concatenated into only one object so
  // that the usual initialization function can be called.
  std::map<std::string, std::vector<double>> rbf_data;
  fill_vectors_from_file(rbf_data, filename, " ");
  number_of_nodes = rbf_data["weight"].size();

  weights.swap(rbf_data["weight"]);
  support_radii.swap(rbf_data["support_radius"]);
  basis_functions.swap(rbf_data["basis_function"]);
  nodes_positions.resize(number_of_nodes);
  Point<dim> temp_point;
  for (unsigned int i = 0; i < number_of_nodes; i++)
    {
      temp_point[0] = rbf_data["node_x"][i];
      temp_point[1] = rbf_data["node_y"][i];
      if constexpr (dim == 3)
        temp_point[2] = rbf_data["node_z"][i];
      nodes_positions[i] = temp_point;
    }

  std::shared_ptr<std::vector<size_t>> temp_nodes_id =
    std::make_shared<std::vector<size_t>>(number_of_nodes);
  std::iota(std::begin(*temp_nodes_id), std::end(*temp_nodes_id), 0);
  iterable_nodes.resize(1);
  iterable_nodes[0] = {Point<dim>{},
                       std::numeric_limits<double>::max(),
                       temp_nodes_id};

  maximal_support_radius =
    *std::max_element(std::begin(support_radii), std::end(support_radii));
  initialize_bounding_box();
  rotate_nodes();
  this->effective_radius = bounding_box->half_lengths.norm();
}

template <int dim>
void
RBFShape<dim>::remove_superfluous_data(DoFHandler<dim> &updated_dof_handler,
                                       const bool mesh_based_precalculations)
{
  if (!mesh_based_precalculations)
    return;

  if (likely_nodes_map.empty())
    this->update_precalculations(updated_dof_handler,
                                 mesh_based_precalculations);

  number_of_nodes                         = useful_rbf_nodes.size();
  size_t current_number_of_kept_rbf_nodes = 0;
  for (size_t i = 0; i < number_of_nodes; i++)
    if (useful_rbf_nodes[i])
      current_number_of_kept_rbf_nodes++;

  std::vector<double>     temp_weights;
  std::vector<Point<dim>> temp_rotated_nodes_positions;
  std::vector<double>     temp_support_radii;
  std::vector<double>     temp_basis_functions;

  // Here we treat all vectors separately to keep the maximum memory low
  std::shared_ptr<std::vector<size_t>> temp_nodes_id =
    std::make_shared<std::vector<size_t>>(current_number_of_kept_rbf_nodes);
  std::iota(std::begin(*temp_nodes_id), std::end(*temp_nodes_id), 0);
  iterable_nodes.resize(1);
  iterable_nodes[0] = {Point<dim>{},
                       std::numeric_limits<double>::max(),
                       temp_nodes_id};
  std::get<2>(iterable_nodes[0])->shrink_to_fit();

  // Weights treatment
  size_t counter = 0;
  temp_weights.resize(current_number_of_kept_rbf_nodes);
  for (size_t i = 0; i < number_of_nodes; i++)
    {
      if (useful_rbf_nodes[i])
        {
          // We put back the value to the real node id
          temp_weights[counter] = weights[i];
          counter++;
        }
    }
  weights.swap(temp_weights);
  temp_weights.clear();

  // Nodes' positions treatment
  counter = 0;
  temp_rotated_nodes_positions.resize(current_number_of_kept_rbf_nodes);
  for (size_t i = 0; i < number_of_nodes; i++)
    {
      if (useful_rbf_nodes[i])
        {
          // We put back the value to the real node id
          temp_rotated_nodes_positions[counter] = rotated_nodes_positions[i];
          counter++;
        }
    }
  rotated_nodes_positions.swap(temp_rotated_nodes_positions);
  temp_rotated_nodes_positions.clear();

  // Support radii treatment
  counter = 0;
  temp_support_radii.resize(current_number_of_kept_rbf_nodes);
  for (size_t i = 0; i < number_of_nodes; i++)
    {
      if (useful_rbf_nodes[i])
        {
          // We put back the value to the real node id
          temp_support_radii[counter] = support_radii[i];
          counter++;
        }
    }
  support_radii.swap(temp_support_radii);
  temp_support_radii.clear();

  // Basis functions treatment
  counter = 0;
  temp_basis_functions.resize(current_number_of_kept_rbf_nodes);
  for (size_t i = 0; i < number_of_nodes; i++)
    {
      if (useful_rbf_nodes[i])
        {
          // We put back the value to the real node id
          temp_basis_functions[counter] = basis_functions[i];
          counter++;
        }
    }
  basis_functions.swap(temp_basis_functions);
  temp_basis_functions.clear();

  this->update_precalculations(updated_dof_handler, mesh_based_precalculations);
}

template <int dim>
void
RBFShape<dim>::rotate_nodes()
{
  rotated_nodes_positions.clear();
  rotated_nodes_positions.resize(weights.size());
  for (unsigned int i = 0; i < weights.size(); ++i)
    rotated_nodes_positions[i] =
      this->reverse_align_and_center(nodes_positions[i]);
  nodes_positions.clear();
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
RBFShape<dim>::swap_iterable_nodes(
  const typename DoFHandler<dim>::active_cell_iterator cell)
{
  // Here we check if the likely nodes have been identified
  auto iterator = likely_nodes_map.find(cell);
  if (iterator != likely_nodes_map.end())
    iterable_nodes.swap(*(likely_nodes_map[cell]));
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
    return std::pow(radius_diff * radius_diff + h_diff * h_diff, 0.5) -
           this->layer_thickening;
  else if (radius_diff <= 0 && h_diff > 0)
    return h_diff - this->layer_thickening;
  else if (radius_diff > 0 && h_diff <= 0)
    return radius_diff - this->layer_thickening;

  return std::max(radius_diff, h_diff) - this->layer_thickening;
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
  double level_set_of_cylinder_hollow = 0;
  double p_radius      = std::pow(centered_point[0] * centered_point[0] +
                               centered_point[1] * centered_point[1],
                             0.5);
  double radius_diff_o = p_radius - (radius + rectangular_base / 2);
  double radius_diff_i = p_radius - (radius - rectangular_base / 2);
  double h_diff_o      = abs(centered_point[2] - height / 2) - height / 2;

  if (radius_diff_o > 0 && h_diff_o > 0)
    level_set_of_cylinder_hollow =
      std::pow(radius_diff_o * radius_diff_o + h_diff_o * h_diff_o, 0.5);
  else if (radius_diff_o > 0 && h_diff_o <= 0)
    level_set_of_cylinder_hollow = radius_diff_o;
  else if (radius_diff_o <= 0 && radius_diff_i > 0 && h_diff_o > 0)
    level_set_of_cylinder_hollow = h_diff_o;
  else if (radius_diff_i <= 0 && h_diff_o > 0)
    level_set_of_cylinder_hollow =
      std::pow(radius_diff_i * radius_diff_i + h_diff_o * h_diff_o, 0.5);
  else if (radius_diff_i <= 0 && h_diff_o <= 0)
    level_set_of_cylinder_hollow = -radius_diff_i;
  else
    level_set_of_cylinder_hollow =
      std::max(std::max(radius_diff_o, h_diff_o), -radius_diff_i);

  return level_set_of_cylinder_hollow - this->layer_thickening;
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

  // Solve the non-linear equation to find the point on the helix with the
  // first guess
  double level_set_tube = 0;
  t                     = t_initial_guess;
  double residual = -(radius * cos(t) - centered_point[0]) * sin(t) * radius +
                    (radius * sin(t) - centered_point[1]) * cos(t) * radius +
                    (pitch / (2 * numbers::PI) * t - centered_point[2]) *
                      pitch / (2 * numbers::PI);

  unsigned int i = 0;

  // Newton iterations.
  while (abs(residual) > 1e-8 * this->effective_radius && i < 20)
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
  while (abs(residual) > 1e-8 * this->effective_radius && i < 20)
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

  return level_set - this->layer_thickening;
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
  double solid_volume = 0;

  return solid_volume;
}

template class Sphere<2>;
template class Sphere<3>;
template class HyperRectangle<2>;
template class HyperRectangle<3>;
template class Superquadric<3>;
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
template class OpenCascadeShape<2>;
template class OpenCascadeShape<3>;
template class Plane<2>;
template class Plane<3>;
// template class Shape<2>;
// template class Shape<3>;
