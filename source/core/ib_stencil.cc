#include <core/ib_stencil.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>



template <int dim>
unsigned int
IBStencil<dim>::nb_points(const unsigned int order)
{
  // The number of points used in the stencil excluding the DOF is equal to the
  // order.
  unsigned int nb_points = order;
  // In the case where the cell is used directly to to find the solution at
  // the IB only one point is needed.
  if (order > 4)
    nb_points = 1;
  return nb_points;
}

template <int dim>
void
IBStencil<dim>::p_base(const unsigned int order)
{
  using numbers::PI;
  // define the sampling point position of the stencil on the reference 1D
  // stencil 0 to 1. 1 being the position of the dof
  reference_points.resize(order + 1);
  for (unsigned int i = 0; i < order + 1; ++i)
    {
      reference_points[i] = 0.5 * (1 - cos(PI * i / order));
    }
}



template <int dim>
std::vector<double>
IBStencil<dim>::coefficients(const unsigned int order,
                             const double       length_ratio)
{
  p_base(order);
  // Initialize the coefficient vector
  std::vector<double> coef(order + 1);

  FullMatrix<double> vandermonde(order + 1, order + 1);
  FullMatrix<double> stencil(order + 1, order + 1);
  FullMatrix<double> inv_vandermonde(order + 1, order + 1);

  Vector<double> rhs(order + 1);
  // Define the Vandermonde matrix
  for (unsigned int i = 0; i < order + 1; ++i)
    {
      for (unsigned int j = 0; j < order + 1; ++j)
        {
          vandermonde[i][j] = std::pow(reference_points[i], j);
        }
      rhs[i] = std::pow(1 + length_ratio, i);
    }
  // Inverte the vandermond matrix
  inv_vandermonde.invert(vandermonde);

  // Multiply each line of the inverted matrix by (1+length_ratio)^i
  for (unsigned int i = 0; i < order + 1; ++i)
    {
      for (unsigned int j = 0; j < order + 1; ++j)
        {
          stencil[i][j] = inv_vandermonde[i][j] * rhs[i];
        }
    }
  // Sum the colones to get the coefficient
  for (unsigned int i = 0; i < order + 1; ++i)
    {
      for (unsigned int j = 0; j < order + 1; ++j)
        {
          coef[order - i] += stencil[j][i];
        }
    }

  // Fill the coefficient vector based on the order.

  if (order > 4)
    {
      // In this case the cell is directly used to find the solution at the IB
      // position. In this case only one point is needed (the position of the
      // point on the IB) and its coefficient is 1.
      coef.resize(1);
      coef[0] = 1;
    }

  return coef;
}

template <int dim>
std::tuple<Point<dim>, std::vector<Point<dim>>>
IBStencil<dim>::points(
  const unsigned int                                   order,
  const double                                         length_ratio,
  IBParticle<dim> &                                    p,
  const Point<dim> &                                   dof_point,
  const typename DoFHandler<dim>::active_cell_iterator cell_guess)
{
  // Create the vector of points used for the stencil based on the order of the
  // stencil. Also return the DOF position or the position of the point on the
  // IB depending if the cell is used directly
  p_base(order);
  Point<dim> point;
  Point<dim> surface_point;
  p.closest_surface_point(dof_point, surface_point, cell_guess);
  std::vector<Point<dim>> interpolation_points;

  if (order > 4)
    {
      // In this case the cell is directly used to find the solution at the IB
      // position. In this case only one point is needed (the position of the
      // point on the IB).
      Tensor<1, dim, double> vect_ib = dof_point - surface_point;
      point                          = surface_point;

      Point<dim, double> interpolation_point_1(dof_point +
                                               vect_ib * 1. / length_ratio);

      interpolation_points.resize(1);
      interpolation_points[0] = interpolation_point_1;
    }
  else
    {
      interpolation_points.resize(order);
      Tensor<1, dim, double> vect_ib = dof_point - surface_point;
      point                          = dof_point;
      for (unsigned int i = 1; i < order + 1; ++i)
        {
          interpolation_points[i - 1] =
            dof_point +
            vect_ib * (1 - reference_points[order - i]) / (length_ratio);
        }
    }
  return {point, interpolation_points};
}
template <int dim>
std::tuple<Point<dim>, std::vector<Point<dim>>>
IBStencil<dim>::points(const unsigned int order,
                       const double       length_ratio,
                       IBParticle<dim> &  p,
                       const Point<dim> & dof_point)
{
  // Create the vector of points used for the stencil based on the order of the
  // stencil. Also return the DOF position or the position of the point on the
  // IB depending if the cell is used directly
  p_base(order);
  Point<dim> point;
  Point<dim> surface_point;
  p.closest_surface_point(dof_point, surface_point);
  std::vector<Point<dim>> interpolation_points;

  if (order > 4)
    {
      // In this case the cell is directly used to find the solution at the IB
      // position. In this case only one point is needed (the position of the
      // point on the IB).
      Tensor<1, dim, double> vect_ib = dof_point - surface_point;
      point                          = surface_point;

      Point<dim, double> interpolation_point_1(dof_point +
                                               vect_ib * 1. / length_ratio);

      interpolation_points.resize(1);
      interpolation_points[0] = interpolation_point_1;
    }
  else
    {
      interpolation_points.resize(order);
      Tensor<1, dim, double> vect_ib = dof_point - surface_point;
      point                          = dof_point;
      for (unsigned int i = 1; i < order + 1; ++i)
        {
          interpolation_points[i - 1] =
            dof_point +
            vect_ib * (1 - reference_points[order - i]) / (length_ratio);
        }
    }
  return {point, interpolation_points};
}

template <int dim>
Point<dim>
IBStencil<dim>::point_for_cell_detection(
  IBParticle<dim> &                                    p,
  const Point<dim> &                                   dof_point,
  const typename DoFHandler<dim>::active_cell_iterator cell_guess)
{
  // Create the vector of points used for the stencil based on the order of the
  // stencil. Also return the DOF position or the position of the point on the
  // IB depending if the cell is used directly

  Point<dim> surface_point;
  p.closest_surface_point(dof_point, surface_point, cell_guess);
  Tensor<1, dim, double> vect_ib = dof_point - surface_point;
  Point<dim>             point   = dof_point + vect_ib * 1.0 / 16;

  return point;
}
template <int dim>
Point<dim>
IBStencil<dim>::point_for_cell_detection(IBParticle<dim> & p,
                                         const Point<dim> &dof_point)
{
  // Create the vector of points used for the stencil based on the order of the
  // stencil. Also return the DOF position or the position of the point on the
  // IB depending if the cell is used directly

  Point<dim> surface_point;
  p.closest_surface_point(dof_point, surface_point);
  Tensor<1, dim, double> vect_ib = dof_point - surface_point;
  Point<dim>             point   = dof_point + vect_ib * 1.0 / 16;

  return point;
}

template <int dim>
double
IBStencil<dim>::ib_velocity(
  IBParticle<dim> &                                    p,
  const Point<dim> &                                   dof_point,
  unsigned int                                         component,
  const typename DoFHandler<dim>::active_cell_iterator cell_guess)
{
  // Return the value of the IB condition for that specific stencil.
  double v_ib = 0;

  Tensor<1, 3, double> radial_vector;
  Point<dim>           closest_point;
  p.closest_surface_point(dof_point, closest_point, cell_guess);
  Point<dim> position_to_surface;
  position_to_surface    = closest_point - p.position;
  double radial_distance = position_to_surface.norm();
  if (dim == 2)
    {
      // have to do that conversion as there is no proper conversion from tensor
      // of dim 2 to 3.
      radial_vector[0] =
        radial_distance * (position_to_surface / radial_distance)[0];
      radial_vector[1] =
        radial_distance * (position_to_surface / radial_distance)[1];
      radial_vector[2]   = 0;
      Tensor<1, 3> v_rot = cross_product_3d(p.omega, radial_vector);
      if (component == 0)
        {
          // vx in 2D
          v_ib = v_rot[0] + p.velocity[0];
        }
      if (component == 1)
        {
          // vy in 2D
          v_ib = v_rot[1] + p.velocity[1];
        }
    }
  if (dim == 3)
    {
      radial_vector[0] =
        radial_distance * (position_to_surface / radial_distance)[0];
      radial_vector[1] =
        radial_distance * (position_to_surface / radial_distance)[1];
      radial_vector[2] =
        radial_distance * (position_to_surface / radial_distance)[2];
      Tensor<1, 3> v_rot = cross_product_3d(p.omega, radial_vector);
      if (component == 0)
        {
          // vx in 3D
          v_ib = v_rot[0] + p.velocity[0];
        }
      if (component == 1)
        {
          // vy in 3D
          v_ib = v_rot[1] + p.velocity[1];
        }
      if (component == 2)
        {
          // vz in 3D
          v_ib = v_rot[2] + p.velocity[2];
        }
    }

  return v_ib;
}
template <int dim>
double
IBStencil<dim>::ib_velocity(IBParticle<dim> & p,
                            const Point<dim> &dof_point,
                            unsigned int      component)
{
  // Return the value of the IB condition for that specific stencil.
  double v_ib = 0;

  Tensor<1, 3, double> radial_vector;
  Point<dim>           closest_point;
  p.closest_surface_point(dof_point, closest_point);
  Point<dim> position_to_surface;
  position_to_surface    = closest_point - p.position;
  double radial_distance = position_to_surface.norm();
  if (dim == 2)
    {
      // have to do that conversion as there is no proper conversion from tensor
      // of dim 2 to 3.
      radial_vector[0] =
        radial_distance * (position_to_surface / radial_distance)[0];
      radial_vector[1] =
        radial_distance * (position_to_surface / radial_distance)[1];
      radial_vector[2]   = 0;
      Tensor<1, 3> v_rot = cross_product_3d(p.omega, radial_vector);
      if (component == 0)
        {
          // vx in 2D
          v_ib = v_rot[0] + p.velocity[0];
        }
      if (component == 1)
        {
          // vy in 2D
          v_ib = v_rot[1] + p.velocity[1];
        }
    }
  if (dim == 3)
    {
      radial_vector[0] =
        radial_distance * (position_to_surface / radial_distance)[0];
      radial_vector[1] =
        radial_distance * (position_to_surface / radial_distance)[1];
      radial_vector[2] =
        radial_distance * (position_to_surface / radial_distance)[2];
      Tensor<1, 3> v_rot = cross_product_3d(p.omega, radial_vector);
      if (component == 0)
        {
          // vx in 3D
          v_ib = v_rot[0] + p.velocity[0];
        }
      if (component == 1)
        {
          // vy in 3D
          v_ib = v_rot[1] + p.velocity[1];
        }
      if (component == 2)
        {
          // vz in 3D
          v_ib = v_rot[2] + p.velocity[2];
        }
    }

  return v_ib;
}

template class IBStencil<2>;
template class IBStencil<3>;
