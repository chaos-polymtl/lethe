#include <core/ib_stencil.h>


template <int dim>
unsigned int
IBStencil<dim>::nb_points(unsigned int order)
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
std::vector<double>
IBStencil<dim>::coefficients(unsigned int order)
{
  // Initialize the coefficient vector
  std::vector<double> coef(order + 1);

  // Fill the coefficient vector based on the order.
  if (order == 1)
    {
      coef[0] = 9;
      coef[1] = -8;
    }
  if (order == 2)
    {
      coef[0] = 153;
      coef[1] = -288;
      coef[2] = 136;
    }
  if (order == 3)
    {
      coef[0] = 2925;
      coef[1] = -8424;
      coef[2] = 8100;
      coef[3] = -2600;
    }
  if (order == 4)
    {
      coef[0] = 58905;
      coef[1] = -228480;
      coef[2] = 332640;
      coef[3] = -215424;
      coef[4] = 52360;
    }
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
IBStencil<dim>::points(unsigned int    order,
                       IBParticle<dim> p,
                       Point<dim>      dof_point)
{
  // Create the vector of points used for the stencil based on the order of the
  // stencil. Also return the DOF position or the position of the point on the
  // IB depending if the cell is used directly
  Point<dim>              point;
  std::vector<Point<dim>> interpolation_points;

  if (order == 1)
    {
      point = dof_point;
      Tensor<1, dim, double> vect_ib =
        (dof_point - p.position -
         p.radius * (dof_point - p.position) / (dof_point - p.position).norm());

      Point<dim, double> interpolation_point_1(dof_point + vect_ib * 1 / 8);

      interpolation_points.resize(1);
      interpolation_points[0] = interpolation_point_1;
    }
  if (order == 2)
    {
      point = dof_point;
      Tensor<1, dim, double> vect_ib =
        (dof_point - p.position -
         p.radius * (dof_point - p.position) / (dof_point - p.position).norm());


      Point<dim, double> interpolation_point_1(dof_point + vect_ib * 1. / 16.);

      Point<dim, double> interpolation_point_2(dof_point + vect_ib * 1. / 8.);

      interpolation_points.resize(2);
      interpolation_points[0] = interpolation_point_1;
      interpolation_points[1] = interpolation_point_2;
    }
  if (order == 3)
    {
      point = dof_point;
      Tensor<1, dim, double> vect_ib =
        (dof_point - p.position -
         p.radius * (dof_point - p.position) / (dof_point - p.position).norm());

      Point<dim, double> interpolation_point_1(dof_point + vect_ib * 1. / 24.);

      Point<dim, double> interpolation_point_2(dof_point + vect_ib * 1 / 12.);

      Point<dim, double> interpolation_point_3(dof_point + vect_ib * 1. / 8.);

      interpolation_points.resize(3);
      interpolation_points[0] = interpolation_point_1;
      interpolation_points[1] = interpolation_point_2;
      interpolation_points[2] = interpolation_point_3;
    }
  if (order == 4)
    {
      point = dof_point;
      Tensor<1, dim, double> vect_ib =
        (dof_point - p.position -
         p.radius * (dof_point - p.position) / (dof_point - p.position).norm());

      Point<dim, double> interpolation_point_1(dof_point + vect_ib * 1. / 32.);

      Point<dim, double> interpolation_point_2(dof_point + vect_ib * 1 / 16.);

      Point<dim, double> interpolation_point_3(dof_point + vect_ib * 3. / 32.);

      Point<dim, double> interpolation_point_4(dof_point + vect_ib * 1. / 8.);

      interpolation_points.resize(4);
      interpolation_points[0] = interpolation_point_1;
      interpolation_points[1] = interpolation_point_2;
      interpolation_points[2] = interpolation_point_3;
      interpolation_points[3] = interpolation_point_4;
    }
  if (order > 4)
    {
      // In this case the cell is directly used to find the solution at the IB
      // position. In this case only one point is needed (the position of the
      // point on the IB).
      Tensor<1, dim, double> vect_ib =
        (dof_point - p.position -
         p.radius * (dof_point - p.position) / (dof_point - p.position).norm());

      point = dof_point - vect_ib;

      Point<dim, double> interpolation_point_1(dof_point + vect_ib * 1. / 8.);

      interpolation_points.resize(1);
      interpolation_points[0] = interpolation_point_1;
    }
  return {point, interpolation_points};
}


template <int dim>
double
IBStencil<dim>::ib_velocity(IBParticle<dim> p,
                            Point<dim>      dof_point,
                            unsigned int    component)
{
  // Return the value of the IB condition for that specific stencil.
  double v_ib = 0;

  Tensor<1, 3, double> radial_vector;
  if (dim == 2)
    {
      // have to do that conversion as there is no proper conversion from tensor
      // of dim 2 to 3.
      radial_vector[0]   = p.radius * ((dof_point - p.position) /
                                     (dof_point - p.position).norm())[0];
      radial_vector[1]   = p.radius * ((dof_point - p.position) /
                                     (dof_point - p.position).norm())[1];
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
      radial_vector[0]   = p.radius * ((dof_point - p.position) /
                                     (dof_point - p.position).norm())[0];
      radial_vector[1]   = p.radius * ((dof_point - p.position) /
                                     (dof_point - p.position).norm())[1];
      radial_vector[2]   = p.radius * ((dof_point - p.position) /
                                     (dof_point - p.position).norm())[2];
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
