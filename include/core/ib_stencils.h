//
// Created by lucka on 2021-05-21.
//


#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <core/ib_particle.h>


#ifndef lethe_ib_stencils_h
#  define lethe_ib_stencils_h



using namespace dealii;


template <int dim>
class IBStencils
{
public:
  /**
  * @brief
  Return the number of points used by the interpolation stencil excluding the
  dof itself for a stencil of a given order.
  *
  * @param order, the stencil order.
  */
  virtual unsigned int
  nb_points(unsigned int order);

  /**
   * @brief
   * Define the points for the IB stencil, based on the order and the particle
   * position as well as the DOF position Depending on the order, the output
   * variable "point" change definition. In the case of stencil orders 1 to 4
   * the variable point returns the position of the DOF directly. In the case of
   * high order stencil, it returns the position of the point that is on the IB.
   * The variable "interpolation points" return the points used to define the
   * cell used for the stencil definition and the locations of the points use in
   * the stencil calculation.
   *
   * @param order, the stencil order.
   * @param p, the IB particle that cuts the cell.
   * @param dof_point, the support point of the DOF.
   */
  virtual std::tuple<Point<dim>, std::vector<Point<dim>>>
  points(unsigned int order, IBParticle<dim> p, Point<dim> dof_point);

  /**
   * @brief
   * Return the coefficient of the stencil based on the order.
   * The coefficients of the IB stencil assume a ratio of length between the
   * distance of the farthest interpolation point and the DOF and the distance
   * between the IB and the DOF of 1/8. The coefficient order goes from the
   * coefficient of the DOF to the coefficient of the farthest interpolation
   * point.
   *
   * @param order, the stencil order.
   */
  virtual std::vector<double>
  coefficients(unsigned int order);

  /**
   * @brief
   * Return the velocity of the IB used in the RHS of the equation
   *
   * @param p, the IB particle that cuts the cell.
   * @param dof_point, the support point of the DOF.
   * @param component, the stencil component of the dof (vx=0,vy=1,vz=2).
   */
  virtual double
  ib_velocity(IBParticle<dim> p, Point<dim> dof_point, unsigned int component);
};


#endif // LETHE_IB_STENCILS_H
