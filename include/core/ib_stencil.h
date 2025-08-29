// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ib_stencil_h
#define lethe_ib_stencil_h

#include <core/ib_particle.h>

#include <vector>



using namespace dealii;

/**
 * @brief Define where the sharp immersed-boundaries are imposed on Eulerian grid.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved.
 */
template <int dim>
class IBStencil
{
public:
  /**
   * @brief Return the number of points used by the interpolation stencil (excluding the
   *  dof itself) for a stencil of a given order.
   *
   * @param[in] order Stencil order.
   *
   * @return Number of points used by the interpolation stencil.
   */
  virtual unsigned int
  number_of_interpolation_support_points(const unsigned int order);

  /**
   * @brief Find the DOF at which we try to apply the IB by interpolation
   * and the location of the support points to do so.
   *
   * Define the points for the IB stencil, based on the order and the particle
   * position as well as the DOF position. Depending on the order, the output
   * variable "point" change definition. In the case of stencil orders 1 to 4
   * the variable point returns the position of the DOF directly. In the case of
   * high order stencil, it returns the position of the point that is on the IB.
   * The variable "interpolation points" return the points used to define the
   * cell used for the stencil definition and the locations of the points used
   * in the stencil calculation.
   *
   * @param[in] order Stencil order.
   *
   * @param[in] length_ratio Length ratio between the 2 side of the stencil.
   *
   * @param[in] p IB particle that cuts the cell.
   *
   * @param[in] dof_point Support point of the DOF.
   *
   * @param[in] cell_guess Guess of the cell containing the evaluation point,
   * which is useful to reduce computation time.
   *
   * @return Tuple containing (point, interpolation_points).
   */
  virtual std::tuple<Point<dim>, std::vector<Point<dim>>>
  support_points_for_interpolation(
    const unsigned int                                    order,
    const double                                          length_ratio,
    IBParticle<dim>                                      &p,
    const Point<dim>                                     &dof_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * @brief See overloaded function.
   */
  virtual std::tuple<Point<dim>, std::vector<Point<dim>>>
  support_points_for_interpolation(const unsigned int order,
                                   const double       length_ratio,
                                   IBParticle<dim>   &p,
                                   const Point<dim>  &dof_point);

  /**
   * @brief Define the point used to define the cell used for the stencil calculation.
   *
   * @param[in] p IB particle that cuts the cell.
   *
   * @param[in] dof_point support point of the DOF.
   *
   * @param[in] cell_guess Guess of the cell containing the evaluation point,
   * which is useful to reduce computation time.
   */
  virtual Point<dim>
  point_for_cell_detection(
    IBParticle<dim>                                      &p,
    const Point<dim>                                     &dof_point,
    const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * @brief See overloaded function.
   */
  virtual Point<dim>
  point_for_cell_detection(IBParticle<dim> &p, const Point<dim> &dof_point);

  /**
   * @brief Return the coefficients (weights) of the stencil.
   *
   * The coefficients of the IB stencil assume a length between the
   * distance of the farthest interpolation point to the DOF, and the distance
   * between the IB and the DOF of 1/8. The coefficients are defined from the
   * coefficient of the DOF to the coefficient of the farthest interpolation
   * point (sorted by increasing distance).
   *
   * @param[in] order Stencil order.
   *
   * @param[in] length_ratio Length ratio between the 2 side of the stencil.
   */
  virtual std::vector<double>
  coefficients(const unsigned int order, const double length_ratio);

  /**
   * @brief Return the velocity of the IB used in the RHS of the equation.
   *
   * @param[in] p IB particle that cuts the cell.
   *
   * @param[in] dof_point Support point of the DOF.
   *
   * @param[in] component Stencil component of the dof (vx=0,vy=1,vz=2).
   *
   * @param[in] cell_guess Guess of the cell containing the evaluation point,
   * which is useful to reduce computation time.
   */
  virtual double
  ib_velocity(IBParticle<dim>                                      &p,
              const Point<dim>                                     &dof_point,
              const unsigned int                                    component,
              const typename DoFHandler<dim>::active_cell_iterator &cell_guess);

  /**
   * @brief See overloaded function.
   */
  virtual double
  ib_velocity(IBParticle<dim>   &p,
              const Point<dim>  &dof_point,
              const unsigned int component);

private:
  /**
   * @brief Return the points on the reference 1D element for the polynomial base.
   *
   * @param[in] order Order of the stencil.
   */
  void
  p_base(unsigned int order);

  /**
   * @brief Points of the reference 1D stencil.
   */
  std::vector<double> reference_points;
};



#endif // LETHE_IB_STENCILS_H
