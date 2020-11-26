/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 -  by the Lethe authors
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

 *
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_per_hills_grid_h
#define lethe_per_hills_grid_h

#include <deal.II/base/config.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_in.h>

#include <sstream>
#include <string>

using namespace dealii;

/**
 * @brief PeriodicHillsGrid. The PeriodicHillsGrid class creates an
 * hyper_rectangle and transforms it to obtain the hill geometry with
 * the hill_geometry function. It also attaches a manifold to the geometry.
 */
template <int dim, int spacedim>
class PeriodicHillsGrid
{
public:
  /**
   * @brief Constructor for the PeriodicHillsGrid.
   *
   * @param grid_arguments. A string with 5 parameters
   *  spacing_y : allows to control the shifting of horizontal line [0 - 1]
   *  alpha : allows to elongate while keeping the same flat region [0.5 - 3]
   *  repetitions_x : number of separation of cells in x before refinement
   *  repetitions_y : number of separation of cells in y before refinement
   *  repetitions_z : number of separation of cells in z before refinement
   */
  PeriodicHillsGrid(const std::string &grid_arguments);
  Point<spacedim> static hill_geometry(const Point<spacedim> &p,
                                       double                 spacing_y,
                                       double                 alpha);

  /**
   * @brief make_grid. The make_grid function generates a hyper rectangle of the size of the domain
   * and then transforms it to the hill geometry. It also constructs the
   * geometry manifold with FunctionManifold and finally sets the manifold.
   *
   * @param triangulation. The triangulation object on which the grid is generated
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
  double      spacing_y;
  int         repetitions_x;
  int         repetitions_y;
  int         repetitions_z;
  double      alpha;
};

/**
 * @brief The push_forward & the pull_back classes create the vector_value functions
 * gradient function and because it inherits from Function<spacedim>. (formula
 * is currently Euler and can be changed. See AutoDerivativeFunction
 * documentation)
 */
template <int dim, int spacedim>
class PeriodicHillsPushForward : public AutoDerivativeFunction<spacedim>
{
public:
  PeriodicHillsPushForward(double spacing_y, double alpha)
    : AutoDerivativeFunction<spacedim>(1e-6, spacedim)
    , spacing_y(spacing_y)
    , alpha(alpha)
  {}

  /**
   * @brief vector_value. This function is used to construct the geometry manifold.
   * It changes the original point (op) of the transformed hyper_rectangle
   * with hill_geometry to a new point (np) of the hill grid with per_hills_grid
   * function.
   *
   * @param p. A point in space
   *
   * @param values. The vector of values which will be calculated at the position p.
   */
  virtual void
  vector_value(const Point<spacedim> &p, Vector<double> &values) const override;

  /**
   * @brief value. The value function does the same thing than vector_value for one component.
   * This implementation is needed to use the gradient function inherited by
   * AutoDerivativeFunction.
   *
   * @param op. An original point in space
   *
   * @param component. The component of the point (x=0, y=1, z=2)
   */
  virtual double
  value(const Point<spacedim> &p, const unsigned int component) const override;

private:
  double spacing_y;
  double alpha;
};


template <int dim, int spacedim>
class PeriodicHillsPullBack : public AutoDerivativeFunction<spacedim>
{
public:
  PeriodicHillsPullBack(double spacing_y, double alpha)
    : AutoDerivativeFunction<spacedim>(1e-6, spacedim)
    , spacing_y(spacing_y)
    , alpha(alpha)
  {}

  /**
   * @brief vector_value. This vector_value function is used to construct the
   * geometry manifold. It changes the new point (np) of the hill grid to the
   * original point (op) of the transformed hyper_rectangle grid. This function
   * is mandatory to use FunctionManifold. First, it finds the minimum value of
   * y depending the x position and then calculates the op with the inverse of
   * the transformation done by hill_geometry function.
   *
   * @param np. A point in space.
   *
   * @param values. The vector of values which will be calculated at the position p.
   */
  virtual void
  vector_value(const Point<spacedim> &np,
               Vector<double> &       values) const override;

  /**
   * @brief value. The value function does the same thing than vector_value for one component.
   * This implementation is needed to use the gradient function inherited by
   * AutoDerivativeFunction.
   *
   * @param p. An original point in space
   *
   * @param component. The component of the point (x=0, y=1, z=2)
   *
   */
  virtual double
  value(const Point<spacedim> &np, const unsigned int component) const override;

private:
  double spacing_y;
  double alpha;
};


#endif
