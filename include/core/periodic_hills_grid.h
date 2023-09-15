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

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <sstream>

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
  PeriodicHillsGrid(const std::string &grid_arguments);
  Point<spacedim> static hill_geometry(const Point<spacedim> &p,
                                       double                 spacing_y,
                                       double                 alpha);
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
  virtual void
  vector_value(const Point<spacedim> &p, Vector<double> &values) const override;
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

  virtual void
  vector_value(const Point<spacedim> &np,
               Vector<double>        &values) const override;

  virtual double
  value(const Point<spacedim> &np, const unsigned int component) const override;

private:
  double spacing_y;
  double alpha;
};

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

template <int dim, int spacedim>
PeriodicHillsGrid<dim, spacedim>::PeriodicHillsGrid(
  const std::string &grid_arguments)
{
  this->grid_arguments = grid_arguments;

  // Separate arguments of the string
  std::vector<std::string> arguments;
  std::stringstream        s_stream(grid_arguments);
  while (s_stream.good())
    {
      std::string substr;
      getline(s_stream, substr, ';');
      arguments.push_back(substr);
    }

  std::vector<double> arguments_double =
    dealii::Utilities::string_to_double(arguments);
  spacing_y     = arguments_double[0];
  alpha         = arguments_double[1];
  repetitions_x = arguments_double[2];
  repetitions_y = arguments_double[3];
  if (dim == 3)
    repetitions_z = arguments_double[4];

  if (abs(alpha - 1) < 1e-6)
    alpha = int(alpha);
}

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
template <int dim, int spacedim>
void
PeriodicHillsPushForward<dim, spacedim>::vector_value(
  const Point<spacedim> &op,
  Vector<double>        &values) const
{
  const Point<spacedim> np =
    PeriodicHillsGrid<dim, spacedim>::hill_geometry(op, spacing_y, alpha);

  values(0) = np[0];
  values(1) = np[1];

  if (spacedim == 3)
    values(2) = np[2];
}

/**
 * @brief value. The value function does the same thing than vector_value for one component.
 * This implementation is needed to use the gradient function inherited by
 * AutoDerivativeFunction.
 *
 * @param op. An original point in space
 *
 * @param component. The component of the point (x=0, y=1, z=2)
 */
template <int dim, int spacedim>
double
PeriodicHillsPushForward<dim, spacedim>::value(
  const Point<spacedim> &op,
  const unsigned int     component) const
{
  const Point<spacedim> np =
    PeriodicHillsGrid<dim, spacedim>::hill_geometry(op, spacing_y, alpha);
  return np[component];
}

/**
 * \brief vector_value. This vector_value function is used to construct the
 * geometry manifold. It changes the new point (np) of the hill grid to the
 * original point (op) of the transformed hyper_rectangle grid. This function
 * is mandatory to use FunctionManifold. First, it finds the minimum value of
 * y depending the x position and then calculates the op with the inverse of
 * the transformation done by hill_geometry function.
 *
 * \param np. A point in space.
 *
 * @param values. The vector of values which will be calculated at the position p.
 */
template <int dim, int spacedim>
void
PeriodicHillsPullBack<dim, spacedim>::vector_value(const Point<spacedim> &np,
                                                   Vector<double> &values) const
{
  double x                  = np[0];
  double max_y              = 3.035;
  double max_x              = 9.0;
  double flat_region_length = 5.142;
  double left_hill          = 1.929;
  double right_hill         = 7.071;
  double min_y;

  // Reversing elongation and shifting of x lines
  if (alpha != 1)
    {
      if (x < left_hill * alpha)
        x = (x / alpha);
      else if (x > alpha * left_hill + flat_region_length)
        x = (x - flat_region_length - alpha * left_hill) / alpha + right_hill;
      else
        x = x - (alpha * left_hill) + left_hill;
    }

  if (alpha > 1)
    {
      if (x < max_x / 2)
        x = (-(1 - 0.5) +
             std::sqrt(std::pow((1 - 0.5), 2) - (4 * (1 / max_x) * -x))) /
            (2 / max_x);
      else if (x > max_x / 2 && x < max_x)
        x = (-(1 + 1.5) + std::sqrt(std::pow((1 + 1.5), 2) -
                                    (4 * (-1 / max_x) * (-0.5 * max_x - x)))) /
            (2 * -1 / max_x);
    }

  // Reversing polynomial transformation and shifting of y lines
  if (spacedim == 2)
    {
      min_y =
        PeriodicHillsGrid<dim, spacedim>::hill_geometry(Point<spacedim>(x, 0),
                                                        spacing_y,
                                                        alpha)[1];
    }
  else if (spacedim == 3)
    {
      min_y = PeriodicHillsGrid<dim, spacedim>::hill_geometry(
        Point<spacedim>(x, 0, np[2]), spacing_y, alpha)[1];
      values(2) = np[2];
    }

  double y = (np[1] - min_y) / (1 - min_y / max_y);

  if (y < max_y / 2 && spacing_y != 0)
    y = (-(1 - 0.5 * spacing_y) + std::sqrt(std::pow((1 - 0.5 * spacing_y), 2) -
                                            (4 * (spacing_y / max_y) * -y))) /
        (2 * spacing_y / max_y);
  else if (y > max_y / 2 && y < max_y && spacing_y != 0)
    y =
      (-(1 + 1.5 * spacing_y) +
       std::sqrt(std::pow((1 + 1.5 * spacing_y), 2) -
                 (4 * (-spacing_y / max_y) * (-0.5 * spacing_y * max_y - y)))) /
      (2 * -spacing_y / max_y);

  values(0) = x;
  values(1) = y;
}

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
template <int dim, int spacedim>
double
PeriodicHillsPullBack<dim, spacedim>::value(const Point<spacedim> &np,
                                            const unsigned int component) const
{
  double x                  = np[0];
  double max_y              = 3.035;
  double max_x              = 9.0;
  double flat_region_length = 5.142;
  double left_hill          = 1.929;
  double right_hill         = 7.071;
  double min_y;

  if (alpha != 1)
    {
      if (x < left_hill * alpha)
        x = (x / alpha);
      else if (x > alpha * left_hill + flat_region_length)
        x = (x - flat_region_length - alpha * left_hill) / alpha + right_hill;
      else
        x = x - (alpha * left_hill) + left_hill;
    }

  if (alpha > 1)
    {
      if (x < max_x / 2)
        x = (-(1 - 0.5) +
             std::sqrt(std::pow((1 - 0.5), 2) - (4 * (1 / max_x) * -x))) /
            (2 / max_x);
      else if (x > max_x / 2 && x <= max_x)
        x = (-(1 + 1.5) + std::sqrt(std::pow((1 + 1.5), 2) -
                                    (4 * (-1 / max_x) * (-0.5 * max_x - x)))) /
            (2 * -1 / max_x);
    }

  if (spacedim == 2)
    min_y =
      PeriodicHillsGrid<dim, spacedim>::hill_geometry(Point<spacedim>(x, 0),
                                                      spacing_y,
                                                      alpha)[1];
  else if (spacedim == 3)
    min_y = PeriodicHillsGrid<dim, spacedim>::hill_geometry(
      Point<spacedim>(x, 0, np[2]), spacing_y, alpha)[1];


  double y = (np[1] - min_y) / (1 - min_y / max_y);

  if (y < max_y / 2 && spacing_y != 0)
    y = (-(1 - 0.5 * spacing_y) + std::sqrt(std::pow((1 - 0.5 * spacing_y), 2) -
                                            (4 * (spacing_y / max_y) * -y))) /
        (2 * spacing_y / max_y);
  else if (y > max_y / 2 && y < max_y && spacing_y != 0)
    y =
      (-(1 + 1.5 * spacing_y) +
       std::sqrt(std::pow((1 + 1.5 * spacing_y), 2) -
                 (4 * (-spacing_y / max_y) * (-0.5 * spacing_y * max_y - y)))) /
      (2 * -spacing_y / max_y);

  Point<spacedim> op;
  if (spacedim == 2)
    op = {x, y};
  else if (spacedim == 3)
    op = {x, y, np[2]};

  return op[component];
}

/**
 * @brief The hill_geometry function calculates all the domain of the geometry with 6
 * polynomials depending the x position. (See Hill Geometry Definition file :
 * https://turbmodels.larc.nasa.gov/Other_LES_Data/2dhill_periodic.html)
 * This code has nondimensionalized geometry, but the coefficients provided
 * need a hill height of 28.
 * This function also does a gradual shifting of the horizontal lines
 * prior to have smaller element on the bottom of the geometry where results
 * are more important.
 *
 * @param p A point in space which will be adapted to the periodic hill geometry
 *
 * @param param Non-linear solver parameters
 *
 */
template <int dim, int spacedim>
Point<spacedim>
PeriodicHillsGrid<dim, spacedim>::hill_geometry(const Point<spacedim> &p,
                                                double spacing_y,
                                                double alpha)
{
  double H = 28; // Height dimension to use with polynomials
  double x = p[0] * H, y = p[1] * H;
  double max_y              = 3.035 * H;
  double max_x              = 9.0 * H;
  double flat_region_length = 5.142 * H;
  double left_hill          = 1.929 * H;
  double right_hill         = 7.071 * H;
  double pos_x_left;
  double pos_x_right;
  double pos_y_bottom;
  double pos_y_top;
  double pos_y;
  double y_0;
  double diff;

  // Gradual spacing and swifting depending on x position if the geometry
  // is elongated
  if (alpha > 1)
    {
      if (x < max_x / 2)
        {
          pos_x_left = x / -max_x + 0.5;
          x -= pos_x_left * x;
        }
      else if (x > max_x / 2 && x < max_x)
        {
          pos_x_right = x / max_x - 0.5;
          x += pos_x_right * (max_x - x);
        }
    }

  // Polynomial coefficients :
  const double a1 = 2.800000000000E+01, b1 = 0.000000000000E+00,
               c1 = 6.775070969851E-03, d1 = -2.124527775800E-03;
  const double a2 = 2.507355893131E+01, b2 = 9.754803562315E-01,
               c2 = -1.016116352781E-01, d2 = 1.889794677828E-03;
  const double a3 = 2.579601052357E+01, b3 = +8.206693007457E-01,
               c3 = -9.055370274339E-02, d3 = 1.626510569859E-03;
  const double a4 = 4.046435022819E+01, b4 = -1.379581654948E+00,
               c4 = 1.945884504128E-02, d4 = -2.070318932190E-04;
  const double a5 = 1.792461334664E+01, b5 = +8.743920332081E-01,
               c5 = -5.567361123058E-02, d5 = 6.277731764683E-04;
  const double a6 = 5.639011190988E+01, b6 = -2.010520359035E+00,
               c6 = 1.644919857549E-02, d6 = 2.674976141766E-05;

  double new_x = (max_x - x); // x for the left side of the geometry

  // Gradual spacing and swifting depending on y position
  if (y < max_y / 2 && (std::abs(spacing_y) > 0))
    {
      pos_y_bottom = y / -max_y + 0.5;
      y -= spacing_y * pos_y_bottom * y;
    }
  else if (y > max_y / 2 && y < max_y && (std::abs(spacing_y) > 0))
    {
      pos_y_top = y / max_y - 0.5;
      y += spacing_y * pos_y_top * (max_y - y);
    }

  pos_y = y / -max_y + 1;

  // Polynomial equations :
  if (x >= 0 && x < 9)
    {
      y += pos_y * (a1 + b1 * x + c1 * std::pow(x, 2) + d1 * std::pow(x, 3));

      // Checking if y is under 28 and correction
      y_0  = (a1 + b1 * x + c1 * std::pow(x, 2) + d1 * std::pow(x, 3));
      diff = y_0 - 28.0;
      if (diff > 0)
        y -= pos_y * diff;
    }
  else if (x >= 9 && x < 14)
    y += pos_y * (a2 + b2 * x + c2 * std::pow(x, 2) + d2 * std::pow(x, 3));
  else if (x >= 14 && x < 20)
    y += pos_y * (a3 + b3 * x + c3 * std::pow(x, 2) + d3 * std::pow(x, 3));
  else if (x >= 20 && x < 30)
    y += pos_y * (a4 + b4 * x + c4 * std::pow(x, 2) + d4 * std::pow(x, 3));
  else if (x >= 30 && x < 40)
    y += pos_y * (a5 + b5 * x + c5 * std::pow(x, 2) + d5 * std::pow(x, 3));
  else if (x >= 40 && x < 54)
    {
      y += pos_y * (a6 + b6 * x + c6 * std::pow(x, 2) + d6 * std::pow(x, 3));
      if (y < 0)
        y = 0;
    }
  else if (x <= 252 && x > 243)
    {
      y += pos_y * (a1 + b1 * new_x + c1 * std::pow(new_x, 2) +
                    d1 * std::pow(new_x, 3));
      y_0  = (a1 + b1 * x + c1 * std::pow(x, 2) + d1 * std::pow(x, 3));
      diff = y_0 - 28.0;
      if (diff > 0)
        y -= pos_y * diff;
    }
  else if (x <= 243 && x > 238)
    y += pos_y *
         (a2 + b2 * new_x + c2 * std::pow(new_x, 2) + d2 * std::pow(new_x, 3));
  else if (x <= 238 && x > 232)
    y += pos_y *
         (a3 + b3 * new_x + c3 * std::pow(new_x, 2) + d3 * std::pow(new_x, 3));
  else if (x <= 232 && x > 222)
    y += pos_y *
         (a4 + b4 * new_x + c4 * std::pow(new_x, 2) + d4 * std::pow(new_x, 3));
  else if (x <= 222 && x > 212)
    y += pos_y *
         (a5 + b5 * new_x + c5 * std::pow(new_x, 2) + d5 * std::pow(new_x, 3));
  else if (x <= 212 && x > 198)
    {
      y += pos_y * (a6 + b6 * new_x + c6 * std::pow(new_x, 2) +
                    d6 * std::pow(new_x, 3));
      if (y < 0)
        y = 0;
    }
  else
    y += 0;

  // Elongation of the geometry with the alpha factor
  // Note : The length of the flat region is always the same length
  if (x < left_hill)
    x = alpha * x;
  else if (x > right_hill)
    x = alpha * (x - right_hill) + flat_region_length + alpha * left_hill;
  else
    x = (x - left_hill) + (alpha * left_hill);

  Point<spacedim> q;
  q[0] = (x / H);
  q[1] = (y / H);

  if (spacedim == 3)
    q[2] = p[2];

  return q;
}

/**
 * @brief make_grid. The make_grid function generates a hyper rectangle of the size of the domain
 * and then transforms it to the hill geometry. It also constructs the
 * geometry manifold with FunctionManifold and finally sets the manifold.
 *
 * @param triangulation. The triangulation object on which the grid is generated
 */
template <int dim, int spacedim>
void
PeriodicHillsGrid<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> &triangulation)
{
  std::vector<unsigned int> repetitions(2);
  repetitions[0] = repetitions_x;
  repetitions[1] = repetitions_y;

  if (dim == 2)
    {
      if (alpha != 1 && (repetitions_x > 1 || repetitions_y > 1))
        throw std::logic_error(
          "When alpha parameter is not 1, repetition parameters should be all set to 1");

      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                repetitions,
                                                Point<dim>(0.0, 0.0),
                                                Point<dim>(9.0, 3.035),
                                                true);
    }
  else if (dim == 3)
    {
      if (alpha != 1 &&
          (repetitions_x > 1 || repetitions_y > 1 || repetitions_z > 1))
        throw std::logic_error(
          "When alpha parameter is not 1, repetition parameters should be all set to 1");
      repetitions.push_back(repetitions_z);
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                repetitions,
                                                Point<dim>(0.0, 0.0, 0.0),
                                                Point<dim>(9.0, 3.035, 4.5),
                                                true);
    }

  // Transformation of the geometry with the hill geometry
  GridTools::transform(
    [this](const Point<spacedim> &p) {
      return this->hill_geometry(p, spacing_y, alpha);
    },
    triangulation);

  // Manifold construction
  static const FunctionManifold<dim, spacedim, spacedim> manifold_func(
    std::make_unique<PeriodicHillsPushForward<dim, spacedim>>(spacing_y, alpha),
    std::make_unique<PeriodicHillsPullBack<dim, spacedim>>(spacing_y, alpha));
  triangulation.set_manifold(1, manifold_func);
  triangulation.set_all_manifold_ids(1);
}

#endif
