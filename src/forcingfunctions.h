/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef LETHE_FORCINGFUNCTIONS_H
#define LETHE_FORCINGFUNCTIONS_H

#include <deal.II/base/function.h>

// Finally, this is as in previous programs:
using namespace dealii;

template <int dim> class MMSSineForcingFunction : public Function<dim> {
public:
  MMSSineForcingFunction() : Function<dim>(3){};
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
};
template <int dim>
void MMSSineForcingFunction<dim>::vector_value(const Point<dim> &p,
                                               Vector<double> &values) const {
  assert(dim == 2);
  const double a = M_PI;
  const double x = p[0];
  const double y = p[1];
  values(0) =
      (2 * a * a * (-sin(a * x) * sin(a * x) + cos(a * x) * (cos(a * x))) *
           sin(a * y) * cos(a * y) -
       4 * a * a * sin(a * x) * sin(a * x) * sin(a * y) * cos(a * y) -
       2.0 * x) *
          (-1.) +
      a * std::pow(sin(a * x), 3.) * std::pow(sin(a * y), 2.) * std::cos(a * x);
  values(1) =
      (2 * a * a * (sin(a * y) * (sin(a * y)) - cos(a * y) * cos(a * y)) *
           sin(a * x) * cos(a * x) +
       4 * a * a * sin(a * x) * sin(a * y) * sin(a * y) * cos(a * x) -
       2.0 * y) *
          (-1) +
      a * std::pow(sin(a * x), 2.) * std::pow(sin(a * y), 3.) * std::cos(a * y);
}
template <int dim> class MMS3DSineForcingFunction : public Function<dim> {
public:
  MMS3DSineForcingFunction() : Function<dim>(3){};
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
};
template <int dim>
void MMS3DSineForcingFunction<dim>::vector_value(const Point<dim> &p,
                                                 Vector<double> &values) const {
  assert(dim == 3);
  const double a = M_PI;

  const double x = p[0];
  const double y = p[1];
  const double z = p[2];
  values(0) = 2 * a * a * (-3 * cos(2 * a * x) + 2.) * sin(a * y) * sin(a * z) *
              cos(a * y) * cos(a * z);
  values(1) = 2 * a * a * (-3 * cos(2 * a * y) + 2) * sin(a * x) * sin(a * z) *
              cos(a * x) * cos(a * z);
  values(2) = 4 * a * a * (3 * cos(2 * a * z) - 2) * sin(a * x) * sin(a * y) *
              cos(a * x) * cos(a * y);

  // Convection terms:
  values(0) += a * (2 * pow(cos(a * y), 2) - pow(cos(a * z), 2)) *
               pow(sin(a * x), 3) * pow(sin(a * y), 2) * pow(sin(a * z), 2) *
               cos(a * x);
  values(1) += a * (2 * pow(cos(a * x), 2) - pow(cos(a * z), 2)) *
               pow(sin(a * x), 2) * pow(sin(a * y), 3) * pow(sin(a * z), 2) *
               cos(a * y);
  values(2) += 2 * a * (pow(cos(a * x), 2) + pow(cos(a * y), 2)) *
               pow(sin(a * x), 2) * pow(sin(a * y), 2) * pow(sin(a * z), 3) *
               cos(a * z);

  // Pressure terms:
  // values(0) += 2*x;
  // values(1) += 2*y;
  // values(2) += 2*z;
}

template <int dim> class NoForce : public Function<dim> {
public:
  NoForce() : Function<dim>(3){};
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
};
template <int dim>
void NoForce<dim>::vector_value(const Point<dim> & /*p*/,
                                Vector<double> &values) const {
  values(0) = 0.;
  values(1) = 0.;
  if (dim == 3)
    values(2) = 0.;
}

#endif
