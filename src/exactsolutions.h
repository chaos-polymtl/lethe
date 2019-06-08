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


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

using namespace dealii;
#ifndef LETHE_EXACTSOLUTIONS_H
#define LETHE_EXACTSOLUTIONS_H

template<int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
    ExactSolutionMMS() : Function<dim>(3) {}
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &values) const;
private:
    double time;
};
template<int dim>
void ExactSolutionMMS<dim>::vector_value(const Point<dim> &p,
                                         Vector<double> &values) const
{
    assert(dim==2);
    const double a = M_PI;
    double x = p[0];
    double y = p[1];
    values(0) = sin(a*x)*sin(a*x)*cos(a*y)*sin(a*y);
    values(1) = -cos(a*x)*sin(a*x)*sin(a*y)*sin(a*y);
}


#endif




