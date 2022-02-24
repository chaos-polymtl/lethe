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
 */

#include "core/bdf.h"

#include <deal.II/lac/vector.h>


Vector<double>
bdf_coefficients(const unsigned int p, const std::vector<double> dt)
{
  // There should be at least p time steps
  assert(dt.size() >= p);

  // Create a time table for the bdf formula
  Vector<double> times(p + 1);
  for (unsigned int i = 0; i < p + 1; ++i)
    {
      times[i] = 0.;
      for (unsigned int j = 0; j < i; ++j)
        times[i] -= dt[j];
    }

  // The alphas are the bdf coefficients
  Vector<double> alpha(p + 1);
  alpha = 0;

  for (unsigned int j = 1; j < p + 1; ++j)
    {
      double factor = 1.;
      for (unsigned int i = 1; i < j; ++i)
        factor *= times[0] - times[i];

      Vector<double> term(p + 1);
      term.equ(factor, delta(p, 0, j, times));
      alpha += term;
    }
  return alpha;
}

Vector<double>
bdf_coefficients(const Parameters::SimulationControl::TimeSteppingMethod method,
                 const std::vector<double>                               dt)
{
  switch (method)
    {
      case (Parameters::SimulationControl::TimeSteppingMethod::bdf1):
        return bdf_coefficients(1, dt);
      case (Parameters::SimulationControl::TimeSteppingMethod::steady_bdf):
        return bdf_coefficients(1, dt);
      case (Parameters::SimulationControl::TimeSteppingMethod::bdf2):
        return bdf_coefficients(2, dt);
      case (Parameters::SimulationControl::TimeSteppingMethod::bdf3):
        return bdf_coefficients(3, dt);
      default:
        throw(std::runtime_error(
          "BDF coefficients were requested without a BDF method"));
        break;
    }
}

Vector<double>
delta(const unsigned int   p,
      const unsigned int   n,
      const unsigned int   j,
      const Vector<double> times)
{
  if (j == 0)
    {
      Vector<double> arr(p + 1);
      arr    = 0.;
      arr[n] = 1;
      return arr;
    }
  else
    {
      Vector<double> delta_1 = delta(p, n, j - 1, times);
      Vector<double> delta_2 = delta(p, n + 1, j - 1, times);
      Vector<double> delta_sol(p + 1);
      for (unsigned int i_del = 0; i_del < p + 1; ++i_del)
        delta_sol[i_del] =
          (delta_1[i_del] - delta_2[i_del]) / (times[n] - times[n + j]);
      return delta_sol;
    }
}
