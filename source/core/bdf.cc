// SPDX-FileCopyrightText: Copyright (c) 2019-2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>

#include <deal.II/lac/vector.h>


Vector<double>
bdf_coefficients(const unsigned int p, const std::vector<double> &time_steps)
{
  // There should be at least p time steps
  assert(time_steps.size() >= p);

  // Create a timetable for the bdf formula
  Vector<double> times(p + 1);
  for (unsigned int i = 0; i < p + 1; ++i)
    {
      times[i] = 0.;
      for (unsigned int j = 0; j < i; ++j)
        times[i] -= time_steps[j];
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
calculate_bdf_coefficients(
  const Parameters::SimulationControl::TimeSteppingMethod method,
  const std::vector<double>                              &time_steps)
{
  switch (method)
    {
      case Parameters::SimulationControl::TimeSteppingMethod::bdf1:
        return bdf_coefficients(1, time_steps);
      case Parameters::SimulationControl::TimeSteppingMethod::steady_bdf:
        return bdf_coefficients(1, time_steps);
      case Parameters::SimulationControl::TimeSteppingMethod::bdf2:
        return bdf_coefficients(2, time_steps);
      case Parameters::SimulationControl::TimeSteppingMethod::bdf3:
        return bdf_coefficients(3, time_steps);
      default:
        throw(std::runtime_error(
          "BDF coefficients were requested without a BDF method"));
    }
}

Vector<double>
delta(const unsigned int    order,
      const unsigned int    n,
      const unsigned int    j,
      const Vector<double> &times)
{
  if (j == 0)
    {
      Vector<double> arr(order + 1);
      arr    = 0.;
      arr[n] = 1;
      return arr;
    }

  // else
  Vector<double> delta_1 = delta(order, n, j - 1, times);
  Vector<double> delta_2 = delta(order, n + 1, j - 1, times);
  Vector<double> delta_sol(order + 1);
  for (unsigned int i_del = 0; i_del < order + 1; ++i_del)
    delta_sol[i_del] =
      (delta_1[i_del] - delta_2[i_del]) / (times[n] - times[n + j]);
  return delta_sol;
}
