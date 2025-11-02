// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters.h>
#include <core/sdirk_stage_data.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/full_matrix.h>

#include <vector>
#include <numbers>

using namespace dealii;

/// Function to construct the Butcher's table for a given SDIRK method
// The coefficients are fixed data during the simulation.
SDIRKTable
sdirk_table(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  SDIRKTable table;

  // When referencing to SDIRK methods the pattern used is sdirkOrderStage. For
  // instance, sdirk22 indicates SDIRK method with order 2 and 2 stages
  if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22)
    {
      // SDIRK22 - Alexander's method (L-stable)
      // Reference: R. Alexander, Diagonally implicit Runge–Kutta methods for
      // stiff ODEs, SIAM Journal on Numerical Analysis, 14(6), 1006–1021, 1977.
      // DOI: https://doi.org/10.1137/0714068

      // alpha is a constant coefficient chosen to ensure the method is stable
      // and consistent
      constexpr double alpha = (2. - std::numbers::sqrt2) / 2.;

      table.A.reinit(2, 2);
      table.A(0, 0) = alpha;
      table.A(1, 0) = 1.0 - alpha;
      table.A(1, 1) = alpha;

      table.c.resize(2);
      table.c[0] = alpha;
      table.c[1] = 1.0;

      table.b.resize(2);
      table.b[0] = 1.0 - alpha;
      table.b[1] = alpha;

      return table;
    }

  else if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33)
    {
      // SDIRK33 - 3-stage, 3rd-order L-stable method (Kennedy-Carpenter)
      // Reference: J. F. G. Kennedy and A. D. Carpenter, "Diagonally Implicit
      // Runge-KuttaMethods for Ordinary Differential Equation"
      // NASA/TM–2016–219173 pp. 77, March 2016. No DOI available URL :
      // https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf

      // gamma is a constant coefficient chosen to ensure the method is stable
      // and consistent
      constexpr double gamma = 0.4358665215;

      table.A.reinit(3, 3);
      table.A(0, 0) = gamma;
      table.A(1, 0) = 0.2820667393;
      table.A(1, 1) = gamma;
      table.A(2, 0) = 1.208496649;
      table.A(2, 1) = -0.6443631707;
      table.A(2, 2) = gamma;

      table.c.resize(3);
      table.c[0] = gamma;
      table.c[1] = 0.7179332608;
      table.c[2] = 1.0;

      table.b.resize(3);
      table.b[0] = 1.208496649;
      table.b[1] = -0.6443631707;
      table.b[2] = gamma;

      return table;
    }

  else if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk43)
    {
      // SDIRK43 - 3-stage, 4th-order method
      // Reference: L. Ferracina, M.N. Spijker, "Strong stability of
      // singly-diagonally-implicit Runge–Kutta methods" Applied Numerical
      // Mathematics 58 (2008) 1675–1686 DOI: doi:10.1016/j.apnum.2007.10.004

      // gamma is a constant coefficient chosen to ensure the method is stable
      // and consistent const double gamma = 0.128886400515;
      constexpr double gamma = 1.068579021301629;

      table.A.reinit(3, 3);
      table.A(0, 0) = gamma;
      table.A(1, 0) = 1.0 / 2.0 - gamma;
      table.A(1, 1) = gamma;
      table.A(2, 0) = 2 * gamma;
      table.A(2, 1) = 1 - 4 * gamma;
      table.A(2, 2) = gamma;

      table.c.resize(3);
      table.c[0] = gamma;
      table.c[1] = 1.0 / 2.0;
      table.c[2] = 1 - gamma;

      table.b.resize(3);
      table.b[0] = 1.0 / (6 * std::pow((2 * gamma - 1), 2.0));
      table.b[1] = 2 * (6 * std::pow(gamma, 2.0) - 6 * gamma + 1) /
                   (3 * std::pow((2 * gamma - 1), 2.0));
      table.b[2] = 1 / (6 * std::pow((2 * gamma - 1), 2.0));

      return table;
    }

  else
    {
      AssertThrow(false, ExcMessage("Unknown SDIRK method in sdirk_table()."));
    }
}


// Function to extract the coefficients of a given stage from a SDIRK Butcher's
// table

SDIRKStageData::SDIRKStageData(const SDIRKTable  &table,
                               const unsigned int stage_i)
{
  AssertThrow(table.c.size() == table.b.size(),
              ExcMessage("b and c vectors must have the same size."));

  AssertThrow(table.A.m() == table.b.size(),
              ExcMessage(
                "Matrix A dimensions must match size of b and c vectors."));

  // The stages take the values 1, ..., nstages

  a_ij.resize(stage_i);

  for (unsigned int j = 0; j < stage_i; ++j)
    // The Butcher's table starts at index 0 whereas the stages take the values
    // 1, ..., nstages
    a_ij[j] = table.A(stage_i - 1, j);

  c_i = table.c[stage_i - 1];
  b_i = table.b[stage_i - 1];
}
