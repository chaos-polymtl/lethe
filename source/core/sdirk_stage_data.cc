// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters.h>
#include <core/sdirk_stage_data.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/full_matrix.h>

#include <numbers>
#include <vector>

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
      // SDIRK33 - 3-stage, 3rd-order, L-stable and stiffly-accurate method
      // (Alexander's method). References:
      // R. Alexander, Diagonally implicit Runge–Kutta methods for stiff ODEs,
      // SIAM Journal on Numerical Analysis, 14(6), 1006–1021, 1977.
      // DOI: https://doi.org/10.1137/0714068
      // C. A. Kennedy and M. H. Carpenter, "Diagonally Implicit Runge-Kutta
      // Methods for Ordinary Differential Equations. A Review"
      // NASA/TM–2016–219173, March 2016. No DOI available URL :
      // https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf

      // gamma is the root of gamma^3 - 3 gamma^2 + 3/2 gamma - 1/6 = 0 that
      // makes the method L-stable. All the other coefficients are expressed
      // in closed form from gamma so that the order conditions are satisfied
      // to machine precision.
      constexpr double gamma = 0.43586652150845899941601945;
      constexpr double c2    = (1.0 + gamma) / 2.0;
      constexpr double b1 = -(6.0 * gamma * gamma - 16.0 * gamma + 1.0) / 4.0;
      constexpr double b2 = (6.0 * gamma * gamma - 20.0 * gamma + 5.0) / 4.0;

      table.A.reinit(3, 3);
      table.A(0, 0) = gamma;
      table.A(1, 0) = c2 - gamma;
      table.A(1, 1) = gamma;
      table.A(2, 0) = b1;
      table.A(2, 1) = b2;
      table.A(2, 2) = gamma;

      table.c.resize(3);
      table.c[0] = gamma;
      table.c[1] = c2;
      table.c[2] = 1.0;

      table.b.resize(3);
      table.b[0] = b1;
      table.b[1] = b2;
      table.b[2] = gamma;

      return table;
    }

  else if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk45)
    {
      // SDIRK45 - 5-stage, 4th-order, L-stable and stiffly-accurate method
      // (SDIRK4 of Hairer and Wanner).
      // Reference: E. Hairer and G. Wanner, "Solving Ordinary Differential
      // Equations II: Stiff and Differential-Algebraic Problems", 2nd edition,
      // Springer, 1996, Section IV.6, Table 6.5.
      // DOI: https://doi.org/10.1007/978-3-642-05221-7
      //
      // The method also has an embedded 3rd-order solution, with weights
      // b_hat = (59/48, -17/96, 225/32, -85/12, 0), which could be used to
      // estimate the time integration error.
      //
      // Some coefficients of the method are large compared to gamma (up to
      // |a_ij / gamma| = 31.25). The error of the non-linear solver at a stage
      // is thus amplified in the following stages, so a tight non-linear
      // tolerance should be used with this method.
      constexpr double gamma = 1.0 / 4.0;

      table.A.reinit(5, 5);
      table.A(0, 0) = gamma;

      table.A(1, 0) = 1.0 / 2.0;
      table.A(1, 1) = gamma;

      table.A(2, 0) = 17.0 / 50.0;
      table.A(2, 1) = -1.0 / 25.0;
      table.A(2, 2) = gamma;

      table.A(3, 0) = 371.0 / 1360.0;
      table.A(3, 1) = -137.0 / 2720.0;
      table.A(3, 2) = 15.0 / 544.0;
      table.A(3, 3) = gamma;

      table.A(4, 0) = 25.0 / 24.0;
      table.A(4, 1) = -49.0 / 48.0;
      table.A(4, 2) = 125.0 / 16.0;
      table.A(4, 3) = -85.0 / 12.0;
      table.A(4, 4) = gamma;

      table.c.resize(5);
      table.c[0] = 1.0 / 4.0;
      table.c[1] = 3.0 / 4.0;
      table.c[2] = 11.0 / 20.0;
      table.c[3] = 1.0 / 2.0;
      table.c[4] = 1.0;

      // The method is stiffly accurate: the weights b are the last row of A
      table.b.resize(5);
      for (unsigned int j = 0; j < 5; ++j)
        table.b[j] = table.A(4, j);

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
