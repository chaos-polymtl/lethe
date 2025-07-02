<<<<<<< HEAD
// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters.h>
#include <core/sdirk_stage_data.h>
#include <core/simulation_control.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>

#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

using namespace dealii;

/// Function to construct the Butcher's table for a given SDIRK method
SDIRKTable
sdirk_table(const Parameters::SimulationControl::TimeSteppingMethod method)
{
  SDIRKTable table;

  // Important note : the nomenclature used for the name of the SDIRK methods
  // are sdirkOrderStage sdirk22 means SDIRK with order 2 and 2 stages, sdirk33
  // means SDIRK with order 3 and 3 stages
  if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk22)
    {
      // SDIRK22 - Alexander's method (L-stable)
      // Reference: R. Alexander, Diagonally implicit Runge–Kutta methods for
      // stiff ODEs, SIAM Journal on Numerical Analysis, 14(6), 1006–1021, 1977.
      // DOI: https://doi.org/10.1137/0714068
      const double alpha = (2. - std::sqrt(2.)) / 2.;

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

  // Important note : the nomenclature used for the name of the SDIRK methods
  // are sdirkOrderStage sdirk22 means SDIRK with order 2 and 2 stages, sdirk33
  // means SDIRK with order 3 and 3 stages
  else if (method == Parameters::SimulationControl::TimeSteppingMethod::sdirk33)
    {
      // SDIRK33 - 3-stage, 3rd-order L-stable method (Kennedy-Carpenter)
      // Reference: J. F. G. Kennedy and A. D. Carpenter, "Diagonally Implicit
      // Runge-KuttaMethods for Ordinary Differential Equation"
      // NASA/TM–2016–219173 pp. 77, March 2016. No DOI available URL :
      // https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf
      const double gamma = 0.4358665215;

      table.A.reinit(3, 3);
      table.A(0, 0) = gamma;
      table.A(1, 0) = 0.2820667393;
      table.A(1, 1) = gamma;
      table.A(2, 0) = 1 - 1.772630128 - gamma;
      table.A(2, 1) = 1.772630128;
      table.A(2, 2) = gamma;

      table.c = {gamma, 0.7179332608, 1.0};
      table.b = {1 - 1.772630128 - gamma, 1.772630128, gamma};

      return table;
    }

  else
    {
      AssertThrow(false, ExcMessage("Unknown SDIRK method in sdirk_table()."));
      return table;
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
=======
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/logstream.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <iomanip>

#include <core/sdirk_stage_data.h>

using namespace dealii;

SDIRKStageData
sdirk_stage_data(const FullMatrix<double> &butcher_table,
                     const std::vector<double> &c,
                     const std::vector<double> &b,
                     const unsigned int         stage_i)
{
  const unsigned int n_stages = b.size();

  if (stage_i >= n_stages)
    throw std::invalid_argument("Invalid stage index.");

  SDIRKStageData data;
  data.a_ij.resize(stage_i + 1);

  for (unsigned int j = 0; j <= stage_i; ++j)
    data.a_ij[j] = butcher_table(stage_i, j);

  data.c_i = c[stage_i];
  data.b_i = b[stage_i];

  return data;
}

>>>>>>> 737d91a8 (sdirk_coefficients function)
