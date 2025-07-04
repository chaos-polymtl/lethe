// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef SDIRK_STAGE_DATA_H
#define SDIRK_STAGE_DATA_H

#include <core/parameters.h>
#include <core/simulation_control.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/full_matrix.h>

#include <string>
#include <vector>

using namespace dealii;

struct SDIRKTable
{
  FullMatrix<double>  A;
  std::vector<double> b;
  std::vector<double> c;
};

/**
 * @brief Construct the Butcher tableau (A, b, c) for a given SDIRK method.
 *
 * This function returns the Runge-Kutta coefficients for a diagonally implicit
 * Runge-Kutta (SDIRK) time integration method, such as SDIRK22 (Alexander's
 * method here) or SDIRK33 (Kennedy-Carpenter method here). These coefficients
 * define the time-stepping scheme and are used in implicit time integration of
 * ODEs or PDEs.
 *
 * The Butcher tableau is defined as follows:
 * - A: the matrix of coefficients a_ij
 * - b: the weights used in the final combination of stages
 * - c: the time locations of each stage
 *
 * @param[in] method The time integration method (e.g.,
 * TimeSteppingMethod::sdirk22, TimeSteppingMethod::sdirk33) selected from the
 * Parameters::SimulationControl::TimeSteppingMethod enum.
 *
 * @return SDIRKTable A structure containing the Butcher tableau (A, b, c) for the requested method.
 */

SDIRKTable
sdirk_table(const Parameters::SimulationControl::TimeSteppingMethod method);



class SDIRKStageData
{
public:
  /**
   * @brief Extract the coefficients of a given stage from a SDIRK Butcher tableau.
   *
   * This function returns the data needed to compute the solution increment
   * at a particular stage of a SDIRK (Singly Diagonally Implicit Runge-Kutta)
   * method. The returned structure contains:
   * - a_ij: the vector of coefficients from stage j to stage i (for j = 0 to
   * i),
   * - c_i: the time coefficient associated with stage i,
   * - b_i: the weight for combining the result of stage i in the final
   * solution.
   *
   * @param[in] table The SDIRK Butcher tableau (structure containing A, b, and
   * c).
   * @param[in] stage_i Index of the stage (starting from 1) for which the data
   * is extracted.
   * @throws dealii::ExcMessage if the tableau is inconsistent (e.g., wrong
   * dimensions).
   */
  SDIRKStageData(const SDIRKTable &table, const unsigned int stage_i);

  std::vector<double> a_ij;
  double              c_i;
  double              b_i;
};

#endif // SDIRK_STAGE_DATA_H
