// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_sdirk_stage_data_h
#define lethe_sdirk_stage_data_h

#include <core/parameters.h>

#include <deal.II/lac/full_matrix.h>

#include <vector>

using namespace dealii;


/**
 * @brief Structure representing a Butcher tableau for a SDIRK method.
 *
 * A singly diagonally implicit Runge-Kutta (SDIRK) method is defined by
 * a Butcher tableau composed of three components:
 *
 * - A: a lower triangular matrix of coefficients a_ij, used to compute each
 * stage.
 * - b: a vector of weights b_i, used to combine the stage values into the final
 * solution.
 * - c: a vector of stage times c_i, indicating at which time fraction (within a
 * time step) each stage is evaluated.
 *
 * These coefficients fully define the behavior of an implicit Runge-Kutta time
 * integration method and are specific to each SDIRK variant (e.g., SDIRK22,
 * SDIRK33).
 */
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

// When referencing to SDIRK methods the pattern used is sdirkOrderStage. For
// instance, sdirk22 indicates SDIRK method with order 2 and 2 stages
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
   * These coefficients are used to compute each stage value \( \boldsymbol{z}_i
   * \) of the Runge-Kutta method according to the formula:
   *
   * \[
   * \boldsymbol{z}_i = \mathcal{F}\left(t_n + c_i \Delta t, \boldsymbol{u}_n +
   * \Delta t \sum_{j=1}^{i} a_{ij} \boldsymbol{z}_j \right)
   * \]
   *
   * Then, the final solution at the next time step is obtained as:
   *
   * \[
   * \boldsymbol{u}_{n+1} = \boldsymbol{u}_n + \Delta t \sum_{i=1}^{s} b_i
   * \boldsymbol{z}_i
   * \]
   *
   * where \( s \) is the total number of stages.
   *
   * Pay attention that the stage indices start from 1, but the indices
   * for the matrix and the vectors start at 0.
   *
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
