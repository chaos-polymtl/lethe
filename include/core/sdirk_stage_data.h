<<<<<<< HEAD
// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
=======
// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2025 The Lethe Authors
>>>>>>> 996dad99 (version 2 après révision de Bruno)
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

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> a5c7f729 (version 3 après révision de Bruno)

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
<<<<<<< HEAD
=======
>>>>>>> 996dad99 (version 2 après révision de Bruno)
=======
>>>>>>> a5c7f729 (version 3 après révision de Bruno)
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

<<<<<<< HEAD
<<<<<<< HEAD
// Important note : the nomenclature used for the name of the SDIRK methods
// are sdirkOrderStage sdirk22 means SDIRK with order 2 and 2 stages, sdirk33
// means SDIRK with order 3 and 3 stages
=======
>>>>>>> 996dad99 (version 2 après révision de Bruno)
=======
// Important note : the nomenclature used for the name of the SDIRK methods
// are sdirkOrderStage sdirk22 means SDIRK with order 2 and 2 stages, sdirk33
// means SDIRK with order 3 and 3 stages
>>>>>>> a5c7f729 (version 3 après révision de Bruno)
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
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> a5c7f729 (version 3 après révision de Bruno)
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
   * Pay attention that the stage indices start from 1, but the but indicices
   * for the matrix and the vectors when coding start from 0.
   *
   *
<<<<<<< HEAD
=======
>>>>>>> 996dad99 (version 2 après révision de Bruno)
=======
>>>>>>> a5c7f729 (version 3 après révision de Bruno)
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
