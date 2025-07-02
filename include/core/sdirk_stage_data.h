// get_sdirk_stage_data.h
#ifndef SDIRK_STAGE_DATA_H
#define SDIRK_STAGE_DATA_H

#include <deal.II/lac/full_matrix.h>
#include <vector>
#include <stdexcept>

using namespace dealii;

struct SDIRKStageData
{
  std::vector<double> a_ij;
  double              c_i;
  double              b_i;
};

/**
 * @brief Extract the coefficients of a given stage from a SDIRK Butcher tableau.
 *
 * This function returns the data needed to compute the solution increment
 * at a particular stage of a SDIRK (Singly Diagonally Implicit Runge-Kutta)
 * method. The returned structure contains:
 * - a_ij: the vector of coefficients from stage j to stage i (for j = 0 to i),
 * - c_i: the time coefficient associated with stage i,
 * - b_i: the weight for combining the result of stage i in the final solution.
 *
 * @param[in] butcher_table Full matrix representing the A coefficients (i.e., the lower triangle of the Butcher tableau).
 * @param[in] c Vector of c_i coefficients (stage times).
 * @param[in] b Vector of b_i coefficients (final weights).
 * @param[in] stage_i Index of the stage for which the data is to be extracted.
 *
 * @return SDIRKStageData A structure containing the relevant coefficients for stage i.
 *
 * @throws std::invalid_argument if the stage index is out of bounds.
 */

SDIRKStageData
sdirk_stage_data(const FullMatrix<double>  &butcher_table,
                     const std::vector<double> &c,
                     const std::vector<double> &b,
                     unsigned int               stage_i);

#endif // SDIRK_STAGE_DATA_H

