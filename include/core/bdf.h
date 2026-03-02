// SPDX-FileCopyrightText: Copyright (c) 2019-2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_bdf_h
#define lethe_bdf_h

#include <core/parameters.h>

#include <deal.II/lac/vector.h>

#include <vector>

using namespace dealii;

/**
 * @brief Calculate the coefficients required for backward differentiation
 * formula (BDF) integration of @p order \f$n\f$.
 *
 * The coefficients are determined through a recursion algorithm.
 * The algorithm is taken from: Hay, Alexander, et al. "hp-Adaptive time
 * integration based on the BDF for viscous flows." Journal of Computational
 * Physics 291 (2015): 151-176.
 *
 * @param[in] method Time stepping method.
 *
 * @param[in] time_steps Vector containing all the time steps. The time steps
 * should be in decreasing order.
 *
 * For example, if the method is a BDF2 (\f$n=2\f$), it uses three values for
 * the time: \f$t\f$, \f$t-\Delta t_1\f$ and \f$t-\Delta t_2\f$. Thus, the @p
 * time_steps vector should contain \f$\Delta t_1\f$ and \f$\Delta t_2\f$.
 *
 * @return Vector containing BDF integration coefficient values.
 *
 * @note At the moment, the highest implemented scheme is @p bdf3 (@p order=3).
 */
Vector<double>
calculate_bdf_coefficients(
  Parameters::SimulationControl::TimeSteppingMethod method,
  const std::vector<double>                        &time_steps);


/**
 * @brief Recursive function for calculating the \f$j^{th}\f$ order divided
 * differences (\f$\delta^j\f$) that is used for calculating backward
 * differentiation formula (BDF) coefficients.
 *
 * The algorithm is taken from: Hay, Alexander, et al. "hp-Adaptive time
 * integration based on the BDF for viscous flows." Journal of Computational
 * Physics 291 (2015): 151-176.
 *
 * @param[in] order Order of the BDF method.
 *
 * @param[in] n Integer required for the recursive function calculation.
 *
 * @param[in] j Order of the divided difference of the value of interest.
 *
 * @param[in] times Vector of times values corresponding to the various
 * iteration in decreasing order.
 *
 * For example, if the method is a BDF2 (\f$n=2\f$), three values of the time
 * must be in the array: \f$t\f$, \f$t-\Delta t_1\f$ and \f$t-\Delta t_2\f$.
 *
 * @return Vector of \f$n\f$ divided differences (\f$\delta^j\f$) used in the
 * calculation of BDF coefficients.
 */
Vector<double>
delta(unsigned int          order,
      unsigned int          n,
      unsigned int          j,
      const Vector<double> &times);

/**
 * @brief Get the maximum number of previous time steps supposed by the BDF
 * schemes implemented in Lethe.
 *
 * @return Value of the maximum number of previous time steps used in BDF
 * schemes implemented in Lethe.
 *
 * @note At the moment, this is hardcoded to 3, but eventually this could be
 * made larger or smaller depending on the methods used.
 */
inline unsigned int
maximum_number_of_previous_solutions()
{
  return 3;
}

/**
 * @brief Get the number of previous time steps for a specific BDF scheme.
 *
 * @param[in] method Time stepping method of the simulation.
 *
 * @return Number of previous time steps used in a given BDF scheme.
 */
inline unsigned int
number_of_previous_solutions(
  const Parameters::SimulationControl::TimeSteppingMethod method)
{
  if (method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    return 1;
  if (method == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    return 2;
  if (method == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    return 3;

  return 0;
}

/**
 * @brief Extrapolate vector of solution to time_vector[0] using previous
 * solution times and previous solutions.
 *
 * @tparam DataType Type of the variable being extrapolated (e.g. double,
 * Tensor<1,dim>)
 *
 * @param[in] time_vector Vector of times. The solution will be extrapolated to
 * time 0. It should be of size @p number_of_previous_solutions+1.
 *
 * @param[in] solution_vector Vector of solution. The solution will be
 * extrapolated to index 0. It should be at least of size @p
 * number_of_previous_solutions.
 *
 * @param[in] number_of_previous_solutions Number of previous solutions used
 * to extrapolate.
 *
 * @param[out] extrapolated_solution Vector of extrapolated solution.
 */

template <typename DataType>
void
bdf_extrapolate(const std::vector<double>                &time_vector,
                const std::vector<std::vector<DataType>> &solution_vector,
                const unsigned int     number_of_previous_solutions,
                std::vector<DataType> &extrapolated_solution)
{
  // If only one previous solution is used, we don't make a Lagrange
  // polynomial but use the previous value
  if (number_of_previous_solutions == 1)
    {
      for (unsigned int q = 0; q < extrapolated_solution.size(); ++q)
        extrapolated_solution[q] = solution_vector[0][q];
      return;
    }

  // Otherwise, we extrapolate with a Lagrange polynomial
  for (unsigned int q = 0; q < extrapolated_solution.size(); ++q)
    {
      // Set extrapolated solution to zero
      extrapolated_solution[q] = 0;

      // For all the previous solutions which we will use to extrapolate
      for (unsigned int p = 0; p < number_of_previous_solutions; ++p)
        {
          // Factor is the weight of the previous solution fixed by the Lagrange
          // polynomial
          double factor = 1;
          for (unsigned int k = 0; k < number_of_previous_solutions; ++k)
            {
              if (p != k)
                {
                  // The time vector also contains the time of the solution to
                  // be extrapolated to, hence the +1
                  factor *= (time_vector[0] - time_vector[k + 1]) /
                            (time_vector[p + 1] - time_vector[k + 1]);
                }
            }
          // Calculate the contribution of this time step to the extrapolated
          // solution
          extrapolated_solution[q] += solution_vector[p][q] * factor;
        }
    }
}


#endif
