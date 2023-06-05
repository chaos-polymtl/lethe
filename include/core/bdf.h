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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef lethe_bdf_h
#define lethe_bdf_h

#include <core/parameters.h>

#include <deal.II/lac/vector.h>

#include <vector>

using namespace dealii;

/**
 * @brief Calculate the coefficients required for BDF integration of order n.
 * The coefficients are determined through a recursion algorithm.
 * The algorithm is taken from : Hay, Alexander, et al. "hp-Adaptive time
 * integration based on the BDF for viscous flows." Journal of Computational
 * Physics 291 (2015): 151-176.
 *
 * @param order The order of the BDF method. The BDF method of order n requires n+1 arrays
 *
 * @param time_steps a vector containing all the time steps. The time steps should be in reverse order.
 * For example, if the method is a BDF2, it uses three values for the time : (t)
 * (t-dt_1) and (t-dt_2). Thus the time step vector should contain dt_1 and
 * dt_2.
 */
Vector<double>
bdf_coefficients(unsigned int order, const std::vector<double> time_steps);


/**
 * @brief Calculate the coefficients required for BDF integration of order n.
 * The coefficients are determined through a recursion algorithm.
 * The algorithm is taken from : Hay, Alexander, et al. "hp-Adaptive time
 * integration based on the BDF for viscous flows." Journal of Computational
 * Physics 291 (2015): 151-176.
 *
 * @param method The time stepping method
 *
 * @param time_steps a vector containing all the time steps. The time steps should be in reverse order.
 * For example, if the method is a BDF2, it uses three values for the time : (t)
 * (t-dt_1) and (t-dt_2). Thus the time step vector should contain dt_1 and
 * dt_2.
 */
Vector<double>
bdf_coefficients(const Parameters::SimulationControl::TimeSteppingMethod method,
                 const std::vector<double> time_steps);


/**
 * @brief Recursion function to calculate the bdf coefficient
 *
 * @param order Order of the bdf method
 *
 * @param n Integer required for the recursive function calculation
 *
 * @param j Integer required for the recursive function calculation
 *
 * @param times Vector of the times corresponding to the various iteration in reverse order.
 * For example in a BDF2 method, three values of the time must be in the array :
 * t, t-dt_1 and t-dt_2.
 */
Vector<double>
delta(const unsigned int   order,
      const unsigned int   n,
      const unsigned int   j,
      const Vector<double> times);



/**
 * @brief Returns the maximum number of previous time steps supposed by the BDF schemes implemented in Lethe.
 *  At the moment this is hardcoded to 3, but eventually this could be made
 * larger or smaller depending on the methods used.
 *
 */
inline unsigned int
maximum_number_of_previous_solutions()
{
  return 3;
}

/**
 * @brief Returns the number of previous time steps for a BDF Scheme
 *
 */
inline unsigned int
number_of_previous_solutions(
  Parameters::SimulationControl::TimeSteppingMethod method)
{
  if (method == Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
      method == Parameters::SimulationControl::TimeSteppingMethod::steady_bdf)
    return 1;
  else if (method == Parameters::SimulationControl::TimeSteppingMethod::bdf2)
    return 2;
  else if (method == Parameters::SimulationControl::TimeSteppingMethod::bdf3)
    return 3;
  else
    return 0;
}

/**
 * @brief Extrapolates vector of solution to time_vector[0] using previous solution times and previous solution
 *
 * @tparam The type of the variable being extrapolated (e.g. double, Tensor<1,dim>)
 * @param time_vector The vector of times. The solution will be extrapolated to time 0. It should be of size number_of_previous_solutions+1
 * @param solution_vector The vector of solutions. The solution will be extrapolated to index 0. It should be at least of size number_of_previous_solutions
 * @param number_of_previous_solution The number of previous solutions to use to extrapolate
 * @param extrapolated_solution The vector of extrapolated solutions
 */

template <typename DataType>
inline void
bdf_extrapolate(const std::vector<double> &               time_vector,
                const std::vector<std::vector<DataType>> &solution_vector,
                const unsigned int     number_of_previous_solution,
                std::vector<DataType> &extrapolated_solution)
{
  // If only one previous solution is used, we don't make a Lagrange
  // polynomial but use the previous value
  if (number_of_previous_solution == 1)
    {
      for (unsigned int q = 0; q < extrapolated_solution.size(); ++q)
        extrapolated_solution[q] = solution_vector[0][q];
      return;
    }

  // Otherwise we extrapolate with a Lagrange polynomial
  for (unsigned int q = 0; q < extrapolated_solution.size(); ++q)
    {
      // Set extrapolated solution to zero
      extrapolated_solution[q] = 0;

      // For all the previous solutions which we will use to extrapolate
      for (unsigned int p = 0; p < number_of_previous_solution; ++p)
        {
          // Factor is the weight of the previous solution fixed by the Lagrange
          // Polynomial
          double factor = 1;
          for (unsigned int k = 0; k < number_of_previous_solution; ++k)
            {
              if (p != k)
                {
                  // The time vector also contains the time of the solution to
                  // be extrapoled to, hence the +1
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
