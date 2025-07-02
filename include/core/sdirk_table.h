#ifndef SDIRK_TABLE_H
#define SDIRK_TABLE_H

#include <deal.II/lac/full_matrix.h>
#include <vector>
#include <string>
#include <stdexcept>

using namespace dealii;
 
struct SDIRKTable
{
  FullMatrix<double> A;
  std::vector<double> b;
  std::vector<double> c;
};


/**
 * @brief Construct the Butcher tableau (A, b, c) for a given SDIRK method.
 *
 * This function returns the Runge-Kutta coefficients for a diagonally implicit
 * Runge-Kutta (SDIRK) time integration method, such as SDIRK2 (Alexander's method here)
 * or SDIRK3 (Kennedy-Carpenter method here). These coefficients define the time-stepping
 * scheme and are used in implicit time integration of ODEs or PDEs.
 *
 * The Butcher tableau is defined as follows:
 * - A: the matrix of coefficients a_ij
 * - b: the weights used in the final combination of stages
 * - c: the time locations of each stage
 *
 * @param[in] method_name A string identifying the SDIRK method (e.g. "SDIRK2", "SDIRK3").
 *
 * @return SDIRKTable A structure containing the Butcher tableau (A, b, c) for the requested method.
 *
 * @throws std::invalid_argument if the method name is not recognized.
 */

inline SDIRKTable
sdirk_table(const std::string &method_name)
{
  SDIRKTable table;

  if (method_name == "SDIRK2")
    {
      // SDIRK2 - Alexander's method (L-stable)
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

  else if (method_name == "SDIRK3")
    {
      // SDIRK3 - 3-stage, 3rd-order L-stable method (Kennedy-Carpenter)
      const double gamma = 0.4358665215;

      table.A.reinit(3, 3);
      table.A(0, 0) = gamma;
      table.A(1, 0) = (1.0 - gamma) / 2.0;
      table.A(1, 1) = gamma;
      table.A(2, 0) = 2.0 * (1.0 - gamma);
      table.A(2, 1) = - (1.0 - 2.0 * gamma);
      table.A(2, 2) = gamma;

      table.c = {gamma, 0.5, 1.0};
      table.b = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};

      return table;
    }

  else
    {
      throw std::invalid_argument("Unknown SDIRK method: " + method_name);
    }
}

#endif // SDIRK_TABLE_H
