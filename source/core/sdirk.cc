#include <deal.II/lac/full_matrix.h>

#include <core/sdirk.h>

// Matrix of coefficients for the SDIRK methods
// The lines store the information required for each step
// Column 0 always refer to outcome of the step that is being calculated
// Column 1 always refer to step n
// Column 2+ refer to intermediary steps
FullMatrix<double>
sdirk_coefficients(unsigned int order, double time_step)
{
  FullMatrix<double> sdirk_coefs;
  double             sdt = 1. / time_step;


  if (order == 2)
    {
      const double alpha = (2. - std::sqrt(2)) / 2.;
      sdirk_coefs.reinit(2, 3);
      sdirk_coefs[0][0] = 1. / alpha * sdt;
      sdirk_coefs[0][1] = -1. / alpha * sdt;
      sdirk_coefs[1][0] = 1. / alpha * sdt;
      sdirk_coefs[1][1] = -(2 * alpha - 1) / alpha / alpha * sdt;
      sdirk_coefs[1][2] = -(1 - alpha) / alpha / alpha * sdt;
    }

  if (order == 3)
    {
      sdirk_coefs.reinit(3, 4);
      sdirk_coefs[0][0] = 2.29428036032357 * sdt;
      sdirk_coefs[0][1] = -2.29428036032357 * sdt;
      sdirk_coefs[1][0] = 2.29428036032357 * sdt;
      sdirk_coefs[1][1] = -0.809559354865318 * sdt;
      sdirk_coefs[1][2] = -1.48472100545825 * sdt;
      sdirk_coefs[2][0] = 2.29428036032357 * sdt;
      sdirk_coefs[2][1] = 2.87009860266089 * sdt;
      sdirk_coefs[2][2] = -8.55612780171952 * sdt;
      sdirk_coefs[2][3] = 3.39174883873506 * sdt;
    }

  return sdirk_coefs;
}
