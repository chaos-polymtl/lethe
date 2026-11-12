// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the moments of the QMOM field functions computed by
 * QMOMKit from a particle size distribution. The moments are M_k = n E[L^k],
 * where n is the number density. With n = 1 + x evaluated at x = 1:
 * - lognormal (0, 0.5): M_k = 2 exp(k^2 / 8);
 * - gamma (shape 2, scale 3): M_k = 2 3^k (k + 1)!;
 * - weibull (shape 2, scale 1): M_k = 2 Gamma(1 + k / 2);
 * - a zero number density gives zero moments whatever the distribution.
 */

#include <core/parameters_multiphysics.h>

#include <solvers/qmom_field_function.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <../tests/tests.h>

#include <iostream>
#include <string>

using namespace dealii;

/**
 * @brief Parse the QMOM functions with one scalar and a particle size
 * distribution.
 *
 * @param[in] type Type of the particle size distribution.
 * @param[in] expressions Number density followed by the parameters of the
 * distribution.
 *
 * @return Parsed functions.
 */
template <int dim>
Parameters::QMOMFunctions<dim>
parse_functions(const std::string &type, const std::string &expressions)
{
  ParameterHandler prm;
  prm.enter_subsection("QMOM");
  Parameters::QMOMFunctions<dim>::declare_parameters(prm);
  prm.leave_subsection();

  prm.parse_input_from_string("subsection QMOM\n"
                              "  subsection distribution\n"
                              "    set type                = " +
                              type +
                              "\n"
                              "    set Function expression = " +
                              expressions +
                              "\n"
                              "  end\n"
                              "end\n");

  Parameters::QMOMFunctions<dim> functions;
  prm.enter_subsection("QMOM");
  functions.parse_parameters(prm);
  prm.leave_subsection();
  return functions;
}

/**
 * @brief Print the moments of a QMOM field function with one scalar and 4
 * moments.
 *
 * @param[in] label Description of the distribution printed in the log.
 * @param[in] type Type of the particle size distribution.
 * @param[in] expressions Number density followed by the parameters of the
 * distribution.
 */
template <int dim>
void
print_moments(const std::string &label,
              const std::string &type,
              const std::string &expressions)
{
  constexpr unsigned int n_scalars = 1;
  constexpr unsigned int n_moments = 4;

  Point<dim> point;
  point[0] = 1.;

  const QMOMFieldFunction<dim> function(parse_functions<dim>(type, expressions),
                                        n_scalars,
                                        n_moments,
                                        "test");
  Vector<double>               values(function.n_components);
  function.vector_value(point, values);

  deallog << label << ":";
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    deallog << " " << values[n_scalars + i_moment];
  deallog << std::endl;
}

template <int dim>
void
test()
{
  print_moments<dim>("Lognormal", "lognormal", "1 + x; 0; 0.5");
  print_moments<dim>("Gamma", "gamma", "1 + x; 2; 3");
  print_moments<dim>("Weibull", "weibull", "1 + x; 2; 1");
  print_moments<dim>("Zero number density", "lognormal", "0; 0; 0.5");
}

int
main()
{
  try
    {
      initlog();
      test<2>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
