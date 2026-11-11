// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests QMOMFieldFunction without a particle size distribution, so
 * that QMOMKit is not required. It checks that:
 * - the components are the scalars followed by the moments, scalar i being
 *   given by expression i, and the time is passed to the expressions;
 * - value() and vector_value() agree;
 * - without scalars, all the components are moments;
 * - a number of scalar expressions that differs from the number of scalars
 *   throws.
 */

#include <core/parameters_multiphysics.h>

#include <solvers/qmom_field_function.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <../tests/tests.h>

#include <exception>
#include <functional>
#include <string>

/**
 * @brief Parse the QMOM functions.
 *
 * @param[in] scalar_expressions Expressions of the scalars.
 *
 * @return Parsed functions, without particle size distribution.
 */
Parameters::QMOMFunctions<2>
parse_functions(const std::string &scalar_expressions)
{
  ParameterHandler prm;
  prm.enter_subsection("QMOM");
  Parameters::QMOMFunctions<2>::declare_parameters(prm);
  prm.leave_subsection();

  prm.parse_input_from_string("subsection QMOM\n"
                              "  subsection scalars\n"
                              "    set Function expression = " +
                              scalar_expressions +
                              "\n"
                              "  end\n"
                              "end\n");

  Parameters::QMOMFunctions<2> functions;
  prm.enter_subsection("QMOM");
  functions.parse_parameters(prm);
  prm.leave_subsection();
  return functions;
}

/**
 * @brief Run a construction of a QMOMFieldFunction and report whether it
 * threw.
 *
 * @param[in] label Description of the scenario printed in the log.
 * @param[in] construct Construction to run.
 */
void
check_throws(const std::string &label, const std::function<void()> &construct)
{
  try
    {
      construct();
      deallog << label << ": did not throw" << std::endl;
    }
  catch (const std::exception &)
    {
      deallog << label << ": threw" << std::endl;
    }
}

void
test()
{
  const Point<2> point(0.25, 0.);

  // 2 scalars and 3 moments: components 0 and 1 are the scalars, components
  // 2 to 4 are the moments, which are zero without a distribution
  {
    const unsigned int   n_scalars = 2;
    QMOMFieldFunction<2> function(parse_functions("1 + x; 2 + t"),
                                  n_scalars,
                                  3,
                                  "test");
    function.set_time(0.5);
    deallog << "2 scalars, 3 moments, t = 0.5, x = 0.25" << std::endl;

    Vector<double> values(function.n_components);
    function.vector_value(point, values);

    bool vector_value_matches_value = true;
    for (unsigned int c = 0; c < function.n_components; ++c)
      {
        if (c < n_scalars)
          deallog << "  scalar " << c << ": ";
        else
          deallog << "  moment " << c - n_scalars << ": ";
        deallog << function.value(point, c) << std::endl;
        vector_value_matches_value &= (values[c] == function.value(point, c));
      }
    deallog << "  vector_value matches value: "
            << (vector_value_matches_value ? "true" : "false") << std::endl;
  }

  // Without scalars, the scalar expression is ignored and all the components
  // are moments
  {
    const QMOMFieldFunction<2> function(parse_functions("0"), 0, 4, "test");
    Vector<double>             values(function.n_components);
    function.vector_value(point, values);
    deallog << "No scalars, 4 moments: " << function.n_components
            << " components, all zero: "
            << (values.l2_norm() == 0. ? "true" : "false") << std::endl;
  }

  check_throws("One expression for 3 scalars", [] {
    const QMOMFieldFunction<2> function(parse_functions("1"), 3, 2, "test");
  });
  check_throws("Two expressions for 3 scalars", [] {
    const QMOMFieldFunction<2> function(parse_functions("1; 2"), 3, 2, "test");
  });
}

int
main()
{
  try
    {
      initlog();
      test();
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
