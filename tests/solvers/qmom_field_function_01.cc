// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the functions that prescribe the fields of the QMOM
 * physics (initial conditions, Dirichlet boundary conditions and analytical
 * solution). It checks that:
 * - the default functions are zero and have a single moment expression, which
 *   applies to every moment;
 * - one expression per moment gives each moment its own function;
 * - the time of the functions is propagated to the time-dependent
 *   expressions;
 * - a number of moment expressions that is neither 1 nor the number of
 *   moments throws.
 */

#include <core/parameters_multiphysics.h>

#include <solvers/qmom.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <../tests/tests.h>

#include <exception>
#include <iostream>
#include <string>

using namespace dealii;

/**
 * @brief Parse the QMOM functions from an input string.
 *
 * @param[in] input Content of the "QMOM" subsection.
 *
 * @return Parsed functions.
 */
template <int dim>
Parameters::QMOMFunctions<dim>
parse_functions(const std::string &input)
{
  ParameterHandler prm;
  prm.enter_subsection("QMOM");
  Parameters::QMOMFunctions<dim>::declare_parameters(prm);
  prm.leave_subsection();

  prm.parse_input_from_string("subsection QMOM\n" + input + "end\n");

  Parameters::QMOMFunctions<dim> functions;
  prm.enter_subsection("QMOM");
  functions.parse_parameters(prm);
  prm.leave_subsection();
  return functions;
}

/**
 * @brief Print all the components of a QMOM field function at a point.
 *
 * @param[in] function Function to evaluate.
 * @param[in] point Evaluation point.
 */
template <int dim>
void
print_components(const QMOMFieldFunction<dim> &function,
                 const Point<dim>             &point)
{
  for (unsigned int component = 0; component < function.n_components;
       ++component)
    deallog << "  component " << component << ": "
            << function.value(point, component) << std::endl;
}

template <int dim>
void
test()
{
  Point<dim> point;
  point[0] = 0.25;

  // Default functions: zero, with a single moment expression applied to the
  // 8 moments of a 4-node quadrature
  {
    const Parameters::QMOMFunctions<dim> functions = parse_functions<dim>("");
    deallog << "Default moment expressions: " << functions.moments->n_components
            << std::endl;
    const QMOMFieldFunction<dim> function(functions, 8, "test");
    deallog << "Default function components: " << function.n_components
            << std::endl;
    print_components(function, point);
  }

  // One time-dependent expression per moment for a 2-node quadrature
  {
    const Parameters::QMOMFunctions<dim> functions = parse_functions<dim>(
      "  subsection ethanol concentration\n"
      "    set Function expression = 1 + x\n"
      "  end\n"
      "  subsection free lipids concentration\n"
      "    set Function expression = 2 + t\n"
      "  end\n"
      "  subsection moments\n"
      "    set Function expression = 10; 20 + t; 30 * x; 40\n"
      "  end\n");
    QMOMFieldFunction<dim> function(functions, 4, "test");
    function.set_time(0.5);
    deallog << "Per-moment function at t = " << function.get_time() << ":"
            << std::endl;
    print_components(function, point);
  }

  // A single expression applied to the 4 moments of a 2-node quadrature
  {
    const Parameters::QMOMFunctions<dim> functions =
      parse_functions<dim>("  subsection moments\n"
                           "    set Function expression = 3 * x\n"
                           "  end\n");
    const QMOMFieldFunction<dim> function(functions, 4, "test");
    deallog << "Single moment expression:" << std::endl;
    print_components(function, point);
  }

  // Three moment expressions for 4 moments
  {
    const Parameters::QMOMFunctions<dim> functions =
      parse_functions<dim>("  subsection moments\n"
                           "    set Function expression = 1; 2; 3\n"
                           "  end\n");
    try
      {
        const QMOMFieldFunction<dim> function(functions, 4, "test");
        deallog << "Wrong number of moment expressions did not throw"
                << std::endl;
      }
    catch (const std::exception &)
      {
        deallog << "Wrong number of moment expressions threw" << std::endl;
      }
  }
}

int
main()
{
  try
    {
      initlog();
      test<2>();
      test<3>();
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
