// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This code tests the functions that prescribe the fields of the QMOM
 * physics (initial conditions, Dirichlet boundary conditions, analytical
 * solution and source terms) without a particle size distribution, so that it
 * does not require QMOMKit. It checks that:
 * - the default functions are zero, including all the moments;
 * - a single scalar expression applies to every scalar;
 * - one expression per scalar gives each scalar its own function, and the
 *   time of the functions is propagated to the time-dependent expressions;
 * - a number of scalar expressions that is neither 1 nor the number of
 *   scalars throws.
 */

#include <core/parameters_multiphysics.h>

#include <solvers/qmom_field_function.h>

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
 * @brief Print all the components of a QMOM field function at a point, with
 * value() and with vector_value().
 *
 * @param[in] function Function to evaluate.
 * @param[in] point Evaluation point.
 */
template <int dim>
void
print_components(const QMOMFieldFunction<dim> &function,
                 const Point<dim>             &point)
{
  Vector<double> values(function.n_components);
  function.vector_value(point, values);
  for (unsigned int component = 0; component < function.n_components;
       ++component)
    deallog << "  component " << component << ": "
            << function.value(point, component) << " " << values[component]
            << std::endl;
}

template <int dim>
void
test()
{
  Point<dim> point;
  point[0] = 0.25;

  // Default functions: zero, with 2 scalars and 4 moments
  {
    const Parameters::QMOMFunctions<dim> functions = parse_functions<dim>("");
    const QMOMFieldFunction<dim>         function(functions, 2, 4, "test");
    deallog << "Default function components: " << function.n_components
            << std::endl;
    print_components(function, point);
  }

  // A single expression applied to the 3 scalars
  {
    const Parameters::QMOMFunctions<dim> functions =
      parse_functions<dim>("  set Function expression = 3 * x\n");
    const QMOMFieldFunction<dim> function(functions, 3, 2, "test");
    deallog << "Single scalar expression:" << std::endl;
    print_components(function, point);
  }

  // One time-dependent expression per scalar
  {
    const Parameters::QMOMFunctions<dim> functions =
      parse_functions<dim>("  set Function expression = 1 + x; 2 + t\n");
    QMOMFieldFunction<dim> function(functions, 2, 2, "test");
    function.set_time(0.5);
    deallog << "Per-scalar function at t = " << function.get_time() << ":"
            << std::endl;
    print_components(function, point);
  }

  // Two scalar expressions for 3 scalars
  {
    const Parameters::QMOMFunctions<dim> functions =
      parse_functions<dim>("  set Function expression = 1; 2\n");
    try
      {
        const QMOMFieldFunction<dim> function(functions, 3, 2, "test");
        deallog << "Wrong number of scalar expressions did not throw"
                << std::endl;
      }
    catch (const std::exception &)
      {
        deallog << "Wrong number of scalar expressions threw" << std::endl;
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
