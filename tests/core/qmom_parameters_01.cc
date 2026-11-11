// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the logic of the parsing of the QMOM subsection and of the
 * functions prescribing the QMOM fields. It checks that:
 * - empty scalar names give generic names, one per scalar;
 * - an odd number of moments uses one more quadrature node;
 * - a number of scalar names that differs from the number of scalars and an
 *   invalid scalar name throw;
 * - the number of expressions of the scalars and of the particle size
 *   distribution is parsed, and a number of distribution expressions that
 *   does not match the distribution type throws.
 */

#include <core/parameters_multiphysics.h>
#include <core/utilities.h>

#include <deal.II/base/parameter_handler.h>

#include <../tests/tests.h>

#include <exception>
#include <functional>
#include <string>

/**
 * @brief Declare and parse a QMOM subsection.
 *
 * @param[in] input Content of the QMOM subsection.
 *
 * @return The parsed QMOM parameters.
 */
Parameters::QMOM
parse_qmom(const std::string &input)
{
  ParameterHandler       prm;
  const Parameters::QMOM declarations;
  declarations.declare_parameters(prm);

  prm.parse_input_from_string("subsection QMOM\n" + input + "end\n");

  Parameters::QMOM qmom;
  qmom.parse_parameters(prm);
  return qmom;
}

/**
 * @brief Print the scalar names and the numbers of moments and nodes.
 *
 * @param[in] qmom QMOM parameters.
 */
void
print_qmom(const Parameters::QMOM &qmom)
{
  deallog << "scalar names:";
  for (const std::string &scalar_name : qmom.scalar_names)
    deallog << " " << scalar_name;
  deallog << std::endl;
  deallog << "number of moments: " << qmom.number_of_moments
          << ", number of nodes: " << qmom.get_number_of_nodes() << std::endl;
}

/**
 * @brief Declare and parse the functions prescribing the QMOM fields.
 *
 * @param[in] input Content of the subsection holding the functions.
 *
 * @return The parsed functions.
 */
Parameters::QMOMFunctions<2>
parse_functions(const std::string &input)
{
  ParameterHandler prm;
  prm.enter_subsection("functions");
  Parameters::QMOMFunctions<2>::declare_parameters(prm);
  prm.leave_subsection();

  prm.parse_input_from_string("subsection functions\n" + input + "end\n");

  Parameters::QMOMFunctions<2> functions;
  prm.enter_subsection("functions");
  functions.parse_parameters(prm);
  prm.leave_subsection();
  return functions;
}

/**
 * @brief Print the distribution type and the number of expressions of the
 * functions prescribing the QMOM fields.
 *
 * @param[in] functions Functions to print.
 */
void
print_functions(const Parameters::QMOMFunctions<2> &functions)
{
  deallog << "scalar expressions: " << functions.scalars->n_components
          << ", distribution type: "
          << enum_to_string(functions.distribution_type)
          << ", distribution expressions: "
          << functions.distribution->n_components << std::endl;
}

/**
 * @brief Run a parsing and report whether it threw.
 *
 * @param[in] label Description of the scenario printed in the log.
 * @param[in] parse Parsing to run.
 */
void
check_throws(const std::string &label, const std::function<void()> &parse)
{
  try
    {
      parse();
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
  deallog << "Generic scalar names, even number of moments" << std::endl;
  print_qmom(parse_qmom("  set number of scalars = 3\n"
                        "  set number of moments = 4\n"));

  deallog << "Given scalar names, odd number of moments" << std::endl;
  print_qmom(parse_qmom("  set number of scalars = 2\n"
                        "  set scalar names      = ethanol, lipids\n"
                        "  set number of moments = 5\n"));

  check_throws("Too few scalar names", [] {
    parse_qmom("  set number of scalars = 2\n"
               "  set scalar names      = ethanol\n");
  });
  check_throws("Invalid scalar name", [] {
    parse_qmom("  set number of scalars = 1\n"
               "  set scalar names      = free lipids\n");
  });

  deallog << "Default functions" << std::endl;
  print_functions(parse_functions(""));

  deallog << "Lognormal distribution" << std::endl;
  print_functions(
    parse_functions("  subsection scalars\n"
                    "    set Function expression = 1; x\n"
                    "  end\n"
                    "  subsection distribution\n"
                    "    set type                = lognormal\n"
                    "    set Function expression = 1e18; log(1e-7); 0.3\n"
                    "  end\n"));

  deallog << "Chi squared distribution" << std::endl;
  print_functions(parse_functions("  subsection distribution\n"
                                  "    set type                = chi_squared\n"
                                  "    set Function expression = 1e18; 3\n"
                                  "  end\n"));

  check_throws("Wrong number of distribution expressions", [] {
    parse_functions("  subsection distribution\n"
                    "    set type                = gamma\n"
                    "    set Function expression = 1e18; 2\n"
                    "  end\n");
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
