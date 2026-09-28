// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the declaration and the parsing of the parameters of the
 * Anderson-Jackson filter.
 *
 * Four scenarios are tested:
 *  1. The subsection is empty                                   -> defaults
 *  2. Every parameter is given a non-default value              -> valid
 *  3. The filter width is zero                                  -> must throw
 *  4. A top-hat kernel with a zero gaussian cutoff, which is
 *     only used by the gaussian kernel                          -> valid
 */

// Deal.II includes
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>

// Lethe
#include <core/parameters.h>
#include <core/parameters_cfd_dem.h>

#include <iomanip>
#include <sstream>
#include <string>

// Tests (with common definitions)
#include <../tests/tests.h>

namespace
{
  /**
   * @brief Format a number in scientific notation.
   *
   * @param[in] value Number to format.
   *
   * @return Formatted number.
   */
  std::string
  format(const double value)
  {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(6) << value;
    return stream.str();
  }

  /**
   * @brief Format a boolean as true or false.
   *
   * @param[in] value Boolean to format.
   *
   * @return Formatted boolean.
   */
  std::string
  format(const bool value)
  {
    return value ? "true" : "false";
  }

  /**
   * @brief Declare and parse an anderson jackson filter subsection.
   *
   * @param[in] entries Entries of the subsection, one "set ... = ..." per
   * line.
   *
   * @return The parsed parameters.
   */
  Parameters::AndersonJacksonFilter
  parse_filter_parameters(const std::string &entries)
  {
    ParameterHandler prm;
    Parameters::AndersonJacksonFilter::declare_parameters(prm);

    prm.parse_input_from_string("subsection anderson jackson filter\n" +
                                entries + "end\n");

    Parameters::AndersonJacksonFilter filter_parameters;
    filter_parameters.parse_parameters(prm);

    return filter_parameters;
  }

  /**
   * @brief Parse the given entries and print the parsed parameters, or report
   * that the parsing threw.
   *
   * @param[in] label Description of the scenario printed in the log.
   *
   * @param[in] entries Entries of the subsection, one "set ... = ..." per
   * line.
   */
  void
  check(const std::string &label, const std::string &entries)
  {
    deallog << "--- " << label << " ---" << std::endl;

    try
      {
        const Parameters::AndersonJacksonFilter p =
          parse_filter_parameters(entries);

        deallog << "  kernel type                    : "
                << (p.kernel_type == Parameters::FilterKernelType::gaussian ?
                      "gaussian" :
                      "top-hat")
                << std::endl;
        deallog << "  filter width                   : "
                << format(p.filter_width) << std::endl;
        deallog << "  gaussian cutoff                : "
                << format(p.gaussian_cutoff) << std::endl;
        deallog << "  output polynomial degree       : " << p.output_degree
                << std::endl;
        deallog << "  quadrature points              : "
                << p.n_quadrature_points << std::endl;
        deallog << "  cut cell subdivisions          : "
                << p.cut_cell_subdivisions << std::endl;
        deallog << "  filter pressure                : "
                << format(p.filter_pressure) << std::endl;
        deallog << "  minimum fluid fraction         : "
                << format(p.minimum_fluid_fraction) << std::endl;
        deallog << "  normalize at domain boundaries : "
                << format(p.normalize_at_domain_boundaries) << std::endl;
        deallog << "  output folder                  : " << p.output_folder
                << std::endl;
        deallog << "  output name                    : " << p.output_name
                << std::endl;
        deallog << "  verbosity                      : "
                << Parameters::to_string(p.verbosity) << std::endl;
      }
    // The message of the exception is not printed since it contains the file
    // name and the line number, which are not portable across builds.
    catch (const ExceptionBase &)
      {
        deallog << "  parsing threw an exception" << std::endl;
      }
  }
} // namespace

int
main()
{
  try
    {
      initlog();

      check("Scenario 1: empty subsection", "");
      check("Scenario 2: non-default values",
            "  set kernel type                    = top-hat\n"
            "  set filter width                   = 0.25\n"
            "  set gaussian cutoff                = 4.0\n"
            "  set output polynomial degree       = 2\n"
            "  set quadrature points              = 3\n"
            "  set cut cell subdivisions          = 8\n"
            "  set filter pressure                = false\n"
            "  set minimum fluid fraction         = 1e-6\n"
            "  set normalize at domain boundaries = true\n"
            "  set output folder                  = ./ajf_output/\n"
            "  set output name                    = ajf\n"
            "  set verbosity                      = extra verbose\n");
      check("Scenario 3: zero filter width", "  set filter width = 0\n");
      check("Scenario 4: top-hat kernel with a zero gaussian cutoff",
            "  set kernel type     = top-hat\n"
            "  set gaussian cutoff = 0\n");
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
