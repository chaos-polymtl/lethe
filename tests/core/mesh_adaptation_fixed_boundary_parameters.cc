// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the parsing of the parameters used to fix the mesh refinement at
 * the boundaries in the mesh adaptation subsection.
 *
 * Enabling 'fix boundary refinement' without populating 'boundaries fixed'
 * used to be silently ignored: the loop that clears the refinement and
 * coarsening flags of the boundary cells never matches any boundary id, so the
 * boundary refinement was progressively coarsened away by the error estimator.
 * The parsing must now throw in that case.
 *
 * Three scenarios are tested:
 *  1. The feature is disabled and no boundary is listed  -> valid
 *  2. The feature is enabled and boundaries are listed   -> valid
 *  3. The feature is enabled and no boundary is listed   -> must throw
 */

// Deal.II includes
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>

// Lethe
#include <core/parameters.h>

#include <string>

// Tests (with common definitions)
#include <../tests/tests.h>

/**
 * @brief Declare and parse a mesh adaptation subsection built from the given
 * values of the two parameters of interest.
 *
 * @param[in] fix_boundary_refinement Value given to the 'fix boundary
 * refinement' parameter.
 *
 * @param[in] boundaries_fixed Value given to the 'boundaries fixed' parameter.
 *
 * @return The parsed mesh adaptation parameters.
 */
Parameters::MeshAdaptation
parse_mesh_adaptation(const std::string &fix_boundary_refinement,
                      const std::string &boundaries_fixed)
{
  ParameterHandler prm;
  Parameters::MeshAdaptation::declare_parameters(prm);

  prm.parse_input_from_string("subsection mesh adaptation\n"
                              "  set type                    = adaptive\n"
                              "  set fix boundary refinement = " +
                              fix_boundary_refinement +
                              "\n"
                              "  set boundaries fixed        = " +
                              boundaries_fixed +
                              "\n"
                              "end\n");

  Parameters::MeshAdaptation mesh_adaptation;
  mesh_adaptation.parse_parameters(prm);

  return mesh_adaptation;
}

/**
 * @brief Parse the given combination of parameters and report whether the
 * parsing succeeded or threw.
 *
 * @param[in] label Description of the scenario printed in the log.
 *
 * @param[in] fix_boundary_refinement Value given to the 'fix boundary
 * refinement' parameter.
 *
 * @param[in] boundaries_fixed Value given to the 'boundaries fixed' parameter.
 */
void
check(const std::string &label,
      const std::string &fix_boundary_refinement,
      const std::string &boundaries_fixed)
{
  deallog << "--- " << label << " ---" << std::endl;

  try
    {
      const Parameters::MeshAdaptation mesh_adaptation =
        parse_mesh_adaptation(fix_boundary_refinement, boundaries_fixed);

      deallog << "  parsing succeeded" << std::endl;
      deallog << "  is boundary refinement fixed : "
              << mesh_adaptation.is_boundary_refinement_fixed << std::endl;
      deallog << "  number of fixed boundaries   : "
              << mesh_adaptation.boundaries_to_fix.size() << std::endl;
      for (const int boundary_id : mesh_adaptation.boundaries_to_fix)
        deallog << "  fixed boundary id            : " << boundary_id
                << std::endl;
    }
  // The message of the exception is not printed since it contains the file
  // name and the line number, which are not portable across builds.
  catch (const ExceptionBase &)
    {
      deallog << "  parsing threw an exception" << std::endl;
    }
}

int
main()
{
  try
    {
      initlog();

      deallog << "Beginning fixed boundary refinement parameter tests"
              << std::endl;

      check("Scenario 1: feature disabled, no boundary listed", "false", "");
      check("Scenario 2: feature enabled, boundaries listed",
            "true",
            "0, 1, 2, 3");
      check("Scenario 3: feature enabled, no boundary listed", "true", "");

      deallog << "OK" << std::endl;
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
