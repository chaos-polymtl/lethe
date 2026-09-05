// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests Parameters::get_size_of_subsections, which pre-parses the size
 * of the variable size subsections of a parameter file before they are
 * declared.
 *
 * The size is the largest value given to a "number" entry anywhere in the
 * file, no matter the subsection it belongs to. Solvers for which a "boundary
 * conditions" subsection is mandatory require that at least one such entry is
 * present. The DEM solvers, which declare their boundary conditions in a "DEM
 * boundary conditions" subsection with a "number of boundary conditions"
 * entry, do not: a file without any "number" entry then yields a size of zero.
 */

// Deal.II
#include <deal.II/base/exceptions.h>

// Lethe
#include <core/parameters.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cstdio>
#include <fstream>
#include <string>

/**
 * @brief Write a parameter file in the working directory of the test.
 *
 * @param[in] file_name Name of the file to write.
 *
 * @param[in] content Content of the file.
 */
void
write_parameter_file(const std::string &file_name, const std::string &content)
{
  std::ofstream file(file_name);
  file << content;
}

void
test()
{
  // A parameter file of the CFD solvers. The size of both subsections is the
  // largest "number" of the file, here the number of boundary conditions.
  const std::string cfd_file_name = "cfd_parameters.prm";
  write_parameter_file(cfd_file_name,
                       "subsection boundary conditions\n"
                       "  set number = 3\n"
                       "end\n"
                       "subsection manifolds\n"
                       "  set number = 1\n"
                       "end\n");

  const Parameters::SizeOfSubsections cfd_sizes =
    Parameters::get_size_of_subsections(cfd_file_name);

  deallog << "CFD parameter file: " << cfd_sizes.boundary_conditions
          << " boundary conditions, " << cfd_sizes.manifolds << " manifolds"
          << std::endl;

  // A parameter file of the DEM solvers. It contains no "number" entry, since
  // "number of boundary conditions" is a different entry, and it declares no
  // manifold.
  const std::string dem_file_name = "dem_parameters.prm";
  write_parameter_file(dem_file_name,
                       "subsection DEM boundary conditions\n"
                       "  set number of boundary conditions = 2\n"
                       "end\n");

  const Parameters::SizeOfSubsections dem_sizes =
    Parameters::get_size_of_subsections(dem_file_name, false);

  deallog << "DEM parameter file: " << dem_sizes.boundary_conditions
          << " boundary conditions, " << dem_sizes.manifolds << " manifolds"
          << std::endl;

  // The same file is refused by the solvers for which the size of at least one
  // subsection must be given.
  try
    {
      Parameters::get_size_of_subsections(dem_file_name);
      deallog << "No exception raised for a mandatory subsection size"
              << std::endl;
    }
  catch (const ExceptionBase &)
    {
      deallog << "Exception raised for a mandatory subsection size"
              << std::endl;
    }

  std::remove(cfd_file_name.c_str());
  std::remove(dem_file_name.c_str());
}

int
main(int argc, char *argv[])
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
