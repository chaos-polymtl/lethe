// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests check_file_exists, which verifies that a file that Lethe must
 * read exists and can be opened before the reader of that file is invoked.
 *
 * The four situations that the function distinguishes are covered: a file that
 * can't be read, a file name that was left empty, a file that does not exist
 * and a directory given where a file is expected. The last case of the test
 * verifies that `get_dimension`, which is the first function to open the
 * parameter file in every application, now reports a parameter file that does
 * not exist instead of raising a generic input/output error.
 */

// Deal.II
#include <deal.II/base/exceptions.h>

// Lethe
#include <core/utilities.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

/**
 * @brief Call check_file_exists and log whether it raised an exception.
 *
 * The message of the exception is deliberately not logged: it contains the
 * name and the line of the source file in which the exception was raised,
 * which would make this test fail whenever utilities.cc is edited.
 *
 * @param[in] case_name Name of the case, used to identify the line of output.
 *
 * @param[in] file_name Path of the file to check.
 */
void
check_and_log(const std::string &case_name, const std::string &file_name)
{
  try
    {
      check_file_exists(file_name, "file used by the test");
      deallog << case_name << ": no exception raised" << std::endl;
    }
  catch (const ExceptionBase &)
    {
      deallog << case_name << ": exception raised" << std::endl;
    }
}

void
test()
{
  // A file that exists and can be read.
  const std::string existing_file_name = "an_existing_file.prm";
  {
    std::ofstream file(existing_file_name);
    file << "set dimension = 2\n";
  }
  check_and_log("Existing file", existing_file_name);

  // A file name that was left empty in the parameter file.
  check_and_log("Empty file name", "");

  // A file that does not exist. This is the case a user hits when the name of
  // the parameter file is mistyped.
  check_and_log("Missing file", "a_file_that_does_not_exist.prm");

  // A directory given where a file is expected. The directory exists, so it is
  // the attempt to open it that must be reported.
  const std::string directory_name = "a_directory";
  std::filesystem::create_directory(directory_name);
  check_and_log("Directory", directory_name);

  // get_dimension is the first function to open the parameter file in every
  // application. It must report a file that does not exist rather than the
  // generic input/output error it used to raise.
  try
    {
      get_dimension("a_file_that_does_not_exist.prm");
      deallog << "Missing parameter file: no exception raised" << std::endl;
    }
  catch (const ExceptionBase &)
    {
      deallog << "Missing parameter file: exception raised" << std::endl;
    }

  // The dimension is still read from a parameter file that exists.
  deallog << "Dimension of an existing parameter file: "
          << get_dimension(existing_file_name) << std::endl;

  std::remove(existing_file_name.c_str());
  std::filesystem::remove(directory_name);
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
