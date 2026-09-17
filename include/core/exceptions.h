// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file exceptions.h
 * @brief Common exception declarations shared across multiple solvers in Lethe.
 *
 * This file defines deal.II-style exception macros used to report errors that
 * are raised from more than one place in Lethe, such as the errors related to
 * the Nitsche immersed boundary method or to the files that Lethe reads,
 * avoiding code repetition.
 */

#ifndef lethe_exceptions_h
#define lethe_exceptions_h

#include <deal.II/base/exceptions.h>

using namespace dealii;

/**
 * @brief Exception raised when a solver does not support Nitsche solid
 * restriction but solids are defined in the parameter file.
 *
 * @param[in] arg1 Number of solids defined in the parameter file.
 * @param[in] arg2 Name of the current solver that does not support Nitsche.
 * @param[in] arg3 Name of the solver that should be used instead.
 */
DeclException3(SolidWarning,
               unsigned int,
               std::string,
               std::string,
               << "'number of solids = " << arg1 << "' but " << arg2
               << " solver does not support nitsche restriction. Use " << arg3
               << " solver instead.");

/**
 * @brief Exception raised when no solid is defined but the solver attempts
 * to assemble a Nitsche restriction.
 *
 * @param[in] arg1 Name of the solver in which the error occurred.
 */
DeclException1(
  NoSolidWarning,
  std::string,
  << "No solid defined: impossible to assemble nitsche restriction in " << arg1
  << " solver. Change the 'number of solids' parameter.");

/**
 * @brief Exception raised when no name was given for a file that Lethe must
 * read.
 *
 * @param[in] arg1 Description of the role of the file.
 */
DeclException1(EmptyFileName,
               std::string,
               << "No file name was given for the " << arg1
               << ". Specify it in the parameter file before running Lethe.");

/**
 * @brief Exception raised when a file that Lethe must read does not exist.
 *
 * @param[in] arg1 Path of the file, as given by the user.
 * @param[in] arg2 Description of the role of the file.
 */
DeclException2(
  FileDoesNotExist,
  std::string,
  std::string,
  << "The " << arg2 << " <" << arg1 << "> does not exist. "
  << "Verify that its name is spelled correctly. Relative paths are "
  << "interpreted from the directory in which the application was launched, "
  << "not from the directory of the parameter file.");

/**
 * @brief Exception raised when a file that Lethe must read exists but cannot
 * be opened for reading.
 *
 * @param[in] arg1 Path of the file, as given by the user.
 * @param[in] arg2 Description of the role of the file.
 */
DeclException2(
  FileIsNotReadable,
  std::string,
  std::string,
  << "The " << arg2 << " <" << arg1 << "> exists but could not be opened for "
  << "reading. Verify that it is a file and not a directory, and that its "
  << "permissions allow Lethe to read it.");

#endif
