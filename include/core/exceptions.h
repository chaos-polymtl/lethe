// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file exceptions.h
 * @brief Common exception declarations shared across multiple solvers in Lethe.
 *
 * This file defines deal.II-style exception macros used to report errors
 * related to the Nitsche immersed boundary method across different solver
 * types, avoiding code repetition.
 */

#ifndef lethe_exceptions_h
#define lethe_exceptions_h

#include <deal.II/base/exceptions.h>

using namespace dealii;

/**
 * @brief Exception raised when a solver does not support Nitsche solid
 * restriction but solids are defined in the parameter file.
 *
 * @param arg1 Number of solids defined in the parameter file.
 * @param arg2 Name of the current solver that does not support Nitsche.
 * @param arg3 Name of the solver that should be used instead.
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
 * @param arg1 Name of the solver in which the error occurred.
 */
DeclException1(
  NoSolidWarning,
  std::string,
  << "No solid defined: impossible to assemble nitsche restriction in " << arg1
  << " solver. Change the 'number of solids' parameter.");

#endif
