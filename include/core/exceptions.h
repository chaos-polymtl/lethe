// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/*
 * This file defines the exceptions common to multiple solvers in Lethe
 * to avoid code repetition.
 */

#ifndef lethe_exceptions_h
#define lethe_exceptions_h

#include <deal.II/base/exceptions.h>

using namespace dealii;

DeclException3(SolidWarning,
               unsigned int,
               std::string,
               std::string,
               << "'number of solids = " << arg1 << "' but " << arg2
               << " solver does not support nitsche restriction. Use " << arg3
               << " solver instead.");

DeclException1(
  NoSolidWarning,
  std::string,
  << "No solid defined: impossible to assemble nitsche restriction in " << arg1
  << " solver. Change the 'number of solids' parameter.");

#endif
