/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

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
