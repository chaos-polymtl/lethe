/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019
 */

#ifndef lethe_pvd_handler_h
#define lethe_pvd_handler_h

#include "parameters.h"

using namespace dealii;

class PVDHandler
{
public:
  std::vector<std::pair<double, std::string>> times_and_names_;
  void
  save(std::string filename);
  void
  read(std::string filename);
  void
  append(double time, std::string pvtu_filename)
  {
    times_and_names_.push_back(
      std::pair<double, std::string>(time, pvtu_filename));
  }
  unsigned
  size()
  {
    return times_and_names_.size();
  }
};

#endif
