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

#include <core/parameters.h>

using namespace dealii;

/**
 * @brief The PVDHandler class manages the storage of the information required to write pvd file.
 * It is particularly necessary for simulations which must be restarted from
 * checkpoint since it manages the continuity of the pvd output.
 */
class PVDHandler
{
public:
  /**
   * @brief save Saves the content of the pvd times_and_names to a file
   *
   * @param filename Name of the file to which the PVDHandler content is save
   */
  void
  save(std::string filename);

  /**
   * @brief read Reads the content of a pvd times_and_names checpoint
   *
   * @param filename Name of the file frin which the PVDHandler content is read
   */
  void
  read(std::string filename);

  void
  append(double time, std::string pvtu_filename)
  {
    times_and_names.push_back(
      std::pair<double, std::string>(time, pvtu_filename));
  }

  // The name of the pvtu files and the time associated with each file
  // is stored in a simple vector of pairs.
  std::vector<std::pair<double, std::string>> times_and_names;
};

#endif
