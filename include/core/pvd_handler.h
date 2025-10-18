// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
  save(const std::string &filename) const;

  /**
   * @brief read Reads the content of a pvd times_and_names checpoint
   *
   * @param filename Name of the file frin which the PVDHandler content is read
   */
  void
  read(const std::string &filename);

  void
  append(double time, const std::string &pvtu_filename)
  {
    times_and_names.emplace_back(time, pvtu_filename);
  }

  // The name of the pvtu files and the time associated with each file
  // is stored in a simple vector of pairs.
  std::vector<std::pair<double, std::string>> times_and_names;
};

#endif
