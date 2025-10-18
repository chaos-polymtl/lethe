// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/pvd_handler.h>

#include <fstream>

using namespace dealii;

void
PVDHandler::save(const std::string &prefix) const
{
  std::string   filename = prefix + ".pvdhandler";
  std::ofstream output(filename.c_str());
  output << times_and_names.size() << std::endl;
  output << "Time File" << std::endl;
  for (const auto &time_name : times_and_names)
    {
      output << time_name.first << " " << time_name.second << std::endl;
    }
}

void
PVDHandler::read(const std::string &prefix)
{
  times_and_names.clear();
  std::string   filename = prefix + ".pvdhandler";
  std::ifstream input(filename.c_str());
  AssertThrow(input, ExcFileNotOpen(filename));

  std::string  buffer;
  unsigned int size;
  input >> size;

  std::getline(input, buffer);
  std::getline(input, buffer);
  for (unsigned int i = 0; i < size; ++i)
    {
      double      time;
      std::string filename_i;
      input >> time;
      input >> filename_i;
      append(time, filename_i);
    }

  if (size != times_and_names.size())
    throw std::runtime_error("Error when reading pvd restart file ");
}
