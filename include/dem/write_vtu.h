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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include "visualization.h"

#ifndef DEM_WRITEVTU_H_
#  define DEM_WRITEVTU_H_

class WriteVTU
{
public:
  WriteVTU();
  void
  write_master_files(const Visualization &           data_out,
                     const std::string &             solution_file_prefix,
                     const std::vector<std::string> &filenames);
  void
  writeVTUFiles(int, float, Visualization data_out);
};

#endif /* CMAKEFILES_WRITEVTU_H_ */
