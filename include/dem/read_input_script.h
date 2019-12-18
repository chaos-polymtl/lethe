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

#ifndef READINPUTSCRIPT_H_
#define READINPUTSCRIPT_H_

class ReadInputScript
{
public:
  float diameter, density, CoR, kn, ethan, kt, ethat, mu, mass;
  float ins_x_min, ins_y_min, ins_z_min, ins_x_max, ins_y_max, ins_z_max;
  float x_min, y_min, z_min, x_max, y_max, z_max;
  float dt;
  int   tInsertion, tFinal, nTotal, nInsert, insertFrequncy, writeFrequency;
  int   numOfParams;
  float g[3];

  ReadInputScript();
};

#endif /* READINPUTSCRIPT_H_ */
