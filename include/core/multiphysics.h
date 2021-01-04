/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

/*
 * This file defines a small enum which is used to identify
 * the various physics solverd within Lethe
 */


#ifndef lethe_multiphysics_h
#define lethe_multiphysics_h

enum class PhysicsID : int
{
  fluid_dynamics = 0,
  heat_transfer  = 1
};

#endif
