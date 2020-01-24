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
 * Author: Bruno Blais, Polytechnique Montreal, 2020-
 */


/* This enum class is reponsible for the handling the specific indexes of the
 * particles properties within the property pool
 *
 */
namespace DEM
{
  // This is the only part where we should use a classical enum because it is
  // used as an integer within the code I think this is a temporary solution for
  // now
  enum PropertiesIndex : int
  {
    v_x     = 7,
    v_y     = 8,
    v_z     = 9,
    acc_x   = 10,
    acc_y   = 11,
    acc_z   = 12,
    force_x = 13,
    force_y = 14,
    force_z = 15,
    mass    = 19
  };
} // namespace DEM
